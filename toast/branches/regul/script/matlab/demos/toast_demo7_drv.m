function toast_demo7

% ======================================================================
% MATLAB-TOAST sample script:
% Test implementation for spectral reconstruction
% Try to recover 2 chromophore concentrations from 2 wavelengths

toastCatchErrors();
disp('MATLAB-TOAST sample script:')
disp('Multispectral reconstruction')
disp('Reconstruct 2 chromophore concentrations from CW intensity')
disp('data at 2 wavelengths (no scatter reconstruction yet)')
disp('-------------------------------------------------------')

% ======================================================================
% Defining file paths and parameters

meshname  = 'circle25_32.msh';          % mesh file
%meshname  = 'ellips_tri10.msh';
qmname    = 'circle25_32x32.qm';        % source-detector file

refind = 1.4;                           % refractive index
c0 = 0.3;                               % speed of light invacuum [mm/ps]
cm = c0/refind;                         % speed of light in medium
grid = [80 80];                       % solution basis: grid dimension
tolGN = 1e-6;                           % Gauss-Newton convergence criterion
tolKrylov = 1e-2;                       % Krylov convergence criterion
itrmax = 100;                           % Gauss-Newton max. iterations

nch = 3;                                % number of chromophores
nlambda = 3;                            % number of wavelengths
extinct = [[0.01 0.02 0.025]; [0.02 0.01 0.01]; [0.015 0.025 0.01]]; % extinction coefficient for both chromophores at both wavelengths
                                        % index1: chromophore, index2: wavelength

% ======================================================================
% Initialisations

rand('state',0);
blen = prod(grid);

% ======================================================================
% Read a TOAST mesh definition from file.
hMesh = toastReadMesh (meshname);
n = toastMeshNodeCount (hMesh);

% ======================================================================
% Generate source and measurement vectors
toastReadQM (hMesh, qmname);
qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);
mvec = toastMvec (hMesh, 'Gaussian', 2);
dmask = toastDataLinkList (hMesh);
nqm = length(dmask);

% ======================================================================
% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis (hMesh, grid);
solmask = toastSolutionMask (hBasis);
solmask2 = [solmask solmask+blen];
slen = length(solmask);

% ======================================================================
% Fixed background parameters
bmus = ones(blen,1)*2;
bref = ones(blen,1)*refind;
mus = toastMapGridToMesh(hBasis,bmus);
ref = toastMapGridToMesh(hBasis,bref);

% ======================================================================
% target chromophore distributions
blobtgt = 5;
for c=1:nch
    img = ones(grid);
    nblob = 3+round(rand()*(blobtgt-3));
    for b=1:nblob
        rad(b) = (0.05+rand()*0.15)*grid(1);
        cont = true;
        while cont
            cont = false;
            cnt(b,:) = rand(1,2).*grid;
            r = norm(cnt(b,:)-grid/2);
            if r+rad(b) >= grid(1)/2, cont=true; end
            for bb=1:b-1
                if norm(cnt(b,:)-cnt(bb,:)) < rad(b)+rad(bb), cont=true; end
            end
        end
        mag(b) = rand()*2+1;
    end
    for j=1:grid(2)
        for i=1:grid(1)
            p = [i j];
            for b=1:nblob
                if norm(p-cnt(b,:)) <= rad(b)
                    img(i,j) = img(i,j) + mag(b);
                end
            end
        end
    end
    img = reshape(img,[],1);
    bC_tgt(c,:) = zeros(prod(grid),1);
    bC_tgt(c,solmask) = img(solmask) * 0.3;
    subplot(2,nch,c);imagesc(reshape(bC_tgt(c,:),grid(1),grid(2)));axis equal tight;title 'C target'; colorbar
end

% ======================================================================
% Generate target forward data
data = [];
for i=1:nlambda
    bmua = GetMua(extinct(:,i),bC_tgt);
    proj = Project (toastMapGridToMesh (hBasis,bmua));
    data = [data; proj];
    m0 = length(proj);
end
m = length(data);

% ======================================================================
% Set up homogeneous initial parameter estimates
for i=1:nch
    C(i,:) = ones(n,1) * 0.5;
end

% ======================================================================
% initial parameter estimates in solution basis
x = [];
for i=1:nch
    bC(i,:) = toastMapMeshToBasis(hBasis,C(i,:));
    sC(i,:) = bC(i,solmask);
    x = [x; sC(i,:)'];
end
p = length(x);

% ======================================================================
% Initial data set f[x0]
proj = ProjectAll(C);

% ======================================================================
% data scaling
for i=1:nlambda
    p = proj((i-1)*m0+1:i*m0);
    d = data((i-1)*m0+1:i*m0);
    nm = norm(p-d);
    sd((i-1)*m0+1:i*m0,1) = ones(m0,1)*nm;
end

% ======================================================================
% initial data error
err0 = objective (proj);           %initial error
err = err0;                        % current error
errp = inf;                        % previous error
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

hReg = 0;
step = 1.0; % initial step length for line search
itr = 0; % iteration counter

% Gauss-Newton loop
while (itr < itrmax) & (err > tolGN*err0) & (errp-err > tolGN)

    errp = err;
    
    % Build the Jacobian
    clear J;
    for i = 1:nlambda
        mua = GetMua(extinct(:,i),C);
        Ji = toastJacobianCW (hMesh, hBasis, qvec, mvec, mua, mus, ref, 'direct');
        for j = 1:nch
            J((i-1)*nqm+1:i*nqm,(j-1)*slen+1:j*slen) = Ji * extinct(j,i);
        end
        clear Ji;
    end
    
    % data normalisation
    for i = 1:m, J(i,:) = J(i,:) / sd(i); end
    
    % Normalisation of Hessian
    for i = 1:size(J,2)
        M(i) = 1 ./ sqrt(sum(J(:,i) .* J(:,i)));
    end
    for i = 1:size(M,2)
        J(:,i) = J(:,i) * M(1,i);
    end
    
    % Gradient of cost function
    r = J' * 2*((data-proj)./sd);

    % Update with implicit Krylov solver
    fprintf (1, 'Entering Krylov solver\n');
    dx = krylov(r);
    
    % Line search
    fprintf (1, 'Entering line search\n');
    [step, err] = toastLineSearch (x, dx, step, err, @objective2);
    
    % Add update to solution
    x = x + dx*step;
    
    % Map parameters back to mesh
    for i=1:nch
        sC(i,:) = x((i-1)*slen+1:i*slen);
        bC(i,solmask) = sC(i,:);
        C(i,:) = toastMapBasisToMesh (hBasis, bC(i,:));
    end

    for c=1:nch
        subplot(2,nch,nch+c), imagesc(reshape(bC(c,:),grid(1),grid(2))), axis equal tight, colorbar;title('C recon');
    end
    drawnow
    
    % update objective function
    itr = itr+1;
    proj = ProjectAll(C);
    err = objective(proj);
    errn(itr) = err;
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);
end


    % ======================================================================
    % forward model: maps mua to boundary measurement data
    
    function p_ = Project (mua)
    smat_ = real(toastSysmat (hMesh, mua, mus, ref, 0)); % system matrix
    p_ = reshape (log(mvec.' * (smat_\qvec)), [], 1);  % source and measurement operators
    p_ = full(p_);
    clear smat_;
    end


    % ======================================================================
    % forward model: maps chromophore concentrations at a given wavelength
    % (defined by extinction coefficients) to boundary measurement data
    
    function p_ = ProjectSingle (extinct, C)
    p_ = Project (GetMua (extinct, C));
    end


    % ======================================================================
    % forward model: maps chromophore concentrations at all wavelengths to
    % boundary measurement data
    
    function p_ = ProjectAll (C)
    for i_ = 1:nlambda
        p_((i_-1)*m0+1:i_*m0,1) = ProjectSingle (extinct(:,i_), C);
    end
    end

    
    % ======================================================================
    % Calculate mua given extinction coefficients and chromophore
    % concentrations
    
    function mua_ = GetMua (extinct, C)
    mua_ = zeros(size(C,2),1);
    for i_ = 1:nch
        mua_ = mua_ + C(i_,:)'*extinct(i_);
    end
    end

    % =====================================================================
    % returns objective function, given model boundary data proj
    
    function of_ = objective(proj)
    of_ = sum(((data-proj)./sd).^2);
    end


    % =====================================================================
    % returns objective function, given model parameters (used as callback
    % function by toastLineSearch)
    
    function of_ = objective2(x)
    if min(x) < 0    
        of_ = 1e8;                               % discourage negative solutions
    else
        for i_ = 1:nch
            C_(i_,:) = toastMapSolToMesh(hBasis,x((i_-1)*slen+1:i_*slen));
        end
        of_ = objective(ProjectAll(C_));        % calc objective function
    end
    end
    

    % =====================================================================
    % Krylov solver subroutine
    
    function dx_ = krylov(r)
    dx_ = gmres (@jtjx, r, 30, tolKrylov, 100);
    end


    % =====================================================================
    % Callback function for matrix-vector product (called by krylov)
    
    function b = jtjx(x)
    b = J' * (J*x);
    end

end