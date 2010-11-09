function recon1

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with Gauss-Newton solver')
disp('using edge prior information from correct target image.')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = '../meshes/ellips_tri10.msh';              % mesh file
qmname    = '../meshes/circle25_32x32.qm';            % source-detector file
fmod_name = '../fwdfem/fmod_ellips_32x32_100MHz.fem'; % data file: log amplitude
farg_name = '../fwdfem/farg_ellips_32x32_100MHz.fem'; % data file: phase

refind = 1.4;                           % refractive index
bx = 100; by = 100;                     % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
tau = 1e-4;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolGN = 1e-6;                           % Gauss-Newton convergence criterion
tolKrylov = 1e-2;                       % Krylov convergence criterion
itrmax = 100;                           % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian matrix
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();
%toastSetVerbosity(1);

% Set up some variables
blen = bx*by;
c0 = 0.3;
cm = c0/refind;

% Read a TOAST mesh definition from file.
hMesh = toastReadMesh (meshname);
toastReadQM (hMesh, qmname);
n = toastMeshNodeCount (hMesh);
dmask = toastDataLinkList (hMesh);

% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;
mus = ones(n,1) * 2;
ref = ones(n,1) * refind;
kap = 1./(3*(mua+mus));

% Read the data
mdata = toastReadRealVector(fmod_name);
pdata = toastReadRealVector(farg_name);

mdata = mdata + mdata.*0.01.*randn(size(mdata));
pdata = pdata + pdata.*0.01.*randn(size(pdata));

data = [mdata;pdata];
m = length(data);

% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis ('LINEAR', hMesh, [bx by]);
solmask = toastSolutionMask (hBasis);
solmask2 = [solmask solmask+blen];

% Generate source vectors
qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);

% Generate measurement vectors
mvec = toastMvec (hMesh, 'Gaussian', 2);

% Initial data set f[x0]
smat = toastSysmat (hMesh, mua, mus, ref, freq);
lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);
lgamma = lgamma(dmask);
proj = [real(lgamma);imag(lgamma)];

% data scaling
sd_lnmod = ones(size(lgamma)) * norm(mdata-real(lgamma)); 
sd_phase = ones(size(lgamma)) * norm(pdata-imag(lgamma));
sd = [sd_lnmod;sd_phase];

% initial parameter estimates in solution basis
bmua = toastMapMeshToBasis (hBasis, mua);
bmus = toastMapMeshToBasis (hBasis, mus);
bmua_itr(1,:) = bmua;
bmus_itr(1,:) = bmus;
bkap = toastMapMeshToBasis (hBasis, kap);
bcmua = bmua*cm;
bckap = bkap*cm;
scmua = bcmua(solmask);
sckap = bckap(solmask);

x = [scmua;sckap];
logx = log(x);
p = length(x);
step = 1.0; % initial step length for line search

% Initialise regularisation
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
%hReg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta);

tgtmua=toastReadNIM('../meshes/tgt_mua_ellips_tri10.nim');
tgtmus=toastReadNIM('../meshes/tgt_mus_ellips_tri10.nim');
tgtkap=1./(3*(tgtmua+tgtmus));
hMesh2 = toastReadMesh('../meshes/ellips_tri10.msh');
hBasis2 = toastSetBasis ('LINEAR',hMesh2,[bx by]);
btgtmua=toastMapMeshToGrid (hBasis2, tgtmua);
btgtkap=toastMapMeshToGrid (hBasis2, tgtkap);
toastDeleteBasis(hBasis2);
toastDeleteMesh(hMesh2);
btgt = [btgtmua; btgtkap];
hReg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta, 'KapRefImage', ...
                   btgt, 'KapRefScale', 1, 'KapRefPMThreshold', 0.1);

% initial data error (=2 due to data scaling)
err0 = privObjective (proj, data, sd, hReg, logx);  %initial error
err = err0;                                         % current error
errp = inf;                                         % previous error
erri(1) = err0;
itr = 1; % iteration counter
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

% Gauss-Newton loop
while (itr <= itrmax) & (err > tolGN*err0) & (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    fprintf (1,'Calculating Jacobian\n');
    J = toastJacobian (hMesh, hBasis, qvec, mvec, mua, mus, ref, freq, 'bicgstab', 1e-16);
    size(J)
    pause 
    % data normalisation
    for i = 1:m
        J(i,:) = J(i,:) / sd(i);
    end
    %J = spdiags(1./sd,0,m,m) * J;

    % parameter normalisation
    for i = 1:p
        J(:,i) = J(:,i) * x(i);
    end
    %J = J * spdiags (x,0,p,p);
    
    % Normalisation of Hessian
    psiHdiag = toastRegulHDiag (hReg, logx);
    for i = 1:size(J,2)
        M(i) = sum(J(:,i) .* J(:,i));
        M(i) = M(i) + psiHdiag(i);
        M(i) = 1 ./ sqrt(M(i));
    end
    for i = 1:p
        J(:,i) = J(:,i) * M(i);
    end
    %J = J * spdiags (M',0,p,p);
    
    % Gradient of cost function
    r = J' * ((data-proj)./sd);
    r = r - toastRegulGradient (hReg, logx) .* M';
    
    if Himplicit == true
        % Update with implicit Krylov solver
        fprintf (1, 'Entering Krylov solver\n');
        dx = toastKrylov (x, J, r, M, 0, hReg, tolKrylov);
    else
        % Update with explicit Hessian
        H = J' * J;
        lambda = 0.01;
        H = H + eye(size(H)).* lambda;
        dx = H \ r;
        clear H;
    end
    
    clear J;
    
    % Line search
    fprintf (1, 'Entering line search\n');
    [step, err] = toastLineSearch (logx, dx, step, err, @objective);
    
    % Add update to solution
    logx = logx + dx*step;
    x = exp(logx);
    
    % Map parameters back to mesh
    scmua = x(1:size(x)/2);
    sckap = x(size(x)/2+1:size(x));
    smua = scmua/cm;
    skap = sckap/cm;
    smus = 1./(3*skap) - smua;
    mua = toastMapSolToMesh (hBasis, smua);
    mus = toastMapSolToMesh (hBasis, smus);
    bmua(solmask) = smua;
    bmus(solmask) = smus;

    subplot(1,2,1), imagesc(reshape(bmua,bx,by),[min(mua) max(mua)]), axis equal, axis tight, colorbar('horiz'); colormap(gray);
    subplot(1,2,2), imagesc(reshape(bmus,bx,by),[min(mus) max(mus)]), axis equal, axis tight, colorbar('horiz'); colormap(gray);
    drawnow
    
    proj = toastProject (hMesh, mua, mus, ref, freq, qvec, mvec);

    % update objective function
    err = privObjective (proj, data, sd, hReg, logx);
    itr = itr+1;
    erri(itr) = err;
    bmua_itr(itr,:) = bmua;
    bmus_itr(itr,:) = bmus;
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);
end

dlmwrite('recon.dat', erri);
disp('recon1: finished')

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (hMesh, hBasis, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = privObjective (proj, data, sd, hReg, x);
    fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end

end
