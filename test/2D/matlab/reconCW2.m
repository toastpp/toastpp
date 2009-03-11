function reconCW2

% Sample code: reconstruction of absorption image from log amplitude data
% using a Gauss-Newton Krylov solver.
% This version:
%  - reconstructs for log mua parameters
%  - does not use regularisation
%  - uses noise-free data

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with Gauss-Newton solver')
disp('Data: CW intensity (log); param: absorption')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = '../meshes/ellips_tri10.msh';         % mesh file
qmname    = '../meshes/circle25_32x32.qm';        % source-detector file
data_name = '../fwdfem/fmod_ellips_32x32_CW.fem'; % data file: log amplitude

refind = 1.4;                        % refractive index
c0 = 0.3;                            % speed of light in vacuum [mm/ps]
cm = c0/refind;                      % speed of light in medium
grid = [100 100];                    % solution basis: grid dimension
tolGN = 1e-6;                        % Gauss-Newton convergence criterion
tolKrylov = 1e-2;                    % Krylov convergence criterion
itrmax = 100;                        % Gauss-Newton max. iterations
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();

% The target image. Only for visual comparison, and not used in any way
% by the reconstruction.
btgt = loadtgt(grid);
subplot(2,2,1); imagesc(reshape(btgt,grid(1),grid(2)),[0.01 0.055]);
axis equal tight; colorbar; title('target');

% Read a TOAST mesh definition from file.
hMesh = toastReadMesh (meshname);        % load the mesh
toastReadQM (hMesh, qmname);             % load the source/detector specs
n = toastMeshNodeCount (hMesh);          % node count

% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;                 % absorption coeff.
mus = ones(n,1) * 2;                     % scattering coeff.
ref = ones(n,1) * refind;                % refractive index

% Read the data
data = toastReadRealVector(data_name);

% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis (hMesh, grid); % basis mapper handle
solmask = toastSolutionMask (hBasis);    % mask for removing exterior pixels

% Generate source vectors
qvec = real(toastQvec (hMesh, 'Neumann', 'Gaussian', 2));

% Generate measurement vectors
mvec = real(toastMvec (hMesh, 'Gaussian', 2));

% Initial data set f[x0]
proj = forwardModel(mua);

% Initial parameter estimates in solution basis
% Note that the solution parameters are given in terms of log(c*mua)
% where c is the speed of light, and mua is the absorption coefficient.
bcmua = toastMapMeshToBasis (hBasis, mua * cm);
scmua = bcmua(solmask);
bmua = zeros(size(bcmua));
dmua = toastMapMeshToSol(hBasis,mua)-btgt(solmask); % initial image residual
x = scmua;
logx = log(x);   % initial solution vector

% initial data error
err0 = objective(proj);
err = err0;                                        % current error
errp = inf;                                        % previous error
erri(1) = err0;
errm(1) = norm(dmua);
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

step = 1.0;      % initial step length for line search
itr = 0;         % iteration counter

% Gauss-Newton loop
while (itr < itrmax) && (err > tolGN*err0) && (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    fprintf (1,'Calculating Jacobian\n');
    J = toastJacobianCW (hMesh, hBasis, qvec, mvec, mua, mus, ref, 'direct');
    
    %J = J * diag(x);      % parameter normalisation
    for i = 1:size(x,1)
        J(:,i) = J(:,i) * x(i,1);
    end
    
    % Normalisation of Hessian
    for i = 1:size(J,2)
        M(i) = sum(J(:,i) .* J(:,i));
        M(i) = 1 ./ sqrt(M(i));
    end
    %J = J * diag(M);
    for i = 1:size(M,2)
        J(:,i) = J(:,i) * M(1,i);
    end
    
    % Gradient of cost function
    r = J' * (data-proj);
    
    % Update with implicit Krylov solver
    fprintf (1, 'Entering Krylov solver\n');
    dx = krylov(r);
    
    % Line search
    fprintf (1, 'Entering line search\n');
    [step, err] = toastLineSearch (logx, dx, step, err, @objective2);
    
    % Add update to solution
    logx = logx + dx*step;
    x = exp(logx);
    
    % Map parameters back to mesh
    scmua = x;
    smua = scmua/cm;
    mua = toastMapSolToMesh (hBasis, smua);
    bmua(solmask) = smua;
    dmua = smua-btgt(solmask); % image residual
    
    proj = forwardModel(mua);
    
    % update objective function
    err = objective(proj);
    itr = itr+1;
    erri(itr+1) = err;
    errm(itr+1) = norm(dmua);
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);

    subplot(2,2,2);imagesc(reshape(bmua,grid(1),grid(2)),[0.01 0.055]);
    axis equal tight, colorbar; title('reconstruction');
    subplot(2,2,3);semilogy([0:itr],erri);axis tight;
    xlabel('iteration');title('posterior')
    subplot(2,2,4);semilogy([0:itr],errm);axis tight;
    xlabel('iteration');title('image residual');
    drawnow
    
end


    % =====================================================================
    % forward model: maps parameters to boundary measurement data
    
    function proj = forwardModel(mua)
    smat = real(toastSysmat (hMesh, mua, mus, ref, 0)); % system matrix
    proj = reshape (log(mvec.' * (smat\qvec)), [], 1);  % source and measurement operators
    end

    
    % =====================================================================
    % returns objective function, given model boundary data proj
    
    function of = objective(proj)
    of = sum((data-proj).^2);
    end


    % =====================================================================
    % returns objective function, given model parameters (used as callback
    % function by toastLineSearch)
    
    function of = objective2(logx)
    x = exp(logx);                          % undo parameter transformation
    smua = x/cm;                            % map from x (=c*mua) to mua
    mua = toastMapSolToMesh (hBasis, smua); % map to mesh basis
    of = objective(forwardModel(mua));      % calc objective function
    end
    

    % =====================================================================
    % Krylov solver subroutine
    
    function dx = krylov(r)
    dx = gmres (@jtjx, r, 30, tolKrylov, 100);
    end
    

    % =====================================================================
    % Callback function for matrix-vector product (called by krylov)
    
    function b = jtjx(x)
    b = J' * (J*x);
    end

end

% =====================================================================
% Load the target image (for comparison only)

function btgt = loadtgt(grid)
hMesh = toastReadMesh('../meshes/ellips_tri10.msh');
hBasis = toastSetBasis(hMesh,grid);
htgt = toastReadNIM('../meshes/tgt_mua_ellips_tri10.nim');
btgt = toastMapMeshToGrid(hBasis,htgt);
toastDeleteBasis(hBasis);
toastDeleteMesh(hMesh);
end
