function toast_tut1

% ======================================================================
% Sample code: reconstruction of absorption image from log amplitude data
% using a Gauss-Newton Krylov solver.
% This version:
%  - reconstructs for linear parameters (no log transformation)
%  - does not use regularisation
%  - uses noise-free data

toastCatchErrors();
disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with Gauss-Newton solver')
disp('Data: CW intensity (log); param: absorption')
disp('-------------------------------------------------------')

% Make sure we are running from the right directory
eval(['cd ' ['''' fileparts(which('toast_tut1')) '''' '/../../../test/2D/matlab']]);

% ======================================================================
% Defining file paths and parameters

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
% The target image. Only for visual comparison, and not used in any way
% by the reconstruction.
btgt = loadtgt(grid);
subplot(2,2,1); imagesc(reshape(btgt,grid(1),grid(2)),[0.01 0.055]);
axis equal tight; colorbar; title('target');

% ======================================================================
% Read a TOAST mesh definition from file.
hMesh = toastReadMesh (meshname);        % load the mesh
n = toastMeshNodeCount (hMesh);          % node count

% ======================================================================
% Generate source and measurement vectors
toastReadQM (hMesh, qmname);             % load the source/detector specs
qvec = real(toastQvec (hMesh, 'Neumann', 'Gaussian', 2));
mvec = real(toastMvec (hMesh, 'Gaussian', 2));

% ======================================================================
% Read the data
data = toastReadRealVector(data_name);

% ======================================================================
% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;                 % absorption coeff.
mus = ones(n,1) * 2;                     % scattering coeff.
ref = ones(n,1) * refind;                % refractive index

% ======================================================================
% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis (hMesh, grid); % basis mapper handle
solmask = toastSolutionMask (hBasis);    % mask for removing exterior pixels

% ======================================================================
% Initial parameter estimates in solution basis
% Note that the solution parameters are given in terms of c*mua
% where c is the speed of light, and mua is the absorption coefficient.
bcmua = toastMapMeshToBasis (hBasis, mua * cm);
x = bcmua(solmask);  % initial solution vector
bmua = zeros(size(bcmua));
dmua = toastMapMeshToSol(hBasis,mua)-btgt(solmask); % initial image residual

% ======================================================================
% Initial data set f[x0]
proj = forwardModel(mua);

% ======================================================================
% initial data error
err0 = objective(proj);
err = err0;                                        % current error
errp = inf;                                        % previous error
erri(1) = err0;
errm(1) = norm(dmua);
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

step = 1.0;      % initial step length for line search
itr = 0;         % iteration counter

% ======================================================================
% Gauss-Newton loop
while (itr < itrmax) && (err > tolGN*err0) && (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    fprintf (1,'Calculating Jacobian\n');
    J = toastJacobianCW (hMesh, hBasis, qvec, mvec, mua, mus, ref, 'direct');
    
    % Normalisation of Hessian
    for i = 1:size(J,2)
        M(i) = 1 ./ sqrt(sum(J(:,i) .* J(:,i)));
    end
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
    [step, err] = toastLineSearch (x, dx, step, err, @objective2);
    
    % Add update to solution
    x = x + dx*step;
    
    % Map parameters back to mesh
    smua = x/cm;
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
    
    function of = objective2(x)
    smua = x/cm;                                % map from x (=c*mua) to mua
    if min(smua) < 0
        of = 1e8;                               % discourage negative solutions
    else
        mua = toastMapSolToMesh (hBasis, smua); % map to mesh basis
        of = objective(forwardModel(mua));      % calc objective function
    end
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
