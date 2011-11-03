function toastReconMultispectralCW(prm)
%toastReconMultispectralCW - Multispectral reconstruction from CW data.
%
% Synopsis: toastReconMultispectralCW(prm)
%    prm: reconstruction parameter structure
%
% Runs a toast reconstruction with the parameters defined in prm.
% prm contains information about measurements, meshes and grids,
% tolerance limits for the forward and inverse solvers, regularisation
% parameters, etc.
%
% This version uses CW intensity data at several wavelengths to
% reconstruct directly for chromophore concentrations, given extinction
% coefficients for all chromophores at all wavelengths.
%
% See also: toastRecon, toastReadParam
  
toastCatchErrors();
disp('---------------------------------------')
disp('Starting multispectral reconstruction from CW data')
disp('---------------------------------------')

% ----------------------------------------------------------------------
% Defining file paths and parameters
global ITR_TERM
ITR_TERM = false;
prm = checkprm(prm);

itrmax = 100;                           % Gauss-Newton max. iterations

nch = prm.nch;                          % number of chromophores
nlambda = prm.nlambda;                  % number of wavelengths
extinct = prm.extinct;                  % extinction coefficients

% ----------------------------------------------------------------------
% Read a TOAST mesh definition from file.
if isfield(prm,'basis') && isfield(prm.basis,'hMesh')
    hMesh = prm.basis.hMesh;
else
    hMesh = toastReadMesh (prm.fwdsolver.meshfile);
    toastReadQM (hMesh, prm.meas.qmfile);
end
n = toastMeshNodeCount (hMesh);

% ----------------------------------------------------------------------
% Generate source and measurement vectors
qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);
mvec = toastMvec (hMesh, 'Gaussian', 2);
dmask = toastDataLinkList (hMesh);
nqm = length(dmask);

% ----------------------------------------------------------------------
% Set up the mapper between FEM and solution bases
if isfield(prm,'basis') && isfield(prm.basis,'hBasis')
    hBasis = prm.basis.hBasis;
else
    hBasis = toastSetBasis (hMesh, prm.basis.bdim);
end
if isfield(prm,'smask')
    solmask = prm.smask;
else
    solmask = toastSolutionMask(hBasis);
end
slen = length(solmask);

% ----------------------------------------------------------------------
% Fixed background parameters
mus = resetprm (prm.initprm.mus, hMesh);
ref = resetprm (prm.initprm.ref, hMesh);
bmus = toastMapMeshToBasis (hBasis, mus);
bref = toastMapMeshToBasis (hBasis, ref);

% ----------------------------------------------------------------------
% Measurement data
data = prm.data;
m = length(data);

% ----------------------------------------------------------------------
% Set up homogeneous initial parameter estimates - generalise!
for i=1:prm.nch
    if isfield(prm,'chinit')
        C(i,:) = ones(n,1) * prm.chinit(i);
    else
        C(i,:) = ones(n,1) * 0.5;
    end
end

% ----------------------------------------------------------------------
% initial parameter estimates in solution basis
x = [];
for i=1:prm.nch
    bC(i,:) = toastMapMeshToBasis(hBasis,C(i,:));
    sC(i,:) = bC(i,solmask);
    x = [x; sC(i,:)'];
end
logx = log(max(x,1e-20));
p = length(x);

% ----------------------------------------------------------------------
% Initial data set f[x0]
proj = ProjectAll(C);

% ----------------------------------------------------------------------
% data scaling
for i=1:nlambda
    pr = proj((i-1)*nqm+1:i*nqm);
    dt = data((i-1)*nqm+1:i*nqm);
    nm = norm(pr-dt);
    sd((i-1)*nqm+1:i*nqm,1) = ones(nqm,1)*nm;
end

% ----------------------------------------------------------------------
% initial data error
err0 = objective (proj);           %initial error
err = err0;                        % current error
errp = inf;                        % previous error
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

hReg = 0; % for now
step = prm.solver.step0; % initial step length for line search
step = 1e2;

itr = 0; % iteration counter

res.of(itr+1) = err0;
res.bC = bC;
feval(prm.callback.iter, prm.callback.context, res);

% ----------------------------------------------------------------------
% Gauss-Newton loop
while (itr < itrmax) && ...
      (err > prm.solver.tol*err0) && ...
      (errp-err > prm.solver.tol) && ...
      (ITR_TERM == false)

    errp = err;
    
    % Build the Jacobian
    clear J;
    for i = 1:nlambda
        mua = GetMua(extinct(:,i),C);
        Ji = toastJacobianCW (hMesh, hBasis, qvec, mvec, mua, mus, ref, prm.linsolver.method, prm.linsolver.tol);
        for j = 1:prm.nch
            J((i-1)*nqm+1:i*nqm,(j-1)*slen+1:j*slen) = Ji * extinct(j,i);
        end
        clear Ji;
    end
    
    % data normalisation
    for i = 1:m, J(i,:) = J(i,:) / sd(i); end
    
    % parameter normalisation
    for i = 1:p, J(:,i) = J(:,i) * x(i); end;

    % Normalisation of Hessian
    for i = 1:size(J,2)
        M(i) = 1 ./ sqrt(sum(J(:,i) .* J(:,i)));
    end
    for i = 1:size(M,2)
        J(:,i) = J(:,i) * M(1,i);
    end
    
    % Gradient of cost function
    dy = ((data-proj)./sd);
    for i=1:size(J,2) % do product by hand to save space
        r(i,1) = J(:,i)' * dy;
    end
    %r = J' * 2*((data-proj)./sd);

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
    for i=1:prm.nch
        sC(i,:) = x((i-1)*slen+1:i*slen);
        bC(i,solmask) = sC(i,:);
        C(i,:) = toastMapBasisToMesh (hBasis, bC(i,:));
    end

    % update objective function
    itr = itr+1;
    proj = ProjectAll(C);
    err = objective(proj);
    
    res.bC = bC;
    res.of(itr+1) = err;
    feval(prm.callback.iter, prm.callback.context, res);
    
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
        proj_ = ProjectSingle (extinct(:,i_), C);
        p_((i_-1)*nqm+1:i_*nqm,1) = proj_(dmask);
    end
    end

    
    % ======================================================================
    % Calculate mua given extinction coefficients and chromophore
    % concentrations
    
    function mua_ = GetMua (extinct, C)
    mua_ = zeros(size(C,2),1);
    for i_ = 1:prm.nch
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
    
    function of_ = objective2(logx)
        x = exp(logx);
        for i_ = 1:prm.nch
            C_(i_,:) = toastMapSolToMesh(hBasis,x((i_-1)*slen+1:i_*slen));
        end
        of_ = objective(ProjectAll(C_));        % calc objective function
    end
    

    % =====================================================================
    % Krylov solver subroutine
    
    function dx_ = krylov(r)
    switch prm.solver.krylov.method
        case 'bicgstab'
            dx_ = bicgstab(@jtjx, r, prm.solver.krylov.tol, prm.solver.krylov.maxit);
        otherwise
            dx_ = gmres (@jtjx, r, 30, prm.solver.krylov.tol, prm.solver.krylov.maxit);
    end
    end


    % =====================================================================
    % Callback function for matrix-vector product (called by krylov)
    
    function b = jtjx(x)
    b = J' * (J*x);
    end

end


% ======================================================================
% check validity of reconstruction parameters

function prm = checkprm(prm)
prm.solver.method = 'LM'; % only supported method for now
if isfield(prm.solver,'tol') == false
    prm.solver.tol = 1e-10;
end
if isfield(prm.linsolver,'tol') == false
    prm.linsolver.tol = 1e-10;
end
if isfield(prm.solver,'krylov') == false || isfield(prm.solver.krylov,'method') == false
    prm.solver.krylov.method = 'gmres';
end
if isfield(prm.solver.krylov,'tol') == false
    prm.solver.krylov.tol = 1e-2;
end
if isfield(prm.solver.krylov,'maxit') == false
    prm.solver.krylov.maxit = 100;
end
end


% =============================================================
% Set up initial parameter distributions

function prm=resetprm(cfg,hmesh)
n = toastMeshNodeCount (hmesh);
switch upper(cfg.reset)
    case 'HOMOG'
        prm = ones(n,1) * cfg.val;
    case 'NIM'
        prm = toastReadNIM(cfg.nim);
        if length(prm) ~= n
            disp('Warning: incompatible size of NIM file')
        end
    case 'NIM_LAST'
        prm = toastReadNIM(cfg.nim,0);
        if length(prm) ~= n
            disp('Warning: incompatible size of NIM file')
        end
    otherwise
        disp('Warning: Unsupported reset method')
end
end
