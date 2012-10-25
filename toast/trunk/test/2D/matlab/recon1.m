function recon1

% Reconstruction of optical parameters in 2D circle

% ======================================================================
% User-defined parameters
% ======================================================================
verbosity = 0;                           % silent mode
meshdir = '../meshes/';                  % example mesh repository
qmname  = [meshdir 'circle25_32x32.qm']; % source-detector definitions
fwdmesh = [meshdir 'ellips_tri10.msh'];  % mesh for data generation
invmesh = [meshdir 'circle25_32.msh'];   % mesh for inverse solver forward model

refind = 1.4;                           % refractive index

grd = [100 100];                        % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
noiselevel = 0.01;                      % additive data noise level
tau = 1e-3;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolGN = 1e-6;                           % Gauss-Newton convergence criterion
tolKrylov = 1e-2;                       % Krylov convergence criterion
itrmax = 5;                             % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian matrix

% ======================================================================
% End user-defined parameters
% ======================================================================

% Make sure the toast paths are available
cwd = pwd;
cd (getenv('TOASTDIR'));
mtoast_install(true);
cd (cwd);

% Initialisations
toastCatchErrors();
toastSetVerbosity(verbosity);
c0 = 0.3;                               % lightspeed in vacuum [mm/ps]
cm = c0/refind;                         % lightspeed in medium

%% Generate target data

hmesh = toastReadMesh (fwdmesh);
toastReadQM (hmesh, qmname);
n = toastMeshNodeCount (hmesh);
dmask = toastDataLinkList (hmesh);
nqm = length(dmask);

mua = toastReadNIM([meshdir 'tgt_mua_ellips_tri10.nim']);
mus = toastReadNIM([meshdir 'tgt_mus_ellips_tri10.nim']);
ref = ones(n,1)*refind;

qvec = toastQvec (hmesh, 'Neumann', 'Gaussian', 2);
mvec = toastMvec (hmesh, 'Gaussian', 2);

smat = toastSysmat (hmesh, mua, mus, ref, freq);
phi = full (smat\qvec);

lgamma = reshape (log(mvec.' * phi), [], 1);
mdata = real(lgamma);  % log amplitude data
pdata = imag(lgamma);  % phase data
data = [mdata;pdata];  % linear data vector
m = length(data);

toastDeleteMesh(hmesh);


%% Inverse solver

% Read a TOAST mesh definition from file.
hmesh = toastReadMesh (invmesh);
toastReadQM (hmesh, qmname);
n = toastMeshNodeCount (hmesh);
dmask = toastDataLinkList (hmesh);

% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;
mus = ones(n,1) * 2;
ref = ones(n,1) * refind;
kap = 1./(3*(mua+mus));

% Set up the mapper between FEM and solution bases
hbasis = toastSetBasis ('LINEAR', hmesh, grd);
blen = prod(grd);
solmask = toastSolutionMask (hbasis);
solmask2 = [solmask solmask+blen];

% Generate source vectors
qvec = toastQvec (hmesh, 'Neumann', 'Gaussian', 2);

% Generate measurement vectors
mvec = toastMvec (hmesh, 'Gaussian', 2);

% Initial data set f[x0]
smat = toastSysmat (hmesh, mua, mus, ref, freq);
lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);
lgamma = lgamma(dmask);
mproj = real(lgamma);
pproj = imag(lgamma);
proj = [mproj;pproj];

% data scaling
msd = ones(size(lgamma)) * norm(mdata-mproj); 
psd = ones(size(lgamma)) * norm(pdata-pproj);
sd = [msd;psd];

% map initial parameter estimates to solution basis
bmua = toastMapMeshToBasis (hbasis, mua);
bmus = toastMapMeshToBasis (hbasis, mus);
bmua_itr(1,:) = bmua;
bmus_itr(1,:) = bmus;
bkap = toastMapMeshToBasis (hbasis, kap);
bcmua = bmua*cm;
bckap = bkap*cm;
scmua = bcmua(solmask);
sckap = bckap(solmask);

% solution vector
x = [scmua;sckap];
logx = log(x);   % transform to log
p = length(x);

% Initialise regularisation
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
hreg = toastRegul ('TV', hbasis, logx, tau, 'Beta', beta);

%tgtmua=toastReadNIM('../meshes/tgt_mua_ellips_tri10.nim');
%tgtmus=toastReadNIM('../meshes/tgt_mus_ellips_tri10.nim');
%tgtkap=1./(3*(tgtmua+tgtmus));
%hMesh2 = toastReadMesh('../meshes/ellips_tri10.msh');
%hBasis2 = toastSetBasis ('LINEAR',hMesh2,grd);
%btgtmua=toastMapMeshToGrid (hBasis2, tgtmua);
%btgtkap=toastMapMeshToGrid (hBasis2, tgtkap);
%toastDeleteBasis(hBasis2);
%toastDeleteMesh(hMesh2);
%btgt = [btgtmua; btgtkap];
%hReg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta, 'KapRefImage', ...
%                   btgt, 'KapRefScale', 1, 'KapRefPMThreshold', 0.1);

% initial data error (=2 due to data scaling)
err0 = privObjective (proj, data, sd, hreg, logx);  %initial error
err = err0;                                         % current error
errp = inf;                                         % previous error
erri(1) = err0;
itr = 1; % iteration counter
step = 1.0; % initial step length for line search

% Gauss-Newton loop
while (itr <= itrmax) & (err > tolGN*err0) & (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    J = toastJacobian (hmesh, hbasis, qvec, mvec, mua, mus, ref, freq, 'bicgstab', 1e-16);

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
    psiHdiag = toastRegulHDiag (hreg, logx);
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
    r = r - toastRegulGradient (hreg, logx) .* M';
    
    if Himplicit == true
        % Update with implicit Krylov solver
        dx = toastKrylov (x, J, r, M, 0, hreg, tolKrylov);
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
    step0 = step;
    [step, err] = toastLineSearch (logx, dx, step0, err, @objective, 'verbose', verbosity>0);
    if errp-err <= tolGN
        dx = r; % try steepest descent
        [step, err] = toastLineSearch (logx, dx, step0, err, @objective, 'verbose', verbosity>0);
    end

    % Add update to solution
    logx = logx + dx*step;
    x = exp(logx);
    
    % Map parameters back to mesh
    scmua = x(1:size(x)/2);
    sckap = x(size(x)/2+1:size(x));
    smua = scmua/cm;
    skap = sckap/cm;
    smus = 1./(3*skap) - smua;
    mua = toastMapSolToMesh (hbasis, smua);
    mus = toastMapSolToMesh (hbasis, smus);
    bmua(solmask) = smua;
    bmus(solmask) = smus;

    proj = toastProject (hmesh, mua, mus, ref, freq, qvec, mvec);

    % update objective function
    err = privObjective (proj, data, sd, hreg, logx);

    itr = itr+1;
    erri(itr) = err;
    bmua_itr(itr,:) = bmua;
    bmus_itr(itr,:) = bmus;
end

% Uncomment this to create a new reference file with the current results
%dlmwrite ('recon1.dat', erri,'precision','%0.12e');

% Compare the results with stored values
erri_ref = dlmread('recon1.dat');
maxerr = max(norm(erri_ref-erri)./norm(erri_ref));
if maxerr > 1e-10
    exit(1); % return with error flag
else
    exit(0);
end

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (hmesh, hbasis, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = privObjective (proj, data, sd, hreg, x);
    end

end
