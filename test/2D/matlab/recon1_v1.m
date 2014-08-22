function recon1

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with Gauss-Newton solver')
disp('using edge prior information from correct target image.')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshdir = '../meshes/';
meshfile_fwd = [meshdir 'ellips_tri10.msh'];          % mesh file for data generation
meshfile_inv = [meshdir 'circle25_32.msh'];           % mesh file for reconstruction
qmfile    = [meshdir 'circle25_32x32.qm'];            % source-detector file
muafile   = [meshdir 'tgt_mua_ellips_tri10.nim'];     % nodal target absorption
musfile   = [meshdir 'tgt_mus_ellips_tri10.nim'];     % nodal target scattering

refind = 1.4;                           % refractive index
grd = [100 100];                        % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
tau = 1e-4;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolGN = 1e-6;                           % Gauss-Newton convergence criterion
tolKrylov = 1e-2;                       % Krylov convergence criterion
itrmax = 100;                           % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian matrix
noiselevel = 0;

qtype  = 'Neumann';                     % source type
qprof  = 'Gaussian';                    % source profile
qwidth = 2;                             % source width
mprof  = 'Gaussian';                    % detector profile
mwidth = 2;                             % detector width

% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();
c0 = 0.3;
cm = c0/refind;

%% Generate target data

% Set up mesh geometry
hmesh_fwd = toastReadMesh(meshfile_fwd);
toastReadQM (hmesh_fwd, qmfile);
dmask = toastDataLinkList (hmesh_fwd);
qvec = toastQvec (hmesh_fwd, qtype, qprof, qwidth);
mvec = toastMvec (hmesh_fwd, mprof, mwidth);
nlen = toastMeshNodeCount (hmesh_fwd);
nqm  = length(dmask);

% Target parameters
mua_tgt = toastReadNIM(muafile);
mus_tgt = toastReadNIM(musfile);
kap_tgt = 1./(3*(mua_tgt+mus_tgt));
ref_tgt = ones(nlen,1) * refind;

% Parameter plotting ranges
mua_range = [0.015 0.055];
mus_range = [1 4.5];

% Solve forward problem
phi = toastFields (hmesh_fwd, 0, qvec, mvec, mua_tgt, mus_tgt, ref_tgt, freq);
data = projection (phi, mvec, dmask);

% Add noise
data = data + data.*noiselevel.*randn(size(data));
m = length(data);

% Split in log amplitude and phase data
lnamp_tgt = data(1:nqm);
phase_tgt = data(nqm+1:end);

% Map target parameters to images for display
hraster_fwd = toastSetBasis (hmesh_fwd, grd);
bmua_tgt = toastMapBasis (hraster_fwd, 'M->B', mua_tgt);
bmus_tgt = toastMapBasis (hraster_fwd, 'M->B', mus_tgt);
bkap_tgt = toastMapBasis (hraster_fwd, 'M->B', kap_tgt);

% Clear objects
toastDeleteBasis (hraster_fwd);
toastDeleteMesh (hmesh_fwd);

%% Solve inverse problem

% Set up mesh geometry
hmesh = toastReadMesh (meshfile_inv);
toastReadQM (hmesh, qmfile);
qvec = toastQvec (hmesh, qtype, qprof, qwidth);
mvec = toastMvec (hmesh, mprof, mwidth);
nlen = toastMeshNodeCount (hmesh);

% Initial parameter estimates
mua = ones(nlen,1) * 0.025;
mus = ones(nlen,1) * 2;
ref = ones(nlen,1) * refind;
kap = 1./(3*(mua+mus));

% Solution basis
hraster = toastSetBasis (hmesh, grd);
solmask = toastSolutionMask (hraster);

% Initial projections
phi = toastFields (hmesh, 0, qvec, mvec, mua, mus, ref, freq);
proj = projection (phi, mvec, dmask);
lnamp = proj(1:nqm);
phase = proj(nqm+1:end);

% Data scaling
sd_lnmod = ones(size(lnamp)) * norm(lnamp_tgt-lnamp); 
sd_phase = ones(size(phase)) * norm(phase_tgt-phase);
sd = [sd_lnmod;sd_phase];

% Map parameter estimates to solution basis
bmua = toastMapBasis (hraster, 'M->B', mua);
bmus = toastMapBasis (hraster, 'M->B', mus);
bmua_itr(1,:) = bmua;
bmus_itr(1,:) = bmus;
bkap = toastMapBasis (hraster, 'M->B', kap);
bcmua = bmua*cm;
bckap = bkap*cm;
scmua = toastMapBasis (hraster, 'B->S', bcmua);
sckap = toastMapBasis (hraster, 'B->S', bckap);
slen = length(scmua);

% Vector of unknowns
x = [scmua;sckap];
logx = log(x);
p = length(logx);

figure(1);
subplot(2,2,1), imagesc(reshape(bmua_tgt,grd), mua_range), axis equal tight off; colorbar; title('\mu_a target');
subplot(2,2,2), imagesc(reshape(bmus_tgt,grd), mus_range), axis equal tight off; colorbar; title('\mu_s target');
subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal tight off; colorbar; title ('\mu_a recon');
subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal tight off; colorbar; title ('\mu_s recon');
drawnow

step = 1.0; % initial step length for line search

% Create regularisation object
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
hReg = toastRegul ('TV', hraster, logx, tau, 'Beta', beta);

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
err0 = privObjective (proj, data, sd, hReg, logx);  %initial error
err = err0;                                         % current error
errp = inf;                                         % previous error
erri(1) = err0;
itr = 1; % iteration counter
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

% Gauss-Newton loop
while (itr <= itrmax) && (err > tolGN*err0) && (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    fprintf (1,'Calculating Jacobian\n');
    J = toastJacobian (hmesh, hraster, qvec, mvec, mua, mus, ref, freq, 'direct');

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
    for i = 1:p
        M(i) = sum(J(:,i) .* J(:,i));
        M(i) = M(i) + psiHdiag(i);
        M(i) = 1 ./ sqrt(M(i));
    end
    for i = 1:p
        J(:,i) = J(:,i) * M(i);
    end
    %J = J * spdiags (M',0,p,p);
    
    % Gradient of cost function
    r = J' * (2*(data-proj)./sd);
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
    step0 = step;
    [step, err] = toastLineSearch (logx, dx, step0, err, @objective);
    if errp-err <= tolGN
        dx = r; % try steepest descent
        [step, err] = toastLineSearch (logx, dx, step0, err, @objective);
    end

    % Add update to solution
    logx = logx + dx*step;
    x = exp(logx);
    
    % Map parameters back to mesh
    scmua = x(1:slen);
    sckap = x(slen+1:end);
    smua = scmua/cm;
    skap = sckap/cm;
    smus = 1./(3*skap) - smua;
    mua  = toastMapBasis (hraster, 'S->M', smua);
    mus  = toastMapBasis (hraster, 'S->M', smus);
    bmua = toastMapBasis (hraster, 'S->B', smua);
    bmus = toastMapBasis (hraster, 'S->B', smus);

    figure(1);
    subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal, axis tight, colorbar('horiz'); colormap(gray);
    subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal, axis tight, colorbar('horiz'); colormap(gray);
    drawnow
    
    proj = toastProject (hmesh, mua, mus, ref, freq, qvec, mvec);

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
    % Projection of photon density fields to boundary data
    function prj_ = projection(phi_,mvec_,dmask_)
    
    gamma_ = mvec_.' * phi_;
    gamma_ = reshape (gamma_, [], 1);
    gamma_ = gamma_(dmask_);
    lgamma_ = log(gamma_);
    lnamp_ = real(lgamma_);
    phase_ = imag(lgamma_);
    prj_ = [lnamp_; phase_];
    end

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (hmesh, hraster, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = privObjective (proj, data, sd, hReg, x);
    fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end

end
