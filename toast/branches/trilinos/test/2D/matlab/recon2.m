function recon2

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with nonlinear conjugate gradient solver')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshdir   = '../meshes/';
meshfile1 = [meshdir 'ellips_tri10.msh'];             % mesh for target data generation
meshfile2 = [meshdir 'circle25_32.msh'];              % mesh for reconstruction
qmfile    = [meshdir 'circle25_32x32.qm'];            % source-detector file
muafile   = [meshdir 'tgt_mua_ellips_tri10.nim'];     % nodal target absorption
musfile   = [meshdir 'tgt_mus_ellips_tri10.nim'];     % nodal target scattering

refind = 1.4;                           % refractive index
freq = 100;                             % modulation frequency [MHz]
noiselevel = 0.01;                      % Additive Gaussian data noise level
grd = [100 100];                        % solution basis: grid dimension

tau = 1e-4;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolCG = 1e-5;                           % Gauss-Newton convergence criterion
resetCG = 10;                           % PCG reset interval
itrmax = 100;                           % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian matrix

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
blen = prod(grd);
c0 = 0.3;
cm = c0/refind;

%% Generate target data

% Set up mesh geometry
hmesh_fwd = toastReadMesh(meshfile1);
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
hmesh = toastReadMesh (meshfile2);
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
x = [scmua; sckap];
logx = log(x);
%p = length(x);

% Create regularisation object
%hreg = toastRegul ('TK1', hraster, logx, tau);
hreg = toastRegul ('TV', hraster, logx, tau, 'Beta', beta);
%hreg = toastRegul ('TV', hraster, logx, tau, 'Beta', beta, 'KapRefImage', ...
%         [bmua_tgt;bkap_tgt], 'KapRefScale', 1, 'KapRefPMThreshold', 0.1);

figure(1);
subplot(2,2,1), imagesc(reshape(bmua_tgt,grd), mua_range), axis equal tight off; colorbar; title('\mu_a target');
subplot(2,2,2), imagesc(reshape(bmus_tgt,grd), mus_range), axis equal tight off; colorbar; title('\mu_s target');
subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal tight off; colorbar; title ('\mu_a recon');
subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal tight off; colorbar; title ('\mu_s recon');
drawnow

% initial data error (=2 due to data scaling)
err0 = privObjective (proj, data, sd, hreg, logx); %initial error
err = err0;                                        % current error
errp = inf;                                       % previous error
erri(1) = err0;

itr = 1;    % iteration counter
step = 1.0; % initial step length for line search
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

% NCG loop
while (itr <= itrmax) && (err > tolCG*err0) && (errp-err > tolCG)

    errp = err;
    
    % Gradient of cost function
    r = -toastGradient (hmesh, hraster, qvec, mvec, mua, mus, ref, freq, ...
                       data, sd, 'direct');
    r = r .* x; % parameter scaling
    
    rr = toastRegulGradient (hreg, logx);
    r = r - rr;
    
    % Plot the gradient
    figure(2);
    rgrd_mua = toastMapBasis (hraster, 'S->B', -rr(1:slen));
    subplot (2,2,1); imagesc (reshape (rgrd_mua, grd)); axis equal tight off; colorbar;title('Gradient of prior, \mu_a');
    rgrd_kap = toastMapBasis (hraster, 'S->B', -rr(slen+1:end));
    subplot (2,2,2); imagesc (reshape (rgrd_kap, grd)); axis equal tight off;colorbar;title('Gradient of prior, \kappa');
    grd_mua  = toastMapBasis (hraster, 'S->B', r(1:slen));
    subplot (2,2,3); imagesc (reshape (grd_mua, grd)); axis equal tight off;colorbar;title('Total gradient, \mu_a');
    grd_kap  = toastMapBasis (hraster, 'S->B', r(slen+1:end));
    subplot (2,2,4); imagesc (reshape (grd_kap, grd)); axis equal tight off;colorbar;title('Total gradient, \kappa');    
    drawnow
    
    if itr > 1
        delta_old = delta_new;
        delta_mid = r' * s;
    end
    
    % Apply PCG preconditioner
    s = r; % dummy for now
    
    if itr == 1
        d = s;
        delta_new = r' * d;
        delta0 = delta_new;
    else
        delta_new = r' * s;
        beta = (delta_new-delta_mid) / delta_old;
        if mod(itr, resetCG) == 0 || beta <= 0
            d = s;
        else
            d = s + d*beta;
        end
    end
    delta_d = d' * d;
    
    % Line search
    fprintf (1, 'Entering line search\n');
    [step, err] = toastLineSearch (logx, d, step, err, @objective);
    
    % Add update to solution
    logx = logx + d*step;
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
    subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal tight off; colorbar; title ('\mu_a recon');
    subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal tight off; colorbar; title ('\mu_s recon');
    drawnow
    
    proj = toastProject (hmesh, mua, mus, ref, freq, qvec, mvec);

    % update objective function
    err = privObjective (proj, data, sd, hreg, logx);
    itr = itr+1;
    erri(itr) = err;
    bmua_itr(itr,:) = bmua;
    bmus_itr(itr,:) = bmus;
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);
end

dlmwrite('recon.dat', erri);
disp('recon2: finished')

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
    [p, p_data, p_prior] = privObjective (proj, data, sd, hreg, x);
    fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end

end
