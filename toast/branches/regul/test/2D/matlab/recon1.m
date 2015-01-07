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
c0 = 0.3;
cm = c0/refind;

%% Generate target data

% Set up mesh geometry
mesh_fwd = toastMesh(meshfile_fwd);
mesh_fwd.ReadQM (qmfile);
dmask = mesh_fwd.DataLinkList;
qvec = mesh_fwd.Qvec (qtype, qprof, qwidth);
mvec = mesh_fwd.Mvec (mprof, mwidth, refind);
nlen = mesh_fwd.NodeCount;
nqm  = length(dmask);

% Target parameters
mua_tgt = toastNim(muafile).Values;
mus_tgt = toastNim(musfile).Values;
kap_tgt = 1./(3*(mua_tgt+mus_tgt));
ref_tgt = ones(nlen,1) * refind;

% Parameter plotting ranges
mua_range = [0.015 0.055];
mus_range = [1 4.5];

% Solve forward problem
phi = toastFields (mesh_fwd, 0, qvec, mua_tgt, mus_tgt, ref_tgt, freq, 'direct');
data = projection (phi, mvec, dmask);

% Add noise
data = data + data.*noiselevel.*randn(size(data));
m = length(data);

% Split in log amplitude and phase data
lnamp_tgt = data(1:nqm);
phase_tgt = data(nqm+1:end);

% Map target parameters to images for display
basis_fwd = toastBasis (mesh_fwd, grd);
bmua_tgt = basis_fwd.Map ('M->B', mua_tgt);
bmus_tgt = basis_fwd.Map ('M->B', mus_tgt);
bkap_tgt = basis_fwd.Map ('M->B', kap_tgt);

% Clear objects
clear basis_fwd;
clear mesh_fwd;

%% Solve inverse problem

% Set up mesh geometry
mesh = toastMesh (meshfile_inv);
mesh.ReadQM (qmfile);
qvec = mesh.Qvec (qtype, qprof, qwidth);
mvec = mesh.Mvec (mprof, mwidth, refind);
nlen = mesh.NodeCount;

% Initial parameter estimates
mua = ones(nlen,1) * 0.025;
mus = ones(nlen,1) * 2;
ref = ones(nlen,1) * refind;
kap = 1./(3*(mua+mus));

% Solution basis
basis = toastBasis (mesh, grd);

% Initial projections
phi = toastFields (mesh, 0, qvec, mua, mus, ref, freq);
proj = projection (phi, mvec, dmask);
lnamp = proj(1:nqm);
phase = proj(nqm+1:end);

% Data scaling
sd_lnmod = ones(size(lnamp)) * norm(lnamp_tgt-lnamp); 
sd_phase = ones(size(phase)) * norm(phase_tgt-phase);
sd = [sd_lnmod;sd_phase];

% Map parameter estimates to solution basis
bmua = basis.Map ('M->B', mua);
bmus = basis.Map ('M->B', mus);
bmua_itr(1,:) = bmua;
bmus_itr(1,:) = bmus;
bkap = basis.Map ('M->B', kap);
bcmua = bmua*cm;
bckap = bkap*cm;
scmua = basis.Map ('B->S', bcmua);
sckap = basis.Map ('B->S', bckap);
slen = length(scmua);

% Vector of unknowns
x = [scmua;sckap];
logx = log(x);
p = length(logx);

figure(1);
subplot(2,2,1), imagesc(reshape(bmua_tgt,grd), mua_range), axis equal tight off; colorbar; colormap('gray'); title('\mu_a target');
subplot(2,2,2), imagesc(reshape(bmus_tgt,grd), mus_range), axis equal tight off; colorbar; colormap('gray'); title('\mu_s target');
subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal tight off; colorbar; colormap('gray'); title ('\mu_a recon');
subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal tight off; colorbar; colormap('gray'); title ('\mu_s recon');
drawnow

step = 1.0; % initial step length for line search

% Create regularisation object
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
reg = toastRegul ('TV', basis, logx, tau, 'Beta', beta);

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
%reg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta, 'KapRefImage', ...
%                   btgt, 'KapRefScale', 1, 'KapRefPMThreshold', 0.1);

% initial data error (=2 due to data scaling)
err0 = objective (proj, data, sd, reg, logx);  %initial error
err = err0;                                    % current error
errp = inf;                                    % previous error
erri(1) = err0;
itr = 1; % iteration counter
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

% Gauss-Newton loop
while (itr <= itrmax) && (err > tolGN*err0) && (errp-err > tolGN)

    errp = err;
    
    % Construct the Jacobian
    fprintf (1,'Calculating Jacobian\n');
    J = toastJacobian (mesh, basis, qvec, mvec, mua, mus, ref, freq, 'direct');

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
    psiHdiag = reg.HDiag (logx);
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
    r = r - reg.Gradient (logx) .* M';
    
    if Himplicit == true
        % Update with implicit Krylov solver
        fprintf (1, 'Entering Krylov solver\n');
        dx = toastKrylov (x, J, r, M, 0, reg, tolKrylov);
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
    [step, err] = toastLineSearch (logx, dx, step0, err, @objective_ls);
    if errp-err <= tolGN
        dx = r; % try steepest descent
        [step, err] = toastLineSearch (logx, dx, step0, err, @objective_ls);
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
    mua  = basis.Map ('S->M', smua);
    mus  = basis.Map ('S->M', smus);
    bmua = basis.Map ('S->B', smua);
    bmus = basis.Map ('S->B', smus);

    figure(1);
    subplot(2,2,3), imagesc(reshape(bmua,grd), mua_range), axis equal tight off; colorbar; colormap(gray); title ('\mu_a recon');
    subplot(2,2,4), imagesc(reshape(bmus,grd), mus_range), axis equal tight off; colorbar; colormap(gray); title ('\mu_s recon');
    drawnow
    
    proj = projection(toastFields(mesh,0,qvec,mua,mus,ref,freq,'direct'),mvec,dmask);

    % update objective function
    err = objective (proj, data, sd, reg, logx);
    itr = itr+1;
    erri(itr) = err;
    bmua_itr(itr,:) = bmua;
    bmus_itr(itr,:) = bmus;
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);
end

dlmwrite('recon.dat', erri);
disp('recon1: finished')

    % =====================================================================
    % Run the forward solver for photon density field phi, given
    % parameters logx
    function phi = fwd(mesh,basis,logx,ref,freq,qvec)
        x_ = exp(logx);
        cm_ = 0.3./ref;
        scmua_ = x_(1:length(x_)/2);
        sckap_ = x_(length(x_)/2+1:end);
        cmua_ = basis.Map('S->M', scmua_);
        ckap_ = basis.Map('S->M', sckap_);
        mua_  = cmua_ ./ cm_;
        kap_  = ckap_ ./ cm_;
        mus_  = 1./(3*kap_) - mua_;
        phi = toastFields(mesh,0,qvec,mua_,mus_,ref,freq,'direct');
    end
    
    % =====================================================================
    % Projection of photon density fields to boundary data
    function prj = projection(phi,mvec,dmask)
        gamma_ = mvec.' * phi;
        gamma_ = reshape (gamma_, [], 1);
        gamma_ = gamma_(dmask);
        lgamma_ = log(gamma_);
        lnamp_ = real(lgamma_);
        phase_ = imag(lgamma_);
        prj = [lnamp_; phase_];
    end

    % =====================================================================
    % objective function
    
    function [p,p_data,p_prior] = objective(proj,data,sd,reg,logx)
        p_data = full(sum(((data-proj)./sd).^2));
        if nargin >= 5
            p_prior = full(reg.Value(logx));
        else
            p_prior = 0;
        end
        p = p_data + p_prior;
    end

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective_ls(logx)
        proj_ = projection(fwd(mesh,basis,logx,ref,freq,qvec),mvec,dmask);
        [p, p_data, p_prior] = objective (proj_, data, sd, reg, logx);
        fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end

end
