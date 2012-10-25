function recon1

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with conjugate gradient solver')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
verbosity = 1;
meshdir = '../../test/2D/meshes/';       % example mesh repository
qmname  = [meshdir 'circle25_32x32.qm']; % source-detector definitions
fwdmesh = [meshdir 'ellips_tri10.msh'];  % mesh for data generation
invmesh = [meshdir 'circle25_32.msh'];   % mesh for inverse solver forward model

refind = 1.4;                           % refractive index

grd = [100 100];                        % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
noiselevel = 0.01;                      % additive data noise level
tau = 1e-3;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolCG = 1e-6;                           % Gauss-Newton convergence criterion
resetCG = 10;                           % PCG reset interval
itrmax = 100;                           % Gauss-Newton max. iterations
cmap = 'gray';

% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();                     % redirect toast library errors
toastSetVerbosity(verbosity);           % output verbosity level
c0 = 0.3;                               % lightspeed in vacuum [mm/ps]
cm = c0/refind;                         % lightspeed in medium

%% Generate target data

hmesh = toastReadMesh (fwdmesh);        % load FEM mesh from file
toastReadQM (hmesh, qmname);            % add source-detector descriptions
n = toastMeshNodeCount (hmesh);         % number of nodes
dmask = toastDataLinkList (hmesh);      % source-detector connectivity

mua = toastReadNIM([meshdir 'tgt_mua_ellips_tri10.nim']); % target absorption
mus = toastReadNIM([meshdir 'tgt_mus_ellips_tri10.nim']); % target scattering
ref = ones(n,1)*refind;   % target refractive index (homogeneous)

qvec = toastQvec (hmesh, 'Neumann', 'Gaussian', 2); % source specification
mvec = toastMvec (hmesh, 'Gaussian', 2);            % detector specification

smat = toastSysmat (hmesh, mua, mus, ref, freq);    % compute FEM system matrix
phi = full (smat\qvec);                             % solve linear FEM problem for photon density

lgamma = reshape (log(mvec.' * phi), [], 1);        % map to measurements
lgamma = lgamma(dmask);                             % remove unused source-detector combinations
mdata = real(lgamma);                               % log amplitude data
pdata = imag(lgamma);                               % phase data

% add some noise
mdata = mdata + mdata.*noiselevel.*randn(size(mdata));
pdata = pdata + pdata.*noiselevel.*randn(size(pdata));
data = [mdata;pdata];                               % linear data vector

% display the target parameter distributions for comparison
muarng = [min(mua)*0.9, max(mua)*1.1];
musrng = [min(mus)*0.9, max(mus)*1.1];
hbasis = toastSetBasis(hmesh,grd);
muatgt_img = reshape (toastMapBasis (hbasis, 'M->B', mua), grd);
mustgt_img = reshape (toastMapBasis (hbasis, 'M->B', mus), grd);
mua_img = [muatgt_img, zeros(size(muatgt_img))];
mus_img = [mustgt_img, zeros(size(mustgt_img))];

hfig = figure;
set(hfig,'Position',[1 1, 800 400]);
subplot(2,2,1);
imagesc (mua_img, muarng);
colormap(cmap); colorbar; axis equal tight off
title ('\mu_a tgt, recon');
set(gca,'Position',[0.01 0.52 0.4 0.4]);
subplot(2,2,3);
imagesc (mus_img, musrng);
colormap(cmap); colorbar; axis equal tight off
title ('\mu_s tgt, recon');
set(gca,'Position',[0.01 0.05 0.4 0.4]);
drawnow

toastDeleteBasis(hbasis);                           % clean up basis mapper
toastDeleteMesh(hmesh);                             % clean up mesh


%% Inverse solver

% Read a TOAST mesh definition from file.
hmesh = toastReadMesh (invmesh);                    % read inverse solver mesh
toastReadQM (hmesh, qmname);                        % add source/detector descriptions
n = toastMeshNodeCount (hmesh);                     % number of nodes
dmask = toastDataLinkList (hmesh);                  % source-detector connectivity

% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;                            % initial mua estimate
mus = ones(n,1) * 2;                                % initial mus estimate
ref = ones(n,1) * refind;                           % refractive index estimate
kap = 1./(3*(mua+mus));                             % diffusion coefficient

% Set up the mapper between FEM and solution bases
hbasis = toastSetBasis ('LINEAR', hmesh, grd);      % maps between mesh and reconstruction basis
solmask = toastSolutionMask (hbasis);               % mask unused voxels

% Generate source vectors
qvec = toastQvec (hmesh, 'Neumann', 'Gaussian', 2); % nodal source vectors

% Generate measurement vectors
mvec = toastMvec (hmesh, 'Gaussian', 2);            % nodal measurement vectors

% Initial data set f[x0]
smat = toastSysmat (hmesh, mua, mus, ref, freq);    % FEM system matrix
lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);% solve for photon density and map to boundary measurements
lgamma = lgamma(dmask);                             % remove unused source-detector combinations
mproj = real(lgamma);                               % log amplitude data
pproj = imag(lgamma);                               % phase data
proj = [mproj;pproj];                               % linear measurement vector

% data scaling
msd = ones(size(lgamma)) * norm(mdata-mproj);       % scale log amp data with data difference
psd = ones(size(lgamma)) * norm(pdata-pproj);       % scale phase data with data difference
sd = [msd;psd];                                     % linear scaling vector

% map initial parameter estimates to solution basis
bmua = toastMapMeshToBasis (hbasis, mua);           % mua mapped to full grid
bmus = toastMapMeshToBasis (hbasis, mus);           % mus mapped to full grid
bkap = toastMapMeshToBasis (hbasis, kap);           % kap mapped to full grid
bcmua = bmua*cm;                                    % scale parameters with speed of light
bckap = bkap*cm;                                    % scale parameters with speed of light
scmua = bcmua(solmask);                             % map to solution basis
sckap = bckap(solmask);                             % map to solution basis

% solution vector
x = [scmua;sckap];                                  % linea solution vector
logx = log(x);                                      % transform to log

% Initialise regularisation
hreg = toastRegul ('TV', hbasis, logx, tau, 'Beta', beta);

% initial data error (=2 due to data scaling)
err0 = privObjective (proj, data, sd, hreg, logx);  %initial error
err = err0;                                         % current error
errp = inf;                                         % previous error
erri(1) = err0;                                     % keep history
itr = 1;                                            % iteration counter
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);
step = 1.0;                                         % initial step length for line search

% Nonlinear conjugate gradient loop
while (itr <= itrmax) && (err > tolCG*err0) && (errp-err > tolCG)

    errp = err;
    
    % Gradient of cost function
    r = -toastGradient (hmesh, hbasis, qvec, mvec, mua, mus, ref, freq, ...
                       data, sd, 'direct');
    r = r .* x; % parameter scaling
    
    rr = toastRegulGradient (hreg, logx);
    r = r - rr;
    
    if itr > 1
        delta_old = delta_new;
        delta_mid = r' * s;
    end
    
    % Apply PCG preconditioner
    s = r; % dummy for now
    
    if itr == 1
        d = s;
        delta_new = r' * d;
    else
        delta_new = r' * s;
        beta = (delta_new - delta_mid) / delta_old;
        if mod (itr, resetCG) == 0 || beta <= 0
            d = s;  % reset CG
        else
            d = s + d*beta;
        end
    end
    
    % Line search
    fprintf (1, 'Entering line search\n');
    step = toastLineSearch (logx, d, step, err, @objective);
    
    % Add update to solution
    logx = logx + d*step;
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

    %figure(1);
    subplot(2,2,1);
    muarec_img = reshape(bmua,grd);
    mua_img(:,size(muarec_img,2)+1:end) = muarec_img;
    imagesc(mua_img, muarng);
    colormap(cmap); colorbar; axis equal tight off
    title ('\mu_a tgt, recon');
    set(gca,'Position',[0.01 0.52 0.4 0.4]);
    subplot(2,2,3);
    musrec_img = reshape(bmus,grd);
    mus_img(:,size(musrec_img,2)+1:end) = musrec_img;
    imagesc(mus_img, musrng); 
    colormap(cmap); colorbar; axis equal tight off
    title ('\mu_s tgt, recon');
    set(gca,'Position',[0.01 0.05 0.4 0.4]);
    
    proj = toastProject (hmesh, mua, mus, ref, freq, qvec, mvec);

    % update objective function
    err = privObjective (proj, data, sd, hreg, logx);
    fprintf (1, '**** GN ITERATION %d, ERROR %f\n\n', itr, err);

    itr = itr+1;
    erri(itr) = err;
    
    % show objective function
    subplot(1,2,2);
    semilogy(erri);
    axis([1 itr 1e-2 2])
    xlabel('iteration');
    ylabel('objective function');
    drawnow
end

disp('recon1: finished')

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (hmesh, hbasis, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = privObjective (proj, data, sd, hreg, x);
    if verbosity > 0
        fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end
    end

end
