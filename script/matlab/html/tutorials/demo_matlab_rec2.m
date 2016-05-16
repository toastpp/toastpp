% Toast-Matlab web example 5:
% A simple reconstruction of absorption and scattering distributions in
% a 2-D problem from frequency domain boundary data
% (c) Martin Schweiger and Simon Arridge
% www.toastplusplus.org

function demo_matlab_rec2

freq = 100;   % modulation frequency [MHz]
tau = 1e-6;   % regularisation parameter
beta = 0.01;  % TV threshold parameter
tolCG = 1e-10; % nonlinear solver convergence criterion
resetCG = 30; % reset interval
itrmax = 200; % nonlinear solver max iteration count

%% Generate the forward data

% Create the mesh
rad = 25;
nsect = 6;
nring = 64;
nbnd = 4;
[vtx,idx,eltp] = mkcircle (rad, nsect, nring, nbnd);
mesh = toastMesh(vtx,idx,eltp);
mesh.Reorder(mesh.Optimise('mmd'));

% Load bitmaps for parameter distributions
global tgt
load circle_targets.mat
bmua = tgt.mua;
bmus = tgt.mus;
grd = tgt.grd;

% Scale to desired parameter ranges
muarng = [0.005,0.04];
musrng = [0.5,4];
figure;
subplot(2,3,1);
imagesc(rot90(bmua),muarng);
axis equal tight off
title('target \mu_a');
subplot(2,3,2);
imagesc(rot90(bmus),musrng);
axis equal tight off
title('target \mu_s');

% Map to mesh basis
basis = toastBasis (mesh,grd);
mua = basis.Map ('B->M',bmua);
mus = basis.Map ('B->M',bmus);

nnd = mesh.NodeCount;
ref_bkg = 1.4;
cm = 0.3/ref_bkg;  % speed of light [mm/ps]
ref = ones(nnd,1) * ref_bkg;

% Create the source and detector positions
nq = 32;
for i=1:nq
    phi_q = 2*pi*(i-1)/nq;
    Q(i,:) = rad * [cos(phi_q) sin(phi_q)];
    phi_m = 2*pi*(i-0.5)/nq;
    M(i,:) = rad * [cos(phi_m) sin(phi_m)];
end
mesh.SetQM(Q,M);

% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2, ref);

% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,freq);
Phi = K\qvec;
Y = mvec.' * Phi;
logY = reshape(log(Y),[],1);
mdata = real(logY);  % log amplitude data
pdata = imag(logY);  % phase data

% add some noise
noiselevel = 0.0;
mdata = mdata + mdata.*noiselevel.*randn(size(mdata));
pdata = pdata + pdata.*noiselevel.*randn(size(pdata));
data = [mdata;pdata];                         % linear data vector

%% Inverse solver

% Create the (coarser) mesh for the inverse solver
% Create the mesh
rad = 25;
nsect = 6;
nring = 64;
nbnd = 2;
[vtx,idx,eltp] = mkcircle (rad, nsect, nring, nbnd);
mesh = toastMesh(vtx,idx,eltp);
mesh.SetQM(Q,M);
n = mesh.NodeCount;

% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.01;                             % initial mua estimate
mus = ones(n,1) * 1;                                % initial mus estimate
ref = ones(n,1) * ref_bkg;                          % refractive index estimate
kap = 1./(3*(mua+mus));                             % diffusion coefficient

% Set up the mapper between FEM and solution bases
grd = [128 128];
basis = toastBasis (mesh, grd, 'LINEAR');           % maps between mesh and reconstruction basis

% Generate source vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);        % nodal source vectors

% Generate measurement vectors
mvec = mesh.Mvec ('Gaussian', 2, ref);              % nodal measurement vectors

% Initial data set f[x0]
K = dotSysmat (mesh, mua, mus, ref, freq);          % FEM system matrix
lgamma = reshape (log(mvec.' * (K\qvec)), [], 1);   % solve for photon density and map to boundary measurements
mproj = real(lgamma);                               % log amplitude data
pproj = imag(lgamma);                               % phase data
proj = [mproj;pproj];                               % linear measurement vector

% data scaling
msd = ones(size(lgamma)) * norm(mdata-mproj);       % scale log amp data with data difference
psd = ones(size(lgamma)) * norm(pdata-pproj);       % scale phase data with data difference
sd = [msd;psd];                                     % linear scaling vector

% map initial parameter estimates to solution basis
bmua = basis.Map ('M->B', mua);                     % mua mapped to full grid
bmus = basis.Map ('M->B', mus);                     % mus mapped to full grid
bkap = basis.Map ('M->B', kap);                     % kap mapped to full grid
bcmua = bmua*cm;                                    % scale parameters with speed of light
bckap = bkap*cm;                                    % scale parameters with speed of light
scmua = basis.Map ('B->S', bcmua);                  % map to solution basis
sckap = basis.Map ('B->S', bckap);                  % map to solution basis

% solution vector
x = [scmua;sckap];                                  % linea solution vector
logx = log(x);                                      % transform to log

% Initialise regularisation
reg = toastRegul ('TV', basis, logx, tau, 'Beta', beta);

% initial data error (=2 due to data scaling)
err0 = toastObjective (proj, data, sd, reg, logx);  %initial error
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
    r = -toastGradient (mesh, basis, qvec, mvec, mua, mus, ref, freq, ...
                       data, sd, 'Method', 'direct');
    r = r .* x;                   % parameter scaling
    r = r - reg.Gradient (logx);  % regularisation contribution
    
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
    fprintf (1, 'Line search:\n');
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
    mua = basis.Map ('S->M', smua);
    mus = basis.Map ('S->M', smus);
    bmua = basis.Map ('S->B', smua);
    bmus = basis.Map ('S->B', smus);

    %figure(1);
    subplot(2,3,4);
    imagesc(rot90(reshape(bmua,grd)),muarng); axis equal tight off
    title ('recon \mu_a');
    subplot(2,3,5);
    imagesc(rot90(reshape(bmus,grd)),musrng); axis equal tight off
    title ('recon \mu_s');
    
    proj = toastProject (mesh, mua, mus, ref, freq, qvec, mvec);

    % update objective function
    err = toastObjective (proj, data, sd, reg, logx);
    fprintf ('GN iteration %d\n', itr);
    fprintf ('--> Objective: %f\n', err);

    itr = itr+1;
    erri(itr) = err;
    
    % show objective function
    subplot(1,3,3);
    semilogy(erri);
    axis tight
    xlabel('iteration');
    ylabel('objective function');
    drawnow
end

disp('recon1: finished')

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    [mua,mus] = dotXToMuaMus (basis, exp(x), ref);
    proj = toastProject (mesh, mua, mus, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = toastObjective (proj, data, sd, reg, x);
    fprintf (1, '--> LH: %f, PR: %f\n', p_data, p_prior);
    end

end
