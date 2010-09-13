function recon2

disp('MATLAB-TOAST sample script:')
disp('2D image reconstruction with nonlinear conjugate gradient solver')
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = '../meshes/circle25_32.msh';              % mesh file
qmname    = '../meshes/circle25_32x32.qm';            % source-detector file
fmod_name = '../fwdfem/fmod_ellips_32x32_100MHz.fem'; % data file: log amplitude
farg_name = '../fwdfem/farg_ellips_32x32_100MHz.fem'; % data file: phase

refind = 1.4;                           % refractive index
bx = 100; by = 100;                     % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
tau = 1e-4;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolCG = 1e-4;                           % Gauss-Newton convergence criterion
resetCG = 10;                           % PCG reset interval
itrmax = 100;                           % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian matrix
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();

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
hBasis = toastSetBasis (hMesh, [bx by]);
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

figure(1);
subplot(1,2,1), imagesc(reshape(bmua,bx,by)), axis equal, axis tight, colorbar('horiz');
subplot(1,2,2), imagesc(reshape(bmus,bx,by)), axis equal, axis tight, colorbar('horiz');
drawnow

x = [scmua;sckap];
logx = log(x);
p = length(x);
step = 1.0; % initial step length for line search

% Initialise regularisation
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
hReg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta);

%tgtmua=toastReadNIM('../meshes/mua_tgt.nim');
%tgtmus=toastReadNIM('../meshes/mus_tgt.nim');
%tgtkap=1./(3*(tgtmua+tgtmus));
%hMesh2 = toastReadMesh('../meshes/ellips_tri10.msh');
%hBasis2 = toastSetBasis (hMesh2,[bx by]);
%btgtmua=toastMapMeshToGrid (hBasis2, tgtmua);
%btgtkap=toastMapMeshToGrid (hBasis2, tgtkap);
%toastDeleteBasis(hBasis2);
%toastDeleteMesh(hMesh2);
%btgt = [btgtmua; btgtkap];
%hReg = toastRegul ('TV', hBasis, logx, tau, 'Beta', beta, 'KapRefImage', ...
%                   btgt, 'KapRefScale', 1, 'KapRefPMThreshold', 0.1);

% initial data error (=2 due to data scaling)
err0 = privObjective (proj, data, sd, hReg, logx); %initial error
err = err0;                                        % current error
errp = inf;                                       % previous error
erri(1) = err0;
itr = 1; % iteration counter
fprintf (1, '\n**** INITIAL ERROR %f\n\n', err);

% NCG loop
while (itr <= itrmax) & (err > tolCG*err0) & (errp-err > tolCG)

    errp = err;
    
    % Gradient of cost function
    r = -toastGradient (hMesh, hBasis, qvec, mvec, mua, mus, ref, freq, ...
                       data, sd, 'bicgstab', 1e-16);
    r = r .* x; % parameter scaling
    
    rr = toastRegulGradient (hReg, logx);
    r = r - rr;
    
    figure(2);
    tmp = zeros(bx*by*2,1);
    tmp(solmask2) = -rr;
    subplot(2,2,1);imagesc(reshape(tmp(1:bx*by),bx,by));axis equal;axis tight;colorbar;title('Gradient of prior, \mu_a');
    subplot(2,2,2);imagesc(reshape(tmp(bx*by+1:bx*by*2),bx,by));axis equal;axis tight;colorbar;title('Gradient of prior, \kappa');
    tmp(solmask2) = r;
    subplot(2,2,3);imagesc(reshape(tmp(1:bx*by),bx,by));axis equal;axis tight;colorbar;title('Total gradient, \mu_a');
    subplot(2,2,4);imagesc(reshape(tmp(bx*by+1:bx*by*2),bx,by));axis equal;axis tight;colorbar;title('Total gradient, \kappa');    
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
    scmua = x(1:size(x)/2);
    sckap = x(size(x)/2+1:size(x));
    smua = scmua/cm;
    skap = sckap/cm;
    smus = 1./(3*skap) - smua;
    mua = toastMapSolToMesh (hBasis, smua);
    mus = toastMapSolToMesh (hBasis, smus);
    bmua(solmask) = smua;
    bmus(solmask) = smus;

    figure(1);
    subplot(1,2,1), imagesc(reshape(bmua,bx,by),[min(mua) max(mua)]), axis equal, axis tight, colorbar('horiz');colormap(gray);
    subplot(1,2,2), imagesc(reshape(bmus,bx,by),[min(mus) max(mus)]), axis equal, axis tight, colorbar('horiz');colormap(gray);
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
disp('recon2: finished')

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (hMesh, hBasis, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = privObjective (proj, data, sd, hReg, x);
    fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end

end
