% Generate the error surface for a parameter perturbation as a function
% of absorption and scatter coefficient of the perturbation.

clear all
close all

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = '../meshes/circle25_32.msh';              % mesh file
qmname    = '../meshes/circle25_32x32.qm';            % source-detector file
refind = 1.4;                           % refractive index
bx = 100; by = 100;                     % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
tau = 1e-4;                             % regularisation parameter
beta = 0.01;                            % TV regularisation parameter
tolGN = 1e-4;                           % Gauss-Newton convergence criterion
tolKrylov = 1e-4;                       % Krylov convergence criterion
itrmax = 100;                           % Gauss-Newton max. iterations
Himplicit = true;                       % Implicit/explicit Hessian
                                        % matrix
mua0 = 0.025;
mus0 = 2;
muamin = 0.005;
muamax = 0.05;
musmin = 0.5;
musmax = 5;
pert_cnt = [(bx*2)/3 (by*2)/3];
pert_rad = bx*0.1;

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
bmua = ones(blen,1) * 0.025;
bmus = ones(blen,1) * 2;
bref = ones(blen,1) * refind;
bkap = 1./(3*(bmua+bmus));

% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis (hMesh, [bx by]);
solmask = toastSolutionMask (hBasis);
solmask2 = [solmask solmask+blen];

mua = toastMapGridToMesh (hBasis, bmua);
mus = toastMapGridToMesh (hBasis, bmus);
ref = ones(n,1) * refind;

bmua = reshape(bmua,bx,by);
bmus = reshape(bmus,bx,by);

% Generate source vectors
qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);

% Generate measurement vectors
mvec = toastMvec (hMesh, 'Gaussian', 2);

% Initial data set f[x0]
smat = toastSysmat (hMesh, mua, mus, ref, freq);
lgamma0 = reshape (log(mvec.' * (smat\qvec)), [], 1);
lgamma0 = lgamma0(dmask);
proj0 = [real(lgamma0);imag(lgamma0)];

% data scaling
sd_lnmod = ones(size(lgamma0)) * norm(real(lgamma0)); 
sd_phase = ones(size(lgamma0)) * norm(imag(lgamma0));
sd = [sd_lnmod;sd_phase];

% Scan over optical parameters of perturbation
h = waitbar(0,'Generating perturbation data');
for i=0:20
    pmua = muamin + i/20*(muamax-muamin);
    for j = 0:20
        pmus = musmin+j/20*(musmax-musmin);
    
        for ii = 1:bx
            for jj = 1:by
                dst = sqrt ((ii-pert_cnt(1))^2 + (jj-pert_cnt(2))^2);
                if (dst <= pert_rad)
                    bmua(ii,jj) = pmua;
                    bmus(ii,jj) = pmus;
                end
            end
        end
    
        mua = toastMapGridToMesh (hBasis, bmua);
        mus = toastMapGridToMesh (hBasis, bmus);
        smat = toastSysmat (hMesh, mua, mus, ref, freq);
        lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);
        lgamma = lgamma(dmask);
        proj = [real(lgamma);imag(lgamma)];

        err(i+1,j+1) = sum(((proj0-proj)./sd).^2);
        erra(i+1,j+1) = sum(((real(lgamma0)-real(lgamma))./sd_lnmod).^2);
        errp(i+1,j+1) = sum(((imag(lgamma0)-imag(lgamma))./sd_phase).^2);
        waitbar((i*20+j)/20^2)
    end
end
close(h)

subplot(1,3,1);contour(erra,60);title('log amplitude');
subplot(1,3,2);contour(errp,60);title('phase');
subplot(1,3,3);contour(err,60);title('combined');
