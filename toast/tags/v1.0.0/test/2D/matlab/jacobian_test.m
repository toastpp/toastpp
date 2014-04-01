% MATLAB-TOAST sample script:
% 2D image reconstruction with nonlinear conjugate gradient solver

clear all
close all

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = '../meshes/circle25_32.msh';              % mesh file
qmname    = '../meshes/circle25_1x1.qm';            % source-detector file

refind = 1.4;                           % refractive index
grd = [64 64];                          % solution basis: grid dimension
freq = 100;                             % modulation frequency [MHz]
test_meshbasis = false;
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();
tol = 1e-12;

disp('jacobian_test');
disp('Test the integrity of toastJacobian by comparing to explicit');
disp('calculation obtained by pixel-wise parameter perturbation.');
% Set up some variables
blen = prod(grd);
c0 = 0.3;
cm = c0/refind;

% Read a TOAST mesh definition from file.
hMesh = toastReadMesh (meshname);
toastReadQM (hMesh, qmname);
n = toastMeshNodeCount (hMesh);
dmask = toastDataLinkList (hMesh);

% Set up the mapper between FEM and solution bases
hBasis = toastSetBasis (hMesh, grd);
solmask = toastSolutionMask (hBasis);
solmask2 = [solmask solmask+blen];
slen = length(solmask);

% Set up homogeneous initial parameter estimates
bmua0 = ones(blen,1) * 0.025;
bmus0 = ones(blen,1) * 2;
mua0 = toastMapBasisToMesh(hBasis,bmua0);
mus0 = toastMapBasisToMesh(hBasis,bmus0);
ref = ones(n,1) * refind;

% Generate source vectors
qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);

% Generate measurement vectors
mvec = toastMvec (hMesh, 'Gaussian', 2);

% Initial data set f[x0]
proj0 = toastProject (hMesh,mua0,mus0,ref,freq,qvec,mvec);

dmua = 1e-5;


%% Test 1: Jacobian in mesh basis
if test_meshbasis == true
    
J = toastJacobian (hMesh, 0, qvec, mvec, mua0, mus0, ref, freq, 'direct');

h = waitbar(0,'Calculating explicit Jacobian');
for i=1:n
    mua = mua0;
    mua(i) = mua0(i)+dmua;
    proj = toastProject (hMesh,mua,mus0,ref,freq,qvec,mvec);
    dy = (proj-proj0)/dmua;
    Je(:,i) = dy / cm; % parameter is c*mua
    waitbar(i/n);
end
delete(h);

mua_d = J(1,1:n);
mua_e = Je(1,1:n);
mn = min([mua_d mua_e]);
mx = max([mua_d mua_e]);

figure(1);
img = reshape(toastMapMeshToBasis(hBasis,mua_d),grd);
subplot(1,3,1);imagesc(img,[mn mx]); axis equal tight; title('direct'); colorbar

img = reshape(toastMapMeshToBasis(hBasis,mua_e),grd);
subplot(1,3,2);imagesc(img,[mn mx]); axis equal tight; title('explicit'); colorbar

ratio = mua_d./mua_e;
img = reshape(toastMapMeshToBasis(hBasis,ratio),grd);
subplot(1,3,3);imagesc(img);axis equal tight;colorbar;title('ratio')

mean(mean(img(8:24,8:24)))
clear J Je
end

%% Test 2: Jacobian in grid basis
J = toastJacobian (hMesh, hBasis, qvec, mvec, mua0, mus0, ref, freq, 'direct');
smat = toastSysmat (hMesh, mua0, mus0, ref, freq);
[L U] = luinc(smat,1e-4);

h = waitbar(0,'Calculating explicit Jacobian');
tic;
for i=1:slen
    idx = solmask(i);
    bmua = bmua0;
    bmua(idx) = bmua0(idx)+dmua;
    mua = toastMapBasisToMesh(hBasis,bmua);
    smat = toastSysmat (hMesh, mua, mus0, ref, freq);
    %[phi,flag] = bicgstab(smat,qvec,tol,1000,L,U);
    phi = smat\qvec;
    cprj = log(mvec.' * phi);
    proj = [real(cprj);imag(cprj)];
    %proj = toastProject (hMesh,mua,mus0,ref,freq,qvec,mvec);
    dy = (proj-proj0)/dmua;
    Je(:,i) = dy / cm;
    waitbar(i/slen);
end
toc
delete(h);

mua_d = J(1,1:slen);
mua_e = Je(1,1:slen);
mn = min([mua_d mua_e]);
mx = max([mua_d mua_e]);

[bbmin bbmax] = toastMeshBB(hMesh);
elsize = prod ((bbmax-bbmin) ./ grd');
mua_d = mua_d * elsize; % scale with element size

figure(2);
img = reshape(toastMapSolToBasis(hBasis,mua_d),grd);
subplot(1,3,1);imagesc(img,[mn mx]); axis equal tight; title('direct'); colorbar

img = reshape(toastMapSolToBasis(hBasis,mua_e),grd);
subplot(1,3,2);imagesc(img,[mn mx]); axis equal tight; title('explicit'); colorbar

ratio = mua_d./mua_e;
img = reshape(toastMapSolToBasis(hBasis,ratio),grd);
subplot(1,3,3);imagesc(img);axis equal tight;colorbar;title('ratio')

mean (ratio)

toastDeleteBasis(hBasis);
toastDeleteMesh(hMesh);