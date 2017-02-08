% Generate frequency DOT boundary data for a simple 2D circle,
% given parameter distributions and source/detector layout


clear all
close all

% Load mesh and optode layout
meshdir = '../../test/2D/meshes/';
hmesh = toastMesh([meshdir 'ellips_tri10.msh']);
hmesh.ReadQM([meshdir 'circle25_32x32.qm']);

% Load parameter distributions
n = hmesh.NodeCount;
mua = toastNim([meshdir 'tgt_mua_ellips_tri10.nim']);
mus = toastNim([meshdir 'tgt_mus_ellips_tri10.nim']);
ref = ones(n,1)*1.4;
freq = 100;

% Construct source and measurement vectors
qvec = hmesh.Qvec('Neumann','Gaussian',2);
mvec = hmesh.Mvec('Gaussian',2,ref);

% Solve FEM linear system
smat = dotSysmat(hmesh,mua,mus,ref,freq);
phi = smat\qvec;
meas = log(mvec.' * phi);

% Plot the log amplitude and phase sinograms
figure;
subplot(1,2,1);
imagesc(real(meas));
axis equal tight
colorbar
xlabel('source #');
ylabel('detector #');
title('log amp.');

subplot(1,2,2);
imagesc(imag(meas));
axis equal tight
colorbar
xlabel('source #');
ylabel('detector #');
title('phase');
