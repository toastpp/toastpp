% Using the high-level convenience method dotFwd to compute a forward
% solution for the DOT problem.

meshpath = '../../test/2D/meshes/';
meshfile = [meshpath 'ellips_tri10.msh'];
qmfile   = [meshpath 'circle25_32x32.qm'];


%% Set up the parameter structure
prm = toastParam;
prm.fwdsolver.meshfile = meshfile;
prm.fwdsolver.method = 'direct';

prm.data.freq = 100;

prm.meas.qmfile = qmfile;
prm.meas.src.type = 'Neumann';
prm.meas.src.prof = 'Gaussian';
prm.meas.src.width = 2;
prm.meas.det.prof = 'Gaussian';
prm.meas.det.width = 2;

prm.initprm.mua.reset = 'nim';
prm.initprm.mua.nim = [meshpath 'tgt_mua_ellips_tri10.nim'];
prm.initprm.mus.reset = 'nim';
prm.initprm.mus.nim = [meshpath 'tgt_mus_ellips_tri10.nim'];
prm.initprm.ref.reset = 'homog';
prm.initprm.ref.val = 1.4;


%% Call the solver
[mdata,pdata] = dotFwd (prm);


%% Display the data
figure;
subplot(1,2,1);

imagesc(reshape(mdata,32,32));
axis equal tight
colorbar
xlabel('source #');
ylabel('detector #');
title('log amp.');

subplot(1,2,2);
imagesc(reshape(pdata,32,32));
axis equal tight
colorbar
xlabel('source #');
ylabel('detector #');
title('phase');
