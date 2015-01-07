clear all
close all

% generate frequency domain data
mesh = toastMesh('../meshes/circle25_32.msh');
mesh.ReadQM('../meshes/circle25_1x32.qm');

n = mesh.NodeCount;
mua = ones(n,1)*0.025;
mus = ones(n,1)*2;
ref = ones(n,1)*1.4;
freq = 1000;   % high enough to cause phase wrap

smat = dotSysmat(mesh,mua,mus,ref,freq);
qvec = mesh.Qvec('Neumann','Gaussian',2);
mvec = mesh.Mvec('Gaussian',2);

phi = smat\qvec;
lphi = log(phi);
lnamp = real(lphi);
phase = imag(lphi);

subplot(1,3,1);
mesh.Display(phase); title('original');

% now unwrap with script
qp = mesh.Qpos;
uphase = unwrap_phase(mesh,phase,qp(1,:));

subplot(1,3,2);
mesh.Display(uphase); title('unwrapped script');

% now unwrap with mex file
uphase2 = mesh.UnwrapPhase (phase,qp(1,:));

subplot(1,3,3);
mesh.Display(uphase2); title('unwrapped mex');

