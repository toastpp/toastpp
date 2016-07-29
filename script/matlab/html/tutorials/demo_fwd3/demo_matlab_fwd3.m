% Toast-Matlab web example 3: 
% Frequency-domain forward solver for a circular 2D inhomogeneous problem.
% (c) Martin Schweiger and Simon Arridge
% www.toastplusplus.org

clear all
close all

freq = 100; % modulation frequency [MHz]

% Create the mesh
rad = 25;
nsect = 6;
nring = 32;
nbnd = 2;
[vtx,idx,eltp] = mkcircle (rad, nsect, nring, nbnd);
mesh = toastMesh(vtx,idx,eltp);

% Load bitmaps for parameter distributions
bmua = imread('demo_matlab_fwd2_mua.png');
bmus = imread('demo_matlab_fwd2_mus.png');

% Scale to desired parameter ranges
mua_bkg = 0.01;
mus_bkg = 1.0;
bmua = double(bmua)./255.*0.02 + mua_bkg;
bmus = double(bmus)./255.*1.0 + mus_bkg;
figure;
subplot(1,2,1);
imagesc(rot90(bmua));
axis equal tight; colorbar
title('\mu_a');
subplot(1,2,2);
imagesc(rot90(bmus));
axis equal tight; colorbar
title('\mu_s');

% Map to mesh basis
grd = size(bmua);
basis = toastBasis (mesh,grd);
mua = basis.Map ('B->M',bmua);
mus = basis.Map ('B->M',bmus);
figure;
mesh.Display(mua);
figure;
mesh.Display(mus);

nnd = mesh.NodeCount;
ref_bkg = 1.4;
ref = ones(nnd,1) * ref_bkg;

% Create the source and detector positions
nq = 16;
for i=1:nq
  phi_q = 2*pi*(i-1)/nq;
  Q(i,:) = rad * [cos(phi_q) sin(phi_q)];
  phi_m = 2*pi*(i-0.5)/nq;
  M(i,:) = rad * [cos(phi_m) sin(phi_m)];
end
mesh.SetQM(Q,M);
hold on
plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
plot(M(:,1),M(:,2),'bs','MarkerFaceColor','b');

% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2, ref);

% Solve the FEM linear system
K = dotSysmat(mesh,mua,mus,ref,freq);
Phi = K\qvec;
Y = mvec.' * Phi;
logY = log(Y);
lnamp = real(logY);
phase = imag(logY);

% For reference, we also generate data for the
% homogeneous problem
mua0 = ones(nnd,1) * mua_bkg;
mus0 = ones(nnd,1) * mus_bkg;
K0 = dotSysmat (mesh,mua0,mus0,ref,freq);
logY0 = log(mvec.' * (K0\qvec));
lnamp0 = real(logY0);
phase0 = imag(logY0);

% Display sinograms
figure
subplot(1,2,1);
imagesc(lnamp);
xlabel('source index q');
ylabel('detector index m');
title('log amplitude');
axis equal tight;
colorbar
subplot(1,2,2);
imagesc(phase);
xlabel('source index q');
ylabel('detector index m');
title('phase');
axis equal tight;
colorbar

figure
subplot(1,2,1);
imagesc(lnamp-lnamp0);
xlabel('source index q');
ylabel('detector index m');
title('log amplitude difference');
axis equal tight;
colorbar
subplot(1,2,2);
imagesc(phase-phase0);
xlabel('source index q');
ylabel('detector index m');
title('phase difference');
axis equal tight;
colorbar

% Display boundary profile
figure
subplot(1,2,1);
hold on
angle = [360/32:360/16:360];
for i=1:size(lnamp,2)
    ywrap = [lnamp(i:end,i); lnamp(1:i-1,i)];
    plot(angle,ywrap,'o-');
end
axis([0 360 -15 -3]);
xlabel('src-det separation [deg]');
ylabel('log intensity');

subplot(1,2,2);
hold on
angle = [360/32:360/16:360];
for i=1:size(phase,2)
    ywrap = [phase(i:end,i); phase(1:i-1,i)];
    plot(angle,ywrap,'o-');
end
axis([0 360 -1.2 0]);
xlabel('src-det separation [deg]');
ylabel('phase');

% Write solver results to file
%data = reshape(log(Y'),[],1);
