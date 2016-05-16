% Toast-Matlab web example 5:
% Non-unique problem: attempt to reconstruct absorption and scattering
% distribution from steady-state intensity data
% (c) Martin Schweiger and Simon Arridge
% www.toastplusplus.org

clear all
close all

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
bmua = double(bmua)./255.*0.02 + 0.01;
bmus = double(bmus)./255.*1.0 + 1.0;
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
figure; mesh.Display(mua);
figure; mesh.Display(mus);

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
mvec = mesh.Mvec ('Gaussian', 2);

% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,0);
Phi = K\qvec;
Y = mvec.' * Phi;

% Display sinogram
figure
imagesc(log(Y));
xlabel('source index q');
ylabel('detector index m');
axis equal tight;
colorbar

% Display boundary profile
figure
hold on
angle = [360/32:360/16:360];
for i=1:size(Y,2)
    ywrap = [Y(i:end,i); Y(1:i-1,i)];
    plot(angle,log(ywrap),'o-');
end
axis([0 360 -14 -2]);
xlabel('angular source-detector separation');
ylabel('log intensity');

% Write solver results to file
data = reshape(log(Y'),[],1);
toastWriteVector('demo_matlab_fwd2.dat', data);

% Show differences to homogeneous results
data_homog = toastReadVector('demo_matlab_fwd1.dat');
logYhomog = reshape(data_homog,nq,nq)';
dlogY = log(Y)-logYhomog;

figure
imagesc(dlogY);
xlabel('source index q');
ylabel('detector index m');
axis equal tight;
colorbar

figure
hold on
for i=1:size(Y,2)
    ywrap = [dlogY(i:end,i); dlogY(1:i-1,i)];
    plot(angle,ywrap,'o-');
end
axis([0 360 -1 0.1]);
xlabel('angular source-detector separation');
ylabel('log intensity perturbation');
