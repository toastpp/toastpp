% Toast-Matlab web example 4: 
% Time-dependent diffuse light propagation forward problem
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

theta = 0.5; % time coupling (0=explicit, 1/2=Crank-Nicholson, 1=implicit)
dt = 2;      % time step interval [ps]
nstep = 3000; % number of time steps
t = [1:nstep]*dt;

% Set up the required FEM matrices
K = real(dotSysmat(mesh, mua, mus, ref, 0)); % stiffness matrix
M = mesh.Massmat;                            % mass matrix
K0 = -(K * (1-theta) - M * 1/dt);            % backward difference matrix
K1 = K * theta + M * 1/dt;                   % forward difference matrix

[L U] = lu(K1);               % LU factorisation

% initial condition
q = qvec/dt;
phi = U\(L\q);
gamma_t = mvec.' * phi;
gamma(1,:) = gamma_t(1,[1,4,8]); % pick out three sample detectors for display

% loop over time steps
h = waitbar(0,'loop over time steps');
for i=2:nstep
    q = K0 * phi;              % new source vector from current field
    phi = U\(L\q);             % new field
    gamma_t = mvec.' * phi;    % project to boundary
    gamma(i,:) = gamma_t(1,[1,4,8]);
    waitbar(i/nstep,h);
end
close(h);

lgamma = log(max(gamma,1e-20));
figure; plot(t,lgamma)
xlabel('time [ps]');
ylabel('log intensity');
legend('detector 1','detector 2','detector 3');