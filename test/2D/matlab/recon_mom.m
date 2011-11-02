% This example computes the first 3 moments of the temporal distribution
% of boundary measurements, and then uses them to reconstruct the optical
% parameters, starting from a homogeneous initial guess.

clear all
close all

meshdir = '../meshes/';
meshname = [meshdir 'ellips_tri10.msh'];
qmname = [meshdir 'circle25_32x32.qm'];
grd = [128 128];
nmom = 2; % highest moment to compute

% create the mesh
hmesh = toastReadMesh(meshname);
toastReadQM (hmesh,qmname);
n = toastMeshNodeCount (hmesh);

% create the basis mapper
hbasis = toastSetBasis (hmesh,grd);

% source and measurement vectors
qvec = toastQvec (hmesh, 'Neumann', 'Gaussian', 2);
mvec = toastMvec (hmesh, 'Gaussian', 2);
nq = size(qvec,2);
nm = size(mvec,2);
nqm = nq*nm;

% initial guess for the parameters
mua = ones(n,1)*0.025;
mus = ones(n,1)*2;
ref = ones(n,1)*1.4;

jac = toastJacobianMoments (hmesh,hbasis,nmom,qvec,mvec,mua,mus,ref,'DIRECT');

figure;
subplot(1,3,1);
imagesc(reshape(toastMapSolToGrid(hbasis,jac(8,:)),grd));axis equal tight
subplot(1,3,2);
imagesc(reshape(toastMapSolToGrid(hbasis,jac(8+nqm,:)),grd));axis equal tight
subplot(1,3,3);
imagesc(reshape(toastMapSolToGrid(hbasis,jac(8+nqm*2,:)),grd));axis equal tight

figure;
for mom=0:nmom
    mjac = jac(mom*nqm+1:(mom+1)*nqm,:);
    jmin = min(min(mjac));
    jmax = max(max(mjac));
    for i=1:nqm
        imagesc(reshape(toastMapSolToGrid(hbasis,mjac(i,:)),grd)); axis equal tight;colorbar
        title(['moment=' num2str(mom) ', meas=' num2str(i)])
        drawnow
        pause(0.1);
    end
end