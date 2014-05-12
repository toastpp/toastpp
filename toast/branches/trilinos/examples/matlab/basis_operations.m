% This script demonstrates the use of the toastBasis object to
% map functions between the unstructured mesh and a regular grid.

clear all
close all
toastSetVerbosity(1);
% Get some feedback

%% Constructing the basis mapper object

hmesh = toastMesh('../../test/2D/meshes/ellips_tri10.msh');
% the unstructured mesh

grd = [128,128];
% the dimensions of the grid

hbasis = toastBasis(hmesh,grd);
% create the basis mapper object

hbasis2 = toastBasis(hmesh,grd,'cubic');
% create a mapper object with a bicubic grid basis

%% Retrieving basis parameters

fprintf('Mesh     DOF: %5d\n', hbasis.nlen);
fprintf('Grid     DOF: %5d\n', hbasis.blen);
fprintf('Solution DOF: %5d\n', hbasis.slen);
fprintf('Grid dimensions: %d %d\n', hbasis.Dims);

%% Mapping functions between basis representations

mua = toastNim('../../test/2D/meshes/tgt_mua_ellips_tri10.nim');
subplot(1,3,1); hmesh.Display(mua,[0 0.05]);
% load a nodal image

bmua = hbasis.Map('M->B', mua);
subplot(1,3,2); imagesc(rot90(reshape(bmua,grd)),[0 0.05]);
axis equal tight; colorbar
% map to grid basis

% Now map the function into a different mesh
hmesh2 = toastMesh('../../test/2D/meshes/circle25_32.msh');
hbasis2 = toastBasis(hmesh2,grd);
mua2 = hbasis2.Map('B->M', bmua);
subplot(1,3,3); hmesh2.Display(mua2,[0 0.05]);
