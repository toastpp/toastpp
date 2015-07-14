% Use gmesh to generate a volume mesh and import into toast

clear all
close all

% Path to gmsh. Edit as required
gmsh_cmd = '~/local/bin/gmsh';

% Generate mesh from geometry definition
system([gmsh_cmd ' -3 layered_sphere.geo -o gmesh.msh']);

% Read gmsh structure into Matlab
gmesh = gmshMeshReader('gmesh.msh');

% rescale mesh (original radius=0.5)
scale = 50;
gmesh.nodes = gmesh.nodes.*scale;

% Convert to toast mesh
nel = gmesh.nElems;
eltp = ones(nel,1)*3;
mesh = toastMesh(gmesh.nodes, gmesh.elems, eltp);

% Save mesh to file
mesh.Write('layered_sphere.msh');

% And display
mesh.Display;

% Display layers
elreg = gmesh.tags(:,1);
grd = [128 128 128];
basis = toastBasis(mesh,grd);
elref = basis.GridElref;

np = prod(grd);
img = zeros(np,1);
for i=1:np
    el = elref(i);
    if el>0
        img(i) = elreg(el);
    end
end
img = reshape(img,grd);
figure; imagesc(img(:,:,round(grd(3)/2))); axis equal tight