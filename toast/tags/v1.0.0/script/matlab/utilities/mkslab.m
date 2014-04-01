function [vtx idx eltp] = mkslab(bb,gdim,varargin)

% Create a toast mesh for a slab geometry consisting of regular
% voxel elements.
%
% Synopsis: [vtx idx eltp] = mkslab (bb,gdim[,p])
%
% bb:    bounding box (mm): 2x3 real matrix containing min and
%        max coordinates of the slab
% gdim:  1x3 integer vector: contains element grid dimensions
% p:     optional matrix of nodal plane positions
% vtx:   per-vertex coodinate list
% idx:   per-element vertex index list
% eltp:  element type list
%
% Note: The grid dimensions refer to the element grid. Node
%     dimensions are element dimensions+1 in each direction.
%     Therefore, gdim=[10 10 10] will create a mesh with 1000
%     elements and 1331 nodes.
%
%     If p is provided, it must be a matrix of size 3 x max(gdim)+1,
%     containing the coordinates of the node planes of the slab, where
%     p(1,:) contains the x-coordinates, p(2,:) the y-coordinates, and
%     p(3,:) the z-coordinates.
%
%     Use the toastMakeMesh function to convert the vtx, idx and eltp
%     lists to a toast mesh.
%
% Usage example:
%     [vtx idx eltp] = mkslab ([[0 0 0];[10 10 10]], [10 10 10]);
%     n = size(vtx,1);
%     mua = ones(n,1)*0.01;
%     kap = ones(n,1)*0.3;
%     ref = ones(n,1)*1.4;
%     hmesh = toastMakeMesh(vtx,idx,eltp,[mua,kap,ref]);
%     toastWriteMesh (hMesh,'slab.msh');

dim = 3;

% Build vertex coordinate matrix

if nargin < 3
    p = zeros(dim,max(gdim)+1);
    for d = 1:dim
        p(d,:) = [[0:gdim(d)] * (bb(2,d)-bb(1,d)) / gdim(d) + bb(1,d), ...
                  zeros(1,max(gdim)-gdim(d))];
    end
else
    p = varargin{1};hmesh = toastMakeMesh(vtx,idx,eltp,[mua,kap,ref]);
end

vdim = gdim+1;
vtx_x = zeros(vdim);
vtx_y = zeros(vdim);
vtx_z = zeros(vdim);

for i=1:size(vtx_x,1), vtx_x(i,:,:) = p(1,i); end
for j=1:size(vtx_y,2), vtx_y(:,j,:) = p(2,j); end
for k=1:size(vtx_z,3), vtx_z(:,:,k) = p(3,k); end

vtx_x = reshape(vtx_x,[],1);
vtx_y = reshape(vtx_y,[],1);
vtx_z = reshape(vtx_z,[],1);

vtx = [vtx_x vtx_y vtx_z];

% Build element index matrixhmesh = toastMakeMesh(vtx,idx,eltp,[mua,kap,ref]);

idx = zeros(prod(gdim),8);
dy = vdim(1);
dz = vdim(1)*vdim(2);
ii = 1;

for k = 1:gdim(3)
    for j = 1:gdim(2)
        for i = 1:gdim(1)
            vidx0 = (i-1) + (j-1)*dy + (k-1)*dz + 1;
            idx(ii,1) = vidx0;
            idx(ii,2) = vidx0+1;
            idx(ii,3) = vidx0+dy;
            idx(ii,4) = vidx0+dy+1;
            idx(ii,5) = vidx0+dz;
            idx(ii,6) = vidx0+dz+1;
            idx(ii,7) = vidx0+dz+dy;
            idx(ii,8) = vidx0+dz+dy+1;
            ii = ii+1;
        end
    end
end

% Build element type array

eltp = ones(prod(gdim),1) * 5;

% Build parameter matrix

%mua = 0.025;
%mus = 2;
%kap = 1./(3*(mua+mus));
%ref = 1.4;

%prm = ones(prod(vdim),3);
%prm(:,1) = mua;
%prm(:,2) = kap;
%prm(:,3) = ref;

% Make mesh

%hMesh = toastMakeMesh(vtx,idx,eltp,prm);
