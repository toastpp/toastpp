function [vtx,idx,eltp] = mkrect(bb,gdim)

% Create a toast mesh for a rectangular geometry consisting of regular
% pixel elements.
%
% Syntax: [vtx,idx,eltp] = mkrect (bb,gdim)
%
% Parameters:
%         bb:  bounding box (mm): 2x2 real matrix containing min and
%              max coordinates of the slab
%         gdim:  1x2 integer vector: contains grid dimensions
%
% Return values:
%         vtx: list of mesh vertices
%         idx: list of element node indices
%         eltp: list of element type indices
%
% Usage example:
%     [vtx,idx,eltp] = mkslab ([[0 0];[10 10]], [10 10]);
%     mesh = toastMesh(vtx,idx,eltp);

dim = 2;

% Build vertex coordinate matrix

p = zeros(dim,max(gdim)+1);
for d = 1:dim
    p(d,:) = [[0:gdim(d)] * (bb(2,d)-bb(1,d)) / gdim(d) + bb(1,d), ...
              zeros(1,max(gdim)-gdim(d))];
end

vdim = gdim+1;
vtx_x = zeros(vdim);
vtx_y = zeros(vdim);

for i=1:size(vtx_x,1), vtx_x(i,:,:) = p(1,i); end
for j=1:size(vtx_y,2), vtx_y(:,j,:) = p(2,j); end

vtx_x = reshape(vtx_x,[],1);
vtx_y = reshape(vtx_y,[],1);

vtx = [vtx_x vtx_y];

% Build element index matrix

idx = zeros(prod(gdim),8);
dy = vdim(1);
ii = 1;

for j = 1:gdim(2)
    for i = 1:gdim(1)
        vidx0 = (i-1) + (j-1)*dy + 1;
        idx(ii,1) = vidx0;
        idx(ii,2) = vidx0+1;
        idx(ii,3) = vidx0+dy;
        idx(ii,4) = vidx0+dy+1;
        ii = ii+1;
    end
end

% Build element type array

eltp = ones(prod(gdim),1) * 14;

% Build parameter matrix

% mua = 0.025;
% mus = 2;
% kap = 1./(3*(mua+mus));
% ref = 1.4;
% 
% prm = ones(prod(vdim),3);
% prm(:,1) = mua;
% prm(:,2) = kap;
% prm(:,3) = ref;

% Make mesh
% hMesh = toastMakeMesh(vtx,idx,eltp,prm);
