% This example creates a simple voxel mesh from a sphere image.
% This mechanism can be used to quickly create toast meshes
% from 3D images (e.g. segmented MRI images).
%
% It employs the mkvoxmesh utility function which expects a
% binary input image (0=outside, 1=inside voxel).
%
% Note that voxel meshes can lead to ragged surfaces which
% may affect the accuracy of forward solvers, in particular
% close to the boundary.

grd = [64 64 64];

img = zeros(grd);
for i=1:grd(1)
    di = i-32.5;
    for j=1:grd(2)
        dj = j-32.5;
        for k=1:grd(3)
            dk = k-32.5;
            dst = sqrt(di^2 + dj^2 + dk^2);
            if dst < 32
                img(i,j,k) = 1;
            end
        end
    end
end

[vtx idx eltp] = mkvoxmesh(img);
n = size(vtx,1);
mua = ones(n,1)*0.01;
kap = ones(n,1)*0.3;
ref = ones(n,1)*1.4;
hmesh = toastMakeMesh(vtx,idx,eltp,[mua,kap,ref]);
toastShowMesh(hmesh);
