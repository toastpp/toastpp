function p = toastMeshDimension
%toastMeshDimension   - Return the spatial dimension of a mesh (2 or 3).
%
% Synopsis: n = toastMeshDimension(hMesh)
%    hMesh: mesh handle
%    n:     mesh dimension (integer)
%
% Returns 2 for a 2-D mesh, 3 for a 3-D mesh. Mixed-dimensional meshes
% (such as surface meshes in 3-D space) are not supported by toast.
