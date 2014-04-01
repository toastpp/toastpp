function p = toastMeshSize
%toastMeshSize        - Return the size of a TOAST mesh.
%
% Synopsis: s = toastMeshSize (hMesh)
%    hMesh: mesh handle
%    s:     mesh size (real)
%
% Returns the area (2-D meshes) or volume (3-D meshes) of a TOAST mesh.
% The mesh size is the sum of all element sizes.
