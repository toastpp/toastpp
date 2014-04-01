function p = toastElementSize
%toastElementSize     - Return an array with mesh element sizes.
%
% Synopsis: s = toastElementSize(hMesh)
%    hMesh: mesh handle
%    s:     array of element sizes (real)
%
% Returns the areas (2-D meshes) or volumes (3-D meshes) of each element
% of a TOAST mesh in a column vector.
% Note that negative element sizes usually indicate an incorrect node order.
