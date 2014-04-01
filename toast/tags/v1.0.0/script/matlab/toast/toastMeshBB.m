function p = toastMeshBB
%toastMeshBB          - Return mesh bounding box.
%
% Synopsis: [bbmin bbmax] = toastMeshBB(hMesh)
%    hMesh: mesh handle
%    bbmin: lower left corner of bounding box
%    bbmax: upper right corner of bounding box
%
% Returns the minimum and maximum corners of the bounding box enclosing
% the mesh. The bounding box is defined by the min and max vertex
% coordinates in each axis.
%
% Note that for isoparametric elements, the mesh bounding box may not
% be correctly represented by node positions. This is not taken into account.
