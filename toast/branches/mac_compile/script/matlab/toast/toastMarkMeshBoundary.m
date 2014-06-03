function p = toastMarkMeshBoundary
%toastMarkMeshBoundary - Mark boundary nodes in a mesh.
%
% Synopsis: bnd = toastMarkMeshBoundary(hMesh,[bnd])
%    hMesh: mesh handle
%    bnd:   list of boundary flags (optional for input and output)
%
% Many mesh operations require the boundary nodes to be identified
% in the mesh. If a mesh does not contain boundary information,
% toastMarkMeshBoundary can be used to supply it. In addition, the
% function can be used to label internal boundaries.
%
% If toastMarkMeshBoundary is called without the second input
% parameter, boundary nodes are identified automatically according
% to their mesh connectivity.
%
% If an explicit boundary label list is provided, this is used
% instead. The length of bnd should be equal to the number of
% nodes (see toastMeshNodeCount), and contain the following entries:
%
%   0 - internal node
%   1 - surface node
%   2 - internal boundary node
%
% In all cases, the list of boundary flags is returned.
