function p = toastMeshData
%toastMeshData        - Return mesh geometry.
%
% Synopsis: [vtx,idx,eltp] = toastMeshData(hMesh)
%    hMesh: mesh handle
%    vtx:   matrix of n d-dimensional vertices (real n x d)
%    idx:   matrix of vertex indices for each element (m x s)
%    eltp:  list of element types (length m)
%
% Returns the node coordinates and element vertex index list for a mesh.
%
% The returned vertex list is a real matrix of dimension n x d, where n
% is the number of nodes in the mesh, and d is the mesh dimension.
%
% The returned index list is an integer matrix of dimension m x s, where
% m is the number of mesh elements, and s is the maximum number of nodes
% per element.
% The index list is 1-based. For elements containing fewer than s nodes,
% the unused entries are set to 0.
%
% For a list of element type indices, see documentation of function
% toastMakeMesh.
