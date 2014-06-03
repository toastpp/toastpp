functionn p = toastElementData
%toastElementData     - Return geometry of a mesh element.
%
% Synopsis: [vtx,idx] = toastElementData(hMesh,el)
%    hMesh: mesh handle
%    el:    element index (>= 1)
%    vtx:   matrix of n d-dimensional element vertices (real n x d)
%    idx:   matrix of vertex indices for each face
%
% Returns the node coordinates and element vertex index lists for each
% face of a mesh element.
%
% The returned vertex list is a real matrix of dimension n x d, where n
% is the number of nodes in the element, and d is the element dimension.
%
% The returned index list is an integer matrix of dimension m x s, where
% m is the number of element faces, and s is the maximum number of nodes
% associated with a face.
% The index list is 1-based. For faces containing fewer than s nodes,
% the unused entries are set to 0.
