function p = toastSurfData
%toastSurfData        - Return mesh surface geometry.
%
% Synopsis: [vtx,idx,perm] = toastSurfData (hMesh)
%    hMesh: mesh handle
%    vtx:   matrix of n d-dimensional surface vertices (real n x d)
%    idx:   matrix of vertex indices for each surface element
%    perm:  permutation array for extracting surface data from a volume array.
%
% Returns the node coordinates and element vertex index list for the surface
% of a TOAST mesh.
%
% The returned vertex list is a real matrix of dimension n x d, where n
% is the number of surface nodes, and d is the mesh dimension.
%
% The returned index list is an integer matrix of dimension m x s, where
% m is the number of surface elements, and s is the maximum number of nodes
% per element.
%
% The permutation array is an n-dimensional integer array that can be used
% to extract surface data from an array of nodal volume data, e.g.
%     surfdata = voldata(perm)
%
% The index list is 1-based. For elements containing fewer than s nodes,
% the unused entries are set to 0.
