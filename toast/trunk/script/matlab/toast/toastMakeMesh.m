function p = toastMakeMesh
%toastMakeMesh        - Create a TOAST mesh from vertex and index lists.
%
% Synopsis: hMesh = toastMakeMesh (vtx, idx, eltp, [prm])
%    vtx:  vertex array
%    idx:  element index list
%    eltp: element type list
%    prm:  optional nodal parameter list
%
% vtx is a double array of dimension n x 2 or n x 3, where n is number of nodes.
% It contains the node coordinates.
%
% idx is a 1-based integer array of dimension m x v, where m is number of
% elements, and v is max. number of nodes in an element.
% It contains the global node indices for all element nodes.
%
% eltp is an integer array of dimension m x 1. It contains the element type
% identifiers for all mesh elements. Currently supported element types are:
%   1: 3-noded triangle (old format, deprecated)
%   3: 4-noded tetrahedron
%   4: 6-noded wedge
%   5: 8-noded regular voxel
%   6: 6-noded triangle
%   7: 10-noded tetrahedron
%   8: 6-noded isoparametric triangle
%   9: 10-noded triangle
%  10: 10-noded isoparametric triangle
%  11: 10-noded isoparametric tetrahedron
%  14: 4-noded regular pixel
%  15: 3-noded triangle
%  16: 3-noded 3-D surface triangle
%  17: 6-noded 3-D surface triangle
%
% prm is an optional double array of dimension n x 1..3, containing nodal
% values of absorption (column 1), scatter (column 2) and refractive index
% (column 3). If prm is not provided, or if it contains less than 3 columns,
% dummy values are inserted for the missing parameters.
%
% Notes:
%   - The element index list must be sorted in the order that TOAST expects.
%     (see element header files for local node order)
%   - toastMakeMesh labels boundary nodes according to their mesh connectivity.
%     To apply different boundary flags, and to add internal boundaries, use
%     the toastMarkMeshBoundary function. 
