function p = toastElmat
%toastElmat           - Return an element matrix for a single element of a mesh
%
% The returned matrix contains the integrals of products of a user-
% specified combination of nodal shape functions and shape function
% derivatives over the element volume or element surface. The dimension
% of the returned matrix depends on the requested integral. The element
% matrices are the building blocks for the FEM system matrices.
% 
%   Syntax
%
%   E = toastElmat (hMesh, elidx, int_type, [extra_param])
%
%   Parameters :
%
%   [In]
%    hMesh    :  handle      :  mesh handle
%    elidx    :  integer     :  element index (>= 1)
%    int_type :  string      :  integral type identifier (see below)
%
%   [Out]
%    E        : real array [dimension is type-dependent]:  element matrix
%
% Notes
%
%   The following integral types are recognised by this function 
%    Fi: shape function for element node i,
%    Di: shape function derivative for element node i,
%    dFi/dxj: partial derivative of shape function for element
%        node i with respect to coordinate j
%     n: number of nodes associated with the element,
%     d: domain dimension.
%
%       'type_string' : integral type      :Dimension of returned  matrix
%
%       'F'           : \int Fi dr         :  n x 1
%       'FF'          : \int Fi Fj dr      :  n x n
%       'FFF'         : \int Fi Fj Fk dr   :  n x n x n
%       'DD'          : \int Di Dj dr      :  n x n
%       'FD'          : \int Fi Dj dr      :  n x n x d
%       'FDD'         : \int Fi Dj Dk dr   :  n x n x n
%       'dd'          : \int dFi/dxj dFk/dxl dr : n x d x n x d
%       'BndF'        : \int Fi ds         :  n x 1
%       'BndFF'       : \int Fi Fj ds      :  n x n
%
% For boundary integrals ('BndF' and 'BndFF'), the integrals are
% performed over all boundary sides of the element (i.e. sides forming
% part of the mesh surface). Alternatively, by specifying a 4th parameter
% (sideidx) to the call to toastElmat, the integral can be performed over
% a single side of the element (whether boundary side or not):
%
%     E = toastElmat(hMesh, elidx, 'BndFF', sideidx)
%
% where sideidx is the local side index (>= 1). 
%
% The returned matrix is of dimension n x n, where n is the number of
% nodes in the element. Nonzero entries are located at positions (i,j)
% where both nodes i and j belong to a boundary side, or to side
% 'sideidx', if applicable.
