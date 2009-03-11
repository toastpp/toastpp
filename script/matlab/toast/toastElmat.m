function p = toastElmat
% toastElmat    - Return an element matrix for a single element of a mesh
% The returned matrix contains the integral of products of shape functions
% or shape function derivatives over the element volume. The dimension of the
% returned matrix depends on the requested integral. The element matrices are
% the building blocks for the FEM system matrices.
% 
%   Syntax
%
%   E = toastElmat (hMesh, elidx, int_type)
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
%       'dd'          : 
%       'BndFF'       :

