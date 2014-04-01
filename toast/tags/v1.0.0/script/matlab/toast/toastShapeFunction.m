function p = toastShapeFunction
%toastShapeFunction   - Return element shape functions at global point
%
% Synopsis: fun = toastShapeFunction(hMesh,idx,glob)
%    hMesh: mesh handle
%    idx:   element index (1-based)
%    glob:  vector of coordinates for global point
%    fun:   array of element shape functions
%
% Returns the values of the shape functions at a given point for all nodes
% of the element 'idx' of mesh 'hMesh'.
% The dimension of the returned vector is equal to the number of nodes in
% the element.
% The order of the shape function array is defined by the local element
% node order.
