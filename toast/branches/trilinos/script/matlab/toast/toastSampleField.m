function p = toastSampleField
%toastSampleField     - Interpolate nodal volume data on a list of points
%
% Synopsis: V = toastSampleField(hMesh,nval,absc)
%    hMesh: mesh handle
%    nval:  array of scalar nodal values (real)
%    absc:  list of n d-dimensional interpolation points (real n x d)
%    V:     list of n interpolation values
%
% Given a mesh hMesh and the nodal basis coefficients nval for the FEM
% approximation of some function f(r), this calculates the interpolated
% values of f at each of the requested sampling points, using the FEM
% shape functions.
%
% For sampling points that are not in the domain of the mesh, the closest
% element is used to provide the shape functions. This should work if the
% point is very close to the mesh boundary, but will produce unexpected
% results if the node is far away from the boundary.
