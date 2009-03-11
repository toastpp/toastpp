function p = toastSetBasis
%toastSetBasis        - Define a mapper between FEM and reconstruction bases.
%
% Synopsis: hBasis = toastSetBasis (hMesh, [bdim], [gdim])
%    hMesh:    mesh handle
%    bdim:     integer vector with reconstruction basis dimensions.
%    gdim:     integer vector with sampling grid dimensions. (optional)
%
% The forward solution is computed on a finite element grid which may be unstructured.
% The inverse solution is normally calculated on an independent (usually regular)
% grid. This function sets up the dimensions of the solution grid, and defines the
% mapping between the bases.
% The sampling grid is an intermediate basis used to sample a function on the
% unstructured FEM mesh, before mapping it into the solution basis. The sampling
% grid has usually a higher resolution than the solution basis. If the sampling
% grid dimensions are not specified, the same dimensions as the solution basis are
% assumed.
