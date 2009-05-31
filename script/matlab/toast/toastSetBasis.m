function p = toastSetBasis
%toastSetBasis        - Define a mapper between FEM and reconstruction bases.
%
% Synopsis: hBasis = toastSetBasis (hMesh, bdim, [gdim], [bb])
%    hMesh:    mesh handle
%    bdim:     integer vector with reconstruction basis dimensions.
%    gdim:     integer vector with sampling grid dimensions. (optional)
%    bb:       grid bounding box, size dim x 2 (optional)
%
% The forward solution is computed on a finite element grid which may be unstructured.
% The inverse solution is normally calculated on an independent (usually regular)
% grid. This function sets up the dimensions of the solution grid, and defines the
% mapping between the bases.
%
% The sampling grid is an intermediate basis used to sample a function on the
% unstructured FEM mesh, before mapping it into the solution basis. The sampling
% grid has usually a higher resolution than the solution basis. If the sampling
% grid dimensions are not specified, the same dimensions as the solution basis are
% assumed.
%
% If bb is specified, it defines the bounding box for the grid basis. It
% should be of size dim x 2, where bb(:,1) contains the lower coordinates
% of the bounding box, and (bb(:,2) contains the upper coordinates.
%
% If bb is not specified, the bounding box is computed from the mesh dimensions.
  
