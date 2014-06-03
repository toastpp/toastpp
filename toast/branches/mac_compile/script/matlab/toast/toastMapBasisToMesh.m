function p = toastMapBasisToMesh
%toastMapBasisToMesh  - Map an image from the solution basis to the mesh.

%
% Synopsis: nim = toastMapBasisToMesh (hBasis, sol)
%    hbasis: basis mapper handle (see toastSetBasis)
%    sol:    solution basis coefficients
%    nim:    nodal image (real vector)
%
% Maps an image defined in the regular grid of the inverse basis into a
% nodal image of an unstructured mesh, using a basis mapper.
% The image can be real or complex-valued.
