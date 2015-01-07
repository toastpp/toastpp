function p = toastMapMeshToSol
%toastMapMeshToSol    - Map a nodal image to the solution basis.
%
% Synopsis: sol = toastMapMeshToSol (hBasis, nim)
%    hbasis: basis mapper handle (see toastSetBasis)
%    nim:    nodal image (real vector)
%    sol:    image in solution basis (real vector)
%
% Maps a nodal image into into the solution basis. Unlike
% toastMapMeshToBasis, it excludes voxels that are outside the domain
% supported by the mesh.
% The image can be real or complex-valued.
%
% See also: toastMapMeshToBasis, toastSolutionMask  
