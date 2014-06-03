function p = toastMapMeshToBasis
%toastMapMeshToBasis  - Map a nodal image to the solution basis.
%
% Synopsis: sol = toastMapMeshToBasis (hBasis, nim)
%    hbasis: basis mapper handle (see toastSetBasis)
%    nim:    nodal image (real vector)
%    sol:    solution basis coefficients of the mapped image
%
% Maps a nodal image into the solution basis defined by the basis mapper.
% The image can be real or complex-valued.
