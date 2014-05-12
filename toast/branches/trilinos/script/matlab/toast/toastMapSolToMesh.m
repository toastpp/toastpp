function p = toastMapSolToMesh
%toastMapSolToMesh    - Map an image from inverse solution basis to mesh basis.
%
% Synopsis: nim = toastMapSolToMesh (hBasis, sol)
%    hbasis: basis mapper handle (see toastSetBasis)
%    sol:    parameter distribution in inverse solution basis (real vector)
%    nim:    nodal parameter distribution (real vector)
%
% Maps an image represented in the inverse solution basis (SOL) to a nodal
% image (NIM).
% The image can be real or complex-valued.
