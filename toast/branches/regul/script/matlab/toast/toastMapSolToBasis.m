function p = toastMapSolToBasis
%toastMapSolToBasis   - Map an image from inverse solution basis to full grid.
%
% Synopsis: img = toastMapSolToBasis (hBasis, sol)
%    hbasis: basis mapper handle (see toastSetBasis)
%    sol:    parameter distribution in inverse solution basis (real vector)
%    img:    solution mapped on full raster grid.
%
% Maps an image represented in the inverse solution basis (SOL) to a regular
% grid that also contains entries for voxels that are not inside the support
% of the mesh. The added external voxels get assigned a value of zero.
% The image can be real or complex-valued.
%
% The full grid basis is useful for display purposes and similar.
