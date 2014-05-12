function p = toastSolutionMask
%toastSolutionMask    - Returns permutation vector for masking exterior.
%
% Synopsis: mask = toastSolutionMask (hBasis)
%    hBasis: basis mapper handle
%    mask:   solution mask permutation vector
%
% For irregular objects, some of the voxels of the regular image basis
% may be entirely outside the domain supported by the mesh. These pixels
% are excluded from the 'Solution' basis. To map between the full grid
% (filling the entire bounding box) and the solution basis, the
% permutation vector returned by toastSolutionMask can be used.
%
% The length of 'mask' is the same as that of the solution vector for a
% single parameter. For solution vectors composed of multiple parameter
% sets, to complete mask vector must be constructed from concatenated and
% appropriately offset single mask vectors.
%
% Example for building a 2-component mask vector:
%
%   mask_single = toastSolutionMask(hBasis);
%   mask_full = [mask_single, mask_single+blen];
%
% where 'blen' is the dimension of the full basis for a single parameter.
%
% Example for mapping between full basis and solution basis:
%
%   sol = basis(mask);  % full->sparse
%   basis(mask) = sol;  % sparse->full
%
% See also: toastMapMeshToSol, toastMapMeshToBasis
