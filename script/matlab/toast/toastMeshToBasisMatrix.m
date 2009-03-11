function p = toastMeshToBasisMatrix
%toastMeshToBasisMatrix - Return mesh->basis transformation matrix.
%
% Synopsis: M = toastMeshToBasisMatrix (hBasis)
%    hbasis: basis mapper handle (see toastSetBasis)
%    M:      transformation matrix (b x h, real, sparse)
%
% Returns the linear transformation matrix that maps a nodal image
% (size n) into a regular grid (size b). The mapping is equivalent to
% the toastMapMeshToBasis function. Therefore, the operation
%
%    M = toastMeshToBasisMatrix (hBasis);
%    bimg = M * himg;
%
% has the same result as
%
%    bimg = toastMapMeshToBasis (hBasis, himg)
%
% (where himg is an image in nodal basis representation, and hBasis is
% a valid basis mapper handle).
