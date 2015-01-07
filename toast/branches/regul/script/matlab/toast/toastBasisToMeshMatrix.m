function p = toastBasisToMeshMatrix
%toastBasisToMeshMatrix - Return basis->mesh transformation matrix.
%
% Synopsis: M = toastBasisToMeshMatrix (hBasis)
%    hbasis: basis mapper handle (see toastSetBasis)
%    M:      transformation matrix (h x b, real, sparse)
%
% Returns the linear transformation matrix that maps a regular grid image
% (size b) into a nodal image (size n). The mapping is equivalent to
% the toastMapBasisToMesh function. Therefore, the operation
%
%    M = toastBasisToMeshMatrix (hBasis);
%    himg = M * bimg;
%
% has the same result as
%
%    himg = toastMapBasisToMesh (hBasis, bimg)
%
% (where bimg is an image in grid basis representation, and hBasis is
% a valid basis mapper handle).
