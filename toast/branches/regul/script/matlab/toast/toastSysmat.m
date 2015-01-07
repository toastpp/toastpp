function p = toastSysmat
%toastSysmat          - Generate an FEM system matrix for frequency OT.
%
% Synopsis: S = toastSysmat (hMesh, mua, mus, ref, freq)
%           S = toastSysmat (hMesh, mua, mus, ref, freq, 'EL')
%           [S,B] = toastSysmat (hMesh, mua, mus, ref, freq)
%           [S,B,alpha] = toastSystmat (hMesh, mua, mus, ref, freq)
%
%    hMesh: mesh handle
%    mua:   nodal absorption coefficient [1/mm] (real column vector)
%    mus:   nodal reduced scattering coefficient [1/mm] (real column vector)
%    ref:   nodal refractive index (real column vector)
%    freq:  modulation frequency [MHz] (real)
%    'EL'   flag to indicate element basis
%    S:     system matrix (complex, sparse)
%    B      boundary part of the system matrix (complex, sparse)
%    alpha: boundary pre-factors, containing refractive index mismatch
%           information
%
% Returns a sparse complex n x n matrix (n: number of nodes in the mesh)
% containing system matrix S. S is the sum of three components: an
% absorption-dependent real matrix of products of shape functions, a
% diffusion-dependent real matrix of products of shape function
% derivatives, and a frequency-dependent imaginary matrix of products of
% shape functions.
%
% If the 'EL' flag is present, the system matrix is calculated on an
% element basis, rather than a node basis. Length 'n' of all parameter
% vectors in that case must be equal to the number of elements.
% Parameters are considered constant over each element, so parameter
% distribution has discontinuities at element boundaries.
%
% S is required for the evaluation of the frequency-domain FEM diffusion
% equation.
%
% If B is provided as an output argument, it is filled with the boundary
% component of the system matrix (but without the refractive index
% mismatch prefactor). The sparsity structure of B is the same as S.
%
% If alpha is provided as an output argument, it is filled with the
% boundary mismatch prefactors. (returned as a vector of size n).
%
% See also: toastMassmat.
