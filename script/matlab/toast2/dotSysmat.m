function smat = dotSysmat (hmesh,mua,mus,ref,freq)
% Generate an FEM system matrix for frequency DOT.
%
% Syntax: S = dotSysmat (mesh, mua, mus, ref, freq)
%         S = dotSysmat (mesh, mua, mus, ref, freq, 'EL')
%         [S,B] = dotSysmat (mesh, mua, mus, ref, freq)
%         [S,B,alpha] = dotSysmat (mesh, mua, mus, ref, freq)
%
% Parameters:
%         mesh [object]:
%             toastMesh object
%         mua [real column vector]:
%             nodal absorption coefficient [1/mm]
%         mus [real column vector]:
%             nodal reduced scattering coefficient [1/mm]
%         ref [real column vector]:
%             nodal refractive index
%         freq [real]:
%             modulation frequency [MHz]
%         'EL'
%             flag to indicate element basis
%
% Return values:
%         S [complex sparse matrix n x n]:
%             system matrix
%         B [complex sparse matrix n x n]:
%             boundary part of the system matrix
%         alpha [real vector n]:
%             boundary pre-factors, containing refractive index mismatch
%             information
%
% Notes:  This is a convenience function for DOT problems, which removes
%         the need for multiple calls to toastMesh.SysmatComponent for
%         assembling the components of the DOT diffusion forward model.
%
%         Returns a sparse complex n x n matrix (n: number of nodes in the
%         mesh) containing system matrix S. S is the sum of three
%         components: an absorption-dependent real matrix of products of
%         shape functions, a diffusion-dependent real matrix of products of
%         shape function derivatives, and a frequency-dependent imaginary
%         matrix of products of shape functions.
%
%         If the 'EL' flag is present, the system matrix is calculated on
%         an element basis, rather than a node basis. Length 'n' of all
%         parameter vectors in that case must be equal to the number of
%         elements. Parameters are considered constant over each element,
%         so parameter distribution has discontinuities at element
%         boundaries.
%
%         S is required for the evaluation of the frequency-domain FEM
%         diffusion equation.
%
%         If B is provided as an output argument, it is filled with the
%         boundary component of the system matrix (but without the
%         refractive index mismatch prefactor). The sparsity structure of B
%         is the same as S.
%
%         If alpha is provided as an output argument, it is filled with the
%         boundary mismatch prefactors. (returned as a vector of size n).
%
% See also:
%         toastMesh.SysmatComponent, toastMesh.Massmat, toastMesh.Elmat
smat = toast(uint32(23),hmesh.handle,double(mua),double(mus),double(ref),freq);
end
