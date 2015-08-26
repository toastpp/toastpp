% toastFields  - Calculate complex photon density fields
%
% Syntax: phi = toastFields(mesh,basis,qvec,mua,mus,ref,freq,method,tol)
%
% Parameters:
%         mesh (toastMesh instance):
%             mesh object
%         basis (toastBasis instance):
%             basis object (set to 0 to return fields in mesh basis)
%         qvec (complex sparse matrix n x nq):
%             matrix of nq source vectors
%         mua (real array n):
%             nodal absorption coefficients [1/mm]
%         mus (real array n):
%             nodal scattering coefficients [1/mm]
%         ref (real array n):
%             nodal refractive index values
%         freq (scalar):
%             modulation frequency [MHz]
%         method (string):
%             solver method (DIRECT|CG|BICGSTAB|GMRES)
%         tol (scalar):
%             solver tolerance
%
% Return values:
%         phi (complex matrix slen x nq or n x nq):
%             photon density fields for all nq sources
%
% Notes:  If the basis parameter is set to 0, the fields are returned in
%         the mesh basis. Otherwise, they are returned in the solution
%         basis. The returned matrix contains the field for source i in
%         column i.
%
%         To compute adjoint fields, pass the matrix of measurement
%         vectors (mvec) instead of qvec.

function [phi,aphi] = toastFields(mesh,basis,qvec,mua,mus,ref,freq,method,tol)

if nargin < 9
    tol = 1e-10;
    if nargin < 8
        method = 'direct';
    end
end

if basis == 0
    bhandle = 0;
else
    bhandle = basis.handle;
end

phi = toastmex(uint32(51),mesh.handle,bhandle,qvec,mua,mus,ref, ...
            freq,method,tol);
