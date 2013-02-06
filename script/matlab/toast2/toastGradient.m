function p = toastGradient(hmesh,hbasis,qvec,mvec,mua,mus,ref,freq,data,sd,solver,tol)
%toastGradient        - Return gradient of objective function.
%
% Synopsis: G = toastGradient(hMesh,hBasis,qvec,mvec,mua,mus,n,f,data,sd,solver,tol)
%    hmesh:  mesh handle
%    hbasis: basis handle
%    qvec:   array of source column-vectors (complex sparse)
%    mvec:   array of measurement column-vectors (complex sparse)
%    mua:    nodal absorption values [1/mm] (real)
%    mus:    nodal scattering values [1/mm] (real)
%    ref:    nodal refractive index values (real)
%    f:      modulation frequency [MHz] (real scalar)
%    data:   measurement data vector (real)
%    sd:     measurement standard-deviation vector (real)
%    solver: DIRECT|CG|BICG|BICGSTAB|GMRES|GAUSSSEIDEL (string)
%    tol:    solver tolerance (optional, iterative solvers only)
%    G:      gradient of objective function (real)
%
% This function runs the forward solver with the nodal optical parameters
% provided, and calculates the derivative of the objective function with
% respect to the optical parameters.
% The returned gradient vector is with respect to the inverse solution basis
% defined by hBasis, rather than the nodal mesh basis.

if nargin < 12
    tol = 1e-10;
    if nargin < 11
        solver = 'DIRECT';
    end
end

p = toast(uint32(52),hmesh.handle,hbasis.handle,qvec,mvec,mua,mus,ref,freq,data,sd,solver,tol);