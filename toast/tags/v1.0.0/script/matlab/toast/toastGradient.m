function p = toastGradient
%toastGradient        - Return gradient of objective function.
%
% Synopsis: G = toastGradient(hMesh,hBasis,qvec,mvec,mua,mus,n,f,data,sd,solver,tol)
%    hMesh:  mesh handle
%    hBasis: basis handle
%    qvec:   array of source column-vectors (complex sparse)
%    mvec:   array of measurement column-vectors (complex sparse)
%    mua:    nodal absorption values [1/mm] (real)
%    mus:    nodal scattering values [1/mm] (real)
%    n:      nodal refractive index values (real)
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
