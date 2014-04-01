function p = toastFields
%toastFields          - Calculate direct and adjoint fields
%
% Synopsis: [dphi aphi] = toastFields (hMesh, hBasis, qvec, mvec, mua,
%                                      mus, ref, freq, solver, tol)
%    hMesh:  mesh handle (see toastReadMesh)
%    hBasis: basis mapper handle (see toastSetBasis)
%    qvec:   Sparse matrix of source vectors (complex column vectors)
%    mvec:   Sparse matrix of measurement vectors (complex column vectors)
%    mua:    nodal absorption coefficient [1/mm] (real column vector)
%    mus:    nodal reduced scattering coefficient [1/mm] (real column vector)
%    ref:    nodal refractive index (real column vector)
%    freq:   modulation frequency [MHz] (real)
%    solver: linear solver (string):
%                 (DIRECT|CG|BICG|BICGSTAB|GMRES|GAUSSSEIDEL)
%    tol:    linear solver tolerance (optional, iterative solvers only)
%
%    dphi:   direct fields (dense complex matrix glen x nQ)
%    aphi:   adjoint fields, optional (dense complex matrix glen x nM)
%
% Calculates the direct and adjoint fields for given optical
% coefficients, and returns them as complex column vectors in dphi and aphi,
% mapped into grid basis.
% Adjoint field matrix (aphi) is optional.
%
% See also: toastJacobian
