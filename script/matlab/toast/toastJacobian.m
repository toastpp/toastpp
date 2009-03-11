function p = toastJacobian
%toastJacobian        - Generate an unscaled frequecny-domain Jacobian matrix
%
% Synopsis: J = toastJacobian (hMesh, hBasis, qvec, mvec, mua, mus, ref,
%                              freq, solver, tol)
%           J = toastJacobian (hMesh, hBasis, dphi, aphi, proj)
%
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
%    dphi:   complex direct fields (nodal basis), size n x nq
%    aphi:   complex adjoint fields (nodal basis), size n x nm
%    proj:   complex boundary projection data, size nqm
%
%    J:      Jacobian (dense real matrix)
%
% Calculates the derivative of the data (log amplitude and phase) with respect
% to the coefficients (absorption and diffusion) of the forward operator.
%
% There are two versions: (a) providing optical coefficients and solver
% parameters, and (b) providing the fields and projections. The first
% version calculates fields and projections on the fly.
%
% J consists of 4 blocks: d lnmod / d mua (top left), d lnmod / d kappa (top
% right), d phase / d mua (bottom left), and d phase / d kappa (bottom right).
% Dimension of J is 2m x 2n, where m is the number of measurements, and n is
% the dimension of the reconstruction basis.
%
% If hBasis is set to 0, the Jacobian is constructed directly in the mesh
% basis. n in that case is equal to the number of nodes.
%
% J is unscaled on return.
