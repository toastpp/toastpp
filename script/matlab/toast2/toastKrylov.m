function p = toastKrylov(x,J,g,M,lambda,hreg,tol)
%toastKrylov          - Solve linear system Hx=g with impicit representation of H.
%
% Synopsis: z = toastKrylov (x, J, g, M, lambda, hReg, tol)
%
%    x:      current solution (real vector, dimension n)
%    J:      Jacobian (dense real matrix, dimension m x n)
%    g:      gradient (dense real vector, dimension n)
%    M:      diagonal of Hessian scaling matrix (dimension n)
%    lambda: LM control parameter (real)
%    hReg:   regularisation handle
%    tol:    tolerance (real): Krylov solver convergence criterion
%
% This solves a linear system with a GMRES solver, using an implicit
% representation of the Hessian matrix H:
%
%    H = J(x)^T J(x) + M psi''(x) M + lambda I
%
% where J is the Jacobian, M is a diagonal scaling matrix, psi is the
% regularisation term, and lambda is the LM control parameter.
if nargin < 7
    tol = 1e-8;
    if nargin < 6
        hreg = 0;
	if nargin < 5
	    lambda = 0;
	end
    end
end

if isobject(hreg)
    hr = hreg.handle;
else
    hr = 0;
end

p = toastmex(uint32(55),x,J,g,M,lambda,hr,tol);
