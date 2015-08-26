function J = toastJacobianCW(mesh,basis,varargin)
% Jacobian for DOT continuous wave problem for mua parameter
%
% Syntax: J = toastJacobianCW (mesh, basis, qvec, mvec, mua, mus, ref,
%                              solver, tol)
%
% Parameters:
%         mesh [toastMesh object]:
%             mesh object
%         basis [toastBasis object]:
%             basis mapper object (or 0 to return Jacobian in mesh basis)
%         qvec [nlen x nq sparse matrix]:
%             array of source column vectors
%         mvec [nlen x nm sparse matrix]:
%             array of measurement operator column vectors
%         mua [nlen vector]:
%             nodal absorption coefficients
%         mus [nlen vector]:
%             nodal scattering coefficients
%         ref [nlen vector]:
%             nodal refractive index coefficients
%         solver [string]:
%             linear solver (DIRECT|CG|BICG|BICGSTAB|GMRES|GAUSSSEIDEL)
%         tol [scalar]:
%             linear solver tolerance (optional, iterative solvers only)
%
% Return values:
%         J: [nqm x slen dense real matrix]:
%             Jacobian matrix
%
% Notes:  Calculates the derivative of the logarithm of the CW amplitude
%         data with respect to the absorption coefficients of the forward
%         operator.
%
%         The returned matrix is of size nqm x slen, where nqm is the number
%         of measurements, and slen is the dimension of the reconstruction
%         basis.

if isobject(basis)
    hb = basis.handle;
else
    hb = 0;
end

if nargin==5
    % This format isn't supported yet
    error('Incorrect argument syntax');
    dphi = varargin{1};
    aphi = varargin{2};
    proj = varargin{3};
    J = toastmex(uint32(54),mesh.handle,hb,dphi,aphi,proj);
else
    qvec = varargin{1};
    mvec = varargin{2};
    mua = varargin{3};
    mus = varargin{4};
    ref = varargin{5};
    solver = 'direct';
    tol = 1e-8;
    if nargin >= 8
        solver = varargin{6};
        if nargin >= 9
            tol = varargin{7};
        end
    end
    J = toastmex(uint32(54),mesh.handle,hb,qvec,mvec,mua,mus,ref,solver,tol);
end

end
