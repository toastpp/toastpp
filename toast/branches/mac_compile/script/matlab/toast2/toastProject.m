function proj = toastProject(hMesh,mua,mus,ref,omega,qvec,mvec,solver,tol)
%toastProject         - Forward operator to calulate boundary data.
%
% Synopsis: proj = toastProject (hMesh,mua,mus,ref,omega,qvec,mvec, ...
%                                solver,tol)
%    hMesh: mesh handle
%    mua:   vector of nodal absorption values [1/mm]
%    mus:   vector of nodal scatter values [1/mm]
%    ref:   vector of nodal refractive index values
%    omega: modulation frequency [MHz]
%    qvec:  column matrix of nodal source distributions
%    mvec:  column matrix of nodal measurement distributions
%    solver: linear solver [DIRECT|CG|BICG|BICGSTAB|GMRES]
%    tol:   linear solver tolerance (ignored for solver DIRECT)
%    proj:  vector of boundary measurements.
%
% Performs forward solutions for all sources, and generates boundary data
% at all detector sites.
% The resulting projection vector is one-dimensional. It consists of two
% parts: the log of the modulation amplitude data, and the phase shift data.
% Each block contains sub-blocks for each source, and each sub-block
% consists of the data for each measurement site.
%
% The returned data vector is real of length 2*nqm, where nqm is the
% length of the permutation vector returned by toastDataLinkList, that
% is, unused source-detector combinations are removed from the result.

% Subsititute missing parameters with defaults
if nargin < 8
    solver = 'DIRECT';
elseif nargin < 9
    tol = 1e-10;
end

% Calculate system matrix
smat = dotSysmat (hMesh, mua, mus, ref, omega);

% Solve linear system for log complex field and apply boundary operator
switch upper(solver)
    case 'DIRECT'
        % solve with backslash operator
        lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);
    otherwise
        % for now, use BiCGSTAB for any iterative solver request
        nq = size(qvec,2);
        nm = size(mvec,2);

        % incomplete LU factorisation of smat
        %ilu_setup.type = 'ilutp';
        %ilu_setup.droptol = 1e-2;
        ilu_setup.type = 'nofill';
        [L U] = ilu(smat,ilu_setup);
        
        for i=1:nq
            [phi,flag] = bicgstab(smat,qvec(:,i),tol,1000,L,U);
            prj = log(mvec.' * phi);
            lgamma((i-1)*nm+1:i*nm,1) = prj;
        end
end

% Strip off unused source-detector combinations
lgamma = lgamma(hMesh.DataLinkList());

% Rearrange data in terms of log amplitude and phase shift blocks
if omega > 0
    proj = [real(lgamma);imag(lgamma)];
else
    proj = lgamma;
end
