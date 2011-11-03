function jac = toastJacobianMoments(hmesh,hbasis,nmom,qvec,mvec,mua,mus,ref,solver,tol)
%toastJacobianMoments - Return Jacobian for moments of temporal distribution
%
% Synopsis: jac = toastJacobianMoments(hmesh,hbasis,qvec,mvec,mua,mus,ref,solver,tol,maxmom)
%
%    hmesh:  mesh handle
%    hbasis: basis mapper handle
%    nmom:   highest moment to compute (integer, >= 0)
%    qvec:   sparse matrix of source vectors
%    mvec:   sparse matrix of measurement vectors
%    mua:    nodal absorption coefficient vector [1/mm]
%    mus:    nodal scattering coefficient vector [1/mm]
%    ref:    nodal refractive index vector
%    solver: linear solver (string)
%    tol:    iterative solver tolerance (scalar)

n = toastMeshNodeCount(hmesh);
nq = size(qvec,2);
nm = size(mvec,2);
nqm = nq*nm; % for now
[bdim gdim] = toastGetBasisSize(hbasis);
blen = prod(bdim);
glen = prod(gdim);
solmask = toastSolutionMask(hbasis);
slen = length(solmask);

smat = real (toastSysmat (hmesh, mua, mus, ref, 0));
mmat = toastMassmat (hmesh);

% 0th moment - direct fields and measurement vector
h_dphi0 = smat\qvec;
g_dphi0 = zeros(glen,nq);
for q=1:nq
    g_dphi0(:,q) = toastMapMeshToGrid(hbasis,h_dphi0(:,q));
end
gamma0 = full(mvec.' * h_dphi0);

% 0th moment - adjoint fields
h_aphi0 = smat\mvec;
g_aphi0 = zeros(glen,nm);
for m=1:nm
    g_aphi0(:,m) = toastMapMeshToGrid(hbasis,h_aphi0(:,m));
end

% buffer for higher moment fields
g_dphi = zeros(nmom,glen,nq);
g_aphi = zeros(nmom,glen,nm);

% buffer for Jacobian
jac = zeros(nqm*(nmom+1),slen);
idx = 1;

% add zero-th moment
for q=1:nq
    dp = g_dphi0(:,q)';
    for m=1:nm
        ap = g_aphi0(:,m)';
        jac(idx,:) = toastMapGridToSol(hbasis,dp .* ap)';
        idx = idx+1;
    end
end

% add higher moments
for mom = 1:nmom
    h_dphi0 = mom*(smat\(mmat*h_dphi0));
    h_aphi0 = mom*(smat\(mmat*h_aphi0));
    gamma = full((mvec.' * h_dphi0) ./ gamma0); % moment boundary measurements
    
    for q=1:nq
        g_dphi(mom,:,q) = toastMapMeshToGrid(hbasis,h_dphi0(:,q));
    end
    for m=1:nm
        g_aphi(mom,:,m) = toastMapMeshToGrid(hbasis,h_aphi0(:,m));
    end
    
    idx0 = 1;
    for q=1:nq
        for m=1:nm
            jac_mm_qm = zeros(1,glen);
            for r=0:mom
                fac = factorial(mom)/(factorial(r)*factorial(mom-r));
                if r==0
                    dp = g_dphi0(:,q)';
                else
                    dp = g_dphi(r,:,q);
                end
                if mom-r==0
                    ap = g_aphi0(:,m)';
                else
                    ap = g_aphi(mom-r,:,m);
                end
                jac_mm_qm = jac_mm_qm + fac * (dp .* ap);
            end
            % convert Mellin transform to moments
            jac(idx,:) = (toastMapGridToSol(hbasis,jac_mm_qm)' - gamma(m,q).*jac(idx0,:)) ./ gamma0(m,q);

            idx = idx+1;
            idx0 = idx0+1;
        end
    end
end