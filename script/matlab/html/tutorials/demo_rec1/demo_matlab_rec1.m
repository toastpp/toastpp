function demo_matlab_rec1

%% Step 1: Compute target data

% read the target image and convert to absorption coefficients [1/mm]
bmua = imread('demo_matlab_fwd2_mua.png');
bmua = double(bmua)./255.*0.02 + 0.01;

% construct the circular mesh
rad = 25;
[vtx,idx,eltp] = mkcircle(rad,6,32,2);
fwdmesh = toastMesh(vtx,idx,eltp);

% map parameter distribution to nodal coefficients
fwdgrd = size(bmua);
fwdbasis = toastBasis(fwdmesh,fwdgrd);
mua = fwdbasis.Map('B->M',bmua);
mus = ones(size(mua))*1;
ref = ones(size(mua))*1.4; 

% construct source and measurement vectors
nq = 16;
for i=1:nq
  phi_q = 2*pi*(i-1)/nq;
  Q(i,:) = rad * [cos(phi_q) sin(phi_q)];
  phi_m = 2*pi*(i-0.5)/nq;
  M(i,:) = rad * [cos(phi_m) sin(phi_m)];
end
fwdmesh.SetQM(Q,M);
qvec = real(fwdmesh.Qvec('Neumann','Gaussian',2));
mvec = real(fwdmesh.Mvec('Gaussian',2,ref)); 

% solve the discretised linear system and project to boundary data
K = dotSysmat(fwdmesh,mua,mus,ref,0);
Phi = K\qvec;
Y = mvec.' * Phi;

% map to log
data = reshape(log(Y), [], 1);

%% Step 2: inverse solver

% construct a coarser mesh for the inverse problem
[vtx,idx,eltp] = mkcircle(rad,6,16,2);
mesh = toastMesh(vtx,idx,eltp);
n = mesh.NodeCount;

grd = [32,32];
basis = toastBasis(mesh,grd);

% set up initial parameter estimates
mua = ones(n,1)*0.01;
mus = ones(n,1)*1;
ref = ones(n,1)*1.4;

% set up source and measurement vectors for inverse mesh
mesh.SetQM(Q,M);
qvec = mesh.Qvec('Neumann','Gaussian',2);
mvec = mesh.Mvec('Gaussian',2,ref);

% solve forward problem for initial data estimate
K = dotSysmat(mesh,mua,mus,ref,0);
proj = reshape(log(mvec.' * (K\qvec)), [], 1);

% standard deviation for data scaling: set to unity
sd = ones(size(proj));

% set up solution vector
c = 0.3./ref; % speed of light in the medium [mm/ps]
x = basis.Map('M->S',mua.*c);
logx = log(x);

% initial objective function
err0 = objective(data,proj);
err = err0;
errp = inf;

% nonlinear conjugate gradient solver loop
tolcg = 1e-6;
resetcg = 20;
step = 1;
itrmax = 100;
itr = 1;

while (itr <= itrmax) && (err > tolcg*err0) && (errp-err > tolcg)

    errp = err;
    
    r = -toastGradient(mesh,basis,qvec,mvec,mua,mus,ref,0,data,sd,'method','direct');
    r = r(1:basis.slen);
    r = r .* x;
    
    if itr > 1
        delta_old = delta_new;
        delta_mid = r' * s;
    end
    s = r;
    if itr == 1
        d = s;
        delta_new = r' * d;
    else
        delta_new = r' * s;
        beta = (delta_new-delta_mid) / delta_old;
        if mod(itr,resetcg) == 0 || beta <= 0
            d = s; % reset NCG
        else
            d = s + d*beta;
        end
    end
    
    step = toastLineSearch (logx, d, step, err, @ls_objective);
    
end

    function of = objective(data,proj)
        of = sum((data-proj).^2);
    end

    function of = ls_objective(logx)
        x_ = exp(logx);
        mua_ = basis.Map('S->M', x_) ./ c;
        K_ = dotSysmat(mesh,mua_,mus,ref,0);
        proj_ = reshape(log(mvec.' * (K_\qvec)), [], 1);
        of = objective(data,proj_);
    end

end