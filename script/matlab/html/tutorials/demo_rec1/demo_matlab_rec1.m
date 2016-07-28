function demo_matlab_rec1

%% Some general parameters
refind = 1.4;   % refractive index of medium;
c0 = 0.3;       % speed of light in vacuum
cm = c0/refind; % speed of light in the medium

%% Step 1: Compute target data

% read the target image and convert to absorption coefficients [1/mm]
bmua = imread('demo_matlab_fwd2_mua.png');
bmua = double(bmua)./255.*0.02 + 0.01;
figure;
subplot(1,2,1); imagesc(bmua,[0.005,0.03]); axis equal tight; colorbar;
title('absorption target');
drawnow

% construct the circular mesh
rad = 25;
[vtx,idx,eltp] = mkcircle(rad,6,64,2);
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
[vtx,idx,eltp] = mkcircle(rad,6,32,2);
mesh = toastMesh(vtx,idx,eltp);
n = mesh.NodeCount;

grd = [64 64];
basis = toastBasis(mesh,grd);

% set up initial parameter estimates
mua = ones(n,1)*0.01;
mus = ones(n,1)*1;
ref = ones(n,1)*refind;

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
x = basis.Map('M->S',mua.*cm);
logx = log(x);

% initial objective function
err0 = full(objective(data,proj));
err = err0;
errp = inf;
fprintf ('initial cost function: %f\n', err);

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
    
    % Line search
    step = toastLineSearch (logx, d, step, err, @ls_objective);
    
    % Add update to solution
    logx = logx + d*step;
    x = exp(logx);
    
    % Map parameters back to mesh
    mua = basis.Map('S->M', x)./cm;
    
    % And display new update
    mua_img = reshape(basis.Map('S->B', x)./cm, grd);
    subplot(1,2,2);
    imagesc(mua_img,[0.005,0.03]); axis equal tight; colorbar
    title ('absorption recon');
    drawnow
    
    % Update measurements and objective function
    K = dotSysmat(mesh,mua,mus,ref,0);
    proj = reshape(log(mvec.' * (K\qvec)), [], 1);
    err = full(objective(data,proj));
    fprintf ('cost function: %f\n', err);
    
    itr = itr+1;
end

    function of = objective(data,proj)
        of = sum((data-proj).^2);
    end

    function of = ls_objective(logx)
        x_ = exp(logx);
        mua_ = basis.Map('S->M', x_) ./ cm;
        K_ = dotSysmat(mesh,mua_,mus,ref,0);
        proj_ = reshape(log(mvec.' * (K_\qvec)), [], 1);
        of = objective(data,proj_);
    end

end