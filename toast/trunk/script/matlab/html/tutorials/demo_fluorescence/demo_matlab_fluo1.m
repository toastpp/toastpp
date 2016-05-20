%% MATLAB-TOAST sample script:
% Reconstruction of fluorochrome concentration from fluorescence
% measurements

%% Optical background parameters
refind = 1.4;     % homogeneous refractive index of the medium
mua_homog = 0.05; % background absorption [1/mm]
mus_homog = 1;    % background scattering [1/mm]

%% Create the mesh
rad = 25;
[vtx,idx,eltp] = mkcircle(rad,6,32,2);
mesh = toastMesh(vtx,idx,eltp);
n = mesh.NodeCount;

%% Source and measurement projection vectors
np = 32;
for i=1:np
    qphi = (i-1)/np * 2*pi;
    qpos(i,:) = rad * [cos(qphi) sin(qphi)];
    mphi = (i-0.5)/np *2*pi;
    mpos(i,:) = rad * [cos(mphi) sin(mphi)];
end
mesh.SetQM(qpos,mpos);
dmask = mesh.DataLinkList;

qvec = mesh.Qvec ('Neumann', 'Gaussian', 0.5);
nQ = size(qvec,2);

mvec = mesh.Mvec ('Gaussian', 0.5, refind);
nM = size(mvec,2);

%% solution basis
grd = [64 64];
blen = prod(grd);
basis = toastBasis (mesh, grd);
nsol = basis.slen;

%% The flu target in grid dimensions: 3 different blobs 
fltgt_g = zeros(grd);
nblobs = 3;
fcontrast = [0.1 0.06 0.1];
blbcx = grd(1)*[0.7 0.75 0.3]; % xcentre of blobs
blbcy = grd(2)*[0.3 0.75 0.7]; % ycentre of blobs
blbrad = grd(1)*[0.1 0.1 0.1]; % radius of blobs

for i = 1:grd(1)
  for j = 1:grd(2)
      for k = 1:nblobs
        if( (i-blbcx(k))^2 + (j-blbcy(k))^2 < blbrad(k)^2)
            fltgt_g(i,j) = fltgt_g(i,j)+fcontrast(k);
        end
      end
  end
end
muatgt_g = mua_homog*ones(grd);
mustgt_g = mus_homog*ones(grd);

% map to solution basis
fltgt_h   = basis.Map ('B->M', reshape(fltgt_g,[],1));
muatgt_h  = basis.Map ('B->M', reshape(muatgt_g,[],1));
mustgt_h  = basis.Map ('B->M', reshape(mustgt_g,[],1));
 
figure;
subplot(2,2,1);imagesc(fltgt_g); colormap(gray);colorbar;axis square;title('target fluoresence');
subplot(2,2,2);imagesc(muatgt_g);colormap(gray);colorbar;axis square;title('Background \mu_a');
subplot(2,2,3);imagesc(mustgt_g);colormap(gray);colorbar;axis square;title('background \mu_s');

%% modulation frequency [MHz]
freq = 0;

% refractive index
ref = ones(n,1) * refind;

%% create the simulated data

% FEM system matrix
smat = dotSysmat (mesh, muatgt_h, mustgt_h, ref, freq);

% excitation photon density
sol1=smat\qvec;
sol1=full(sol1); % mappings dont work with sparse matrices

% fluoresence photon density
sol_fl=smat\(sol1.*(fltgt_h*ones(1,nQ)));

plot_sol1eta=1;
if plot_sol1eta==1
    dr=(sol1.*(fltgt_h*ones(1,nQ)));
    figure;
    for jj=1:nQ
        subplot(6,6,jj);imagesc(reshape(basis.Map ('M->B',dr(:,jj)),grd)),colormap(gray);axis square;
    end
end

% excitation boundary data
lgamma = reshape ((mvec.' * sol1), (size(mvec,2))*size(sol1,2), 1);
lgamma = lgamma(dmask);
xdata = [real(lgamma)];

% fluorescence boundary data
lgamma_fl = reshape ((mvec.' * sol_fl), (size(mvec,2))*size(sol1,2), 1);
lgamma_fl = lgamma_fl(dmask);
fdata = [real(lgamma_fl)];

R = sparse(diag(1./(xdata))); % scaling matrix: inverse of excitation data
figure;
subplot(2,2,1);imagesc(reshape(xdata,nM,nQ));colormap(gray);colorbar;axis square;title('Excitation Data RtN map');
subplot(2,2,2);imagesc(reshape(fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Fluoresence Data RtN map');
subplot(2,2,3);imagesc(reshape(R*fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Ratio of Fluoresence/Excitation');


%% inverse solution: reconstruct chromophore concentration

r = R*fdata;  % rescale fluorescence data by excitation.
asol1 = smat\mvec; % adjoint

%% Jacobian
Jr = zeros(nM*nQ,nsol);
for i = 1:nQ
    for j = 1:nM
        tmp_h = full(sol1(:,i).*asol1(:,j));
        tmp_s = basis.Map ('M->S',tmp_h);
        Jr(( i-1)*nM +j,:) = tmp_s;
    end
end
Jr = R*Jr; % rescale Jacobian

lambda = 1e-4*(trace(Jr' * Jr)); % regularisation proportional to trace
disp('solving using Conjugate Gradients');
s_tol=1e-6;
s_iter=200;
[x,CGflag,relres,niter] = cgs(Jr'*Jr + lambda*speye(nsol), Jr'*r,s_tol,s_iter);
disp(['CG iterations: ',num2str(niter)]);

xim = basis.Map('S->B', x);
 
figure;
subplot(1,2,1);imagesc(fltgt_g);colorbar;axis square;title('target fluoresence');
subplot(1,2,2);imagesc(reshape(xim,grd));colormap(gray);colorbar;axis square;title('reconstructed fluoresence');
