%% MATLAB-TOAST sample script:

% note that the script contains parameters depending on grid size
% at various places
clear all
close all

% refractive index and speed of light
refind = 1.4;
c0 = 0.3;
cm = c0/refind;

%% ======================================================================
% Read mesh and qm file

% mesh file
PLOT_ON=0;

meshname = [getenv('TOASTDIR') '/test/2D/meshes/circle25_32.msh'];
qmname = [getenv('TOASTDIR') '/test/2D/meshes/circle25_32x32.qm'];

mesh = toastMesh(meshname);
mesh.ReadQM (qmname);
n = mesh.NodeCount;
dmask = mesh.DataLinkList;


%% solution basis
bx = 64; by = 64;
blen = bx*by;

% Set up the mapper between FEM and solution bases
basis = toastBasis (mesh, [bx by]);
nsol = basis.slen;
%solmask = toastSolutionMask (hBasis);
%nsol = size(solmask,2);

%% Generate source vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 0.5);
nQ = size(qvec,2);
nh = size(qvec,1); % number of nodes
% Generate measurement vectors
mvec = mesh.Mvec ('Gaussian', 0.5, refind);
nM = size(mvec,2);

%+The flu target in grid dimensions 3 different blobs 

fltgt_g = zeros(bx,by);
nblobs = 3;
fcontrast = [0.1 0.06 0.1];
%fcontrast = [1,1,1];
blbcx = bx*[0.7 0.75 0.3];  % xcentre of blobs
blbcy = by*[0.3 0.75 0.7];  % ycentre of blobs
blbrad = bx*[0.1 0.1 0.1]; % radius of blobs

for i = 1:bx
  for j = 1:by
      for k = 1:nblobs
        if( (i-blbcx(k))^2 + (j-blbcy(k))^2 < blbrad(k)^2)
            fltgt_g(i,j) = fltgt_g(i,j)+fcontrast(k);
        end
      end
  end
end
muatgt_g=0.05*ones(bx,by);
mustgt_g=1*ones(bx,by);


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


%projection data f[mua,kap]
smat = dotSysmat (mesh, muatgt_h, mustgt_h, ref, freq);

% excitation 
sol1=smat\qvec;
sol1=full(sol1); % mappings dont work with sparse matrices
% fluoresence 
sol_fl=smat\(sol1.*(fltgt_h*ones(1,nQ)));


plot_sol1eta=1;

if plot_sol1eta==1
     dr=(sol1.*(fltgt_h*ones(1,nQ)));
       figure;
    for jj=1:nQ
        subplot(6,6,jj);imagesc(reshape(basis.Map ('M->B',dr(:,jj)),bx,by)),colormap(gray);axis square;
    end

end


% projecting with the measurements vector 

lgamma = reshape ((mvec.' * sol1), (size(mvec,2))*size(sol1,2), 1);
lgamma = lgamma(dmask);
xdata = [real(lgamma)];

lgamma_fl = reshape ((mvec.' * sol_fl), (size(mvec,2))*size(sol1,2), 1);
lgamma_fl = lgamma_fl(dmask);
fdata = [real(lgamma_fl)];
R = sparse(diag(1./(xdata))); % must use this if comparing to TOAST
toastWriteVector('./data/fl_sim_mod_mat.fem',fdata)
figure;
subplot(2,2,1);imagesc(reshape(xdata,nM,nQ));colormap(gray);colorbar;axis square;title('Excitation Data RtN map');
subplot(2,2,2);imagesc(reshape(fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Fluoresence Data RtN map');
subplot(2,2,3);imagesc(reshape(R*fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Ratio of Fluoresence/Excitation');


%% solution

r = R*fdata;  % rescale fluorescence data by excitation.
asol1 = smat\mvec; % adjoint
%% Jacobian
Jnew = zeros(nM*nQ,nsol); % should use 'size' to make general
for i = 1:nQ
  for j = 1:nM
    tmp_h = full(sol1(:,i).*asol1(:,j));
    tmp_s = basis.Map ('M->S',tmp_h);
    Jnew(( i-1)*nM +j,:) = tmp_s;
  end
end
Jr = R*Jnew;
clear Jnew;
% column scaling if required
c = ones(nsol,1);
for k = 1:nsol
    c(k) = norm(Jr(:,k));
end
S = sparse(diag(1./c));
S = 0.001*speye(nsol); % reset to identity (i.e. no image scaling)
Jr = Jr *S;
Algorithm=2;

lambda = 1e-4*(trace(Jr' * Jr)); % regularisation proportional to trace
if(Algorithm==1)
    disp('solving using Backslash');
    tr = [r; ones(nsol,1)];
    tJ = [Jr ; lambda*eye(nsol)];
    x = tJ\tr;
else
    disp('solving using Conjugate Gradients');
    s_tol=1e-6;
    s_iter=200;
    [x,CGflag,relres,niter] = cgs(Jr'*Jr + lambda*speye(nsol), Jr'*r,s_tol,s_iter);
    disp(['CG iterations: ',num2str(niter)]);
end
x = S*x; %unscale
xim = basis.Map('S->B', x);
sim = basis.Map('S->B', diag(S));
 
figure;
subplot(2,2,1);imagesc(fltgt_g);colorbar;axis square;title('target fluoresence');
subplot(2,2,2);imagesc(reshape(sim,bx,by));colormap(gray);colorbar;axis square;title('Image recscaling');
subplot(2,2,3);imagesc(reshape(xim,bx,by));colormap(gray);colorbar;axis square;title('reconstructed fluoresence');
