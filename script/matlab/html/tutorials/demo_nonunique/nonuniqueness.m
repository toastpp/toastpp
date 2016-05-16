% This example demonstrates a non-uniqueness condition in DOT:
% Transillumination amplitude data from a steady-state measurement at a
% single wavelength are not sufficient for reconstructing both absorption
% and scattering distributions.
%
% This is demonstrated by generating data from a model with homogeneous
% absorption and scattering containing an inclusion.
% Using these data, an absorption reconstruction is performed under the
% assumption that scattering is constant.
% The reconstruction succeeds in finding an absorption distribution that
% generates data matching the mus-perturbed target data.

function nonuniqueness

close all

% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
mus_bkg = 1;    % background scattering [1/mm];
rad = 25;       % mesh radius [mm]
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion

% load mesh
mesh = toastMesh('circle_blob.msh','gmsh');
ne = mesh.ElementCount;
nv = mesh.NodeCount;
regidx = mesh.Region;
regno = unique(regidx);
blobel = find(regidx == regno(2));

% define source and detector locations
nopt = 32;
for i=1:nopt
    phiq = (i-1)/32*2*pi;
    qpos(i,:) = rad*[cos(phiq), sin(phiq)];
    phim = (i-0.5)/32*2*pi;
    mpos(i,:) = rad*[cos(phim), sin(phim)];
end
mesh.SetQM(qpos,mpos);
qvec = real(mesh.Qvec('Neumann','Gaussian',2));
mvec = real(mesh.Mvec('Gaussian',2,refind));

% assign elementwise optical coefficients - mus perturbation
ref = ones(ne,1)*refind;
mua = ones(ne,1)*mua_bkg;
mus = ones(ne,1)*mus_bkg;
mus(blobel) = mus_bkg*2;

figure('position',[0,0,640,420]);
subplot(2,3,1); mesh.Display(mua, 'range',[0.005,0.025]); axis off; title('\mu_a target');
subplot(2,3,2); mesh.Display(mus, 'range',[0.8,2.2]);     axis off; title('\mu_s target');

% solve FEM linear system
smat = dotSysmat(mesh, mua, mus, ref, 'EL');
data = log(mvec' * (smat\qvec));

% for reference, also solve the homogeneous problem
mus = ones(ne,1)*mus_bkg;
smat = dotSysmat(mesh, mua, mus, ref, 'EL');
data_homog = log(mvec' * (smat\qvec));

subplot(2,3,3);
imagesc(data-data_homog, [-0.26,0.015]); axis equal tight; colorbar
title('target data');

% now reconstruct for mua, assuming (wrongly) that mus=mus_bkg
mua = ones(nv,1)*mua_bkg;
mus = ones(nv,1)*mus_bkg;
ref = ones(nv,1)*refind;
smat = dotSysmat(mesh, mua, mus, ref);
proj = reshape(log(mvec' * (smat\qvec)), [], 1);
sd = ones(size(proj));
data = reshape(data, [], 1);
subplot(2,3,4); mesh.Display(mua, 'range',[0.005,0.025]); axis off; title('\mu_a recon');
subplot(2,3,5); mesh.Display(mus, 'range',[0.8,2.2]); axis off; title('\mu_s recon');
subplot(2,3,6); imagesc(reshape(proj,nopt,nopt)-data_homog, [-0.26,0.015]); axis equal tight; colorbar
title('recon data');
drawnow

grd = [32,32];
basis = toastBasis(mesh, grd);

bmua = basis.Map('M->B', mua);
bcmua = bmua*cm;
scmua = basis.Map('B->S', bcmua);
x = scmua;
logx = log(x);
slen = length(x);

regul = toastRegul('TV', basis, logx, 1e-4, 'Beta', 0.01);

err0 = toastObjective(proj, data, sd, regul, logx);
err = err0;
errp = inf;
itr = 1;
step = 1.0;
fprintf('Iteration %d, objective %f\n', 0, err);

while (itr <= itrmax) && (err > tolCG*err0) && (errp-err > tolCG)

    errp = err;
    r = -toastGradient(mesh, basis, qvec, mvec, mua, mus, ref, 0, ...
        data, sd, 'method', 'cg', 'tolerance', 1e-12);
    r = r(1:slen); % drop mus gradient
    r = r .* x;    % parameter scaling
    r = r - regul.Gradient(logx);
    
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
        beta = (delta_new - delta_mid) / delta_old;
        if mod(itr, 20) == 0 || beta <= 0
            d = s;
        else
            d = s + d*beta;
        end
    end
    step = toastLineSearch(logx, d, step, err, @objective);
    logx = logx + d*step;
    mua = basis.Map('S->M',exp(logx)/cm);
    subplot(2,3,4); mesh.Display(mua, 'range',[0.005,0.025]); axis off; title('\mu_a recon');
    
    proj = reshape(log(mvec' * (dotSysmat(mesh, mua, mus, ref)\qvec)), [], 1);
    subplot(2,3,6); imagesc(reshape(proj,nopt,nopt)-data_homog, [-0.26,0.015]); axis equal tight; colorbar
    title('recon data');
    err = toastObjective(proj, data, sd, regul, logx);
    fprintf('Iteration %d, objective %f\n', itr, err);
    itr = itr+1;
    drawnow
end

    function p = objective(logx)
        
        mua_ = basis.Map('S->M',exp(logx))/cm;
        proj_ = reshape(log(mvec' * (dotSysmat(mesh, mua_, mus, ref)\qvec)), [], 1);
        p = toastObjective(proj_, data, sd, regul, logx);
        
    end
end


