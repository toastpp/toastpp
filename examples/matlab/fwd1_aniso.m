% Simple 2D example for anisotropic diffusion using a diffusion tensor

function fwd1_aniso
close all

meshdir = '../../test/2D/meshes/';
hmesh = toastMesh([meshdir 'circle25_32.msh']);
hmesh.ReadQM([meshdir 'circle25_32x32.qm']);

refind = 1.4;
c0 = 0.3;
cm = c0/refind;

n = hmesh.NodeCount;
ne = hmesh.ElementCount;

mua_homog = 0.01;
mus_homog = 1;
kap_homog = 1/(3*(mua_homog+mus_homog));

freq = 100;

qvec = hmesh.Qvec('Neumann','Gaussian',2);
mvec = hmesh.Mvec('Gaussian',2,refind);

solve_aniso(1, 0, 1, 'iso');
solve_aniso(2, 0, 4, 'aniso 0 lo');
solve_aniso(2, pi/4, 5, 'aniso 45 lo');
solve_aniso(2, pi/2, 6, 'aniso 90 lo');
solve_aniso(4, 0, 7, 'aniso 0 hi');
solve_aniso(4, pi/4, 8, 'aniso 45 hi');
solve_aniso(4, pi/2, 9, 'aniso 90 hi');

figure(4);
subplot(1,2,1);
legend('iso', 'aniso 0 lo', 'aniso 45 lo', 'aniso 90 lo', 'aniso 0 hi', 'aniso 45 hi', 'aniso 90 hi');
title('log amplitude');

subplot(1,2,2);
legend('iso', 'aniso 0 lo', 'aniso 45 lo', 'aniso 90 lo', 'aniso 0 hi', 'aniso 45 hi', 'aniso 90 hi');
title('phase');

    function solve_aniso(strength, dir, idx, lbl)
    phi_ = fwd_aniso(hmesh, mua_homog, kap_homog*strength, kap_homog/strength, dir, refind, freq);
    meas_ = log(mvec.' * phi_);
    show_photon_density(phi_, idx, lbl);
    show_sinogram(meas_, idx, lbl);
    show_meas(meas_);
    end

    function show_photon_density(phi, i, t)
    figure(1);
    subplot(3,3,i);
    hmesh.Display(real(log(phi(:,1))),'range',[-12,1]);
    title(t);
    axis off
    end

    function show_sinogram(meas, i, t)
    figure(2);
    subplot(3,3,i);
    imagesc(real(meas),[-14,-2]);
    axis equal tight
    colorbar
    xlabel('source #');
    ylabel('detector #');
    title(t);
    figure(3);
    subplot(3,3,i);
    imagesc(imag(meas),[-1.2,0]);
    axis equal tight
    colorbar
    xlabel('source #');
    ylabel('detector #');
    title(t);
    end

    function show_meas(meas)
    figure(4);
    subplot(1,2,1);
    hold on
    plot(real(meas(:,1)));
    subplot(1,2,2);
    hold on
    plot(imag(meas(:,1)));
    end
end

function phi = fwd_aniso(mesh, mua, kap_hi, kap_lo, theta, refind, freq)
cm_ = 0.3/refind;
omega_ = freq*1e6*1e-12*2*pi;
n_ = mesh.NodeCount;
ne_ = mesh.ElementCount;
R_ = [cos(theta), -sin(theta); sin(theta), cos(theta)];
D_ = [kap_hi, 0; 0, kap_lo];
kapT_ = R_*D_*R_';
ref_ = ones(n_,1) * refind;
cmua_ = ones(n_,1) * (mua*cm_);
zeta_ = ones(n_,1) * (cm_./(2*toastDotBndterm(refind,'Keijzer')));
ckap_xx_ = ones(n_,1) * (kapT_(1,1)*cm_);
ckap_xy_ = ones(n_,1) * (kapT_(1,2)*cm_);
ckap_yx_ = ones(n_,1) * (kapT_(2,1)*cm_);
ckap_yy_ = ones(n_,1) * (kapT_(2,2)*cm_);

qvec_ = mesh.Qvec('Neumann','Gaussian',2);

smat_ = sparse(n_, n_);
for i_=1:ne_
    el_ = mesh.Element(i_);
    int_kapdd_xx = el_.Mat('Pdd', ckap_xx_);  int_kapdd_xx = squeeze(int_kapdd_xx(:,1,:,1));
    int_kapdd_xy = el_.Mat('Pdd', ckap_xy_);  int_kapdd_xy = squeeze(int_kapdd_xy(:,1,:,2));
    int_kapdd_yx = el_.Mat('Pdd', ckap_yx_);  int_kapdd_yx = squeeze(int_kapdd_yx(:,2,:,1));
    int_kapdd_yy = el_.Mat('Pdd', ckap_yy_);  int_kapdd_yy = squeeze(int_kapdd_yy(:,2,:,2));
    int_kapDD = int_kapdd_xx + int_kapdd_xy + int_kapdd_yx + int_kapdd_yy;
    int_muaFF = el_.Mat('PFF', cmua_);
    int_omegaFF = omega_.*el_.Mat('FF');
    bint_zetaFF = el_.Mat('BndPFF', zeta_);
    elidx_ = el_.Dof;
    smat_(elidx_, elidx_) = smat_(elidx_, elidx_) + ...
        int_kapDD + int_muaFF + bint_zetaFF + 1i*int_omegaFF;
end
phi = smat_\qvec_;
end


