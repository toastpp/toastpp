function tmp
close all

meshdir = '../../test/3D/meshes/';
hmesh = toastMesh([meshdir 'cyl3.msh']);
hmesh.ReadQM([meshdir 'cyl_3ring.qm']);

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

grd = [64 64 64];
basis = toastBasis(hmesh,grd);

mag = [1, 4, 4];
az = [0, 0, pi/4];
el = [0, 0, pi/5];

for pass = 1:3
    
lambda = [kap_homog*mag(pass), kap_homog/mag(pass), kap_homog/mag(pass)];
phi = fwd_aniso(hmesh, mua_homog, lambda, az(pass), el(pass), refind, freq);

imphi = basis.Map('M->B', real(log(phi(:,1))));
imphi = reshape(imphi, grd);

elref = basis.GridElref;
imphi(find(elref==0)) = -50;

%figure; imagesc(imphi(:,:,32), [-12, 0]); axis equal tight; colorbar

figure;
p = patch(isosurface(imphi, -3));
isonormals(imphi, p);
p.FaceColor = 'red';
p.EdgeColor = 'none';

p = patch(isosurface(imphi, -5));
isonormals(imphi, p);
p.FaceColor = 'red';
p.FaceAlpha = 0.4;
p.EdgeColor = 'none';

p = patch(isosurface(imphi, -7));
isonormals(imphi, p);
p.FaceColor = 'red';
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';

p = patch(isosurface(imphi, -40));
isonormals(imphi, p);
p.FaceColor = 'green';
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';

daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud

end

function phi = fwd_aniso(mesh, mua, lambda, az, el, refind, freq)
R1 = [cos(el), 0, sin(el); 0, 1, 0; -sin(el), 0, cos(el)];
R2 = [cos(az), sin(az), 0; -sin(az), cos(az), 0; 0, 0, 1];
R = R2*R1;
kapT = R * diag(lambda) * R';
ckap_xx = ones(n,1) * (kapT(1,1)*cm);
ckap_xy = ones(n,1) * (kapT(1,2)*cm);
ckap_xz = ones(n,1) * (kapT(1,3)*cm);
ckap_yx = ones(n,1) * (kapT(2,1)*cm);
ckap_yy = ones(n,1) * (kapT(2,2)*cm);
ckap_yz = ones(n,1) * (kapT(2,3)*cm);
ckap_zx = ones(n,1) * (kapT(3,1)*cm);
ckap_zy = ones(n,1) * (kapT(3,2)*cm);
ckap_zz = ones(n,1) * (kapT(3,3)*cm);
cmua = ones(n,1) * (mua*cm);
ref = ones(n,1) * refind;
zeta = ones(n,1) * (cm./(2*toastDotBndterm(refind,'Keijzer')));
omega = freq*1e6*1e-12*2*pi;
smat = sparse(n,n);
for i=1:ne
    e = mesh.Element(i);
    int_kapdd_xx = e.Mat('Pdd', ckap_xx); int_kapdd_xx = squeeze(int_kapdd_xx(:,1,:,1));
    int_kapdd_xy = e.Mat('Pdd', ckap_xy); int_kapdd_xy = squeeze(int_kapdd_xy(:,1,:,2));
    int_kapdd_xz = e.Mat('Pdd', ckap_xz); int_kapdd_xz = squeeze(int_kapdd_xz(:,1,:,3));
    int_kapdd_yx = e.Mat('Pdd', ckap_yx); int_kapdd_yx = squeeze(int_kapdd_yx(:,2,:,1));
    int_kapdd_yy = e.Mat('Pdd', ckap_yy); int_kapdd_yy = squeeze(int_kapdd_yy(:,2,:,2));
    int_kapdd_yz = e.Mat('Pdd', ckap_yz); int_kapdd_yz = squeeze(int_kapdd_yz(:,2,:,3));
    int_kapdd_zx = e.Mat('Pdd', ckap_zx); int_kapdd_zx = squeeze(int_kapdd_zx(:,3,:,1));
    int_kapdd_zy = e.Mat('Pdd', ckap_zy); int_kapdd_zy = squeeze(int_kapdd_zy(:,3,:,2));
    int_kapdd_zz = e.Mat('Pdd', ckap_zz); int_kapdd_zz = squeeze(int_kapdd_zz(:,3,:,3));

    int_kapDD = int_kapdd_xx + int_kapdd_xy + int_kapdd_xz + ...
        int_kapdd_yx + int_kapdd_yy + int_kapdd_yz + ...
        int_kapdd_zx + int_kapdd_zy + int_kapdd_zz;
    
    int_muaFF = e.Mat('PFF', cmua);
    int_omegaFF = omega.*e.Mat('FF');
    bint_zetaFF = e.Mat('BndPFF', zeta);
    elidx = e.Dof;
    smat(elidx, elidx) = smat(elidx, elidx) + ...
        int_kapDD + int_muaFF + bint_zetaFF + 1i*int_omegaFF;
end
phi = smat\qvec;
end

end
