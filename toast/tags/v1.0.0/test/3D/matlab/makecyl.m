clear all
close all

meshname = '../meshes/cyl2.msh';
nx = 32; ny = 32; nz = 32;

hMesh = toastReadMesh(meshname);
hBasis = toastSetBasis(hMesh,[nx ny nz]);

mode = 1; % 'truth' prior
[mua mus] = cylphantom(nx,ny,nz,mode);

ikap = 3.0*(mua+mus+eps);
z = find(ikap == 0);
ikap(z) = eps; % avoid division by zero errors
kap = 1./ikap;

hmua = toastMapGridToMesh(hBasis,mua);
hmus = toastMapGridToMesh(hBasis,mus);

% write out the target images as NIM files
toastWriteNIM('mua_tgt.nim',meshname,hmua);
toastWriteNIM('mus_tgt.nim',meshname,hmus);

% wite out the target images in raw raster format
% (used by the regulariser)
toastWriteRealVector('mua_prior_truth_80x80x80.dat',mua);
toastWriteRealVector('mus_prior_truth_80x80x80.dat',mus);
toastWriteRealVector('kap_prior_truth_80x80x80.dat',kap);


mode = 2; % 'partial' prior
[mua mus] = cylphantom(nx,ny,nz,mode);
ikap = 3.0*(mua+mus+eps);
ikap(z) = eps; % avoid division by zero errors
kap = 1./ikap;
toastWriteRealVector('mua_prior_part_80x80x80.dat',mua);
toastWriteRealVector('mus_prior_part_80x80x80.dat',mus);
toastWriteRealVector('kap_prior_part_80x80x80.dat',kap);


mode = 3; % 'sum' prior
[mua mus] = cylphantom(nx,ny,nz,mode);
ikap = 3.0*(mua+mus+eps);
ikap(z) = eps; % avoid division by zero errors
kap = 1./ikap;
toastWriteRealVector('mua_prior_sum_80x80x80.dat',mua);
toastWriteRealVector('mus_prior_sum_80x80x80.dat',mus);
toastWriteRealVector('kap_prior_sum_80x80x80.dat',kap);
