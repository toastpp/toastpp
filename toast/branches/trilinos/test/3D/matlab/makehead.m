clear all
close all

meshname = '../meshes/vox32_2blobs.msh';
nx = 25; ny = 30; nz = 20;

hMesh = toastReadMesh(meshname);
hBasis = toastSetBasis(hMesh,[nx ny nz]);

mode = 1; % 'truth' prior
[bmua bmus] = headphantom(hMesh,nx,ny,nz,mode);

hmua = toastMapGridToMesh(hBasis,bmua);
hmus = toastMapGridToMesh(hBasis,bmus);

% write out the target images as NIM files
toastWriteNIM('mua_tgt.nim',meshname,hmua);
toastWriteNIM('mus_tgt.nim',meshname,hmus);

save head_vox32.mat nx ny nz bmua bmus
