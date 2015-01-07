bdim = [128 128];
hMesh = toastReadMesh ('../meshes/ellips_tri10.msh');
hBasis = toastSetBasis (hMesh,bdim);

tgtmua = toastReadNIM ('../meshes/tgt_mua_ellips_tri10.nim');
tgtmus = toastReadNIM ('../meshes/tgt_mus_ellips_tri10.nim');
tgtkap = 1./(3*(tgtmua+tgtmus));

btgtmua = toastMapMeshToGrid (hBasis,tgtmua);
btgtkap = toastMapMeshToGrid (hBasis,tgtkap);
btgt = [btgtmua; btgtkap];
toastWriteRealVector('prior_true_128x128.dat',btgt);

toastDeleteBasis (hBasis);
toastDeleteMesh (hMesh);
