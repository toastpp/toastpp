hMesh = toastReadMesh ('neonatebrain_quad.tmsh');
phase = toastReadNIM ('phi_phase_100.nim');
lnamp = toastReadNIM ('phi_lnmod_100.nim');

[vtx,idx,perm] = toastSurfData(hMesh);
nvtx = size(vtx,1);

sphase = phase(perm);
slnamp = lnamp(perm);

col = zeros (nvtx);
cmin = min(sphase);
cmax = max(sphase);
col = fix((sphase-cmin)/(cmax-cmin)*length(colormap))+1;

% cut 6-noded triangles into 4 3-noded ones
ntri = size(idx,1);
idx3 = zeros(ntri*4,3);
for i=1:ntri
    idxi = idx(i,:);
    ofs = (i-1)*4+1;
    idx3(ofs,:)   = idxi([1 4 6]);
    idx3(ofs+1,:) = idxi([4 2 5]);
    idx3(ofs+2,:) = idxi([4 5 6]);
    idx3(ofs+3,:) = idxi([6 5 3]);
end

patch('Faces',idx3,'Vertices',vtx,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','none')
axis equal

% read bem mesh geometry, surface 000
sphere000
bem000_sphase = toastSampleField (hMesh, phase, p);
nvtx = size(bem000_sphase,1);
col = zeros (nvtx);
cmin = min(bem000_sphase);
cmax = max(bem000_sphase);
col = fix((bem000_sphase-cmin)/(cmax-cmin)*length(colormap))+1;
figure

% cut 6-noded triangles into 4 3-noded ones
ntri = size(node,1);
node3 = zeros(ntri*4,3);
for i=1:ntri
    nodei = node(i,:);
    ofs = (i-1)*4+1;
    node3(ofs,:) = nodei([1 2 6]);
    node3(ofs+1,:) = nodei([2 3 4]);
    node3(ofs+2,:) = nodei([2 4 6]);
    node3(ofs+3,:) = nodei([4 5 6]);
end

patch('Faces',node3+1,'Vertices',p,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','none')
axis equal

% read bem mesh geometry, surface 001
clear p
clear node
sphere001
bem001_sphase = toastSampleField (hMesh, phase, p);
nvtx = size(bem001_sphase,1);
col = zeros (nvtx);
%cmin = min(bem001_sphase);
%cmax = max(bem001_sphase);
col = fix((bem001_sphase-cmin)/(cmax-cmin)*length(colormap))+1;
figure

% cut 6-noded triangles into 4 3-noded ones
ntri = size(node,1);
node3 = zeros(ntri*4,3);
for i=1:ntri
    nodei = node(i,:);
    ofs = (i-1)*4+1;
    node3(ofs,:) = nodei([1 2 6]);
    node3(ofs+1,:) = nodei([2 3 4]);
    node3(ofs+2,:) = nodei([2 4 6]);
    node3(ofs+3,:) = nodei([4 5 6]);
end

patch('Faces',node3+1,'Vertices',p,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','none')
axis equal

% read bem mesh geometry, surface 002
clear p
clear node
sphere002
bem002_sphase = toastSampleField (hMesh, phase, p);
nvtx = size(bem002_sphase,1);
col = zeros (nvtx);
%cmin = min(bem002_sphase);
%cmax = max(bem002_sphase);
col = fix((bem002_sphase-cmin)/(cmax-cmin)*length(colormap))+1;
figure

% cut 6-noded triangles into 4 3-noded ones
ntri = size(node,1);
node3 = zeros(ntri*4,3);
for i=1:ntri
    nodei = node(i,:);
    ofs = (i-1)*4+1;
    node3(ofs,:) = nodei([1 2 6]);
    node3(ofs+1,:) = nodei([2 3 4]);
    node3(ofs+2,:) = nodei([2 4 6]);
    node3(ofs+3,:) = nodei([4 5 6]);
end

patch('Faces',node3+1,'Vertices',p,'FaceVertexCData',col,'FaceColor','interp','EdgeColor','none')
axis equal

toastDeleteMesh (hMesh);
