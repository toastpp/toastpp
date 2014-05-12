% Shows the surface of a 3-D mesh and the source and detector locations.

clear all
close all

meshname = input ('Mesh name: ','s');
qmname = input ('QM name: ','s');

hMesh = toastReadMesh (meshname);
toastReadQM (hMesh, qmname);

[vtx,idx,perm] = toastSurfData (hMesh);
nvtx = size(vtx,1);

% show the mesh surface
col = zeros(nvtx,3); col(:,2)=1;
patch('Vertices',vtx,'Faces',idx,'FaceVertexCData',col,'FaceColor','interp')
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
hold on

% show source locations
qp = toastQPos (hMesh);
plot3 (qp(:,1), qp(:,2), qp(:,3), '*r')

% show detector locations
mp = toastMPos (hMesh);
plot3 (mp(:,1), mp(:,2), mp(:,3), '*b')
