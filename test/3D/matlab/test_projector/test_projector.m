%% Projector class test
% Note: this requires gmsh installed and on the search path

clear all
close all

if ~exist('cube_hole.msh','file')
    system('gmsh -3 cube_hole.geo');
    mesh=toastMesh('cube_hole.msh','gmsh');
    mesh.Write('cube_hole.msh');
end
mesh=toastMesh('cube_hole.msh');

if ~exist('cube.qm','file')
    fid = fopen('cube.qm','wt');
    fprintf(fid, 'QM file\nDimension 3\n\n');
    fprintf(fid, 'SourceList 1 external\n');
    fprintf(fid, '100 -0.000000 -0.000000 -1.000000 0.000000 -0.000000\n\n');
    fprintf(fid, 'MeasurementList 1 external image\n');
    tht = 30*pi/180;
    rad = 30;
    fprintf(fid, '%f %f %f %f %f %f\n\n', -rad*sin(tht), -rad*cos(tht), 0, sin(tht), cos(tht), 0);
    fprintf(fid, 'LinkList\n1: 0\n');
    fclose(fid);
end

mesh.ReadQM('cube.qm');
n = mesh.NodeCount;

projgrid = [128 128];
tic;
%hprojlist = toastMakeProjectorList(mesh,projgrid,[0 0 200],'ORTHO','pixelsize',0.25,'shift',[0 0],'driver','software');
hprojlist = toastMakeProjectorList(mesh,projgrid,[0 0 40],'PINHOLE','pixelsize',0.15,'flen',40,'shift',[0 0],'driver','software');
toc

% do a forward simulation to have something interesting to project
mesh.SetQM([-10 0 0],[10 0 0]);
mua = ones(n,1)*0.01;
mus = ones(n,1)*1;
ref = ones(n,1)*1.4;
K = dotSysmat(mesh,mua,mus,ref);
qvec = mesh.Qvec('Neumann','Gaussian',2);
phi = K\qvec;
lphi = log(phi);

%phi = ones(n,1);
figure; h=mesh.Display(lphi,'range',[-13 -1],'showgrid',true);
lighting none
view(-150,0);
campos([-42.9889 74.4590 0]);
camva(16.2);
set(gca,'Projection','Perspective');

ex = toastProjectToImage(hprojlist(1),lphi);   % FLUORESCENCE DATA CCD Camera
proj_ex = reshape(ex,projgrid);
figure; imagesc(proj_ex,[-13 -1]); axis equal tight; colorbar
