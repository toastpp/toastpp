function plottoastmesh2d(Filein)

[ntype,Vert,nreg,Tri,mua,diff,rind,pt] = readtoastmesh2d(Filein);
% plot
figure(1);
hold on;
plotmesh2d(Vert,Tri);
hold off;
