function [NoV, NoTr,Vertex,Triangle] = liorplotquad( FileName )
% Lior & Abdel
% read and plot the .dat format
% liorplot( FileName ); 


fid = fopen(FileName) ;
tline = fgetl(fid);
NoV = str2num(tline);
Triangle = [];
Vertex = [];
for i=1:NoV
    tline = fgetl(fid);
    Vertex(i,:) = str2num(tline);
end
tline = fgetl(fid);
NoTr = str2num(tline);
for i=1:NoTr
    tline = fgetl(fid);
    Triangle(i,:) = str2num(tline); 
end
fclose(fid);

vtx=Vertex;
srf=Triangle(:,1:2:6)+1;


fig_title='Voila !';
f1=figure('Name',fig_title);

ta=vtx(srf(:,1),:);
tb=vtx(srf(:,2),:);
tc=vtx(srf(:,3),:);
cx=cross(tb-ta,tc-tb);
cc=(ta+tb+tc)/3;
quiver3(cc(:,1),cc(:,2),cc(:,3),cx(:,1),cx(:,2),cx(:,3),2);
hold on

%trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',1);
colormap([0 0 0]);
daspect([1 1 1]);

hold on;
axis image;
col = rand(1,3)/2 + [ 0.2 0.2 0.2];
set(gcf,'Colormap',col);
grid off
view(3);
axis off
