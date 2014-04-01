function [NoV, NoTr,Vertex,Triangle] = joaoplot( FileName )
% Abdel
% read and plot the .dat format from joao
% joaoplot( FileName ); 

fid = fopen(FileName);
tline = fgetl(fid);
tline = fgetl(fid);
a=str2num(tline);
NoV=a(1);
NoTr=a(2);
Triangle = [];
Vertex = [];
for i=1:NoV
    tline = fgetl(fid);
    v=str2num(tline);
    Vertex(i,1:3) =v(1:3);
end
for i=1:NoTr
    tline = fgetl(fid);
     v=str2num(tline);
    Triangle(i,:) = str2num(tline); 
end
fclose(fid);

vtx=Vertex;
srf=Triangle;


fig_title='Voila !';
f1=figure('Name',fig_title);

ta=vtx(srf(:,1),:);
tb=vtx(srf(:,2),:);
tc=vtx(srf(:,3),:);
cx=cross(tb-ta,tc-tb);
cc=(ta+tb+tc)/3;
quiver3(cc(:,1),cc(:,2),cc(:,3),cx(:,1),cx(:,2),cx(:,3),2);
hold on

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