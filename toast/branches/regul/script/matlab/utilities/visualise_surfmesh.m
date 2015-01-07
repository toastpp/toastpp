function visualise_surfmesh( NoV,Vertex,NoF,Tri)

x0=Vertex(:,1);
y0=Vertex(:,2);
z0=Vertex(:,3);

NoV
NoF

% visualise

vtx=Vertex;
srf=Tri(:,1:3);


fig_title='Surface Mesh';
f1=figure('Name',fig_title);

ta=vtx(srf(:,1),:);
tb=vtx(srf(:,2),:);
tc=vtx(srf(:,3),:);
cx=cross(tc-tb,tb-ta);
for j = 1:length(cx(:,1))
    cxlon(j)=norm(cx(j,:));
    cx(j,:) = 5*cx(j,:)/cxlon(j);
end
%    size(cx);
%    size(cxlon);

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
