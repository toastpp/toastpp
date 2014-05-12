
function [nodes,elem] = readmesh(lou);

% read mesh 

% the mesh files 
fid = fopen(lou) ;% opens the file for reading: external mesh

%external mesh
NoV=fscanf(fid,'%i',1)

x0=zeros(NoV,1);
y0=zeros(NoV,1);
z0=zeros(NoV,1);

%Reading from the file 
% x,y,z
for j=1:(NoV),
    fscanf(fid,'%c',[2]);
    
    a=fscanf(fid,'%g',[3]);
    x0(j,1)=a(1);
    y0(j,1)=a(2);
    z0(j,1)=a(3);   
    fscanf(fid,'%c',[1]);
    
    % fscanf(fid,'%c',[2]);
end;

nodes(:,1) = x0;
nodes(:,2) = y0;
nodes(:,3) = z0;

NoE=fscanf(fid,'%i',1)

elem=zeros(NoE,6);

for j=1:(NoE),
    fscanf(fid,'%c',[2]);
    a=fscanf(fid,'%g',[6]);
    elem(j,1:6)=a(1:6)';         
    fscanf(fid,'%c',[1]);
end;
fclose(fid); % closes the file




plot3(x0,y0,z0,'g*');
hold on;

     % picture of 6 nodes triangles

    % Plot the numbers
    
% for d=1:6
%     ii1= int2str(elem(1,d));
%     text(x0(elem(1,d)+1),y0(elem(1,d)+1),z0(elem(1,d)+1),ii1);
% end; 
% for d=1:NoV
%     ii1= int2str(d-1);
%     text(x0(d),y0(d),z0(d),ii1);
% end;   


%elem(1,:)
% hold on

for i=1:NoE%drawing the triangles
    
   ie=elem(i,1)+1;
   je=elem(i,2)+1;
   ke=elem(i,3)+1;
   le=elem(i,4)+1;
   me=elem(i,5)+1;
   ne=elem(i,6)+1;
    
   xi=x0(ie);
   yi=y0(ie);
   zi=z0(ie);
    
   xj=x0(je);
   yj=y0(je);
   zj=z0(je);
    
   xk=x0(ke);
   yk=y0(ke);
   zk=z0(ke);
    
   xl=x0(le);
   yl=y0(le);
   zl=z0(le);
    
   xm=x0(me);
   ym=y0(me);
   zm=z0(me);
    
   xn=x0(ne);
   yn=y0(ne);
   zn=z0(ne);
    
   px=[xi,xj,xk,xl,xm,xn,xi];
   py=[yi,yj,yk,yl,ym,yn,yi];
   pz=[zi,zj,zk,zl,zm,zn,zi];
    
  fill3(px,py,pz,[0 0 1]);
   % plot3(px,py,pz);                   
   hold on
    
end

box on 
axis equal



