
%plot_all_mesh.m

% Picture both the brain and the head mesh 
figure(1);
%first the brain mesh 

% sk4_lod_6n_anti
% sc4_lod_6n_anti
sphere000
bem3D

C=x_amp0;
C2=log(C);
%C2=x_pha0;

n_elements=n_e0;
n_nodes=n_n0;

for i=1:n_nodes
    
    x(i)=p(i,1);
    y(i)=p(i,2);
    z(i)=p(i,3);
end
% plot3(x,y,z,'b.');
scale=1.0;

for i=1:n_elements
    
    
    ie=node(i,1)+1;
    je=node(i,2)+1;
    ke=node(i,3)+1;
    le=node(i,4)+1;
    me=node(i,5)+1;
    ne=node(i,6)+1;
    
    xi=x(ie)*scale;
    yi=y(ie)*scale;
    zi=z(ie)*scale;
    pal(1)=C2(ie);
    pal(7)=pal(1);
    
    xj=x(je)*scale;
    yj=y(je)*scale;
    zj=z(je)*scale;
    pal(2)=C2(je);
    
    xk=x(ke)*scale;
    yk=y(ke)*scale;
    zk=z(ke)*scale;
    pal(3)=C2(ke);
    
    xl=x(le)*scale;
    yl=y(le)*scale;
    zl=z(le)*scale;
    pal(4)=C2(le);
    
    xm=x(me)*scale;
    ym=y(me)*scale;
    zm=z(me)*scale;
    pal(5)=C2(me);
    
    xn=x(ne)*scale;
    yn=y(ne)*scale;
    zn=z(ne)*scale;
    pal(6)=C2(ne);
    
    px=[xi,xj,xk,xl, xm, xn, xi];
    py=[yi,yj,yk,yl, ym, yn, yi];
    pz=[zi,zj,zk,zl, zm, zn, zi];
    
    fill3(px,py,pz,pal);
%     plot3(px,py,pz,'k.-');
    hold on
end
plot3(x,y,z,'k.');
XMIN=-45.;
XMAX=45.;
YMIN=-55.;
YMAX=55.;
ZMIN=-65.;
ZMAX=65.;
 box on 
  axis equal
 hold on 
 b1=25;
 b2=2;
%  beta=-.4;
% caxis([-b1 -b2]);
 axis([XMIN XMAX YMIN YMAX ZMIN ZMAX]); 
  colorbar('vert');
%  colormap(gray);
% brighten(beta);
hold on

