%
% plot function C2 on mesh (n_nodes, p)
%

function plotmeshsolnofig(n_nodes, p, n_elements, node, C2, cmin, cmax) 
scale = 1.0;
%figure;
for i=1:n_nodes
    
    x(i)=p(i,1);
    y(i)=p(i,2);
    z(i)=p(i,3);
end

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
%plot3(x,y,z,'k.');
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
caxis([cmin cmax]);
 axis([XMIN XMAX YMIN YMAX ZMIN ZMAX]); 
  colorbar('vert');

hold on

