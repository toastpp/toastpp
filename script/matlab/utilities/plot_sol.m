
function plot_sol(elem,n_elem,nodes,n_nodes,x_sol,source)

% draw the solution on the surface.
% picture of 6 nodes triangles 

for i1=1:n_nodes% turning ampitude to log scale
    x(i1)=nodes(i1,1);
    y(i1)=nodes(i1,2);
    z(i1)=nodes(i1,3);
    x_amp(i1)=log(abs(x_sol(i1)));
  %  x_pha(i1)=angle(x_sol(i1));  
  %  x_amp(i1)=log(abs(real(x_sol(i1))));
  %  x_pha(i1)=log(abs(imag(x_sol(i1))));
x_pha(i1)=angle(abs(real(x_sol(i1)))+i*imag(x_sol(i1)));
end


figure;

for j=1:n_elem %drawing the quadrilaterals
  
  hold on;  

  ie=elem(j,1)+1;
  je=elem(j,2)+1;
  ke=elem(j,3)+1;
  le=elem(j,4)+1;
  me=elem(j,5)+1;
  ne=elem(j,6)+1;
  
  xi=x(ie);
  yi=y(ie);
  zi=z(ie);
  pil=x_amp(ie);
  pilh=x_pha(ie);
  
  xj=x(je);
  yj=y(je);
  zj=z(je);
  pjl=x_amp(je);
  pjlh=x_pha(je);
 
  xk=x(ke);
  yk=y(ke);
  zk=z(ke);
  pkl=x_amp(ke);
  pklh=x_pha(ke);
  
  xl=x(le);
  yl=y(le);
  zl=z(le);
  pll=x_amp(le);
  pllh=x_pha(le);
  
  xm=x(me);
  ym=y(me);
  zm=z(me);
  pml=x_amp(me);
  pmlh=x_pha(me);

  xn=x(ne);
  yn=y(ne);
  zn=z(ne);
  pnl=x_amp(ne);
  pnlh=x_pha(ne);
 
  px=[xi,xj,xn,xi];
  py=[yi,yj,yn,yi];
  pz=[zi,zj,zn,zi];
  pal= [pil,pjl,pnl,pil];
  palh=[pilh,pjlh,pnlh,pilh];
  fill3(px,py,pz,pal);
    
  px=[xj,xl,xn,xj];
  py=[yj,yl,yn,yj];
  pz=[zj,zl,zn,zj];
  pal= [pjl,pll,pnl,pjl];
  palh=[pjlh,pllh,pnlh,pjlh];
  fill3(px,py,pz,pal);
  
  px=[xj,xk,xl,xj];
  py=[yj,yk,yl,yj];
  pz=[zj,zk,zl,zj];
  pal= [pjl,pkl,pll,pjl];
  palh=[pjlh,pklh,pllh,pjlh];
  fill3(px,py,pz,pal);
  
  px=[xn,xl,xm,xn];
  py=[yn,yl,ym,yn];
  pz=[zn,zl,zm,zn];
  pal= [pnl,pll,pml,pnl];
  palh=[pnlh,pllh,pmlh,pnlh];
  fill3(px,py,pz,pal);

     
hold on;
 
end

  plot3(source(1),source(2),source(3),'*r');

 figure;

for j=1:n_elem %drawing the quadrilaterals
  
  hold on;


  ie=elem(j,1)+1;
  je=elem(j,2)+1;
  ke=elem(j,3)+1;
  le=elem(j,4)+1;
  me=elem(j,5)+1;
  ne=elem(j,6)+1;
  
  xi=x(ie);
  yi=y(ie);
  zi=z(ie);
  pil=x_amp(ie);
  pilh=x_pha(ie);
  
  xj=x(je);
  yj=y(je);
  zj=z(je);
  pjl=x_amp(je);
  pjlh=x_pha(je);
 
  xk=x(ke);
  yk=y(ke);
  zk=z(ke);
  pkl=x_amp(ke);
  pklh=x_pha(ke);
  
  xl=x(le);
  yl=y(le);
  zl=z(le);
  pll=x_amp(le);
  pllh=x_pha(le);
  
  xm=x(me);
  ym=y(me);
  zm=z(me);
  pml=x_amp(me);
  pmlh=x_pha(me);

  xn=x(ne);
  yn=y(ne);
  zn=z(ne);
  pnl=x_amp(ne);
  pnlh=x_pha(ne);
 
  px=[xi,xj,xn,xi];
  py=[yi,yj,yn,yi];
  pz=[zi,zj,zn,zi];
  pal= [pil,pjl,pnl,pil];
  palh=[pilh,pjlh,pnlh,pilh];
  fill3(px,py,pz,palh);
    
  px=[xj,xl,xn,xj];
  py=[yj,yl,yn,yj];
  pz=[zj,zl,zn,zj];
  pal= [pjl,pll,pnl,pjl];
  palh=[pjlh,pllh,pnlh,pjlh];
  fill3(px,py,pz,palh);
  
  px=[xj,xk,xl,xj];
  py=[yj,yk,yl,yj];
  pz=[zj,zk,zl,zj];
  pal= [pjl,pkl,pll,pjl];
  palh=[pjlh,pklh,pllh,pjlh];
  fill3(px,py,pz,palh);
  
  px=[xn,xl,xm,xn];
  py=[yn,yl,ym,yn];
  pz=[zn,zl,zm,zn];
  pal= [pnl,pll,pml,pnl];
  palh=[pnlh,pllh,pmlh,pnlh];
  fill3(px,py,pz,palh);
  
  hold on
 
end
  plot3(source(1),source(2),source(3),'*r');
 
