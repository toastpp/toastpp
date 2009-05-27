
function plot_surf_sol(elem,n_elem,nodes,n_nodes,x_sol)

% draw the solution on the surface.
% picture of 3 nodes triangles 

for i1=1:n_nodes
    x(i1)=nodes(i1,1);
    y(i1)=nodes(i1,2);
    z(i1)=nodes(i1,3);
end


figure;

for j=1:n_elem %drawing the quadrilaterals
  
  hold on;  

  ie=elem(j,1);
  je=elem(j,2);
  ke=elem(j,3);
  
  xi=x(ie);
  yi=y(ie);
  zi=z(ie);
  pil=x_sol(ie);
  
  xj=x(je);
  yj=y(je);
  zj=z(je);
  pjl=x_sol(je);
 
  xk=x(ke);
  yk=y(ke);
  zk=z(ke);
  pkl=x_sol(ke);  
  
  px=[xi,xj,xk,xi];
  py=[yi,yj,yk,yi];
  pz=[zi,zj,zk,zi];
  pal= [pil,pjl,pkl,pil];
  fill3(px,py,pz,pal);
     
hold on;
 
end
