function maketworeg_elmesh(Fileout,ntype,Vert,nreg,Tet,mua,diff,rind)
%
% take the parameters of a mesh and output two region "element" mesh.
%
NoV = size(Vert,1);
NoF = size(Tet,1);
%
% set some element wise parameters
%
radsq = 9.8*9.8;
for k = 1:size(Tet,1) 
  distsq = [0,0,0,0];
   for j = 1:4
     distsq(j) = Vert(Tet(k,j),:)*Vert(Tet(k,j),:)';
   end
%   distsq
   if( distsq(1) < radsq && distsq(2) < radsq && distsq(3) < radsq && distsq(4) < radsq ) 
     ereg(k) = 1;
     emua(k) = 0.015;
     ediff(k) = 1/(3*(0.8+0.015));
     erind(k) = 1.2;
   else
     ereg(k) = 0;
     emua(k) = 0.01;
     ediff(k) = 1/(3*(1+0.01));
     erind(k) = 1.4;
  end
end

writetoastelmesh3d(Fileout,ntype,Vert,ereg,Tet,emua,ediff,erind);
