function writetoastmesh2d(Fileout,ntype,Vertex,nreg,Tri,mua,diff,rind)

NoV = size(Vertex,1);
NoF = size(Tri,1);
% assume it is "nodal" mesh : i.e. parameters defined by node
% output TOAST format
fid = fopen(Fileout,'w');
fprintf(fid,'MeshData 5.0\n\nNodeList %d 1\n',NoV);
for i=1:NoV
  fprintf(fid,'%c[%f %f]R%d\n',ntype{i},Vertex(i,1),Vertex(i,2),nreg(i));
end
fprintf(fid,'\nElementList %d\n',NoF);
for i=1:NoF
fprintf(fid,'o %6d %6d %6d\n',Tri(i,1),Tri(i,2),Tri(i,3));
end
fprintf(fid,'\n[ParameterList]N\nSize %d\nParam1 MUA\nParam2 KAPPA\nParam3 N\nData\n',NoV);
for i=1:NoV
fprintf(fid,'%f %f %f\n',mua(i),diff(i),rind(i));
end
fclose(fid);
