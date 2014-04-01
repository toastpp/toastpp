function writetoastelmesh3d(Fileout,ntype,Vertex,ereg,Tet,mua,diff,rind)

NoV = size(Vertex,1);
NoF = size(Tet,1);
% This assumes it is "element" mesh : i.e. parameters defined by element
% output TOAST format
fid = fopen(Fileout,'w');
fprintf(fid,'MeshData 5.0\n\nNodeList %d 1\n',NoV);
for i=1:NoV
  fprintf(fid,'%c[%f %f %f]\n',ntype{i},Vertex(i,1),Vertex(i,2),Vertex(i,3));
end
fprintf(fid,'\nElementList %d\n',NoF);
for i=1:NoF
fprintf(fid,'c %6d %6d %6d %6d R%d\n',Tet(i,1),Tet(i,2),Tet(i,3),Tet(i,4),ereg(i));
end
fprintf(fid,'\n[ParameterList]E\nSize %d\nParam1 MUA\nParam2 KAPPA\nParam3 N\nData\n',NoF);
for i=1:NoF
fprintf(fid,'%f %f %f\n',mua(i),diff(i),rind(i));
end
fclose(fid);
