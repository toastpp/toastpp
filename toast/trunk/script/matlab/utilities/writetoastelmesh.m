function writetoastelmesh(Fileout,ntype,Vertex,ereg,Tri,mua,diff,rind)

NoV = size(Vertex,1);
NoF = size(Tri,1);
% This assumes it is "element" mesh : i.e. parameters defined by element
% output TOAST format
fid = fopen(Fileout,'w');
fprintf(fid,'MeshData 5.0\n\nNodeList %d 1\n',NoV);
for i=1:NoV
  fprintf(fid,'%c[%f %f]\n',ntype{i},Vertex(i,1),Vertex(i,2));
end
fprintf(fid,'\nElementList %d\n',NoF);
for i=1:NoF
fprintf(fid,'o %6d %6d %6d R%d\n',Tri(i,1),Tri(i,2),Tri(i,3),ereg(i));
end
fprintf(fid,'\n[ParameterList]E\nSize %d\nParam1 MUA\nParam2 KAPPA\nParam3 N\nData\n',NoF);
for i=1:NoF
fprintf(fid,'%f %f %f\n',mua(i),diff(i),rind(i));
end
fclose(fid);
