function jan2mesh(NoV,NoF,FileName,Fileout)
%
% convert mesh to 'TOAST' format
%
fid = fopen(FileName) ;
Tri = [];
Vertex = [];
%tline = fgetl(fid);
%NoV = str2num(tline);
for i=1:NoV
    tline = fgetl(fid);
    Vertex(i,:) = str2num(tline);
end
%tline = fgetl(fid);
%NoF = str2num(tline);
for i=1:NoF
    tline = fgetl(fid);
    Tri(i,:) = str2num(tline); 
end
fclose(fid);
Tri = Tri+1;
disp('number of vertices');NoV
disp('number of faces');NoF
% output TOAST format
fid = fopen(Fileout,'w');
fprintf(fid,'MeshData 5.0\n\nNodeList %d 1\n',NoV);
for i=1:NoV
  fprintf(fid,'B[%f %f %f]R0\n',Vertex(i,1),Vertex(i,2),Vertex(i,3));
end
fprintf(fid,'\nElementList %d\n',NoF);
for i=1:NoF
fprintf(fid,'q %6d %6d %6d %6d %6d %6d\n',Tri(i,1),Tri(i,3),Tri(i,5),Tri(i,2),Tri(i,4),Tri(i,6));
end
fprintf(fid,'\n[ParameterList]\nSize %d\nParam1 MUA\nParam2 KAPPA\nParam3 N\nData\n',NoV);
for i=1:NoV
fprintf(fid,'0.01 0.330033 1\n');
end
fclose(fid);
