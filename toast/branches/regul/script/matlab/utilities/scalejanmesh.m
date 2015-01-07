function scalejanmesh(NoV,NoF,FileName,Fileout,scale)

% scale mesh
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
%Tri = Tri+1;
disp('number of vertices');NoV
disp('number of faces');NoF
% output
fid = fopen(Fileout,'w');
for i=1:NoV
  fprintf(fid,'[%12.6f %12.6f %12.6f]\n',scale*Vertex(i,1),scale*Vertex(i,2),scale*Vertex(i,3));
end
for i=1:NoF
fprintf(fid,'[%6d %6d %6d %6d %6d %6d]\n',Tri(i,1),Tri(i,2),Tri(i,3),Tri(i,4),Tri(i,5),Tri(i,6));
end
fclose(fid);
