function lin2quadmesh(Filein, Fileout)
%
% read triangular BEM mesh and convert to quadratic.
%
%[LNoV, LVertex, LNoF, node] = read_linear_mesh

fid = fopen(Filein) ;
Tri = [];
Vertex = [];
tline = fgetl(fid);
NoV = str2num(tline);
for i=1:NoV
    tline = fgetl(fid);
    Vertex(i,:) = str2num(tline);
end
tline = fgetl(fid);
NoF = str2num(tline);
for i=1:NoF
    tline = fgetl(fid);
    Tri(i,:) = str2num(tline); 
end
fclose(fid);
disp('number of vertices');NoV
disp('number of faces');NoF

tri2quadmesh(NoV, Vertex, NoF, Tri, Fileout);

