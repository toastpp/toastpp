function [ntype,Vert,ereg,Tri,mua,diff,rind,ptype] = readtoastmesh2d(Filein)

% read mesh 

% the mesh files 
%fid = fopen(lou) ;% opens the file for reading: external mesh
% output TOAST format
fid = fopen(Filein,'r');

% read the vertices
fgets(fid);
fgets(fid);
fscanf(fid,'%s',1);
NoV = fscanf(fid,'%d',1);
Vert=zeros(NoV,2);
ntype=cell(NoV,1);
nreg =zeros(NoV,1);

fgetl(fid); % skip end of line
for i=1:NoV
     str = fgetl(fid);
%     display(str);
     [a1,a2,c1,c2,b1] = strread(str,'%c%c%f %f%c');
     ntype{i} = a1;
     Vert(i,1) = c1;
     Vert(i,2) = c2;
end

% read the elements

fgets(fid);
fscanf(fid,'%s',1);
NoF = fscanf(fid,'%d',1);
Tri=zeros(NoF,3);
etype=cell(NoF,1);

fgetl(fid); % skip end of line
for i=1:NoF
     str = fgetl(fid);
%     display(str);
     [e v1 v2 v3 b1 d] = strread(str,'%c %d %d %d %c%d');
     etype{i} = e;
     Tri(i,1) = v1;
     Tri(i,2) = v2;
     Tri(i,3) = v3;
     ereg(i) = d;
end

% read parameters

mua = zeros(NoV,1);
diff = zeros(NoV,1);
rind = zeros(NoV,1);

fgets(fid);
str=fgets(fid);
if(str(16)=='E') 
    ptype = 1;
else
    ptype = 0;
end
str = fgets(fid);
[a1 NoP] = strread(str,'%s %d');
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);
for i=1:NoP
     str = fgetl(fid);
%     display(str);
     [a d r] = strread(str,'%f %f %f');
     mua(i) = a;
     diff(i) = d;
     rind(i) = r;
end

fclose(fid);
