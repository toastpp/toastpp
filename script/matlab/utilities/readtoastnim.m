function [vec,meshname] = readtoastnim(Filein,n)
fid = fopen(Filein,'r');
vec = [];
str = fgets(fid);
str = fscanf(fid,'%c',6);
meshname = fgets(fid);
disp(meshname);
str = fgets(fid);
str = fscanf(fid,'%c',11);
npix =  fscanf(fid,'%d')
str = fgets(fid);
for k = 0:n
    str = fscanf(fid,'%c',5);
    nim = fscanf(fid,'%d',1);
    disp('nim : '); disp(nim);
    vec = fscanf(fid,'%f');
    disp(size(vec));
end;

fclose(fid);
