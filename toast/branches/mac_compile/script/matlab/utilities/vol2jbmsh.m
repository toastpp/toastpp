function [NoV, verts, NoF, els] = vol2jbmsh(filename)
% 
fid = fopen(filename);
for i = 1:6 % 6 useless lines
 tline = fgetl(fid);
end
tline = fgetl(fid);
NoF = str2num(tline);
NoF
els = zeros(NoF,3);
% read 3 node triangular elements
for i=1:NoF
    tline = fgetl(fid);
%    [surfnr    bcnr   domin  domout      np      p1      p2      p3 p4 p5 p6] = str2num(tline);
    narray = str2num(tline);
    els(i,:) = narray(6:8);
end
for i = 1:9 % 9 useless lines
 tline = fgetl(fid);
end
tline = fgetl(fid);
NoE = str2num(tline);
NoE
for i = 1:2*NoE % I dont even know what these are!
 tline = fgetl(fid);
end
for i = 1:4 % 4 useless lines
 tline = fgetl(fid);
end
tline = fgetl(fid);
NoV = str2num(tline);
NoV
verts = zeros(NoV,3);
% read vertices
for i=1:NoV
    tline = fgetl(fid);
    verts(i,:) = str2num(tline);
end

fclose(fid);





