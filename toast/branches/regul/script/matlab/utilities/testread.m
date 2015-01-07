function [ntype,Vert,nreg,Tri,mua,diff,rind] = testread(Filein)

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
     [a1,a2,c1,c2,b1,b2,d] = strread(str,'%c%c%f %f%c%c%d');
     ntype{i} = a1;
     Vert(i,1) = c1;
     Vert(i,2) = c2;
     nreg(i) = d;
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
     [e v1 v2 v3] = strread(str,'%c %d %d %d');
     etype{i} = e;
     Tri(i,1) = v1;
     Tri(i,2) = v2;
     Tri(i,3) = v3;
end

% read parameters

mua = zeros(NoV,1);
diff = zeros(NoV,1);
rind = zeros(NoV,1);

fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);
for i=1:NoV
     str = fgetl(fid);
%     display(str);
     [a d r] = strread(str,'%f %f %f');
     mua(i) = a;
     diff(i) = d;
     rind(i) = r;
end

fclose(fid);

% plot
figure(1);
hold on;
for k = 1:size(Tri,1)
   for j = 1:3
     x1 = Vert(Tri(k,j),1);
     y1 = Vert(Tri(k,j),2);
     x2 = Vert(Tri(k,mod(j,3)+1),1);
     y2 = Vert(Tri(k,mod(j,3)+1),2);
     plot([x1,x2],[y1,y2],'g');
  end
end;
plot(Vert(:,1),Vert(:,2),'+');
hold off;
