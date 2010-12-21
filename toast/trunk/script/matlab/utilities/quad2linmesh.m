
function [QVertex,LTri] = quad2linmesh(QNoV, QNoF, Filein, Fileout)
%
% read quadratic BEM mesh and convert to lin by splitting into 4.
%

% the mesh files 
fid = fopen(Filein) ;% opens the file for reading: external mesh
QTri = [];	% holds original quadratic element indices
QVertex = [];	% holds original quadratic nodes
for i = 1:QNoV
	tline = fgetl(fid);
	QVertex(i,:) = str2num(tline);
end
for i = 1:QNoF
	tline = fgetl(fid);
	QTri(i,:) = str2num(tline);
end
min(min(QTri))

LNoV = QNoV
LNoF = 4*QNoF
LTri = zeros(LNoF,3);

for i = 1:QNoF
	LTri(4*i-3,:) = [QTri(i,1) QTri(i,2) QTri(i,6)];
	LTri(4*i-2,:) = [QTri(i,2) QTri(i,4) QTri(i,6)];
	LTri(4*i-1,:) = [QTri(i,2) QTri(i,3) QTri(i,4)];
	LTri(4*i,:)   = [QTri(i,6) QTri(i,4) QTri(i,5)];
end

fis = fopen(Fileout,'w');

for i=1:LNoV
    fprintf(fis,'v %3.5f %3.5f %3.5f\n',QVertex(i,1),QVertex(i,2), QVertex(i,3));
end;

LTri = LTri + 1;
for i=1:LNoF
    fprintf(fis,'f %d %d %d \n',LTri(i,1),LTri(i,2),LTri(i,3));
end;

fclose (fis);


