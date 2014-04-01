
function toastquad2linmesh(Filein, Fileout)
%
% read quadratic BEM mesh and convert to lin by splitting into 4.
%

% the mesh files 

[ntype,QVertex,nreg,QTri,mua,diff,rind,ptype] = readtoastquadmeshsurf3d(Filein);
QNoV = size(QVertex,1);
QNoF = size(QTri,1);
disp(['minimum vertex ' num2str(min(min(QTri)))]);

LNoV = QNoV
LNoF = 4*QNoF
LTri = zeros(LNoF,3);

for i = 1:QNoF
	LTri(4*i-3,:) = [QTri(i,1) QTri(i,4) QTri(i,6)];
	LTri(4*i-2,:) = [QTri(i,2) QTri(i,5) QTri(i,4)];
	LTri(4*i-1,:) = [QTri(i,3) QTri(i,6) QTri(i,5)];
	LTri(4*i,:)   = [QTri(i,6) QTri(i,4) QTri(i,5)];
end

%LTri = LTri + 1;

fid = fopen(Fileout,'w');
fprintf(fid,'MeshData 5.0\n\nNodeList %d 1\n',LNoV);
for i=1:LNoV
    fprintf(fid,'B[%f %f %f]R0\n',QVertex(i,1),QVertex(i,2),QVertex(i,3));
end
fprintf(fid,'\nElementList %d\n',LNoF);
for i=1:LNoF
    fprintf(fid,'p %6d %6d %6d\n',LTri(i,1),LTri(i,2),LTri(i,3));
end
fprintf(fid,'\n[ParameterList]\nSize %d\nParam1 MUA\nParam2 KAPPA\nParam3 N\nData\n',LNoV);
for i=1:LNoV
    fprintf(fid,'%f %f %f\n',0.01, 0.330033, 1.4);
end
fprintf(fid,'\n\n');
fclose(fid);
