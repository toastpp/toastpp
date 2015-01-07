function refine_quadmesh(QNoV, QNoF, Filename, Fileout)
%
% read a quadratic BEM mesh, and split each element into 4.
%

% the mesh files 
fid = fopen(Filename) ;% opens the file for reading: external mesh
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

LNoV = QNoV;
LNoF = 4*QNoF
LTri = zeros(LNoF,3);

for i = 1:QNoF
	LTri(4*i-3,:) = [QTri(i,1) QTri(i,2) QTri(i,6)];
	LTri(4*i-2,:) = [QTri(i,2) QTri(i,4) QTri(i,6)];
	LTri(4*i-1,:) = [QTri(i,2) QTri(i,3) QTri(i,4)];
	LTri(4*i,:)   = [QTri(i,6) QTri(i,4) QTri(i,5)];
end
%QTri
%LTri
%---------------------------------------------------------------------
% CONVERT 3 node -> 6 node
%
%

x0=QVertex(:,1);
y0=QVertex(:,2);
z0=QVertex(:,3);
node = LTri+1;

d=1;
for iq=1:LNoF%drawing the triangles
    
    ie=node(iq,1);
    je=node(iq,2);
    ke=node(iq,3);
    xi=x0(ie);
    yi=y0(ie);
    zi=z0(ie);
    
    xij=(x0(je)+x0(ie))/2;
    yij=(y0(je)+y0(ie))/2;
    zij=(z0(je)+z0(ie))/2;
    new_nodes(d,:)=[xij yij zij ie je];
    d=d+1;
    xj=x0(je);
    yj=y0(je);
    zj=z0(je);
    
    xjk=(x0(je)+x0(ke))/2;
    yjk=(y0(je)+y0(ke))/2;
    zjk=(z0(je)+z0(ke))/2;
    new_nodes(d,:)=[xjk yjk zjk je ke];
    d=d+1;
    
    xk=x0(ke);
    yk=y0(ke);
    zk=z0(ke);
    
    xki=(x0(ke)+x0(ie))/2;
    yki=(y0(ke)+y0(ie))/2;
    zki=(z0(ke)+z0(ie))/2;
    new_nodes(d,:)=[xki yki zki ke ie];
    d=d+1;    
end

%add_nodes=zeros(size(new_nodes,1)/2,5);
idx=1;
for is=1:size(new_nodes,1)
    for isx=is+1:size(new_nodes,1)
        if ((new_nodes(is,4)==new_nodes(isx,5)) & (new_nodes(is,5)==new_nodes(isx,4)))
            add_nodes(idx,:)=new_nodes(is,:);
            idx=idx+1;
        end
    end
end
idx-1
%add_nodes
%clear new_nodes;
new_elem=[node(:,1) zeros(LNoF,1) node(:,2) zeros(LNoF,1) node(:,3) zeros(LNoF,1)];

for is=1:size(add_nodes,1)
    for it=1:LNoF,
        if (new_elem(it,1)==add_nodes(is,4)) & (new_elem(it,3)==add_nodes(is,5)) & (new_elem(it,2)==0)
            new_elem(it,2)=is+LNoV;
        end
        if (new_elem(it,3)==add_nodes(is,4)) &(new_elem(it,5)==add_nodes(is,5)) &(new_elem(it,4)==0)
            new_elem(it,4)=is+LNoV;
        end
        if (new_elem(it,5)==add_nodes(is,4)) & (new_elem(it,1)==add_nodes(is,5)) & (new_elem(it,6)==0)
            new_elem(it,6)=is+LNoV;
        end
        if (new_elem(it,1)==add_nodes(is,5)) & (new_elem(it,3)==add_nodes(is,4)) & (new_elem(it,2)==0)
            new_elem(it,2)=is+LNoV;
        end
        if (new_elem(it,3)==add_nodes(is,5)) &(new_elem(it,5)==add_nodes(is,4)) &(new_elem(it,4)==0)
            new_elem(it,4)=is+LNoV;
        end
        if (new_elem(it,5)==add_nodes(is,5)) & (new_elem(it,1)==add_nodes(is,4)) & (new_elem(it,6)==0)
            new_elem(it,6)=is+LNoV;
        end
    end
end
min(min(new_elem))

new_nodes=[QVertex ; [add_nodes(:,1) add_nodes(:,2) add_nodes(:,3)]];

fis = fopen(Fileout,'w');
QNoV2 = size(new_nodes,1)
%fprintf(fis,'%4.0f\r',size(new_nodes,1));
for i=1:QNoV2
    fprintf(fis,'[%3.5f %3.5f %3.5f]\r',new_nodes(i,1),new_nodes(i,2), new_nodes(i,3));
end;

new_elem=new_elem-1;
%new_elem

QNoF2 = size(new_elem,1)
%fprintf(fis,'%4.0f\r',size(new_elem,1));

for i=1:QNoF2
    fprintf(fis,'[%d %d %d %d %d %d]\n',new_elem(i,1),new_elem(i,2),new_elem(i,3),new_elem(i,4),new_elem(i,5),new_elem(i,6));
end;

fclose (fis);


