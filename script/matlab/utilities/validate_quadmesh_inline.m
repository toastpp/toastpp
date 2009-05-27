function validate_quadmesh_inline( NoV,Vertex,NoF,Tri)

% validate a mesh by seeing if any two polygons have an edge 
% in common, oriented in the same direction.

x0=Vertex(:,1);
y0=Vertex(:,2);
z0=Vertex(:,3);

NoV
NoF
NoE = 3*NoF
NoV + NoF
NoE + 2
edges =zeros(NoE,4);

d=1;
for iq=1:NoF  % find all edges
    
    ie=Tri(iq,1)+1;
    je=Tri(iq,3)+1;
    ke=Tri(iq,5)+1;
%    le=Tri(iq,4)+1;
%    me=Tri(iq,5)+1;
%    ne=Tri(iq,6)+1;
    
%[ie je ke le me ne];
    edges(d,:)=[ ie je iq -1];
    d=d+1;
    edges(d,:)=[ je ke iq -1];
    d=d+1;
    edges(d,:)=[ ke ie iq -1];
    d=d+1;

end

pairidx=1;
flipidx=1;
nopair = 0;
for is=1:NoE
    if(edges(is,4) > -1) % paired edge already found
        continue;        % skip to next
    end
    nopair = nopair +1;
    for isx=is+1:NoE
        if ((edges(is,1)==edges(isx,2)) & (edges(is,2)==edges(isx,1)))
            pairidx=pairidx+1;
            nopair = nopair -1;
            edges(is,4) = isx; % mark corresponding edges
            edges(isx,4) = is;
	    break;
	end            
        if ((edges(is,1)==edges(isx,1)) & (edges(is,2)==edges(isx,2)))
            pairidx=pairidx+1;
            flip(flipidx) = edges(is,3);
            flipidx=flipidx+1;
            nopair = nopair-1;
            edges(is,4) = isx; % mark corresponding edges
            edges(isx,4) = is;
	    break;
        end
    end
end
pairidx
flipidx
nopair
if (flipidx > 1)
	flip
end

% visualise

vtx=Vertex;
srf=Tri(:,1:2:6)+1;


fig_title='Voila !';
f1=figure('Name',fig_title);

ta=vtx(srf(:,1),:);
tb=vtx(srf(:,2),:);
tc=vtx(srf(:,3),:);
cx=cross(tc-tb,tb-ta);
for j = 1:length(cx(:,1))
    cxlon(j)=norm(cx(j,:));
    cx(j,:) = 5*cx(j,:)/cxlon(j);
end
%    size(cx);
%    size(cxlon);

cc=(ta+tb+tc)/3;
quiver3(cc(:,1),cc(:,2),cc(:,3),cx(:,1),cx(:,2),cx(:,3),2);
hold on

%trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',1);
colormap([0 0 0]);
daspect([1 1 1]);

hold on;
axis image;
col = rand(1,3)/2 + [ 0.2 0.2 0.2];
set(gcf,'Colormap',col);
grid off
view(3);
axis off
