close all;clear all;
load neonate_vol_srf_46;
p=vtx;
t=tri;
clear vtx; clear tri;
%simpplot(p,t,'p(:,1)>-10');
brcsf=union(br,csf);
sk;
sc;

face=setdiff(srf_sc,srf_sk,'rows'); % note that there are a repetition in Lior format 
[b,m,n] = unique(face);


% read the fem file
file='neonatebrain_opt_100_fa.fem';
fid=fopen(file,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    valfem = str2num(tline);
    disp(valfem)
end
fclose(fid);


cdata= 0 * p;
cdata(b,1)=valfem;

figure
trimesh(face,p(:,1),p(:,2),p(:,3),'CData',cdata(:,1),'facecolor','interp','Edgecolor','none','FaceLighting','phong','FaceAlpha',1);
axis off;
daspect([1 1 1]);
material dull
axis([min(p(face,1)) max(p(face,1)) min(p(face,2)) max(p(face,2)) min(p(face,3)) max(p(face,3))])
view([50 10])
drawnow

