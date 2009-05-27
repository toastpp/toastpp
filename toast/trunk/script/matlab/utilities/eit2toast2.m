close all;clear all;
load neonate_vol_srf_46;
p=vtx;
t=tri;
clear vtx; clear tri;
%simpplot(p,t,'p(:,1)>-10');
brcsf=union(br,csf);
sk;
sc;

NoV=length(p);
NoTet=length(t);
mua=0.01*ones(NoV,1);
mus=ones(NoV,1);
N=1.4 * ones(NoV,1);
region=zeros(NoV,1);
% external boundary
for i=1:NoV
    bnd(i)='N';
end
bnd(srf_sc(:))='B';
bnd(srf_sk(:))='N';
% regions and optical parameters
clear a; a=(t(sc,:));
mus(a(:))=0.8;mua(a(:))=0.0149;region(a(:))=0;
clear a; a=(t(sk,:));
mus(a(:))=1.0;mua(a(:))=0.01;region(a(:))=1;
clear a; a=(t(brcsf,:));
mus(a(:))=1.25;mua(a(:))=0.0178;region(a(:))=2;
clear a;



FileNameOut=['neonatebrain','.tmsh'];
fid1 = fopen(FileNameOut,'w');
fprintf(fid1,'MeshData 5.0\n\n');
fprintf(fid1,'NodeList %d 1\n',NoV);
for i=1:NoV
    fprintf(fid1,'%s[%g %g %g]R%d\n',char(bnd(i)),p(i,1),p(i,2),p(i,3),region(i));
end
fprintf(fid1,'\nElementList %d\n',NoTet);
for i=1:NoTet
    p1=t(i,1);
    p2=t(i,2);
    p3=t(i,3);
    p4=t(i,4);

    px=cross(p(p2,:)-p(p1,:),p(p3,:)-p(p2,:));
    pc=((p(p1,:)+p(p2,:)+p(p3,:))/3);
    cc=dot(px,p(p4,:)-pc);

    if(cc >= 0)
        fprintf(fid1,'c %d %d %d %d\n',t(i,1),t(i,2),t(i,3),t(i,4));
    else
        fprintf(fid1,'c %d %d %d %d\n',t(i,1),t(i,3),t(i,2),t(i,4));
    end
end
fprintf(fid1,'\n[ParameterList]\n');
fprintf(fid1,'Size %d\n',NoV);
fprintf(fid1,'Param1 MUA\n');
fprintf(fid1,'Param2 MUS\n');
fprintf(fid1,'Param3 N\n');
fprintf(fid1,'Data\n');
for i=1:NoV
    fprintf(fid1,'%g %g %g\n',mua(i),mus(i),N(i));
end
fclose(fid1);

%%%%%%% qm file for one source
ind=find(bnd=='B');
nind=length(ind);
FileNameqm=['neonatebrain','.qm'];
fid2 = fopen(FileNameqm,'w');
fprintf(fid2,'QM file 3D\n');
fprintf(fid2,'Dimension 3\n\n');
fprintf(fid2,'SourceList 1\n');
fprintf(fid2,'%g %g %g\n\n',p(ind(4),1),p(ind(4),2),p(ind(4),3));
fprintf(fid2,'MeasurementList %d\n',nind);
for i=1:nind
fprintf(fid2,'%g %g %g\n',p(ind(i),1),p(ind(i),2),p(ind(i),3));
end
fprintf(fid2,'\nLinkList\n');
fprintf(fid2,'%d: ',nind);
for i=1:nind
fprintf(fid2,'%d ',i-1);
end
fprintf(fid2,'\n');
