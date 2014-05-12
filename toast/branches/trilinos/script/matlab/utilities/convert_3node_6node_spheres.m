
% draw the solution on the head.

% Picture all the meshes Created from the 3D BEM inverse
% Created by Athanasios Feb 2003

clear all;

% the mesh files 
fid = fopen('sphere20.dat') ;% opens the file for reading: external mesh
exitname= 'sphere20_6n.dat';
%external mesh
NoV=fscanf(fid,'%i',1);

x0=zeros(NoV,1);
y0=zeros(NoV,1);
z0=zeros(NoV,1);

%Reading from the file 
% x,y,z
for j=1:(NoV),
    fscanf(fid,'%c',[4]);
    
    a=fscanf(fid,'%g',[3]);
    x0(j,1)=a(1);
    y0(j,1)=a(2);
    z0(j,1)=a(3);   
    fscanf(fid,'%c',[1]);
    
    % fscanf(fid,'%c',[2]);
end;

NoE=fscanf(fid,'%i',1);

node=zeros(NoE,3);

for j=1:(NoE),
    fscanf(fid,'%c',[4]);
    a=fscanf(fid,'%g',[3]);
    node(j,1:3)=a(1:3)';         
    fscanf(fid,'%c',[1]);
end;
fclose(fid); % closes the file

min(min(node))
%figure;

plot3(x0,y0,z0,'g.');
hold on;

Vertex(:,1)=x0;
Vertex(:,2)=y0;
Vertex(:,3)=z0;

% picture of 6 nodes triangles

% Plot the numbers
% ii1= int2str(0:NoV);
%for h=1:101,
%    text(x0(h),y0(h),z0(h),int2str(h-1));
%    hold on;
%end;    

d=1;
for iq=1:NoE%drawing the triangles
    
    ie=node(iq,1);%+1;
    je=node(iq,2);%+1;
    ke=node(iq,3);%+1;
    % le=node(i,4)+1;
    % me=node(i,5)+1;
    % ne=node(i,6)+1;
    
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

    px=[xi,xij,xj,xjk,xk,xki,xi];%,xl,xm,xn,xi];
    py=[yi,yij,yj,yjk,yk,yki,yi];%,yl,ym,yn,yi];
    pz=[zi,zij,zj,zjk,zk,zki,zi];%,zl,zm,zn,zi];
    
%     px=[xij,xjk,xki,xij];%,xl,xm,xn,xi];
%     py=[yij,yjk,yki,yij];%,yl,ym,yn,yi];
%     pz=[zij,zjk,zki,zij];%,zl,zm,zn,zi];
    
    plot3(xij,yij,zij,'r*');
    plot3(xjk,yjk,zjk,'r*');
    plot3(xki,yki,zki,'r*');
    
    fill3(px,py,pz,[1 0 1]);
    %    plot3(px,py,pz);                   
    hold on
    
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
idx
%clear new_nodes;
new_elem=[node(:,1) zeros(NoE,1) node(:,2) zeros(NoE,1) node(:,3) zeros(NoE,1)];



for is=1:size(add_nodes,1)
    for it=1:NoE,
        if (new_elem(it,1)==add_nodes(is,4)) & (new_elem(it,3)==add_nodes(is,5)) & (new_elem(it,2)==0)
            new_elem(it,2)=is+NoV;
        end
        if (new_elem(it,3)==add_nodes(is,4)) &(new_elem(it,5)==add_nodes(is,5)) &(new_elem(it,4)==0)
            new_elem(it,4)=is+NoV;
        end
        if (new_elem(it,5)==add_nodes(is,4)) & (new_elem(it,1)==add_nodes(is,5)) & (new_elem(it,6)==0)
            new_elem(it,6)=is+NoV;
        end
        if (new_elem(it,1)==add_nodes(is,5)) & (new_elem(it,3)==add_nodes(is,4)) & (new_elem(it,2)==0)
            new_elem(it,2)=is+NoV;
        end
        if (new_elem(it,3)==add_nodes(is,5)) &(new_elem(it,5)==add_nodes(is,4)) &(new_elem(it,4)==0)
            new_elem(it,4)=is+NoV;
        end
        if (new_elem(it,5)==add_nodes(is,5)) & (new_elem(it,1)==add_nodes(is,4)) & (new_elem(it,6)==0)
            new_elem(it,6)=is+NoV;
        end
    end
end
min(min(new_elem))

 new_nodes=[Vertex ; [add_nodes(:,1) add_nodes(:,2) add_nodes(:,3)]]

fis = fopen(exitname,'w');


fprintf(fis,'%4.0f\r',size(new_nodes,1));
for i=1:size(new_nodes,1),
    fprintf(fis,'[%3.5f %3.5f %3.5f]\r',new_nodes(i,1),new_nodes(i,2), new_nodes(i,3));
     fprintf(1,'[%3.5f %3.5f %3.5f]\r',new_nodes(i,1),new_nodes(i,2), new_nodes(i,3));

end;


new_elem=new_elem-1;

fprintf(fis,'%4.0f\r',size(new_elem,1));

for i=1:size(new_elem,1),
    fprintf(fis,'[%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f]\r',new_elem(i,1),new_elem(i,2),new_elem(i,3),new_elem(i,4),new_elem(i,5),new_elem(i,6));
    fprintf(1,'[%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f]\r',new_elem(i,1),new_elem(i,2),new_elem(i,3),new_elem(i,4),new_elem(i,5),new_elem(i,6));
   
end;

fclose (fis);