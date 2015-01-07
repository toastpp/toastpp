function  plotregionmesh2d(Vert,Tri,nt,nr)
%
% plot
hold on
for k = 1:size(Tri,1)
   for j = 1:3
     x1 = Vert(Tri(k,j),1);
     y1 = Vert(Tri(k,j),2);
     x2 = Vert(Tri(k,mod(j,3)+1),1);
     y2 = Vert(Tri(k,mod(j,3)+1),2);
     plot([x1,x2],[y1,y2],'g');
  end
end;
NoV = size(Vert,1);
% how many regions
runiq = unique(nr);
NoR = size(runiq,1);
syms = {'+b','+r','+c','+k'};
for k = 0:NoR-1
   plot(Vert(find(nr==k),1),Vert(find(nr==k),2),syms{k+1});
end
% boundary nodes
%plot(Vert(find(nt=='B'),1),Vert(find(nt=='B'),2),'+m');
hold off
