function  plotmesh2d(Vert,Tri)

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
plot(Vert(:,1),Vert(:,2),'+');
hold off
