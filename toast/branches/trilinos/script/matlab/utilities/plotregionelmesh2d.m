function  plotregionelmesh2d(Vert,Tri,nt,er)
%
% plot
syms = {'*b','*r','*c','*k'};
hold on
for k = 1:size(Tri,1)
   xmid = 0.0;ymid = 0.0;
   for j = 1:3
     x1 = Vert(Tri(k,j),1);
     y1 = Vert(Tri(k,j),2);
     x2 = Vert(Tri(k,mod(j,3)+1),1);
     y2 = Vert(Tri(k,mod(j,3)+1),2);
     plot([x1,x2],[y1,y2],'g');
     xmid = xmid + x1/3.0;
     ymid = ymid + y1/3.0;
  end
  plot(xmid,ymid,syms{er(k)+1});   %region label
end;

for k= 1:size(Vert,1)
    if(nt{k}=='I')
        plot(Vert(k,1),Vert(k,2),'+m');
    end
end
hold off
