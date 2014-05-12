function  plotsolemesh2d(Vert,Tri,sol,eflag,vflag)

cols = cell(6,1);
cols = {'m','c','r','g','b','y'};
syms = cell(3,1);
syms = {'N','B','I'};
% plot
hold on
mins = min(sol);
maxs = max(sol);
%sols = (sol - mins)./(maxs - mins);
sols = sol;
for k = 1:size(Tri,1)
   for j = 1:3
     x1 = Vert(Tri(k,j),1);
     y1 = Vert(Tri(k,j),2);
     x2 = Vert(Tri(k,mod(j,3)+1),1);
     y2 = Vert(Tri(k,mod(j,3)+1),2);
     if eflag
        plot([x1,x2],[y1,y2],'g');
     end
   end
   ind = Tri(k,:); % [Tri(k,1) Tri(k,2) Tri(k,3)]; 
   V = Vert(ind,:);
   FX = V(:,1); %[Vert(Tri(k,1),1) Vert(Tri(k,2),1) Vert(Tri(k,3),1) ]  ;
   FY = V(:,2); %[Vert(Tri(k,1),2) Vert(Tri(k,2),2) Vert(Tri(k,3),2) ]  ;
%   pal = [sols(Tri(k,1)) sols(Tri(k,3)) sols(Tri(k,3)) ];
    pal = sols(ind);
%   fill(FX,FY,pal);
   patch(FX,FY,pal);
end;
for j=1:size(Vert,1)
%    opt = [cols{6-nr(j)} '+'];
%    if nt(j)=='B'
%        opt(2) = '.';
%    else if nt(j)=='I'
%            opt(2) = '*';
%        end
%    end
if vflag
    plot(Vert(:,1),Vert(:,2),'r+');
end%[ntype,Vert,nreg,Tri,nereg, mua,diff,rind,ptype]=readtoastemesh2d(meshname);
end
hold off
