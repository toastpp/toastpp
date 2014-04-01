function plot2_sol(elem,n_elem,nodes,n_nodes,x_sol,source,col)





% draw the solution on the surface.
% picture of 6 nodes triangles 

for i1=1:n_nodes% turning ampitude to log scale
    
    %  xyzmax=norm(nodes(i1,:));%max(abs(nodes(i1,:)));
    %  norm(nodes(i1,:))
    
    x(i1)=nodes(i1,1);
    y(i1)=nodes(i1,2);
    z(i1)=nodes(i1,3);
    
    theta(i1)= acos(z(i1)/norm(nodes(i1,:)));
    
    
    x_amp(i1)=abs(x_sol(i1));
  %  x_pha(i1)=angle(x_sol(i1));  
  x_pha(i1)=angle(abs(real(x_sol(i1)))+i*imag(x_sol(i1)));
    
    

end



for i2=1:n_nodes,
    for i3=1:n_nodes-1,
        if(theta(i3)>theta(i3+1))
            t_amp=x_amp(i3);
            x_amp(i3)=x_amp(i3+1);
            x_amp(i3+1)=t_amp;
            
            t_pha=x_pha(i3);
            x_pha(i3)=x_pha(i3+1);
            x_pha(i3+1)=t_pha;
            
            t_thet=theta(i3);
            theta(i3)=theta(i3+1);
            theta(i3+1)=t_thet;
        end
    end
end
            

figure(9);
hold on
plot(180*(theta/pi),log(x_amp),col);

figure(10);
hold on;
plot(180*(theta/pi),(x_pha),col);


