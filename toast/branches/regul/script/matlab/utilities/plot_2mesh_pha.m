%plot_2_mesh_pha.m

% Picture both the brain and the head mesh 

sphere000;
bem3D;
cmax = max(x_pha0);
plotmeshsol(n_n0, p, n_e0, node, x_pha0, 0, cmax) ;
sphere001;
plotmeshsol(n_n1, p, n_e1, node, x_pha1, 0, cmax) ;

