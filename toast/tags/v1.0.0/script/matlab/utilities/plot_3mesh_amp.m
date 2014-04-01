%plot_2_mesh_pha.m

% Picture both the brain and the head mesh 

sphere000;
bem3D;
cmin = min(log10(x_amp0));
plotmeshsol(n_n0, p, n_e0, node, log10(x_amp0), cmin, 0) ;
sphere001;
plotmeshsol(n_n1, p, n_e1, node, log10(x_amp1), cmin, 0) ;
sphere002;
plotmeshsol(n_n2, p, n_e2, node, log10(x_amp2), cmin, 0) ;

