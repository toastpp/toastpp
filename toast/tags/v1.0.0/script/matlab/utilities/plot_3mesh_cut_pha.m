%plot_3mesh_cut_pha.m

% Picture both the brain and the head mesh 


bem3D;
cmin0 = min(x_pha0)
cmax0 = max(x_pha0)
cmax = 2;
sphere002;
plotmeshsol(n_n2, p, n_e2, node, x_pha2, 0, cmax) ;
sphere001;
cmin1=min(x_pha1)
cmin1=max(x_pha1)
hold on;
%plotmeshsolcut(n_n1, p, n_e1, node, x_pha1, 0, cmax,-1,0,0,0) ;
sphere000;
cmin2=min(x_pha2)
cmax2=max(x_pha2)
hold on;
%plotmeshsolcut(n_n0, p, n_e0, node, x_pha0, 0, cmax,-1,0,0,0) ;
hold off;
