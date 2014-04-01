function [img_mua,img_mus] = make_cyl_targets(grd,display)

mua_bkg = 0.01;
mus_bkg = 1;

img_mua = ones(grd) .* mua_bkg;
img_mus = ones(grd) .* mus_bkg;

%% Create inclusions

ellipsoid.r = [16,16,7];
ellipsoid.p = [0,0,-13];
ellipsoid.phi = -30*pi/180;
ellipsoid.mua = 0.02;
ellipsoid.mus = mus_bkg;

sphere1.r = 11;
sphere1.p = [0,5,10];
sphere1.mua = mua_bkg;
sphere1.mus = 2;

sphere2.r = 5;
sphere2.p = [-8,-13,15];
sphere2.mua = 0.02;
sphere2.mus = 0.5;

sphere3.r = 5;
sphere3.p = [8,-13,15];
sphere3.mua = 0.005;
sphere3.mus = 2;

sphere4.r = 3;
sphere4.p = [0,17,-6];
sphere4.mua = 0.04;
sphere4.mus = 0.25;

sphere5.r = 3;
sphere5.p = [0,17,-14];
sphere5.mua = 0.0025;
sphere5.mus = 4;

sphere6.r = 4;
sphere6.p = [0,8,14];
sphere6.mua = 0.02;
sphere6.mus = mus_bkg;

for k=1:grd(3)
    z = ((k-0.5)/grd(3)-0.5)*50;
    for j=1:grd(2)
        y = ((j-0.5)/grd(2)-0.5)*50;
        for i=1:grd(1)
            x = ((i-0.5)/grd(1)-0.5)*50;
            
            % create mua ellipsoid object
            yt = y*cos(ellipsoid.phi) + z*sin(ellipsoid.phi);
            zt = -y*sin(ellipsoid.phi) + z*cos(ellipsoid.phi);
            dr = [x,yt,zt]-ellipsoid.p;
            if (dr(1)/ellipsoid.r(1))^2 + (dr(2)/ellipsoid.r(2))^2 + (dr(3)/ellipsoid.r(3))^2 < 1
                img_mua(i,j,k) = ellipsoid.mua;
                img_mus(i,j,k) = ellipsoid.mus;
            end
            
            dr = [x,y,z]-sphere1.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere1.r^2
                img_mua(i,j,k) = sphere1.mua;
                img_mus(i,j,k) = sphere1.mus;
            end
            
            dr = [x,y,z]-sphere2.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere2.r^2
                img_mua(i,j,k) = sphere2.mua;
                img_mus(i,j,k) = sphere2.mus;
            end
            
            dr = [x,y,z]-sphere3.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere3.r^2
                img_mua(i,j,k) = sphere3.mua;
                img_mus(i,j,k) = sphere3.mus;
            end
            
            dr = [x,y,z]-sphere4.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere4.r^2
                img_mua(i,j,k) = sphere4.mua;
                img_mus(i,j,k) = sphere4.mus;
            end
            
            dr = [x,y,z]-sphere5.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere5.r^2
                img_mua(i,j,k) = sphere5.mua;
                img_mus(i,j,k) = sphere5.mus;
            end
            
            dr = [x,y,z]-sphere6.p;
            if dr(1)^2 + dr(2)^2 + dr(3)^2 < sphere6.r^2
                img_mua(i,j,k) = sphere6.mua;
                img_mus(i,j,k) = sphere6.mus;
            end
            
        end
    end
end

%% display isosurfaces

if nargin > 1 && display == true
    
    img_mua = smooth3(img_mua);
    img_mus = smooth3(img_mus);
    
    cyl = zeros(grd+4);
    for k=1:size(cyl,3)
        z = (k-(size(cyl,3)+1)/2) * 50/grd(3);
        for j=1:size(cyl,2)
            y = (j-(size(cyl,2)+1)/2) * 50/grd(2);
            for i=1:size(cyl,1)
                x = (i-(size(cyl,1)+1)/2) * 50/grd(1);
                if x^2 + y^2 < 25^2
                    cyl(i,j,k) = 1;
                end
            end
        end
    end
    figure;
    %set(gcf,'Renderer','zbuffer');
    
    p = patch(isosurface(cyl,0.5));
    isonormals(cyl,p);
    set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.3);
    daspect([1,1,1])
    view(3); axis tight
    camlight
    lighting gouraud
    
    p = patch(isosurface(img_mua,0.015));
    isonormals(img_mua,p);
    set(p,'FaceColor','red','EdgeColor','none');
    daspect([1,1,1])
    
    hold on
    p = patch(isosurface(img_mus,1.5));
    isonormals(img_mus,p);
    set(p,'FaceColor','blue','EdgeColor','none');
end

