function [mua mus] = cylphantom(nx,ny,nz,mode)
    mua = ones(nx,ny,nz) * 0.01;
    mus = ones(nx,ny,nz) * 1;
    
    % geometry
    cnt1 = [0.5 0.5 0.3];
    rad1 = [0.3 0.15 0.3];
    phi1 = 0;
    theta1 = pi/3;
    
    cnt2 = [0.5 0.3 0.7];
    rad2 = 0.15;

    cnt3 = [0.35 0.7 0.75];
    rad3 = 0.07;

    cnt4 = [0.65 0.7 0.75];
    rad4 = 0.07;

    cnt5 = [0.7 0.75 0.2];
    rad5 = 0.05;

    cnt6 = [0.3 0.75 0.2];
    rad6 = 0.05;

    switch mode
        case 1 % 'truth' prior
            v1 = 0.02;
            mua = DrawEllipsoid(mua,cnt1,rad1,phi1,theta1,v1,nx,ny,nz);
            v2 = 2;
            mus = DrawSphere(mus,cnt2,rad2,v2,nx,ny,nz);
            v3 = 0.04;
            mua = DrawSphere(mua,cnt3,rad3,v3,nx,ny,nz);
            v4 = 4;
            mus = DrawSphere(mus,cnt4,rad4,v4,nx,ny,nz);
            v5 = 0.04;
            mua = DrawSphere(mua,cnt5,rad5,v5,nx,ny,nz);
            v6 = 4;
            mus = DrawSphere(mus,cnt6,rad6,v6,nx,ny,nz);
        case 2 % 'partial' prior
            v1 = 0.02;
            mua = DrawEllipsoid(mua,cnt1,rad1,phi1,theta1,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt2,rad2,v1,nx,ny,nz);
            v2 = 2;
            mus = DrawEllipsoid(mus,cnt1,rad1,phi1,theta1,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt2,rad2,v2,nx,ny,nz);
        case 3 % 'sum' prior
            v1 = 0.02;
            mua = DrawEllipsoid(mua,cnt1,rad1,phi1,theta1,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt2,rad2,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt3,rad3,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt4,rad4,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt5,rad5,v1,nx,ny,nz);
            mua = DrawSphere(mua,cnt6,rad6,v1,nx,ny,nz);
            v2 = 2;
            mus = DrawEllipsoid(mus,cnt1,rad1,phi1,theta1,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt2,rad2,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt3,rad3,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt4,rad4,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt5,rad5,v2,nx,ny,nz);
            mus = DrawSphere(mus,cnt6,rad6,v2,nx,ny,nz);
    end
    mua_pad = zeros(nx+2,ny+2,nz);
    mua_pad(2:nx+1,2:ny+1,:) = mua(:,:,:);
    mus_pad = zeros(nx+2,ny+2,nz);
    mus_pad(2:nx+1,2:ny+1,:) = mus(:,:,:);

    [f,v] = isosurface(mua_pad,0.015);
    p = patch('Faces',f,'Vertices',v);
    isonormals(mua_pad,p);
    set(p,'FaceColor','red','EdgeColor','none');
    
    [f,v] = isosurface(mus_pad,1.5);
    p = patch('Faces',f,'Vertices',v);
    isonormals(mus_pad,p);
    set(p,'FaceColor','blue','EdgeColor','none');
    
    cylsurf = zeros(nx,ny,nz);
    radi = [nx/2 ny/2];
    cnti = [(nx+1)/2 (ny+1)/2];
    for j = 1:ny
        dr(2) = j-cnti(2);
        for i = 1:nx
            dr(1) = i-cnti(1);
            if norm(dr./radi) <= 1
                cylsurf(i,j,:) = 1;
            end
        end
    end
    cylsurf_pad = zeros(nx+2,ny+2,nz);
    cylsurf_pad(2:nx+1,2:ny+1,:) = cylsurf(:,:,:);
    
    [f,v] = isosurface(cylsurf_pad,0.5);
    p = patch('Faces',f,'Vertices',v);
    isonormals(cylsurf_pad,p);
    set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.5);
        
    % show the locations of the cross sections
    %p = patch([0 nx nx 0],[0 0 ny ny],[nz*0.2 nz*0.2 nz*0.2 nz*0.2],[0.5 0.5 0.5 0.5]);
    %set(p,'FaceColor','black','EdgeColor','none','FaceAlpha',0.2);

    %p = patch([0 nx nx 0],[0 0 ny ny],[nz*0.75 nz*0.75 nz*0.75 nz*0.75],[0.5 0.5 0.5 0.5]);
    %set(p,'FaceColor','black','EdgeColor','none','FaceAlpha',0.2);

    %p = patch([0 nx nx 0],[ny*0.5 ny*0.5 ny*0.5 ny*0.5],[0 0 nz nz],[0.5 0.5 0.5 0.5]);
    %set(p,'FaceColor','black','EdgeColor','none','FaceAlpha',0.2);

    daspect([1 1 1])
    view(52,40);
    axis equal;axis tight
    camlight
    lighting gouraud
    xlabel('x');
    ylabel('y');
    zlabel('z');

end

function p = DrawSphere(pi,cnt,rad,v,nx,ny,nz)
    p = pi;
    cnti = [cnt(1)*nx cnt(2)*ny cnt(3)*nz];
    radi = rad*nx;
    for k = 1:nz
        dz = k-cnti(3);
        for j = 1:ny
            dy = j-cnti(2);
            for i = 1:nx
                dx = i-cnti(1);
                dst = sqrt(dx^2 + dy^2 + dz^2);
                if dst <= radi
                    p(i,j,k) = v;
                end
            end
        end
    end
end

function p = DrawEllipsoid(pi,cnt,rad,phi,theta,v,nx,ny,nz)
    p = pi;
    cnti = [cnt(1)*nx cnt(2)*ny cnt(3)*nz];
    radi = [rad(1)*nx rad(2)*ny rad(3)*nz];
    rot1 = [[cos(phi) -sin(phi) 0]; [sin(phi) cos(phi) 0]; [0 0 1]];
    rot2 = [[1 0 0]; [0 cos(theta) -sin(theta)]; [0 sin(theta) cos(theta)]];
    rot = rot1 * rot2;
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                s = [i-cnti(1) j-cnti(2) k-cnti(3)];
                st = rot * s';
                if norm (st./radi') <= 1
                    p(i,j,k) = v;
                end
            end
        end
    end
end
                