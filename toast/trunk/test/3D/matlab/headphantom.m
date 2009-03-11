function [mua mus] = headphantom(hMesh,nx,ny,nz,mode)
    mua = ones(nx,ny,nz) * 0.01;
    mus = ones(nx,ny,nz) * 1;
    
    % geometry
    cnt1 = [0.4 0.7 0.4];
    rad1 = 0.15;
    
    cnt2 = [0.7 0.8 0.6];
    rad2 = 0.12;

    switch mode
        case 1 % 'truth' prior
            v1 = 0.02;
            mua = DrawSphere(mua,cnt1,rad1,v1,nx,ny,nz);
            v2 = 2;
            mus = DrawSphere(mus,cnt2,rad2,v2,nx,ny,nz);
    end

    [vtx,idx,perm] = toastSurfData(hMesh);
    [pmin pmax] = toastMeshBB(hMesh);
    msize=pmax-pmin;
    scale = norm(msize)/100;
    nvtx = size(vtx,1);
    col = zeros(nvtx,3); col(:,2)=1;
    patch('Vertices',vtx,'Faces',idx,'FaceVertexCData',col,'FaceColor', ...
            'interp','EdgeColor','none','FaceAlpha',0.3);

    [f,v] = isosurface(mua,0.015);
    v = [v(:,2) v(:,1) v(:,3)];
    for i=1:size(v,1)
        v(i,1) = v(i,1)/nx*(pmax(1)-pmin(1)) + pmin(1);
        v(i,2) = v(i,2)/ny*(pmax(2)-pmin(2)) + pmin(2);
        v(i,3) = v(i,3)/nz*(pmax(3)-pmin(3)) + pmin(3);
    end
    p = patch('Faces',f,'Vertices',v);
    isonormals(mua,p);
    set(p,'FaceColor','red','EdgeColor','none');
    
    [f,v] = isosurface(mus,1.5);
    v = [v(:,2) v(:,1) v(:,3)];
    for i=1:size(v,1)
        v(i,1) = v(i,1)/nx*(pmax(1)-pmin(1)) + pmin(1);
        v(i,2) = v(i,2)/ny*(pmax(2)-pmin(2)) + pmin(2);
        v(i,3) = v(i,3)/nz*(pmax(3)-pmin(3)) + pmin(3);
    end
    p = patch('Faces',f,'Vertices',v);
    isonormals(mus,p);
    set(p,'FaceColor','blue','EdgeColor','none');
    
    axis equal;

    daspect([1 1 1])
    view(-18,54);
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
    radi = [rad*nx rad*nx rad*nx];
    for k = 1:nz
        dz = k-cnti(3);
        for j = 1:ny
            dy = j-cnti(2);
            for i = 1:nx
                dx = i-cnti(1);
                dst = dx^2/radi(1)^2 + dy^2/radi(2)^2 + dz^2/radi(3)^2;
                if dst <= 1
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
                