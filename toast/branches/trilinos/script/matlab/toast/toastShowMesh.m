function p = toastShowMesh (hMesh, vtxdata, range)
%toastShowMesh        - Render a mesh surface in 3-D view.
%
% Synopsis: toastShowMesh (hMesh)
%           toastShowMesh (hMesh, vtxdata)
%           toastShowMesh (hMesh, vtxdata, [dmin dmax])
%    hMesh:   mesh handle
%    vtxdata: optional nodal parameter array
%    [dmin dmax]: data range
%
% Displays a 2-D mesh or the outer surface of a 3-D toast mesh. If the
% nodal parameter array is provided, its surface values are used to colour
% the mesh surface.

dim = toastMeshDimension(hMesh);

if dim == 3
    [vtx,idx,perm] = toastSurfData (hMesh);
    if size(vtx,1)==0 % try as surface mesh
        [vtx,idx,eltp] = toastMeshData(hMesh);
        eltp = unique(eltp);
        for i=1:length(eltp)
            if eltp(i) ~= 16 && eltp(i) ~= 17
                error('Invalid mesh');
            end
        end
        perm = [1:size(vtx,1)];
    end
else
    [vtx,idx] = toastMeshData (hMesh);
end

nel = size(idx,1);
nnd = size(vtx,1);

if nargin > 1
    showsurfdata = true;
    if dim == 3
        surfdata = vtxdata(perm);
    else
        surfdata = vtxdata;
    end
    if (nargin > 2)
        smin = range(1);
        smax = range(2);
    else
        smin = min(surfdata);
        smax = max(surfdata);
    end
    if smin < smax
        for i=1:length(surfdata)
            if surfdata(i) > smax
	        surfdata(i) = smax;
	    elseif surfdata(i) < smin
                surfdata(i) = smin;
            end
        end
        surfdata = (surfdata-smin)/(smax-smin)*255;
    else
        surfdata = ones(nnd,1);
    end
else
    showsurfdata = false;
    surfdata = ones(nnd,1);
end

if dim==3

    if showsurfdata == true
        patch('Vertices', vtx, 'Faces', idx, 'FaceLighting', 'Phong','CData',surfdata);
        view([1 1 1]);
        shading interp
        lightangle(45,30)
        set(gcf,'Renderer','zbuffer')
        set(findobj(gca,'type','surface'),...
            'FaceLighting','phong',...
            'AmbientStrength',.2,'DiffuseStrength',.8,...
            'SpecularStrength',.3,'SpecularExponent',5,...
            'BackFaceLighting','unlit')
    else
        patch('Vertices', vtx, 'Faces', idx, 'FaceColor', 'white', 'EdgeColor', 'black');
        view([1 1 1]);
        set(gcf,'Renderer','zbuffer')
    end

else

    if size(surfdata,1) == size(vtx,1)
        patch('Vertices', vtx, 'Faces', idx, 'CData',surfdata, 'EdgeColor',[0 0 0],'FaceColor','interp');
    else
        patch('Vertices', vtx, 'Faces', idx, 'FaceVertexCData', surfdata, 'FaceColor', 'flat');
    end
    view([0 0 1]);

end


axis equal
axis tight
