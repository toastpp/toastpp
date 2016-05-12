function h = toastShowMesh (mesh, vtxdata, varargin)
%toastShowMesh        - Render a mesh surface in 3-D view.
%
% Synopsis: toastShowMesh (mesh)
%           toastShowMesh (mesh, vtxdata)
%           toastShowMesh (mesh, vtxdata, [property, value, ...])
%    mesh:        mesh object
%    vtxdata:     optional nodal parameter array
%    [dmin dmax]: data range
%
% Displays a 2-D mesh or the outer surface of a 3-D toast mesh. If the
% nodal parameter array is provided, its surface values are used to colour
% the mesh surface.

dim = mesh.Dimension;

[vtx,idx,eltp] = mesh.Data;
eltps = unique(eltp);

if dim == 3
    % Test for surface mesh
    is_surfmesh = true;
    eltps = unique(eltp);
    for i=1:length(eltps)
        if eltps(i) ~= 16 && eltps(i) ~= 17
            is_surfmesh = false;
            break;
        end
    end
    if ~is_surfmesh
        [vtx,idx,perm] = mesh.SurfData;
    else
        perm = [1:size(vtx,1)];
    end
end

nel = size(idx,1);
nnd = size(vtx,1);
showgrid = true;
scnlighting = true;
showsurfdata = false;
showcolorbar = false;

if nargin > 1
    showcolorbar = true;
    vtxdata = double(vtxdata);
    showsurfdata = true;
    showgrid = false;
    if dim == 3
        surfdata = vtxdata(perm);
    else
        surfdata = vtxdata;
    end
    smin = min(surfdata);
    smax = max(surfdata);
    optarg = 1;
    while nargin-2 >= optarg
        label = varargin{optarg};
        if ischar(label) && nargin-1 >= optarg
            value = varargin{optarg+1};
            optarg = optarg+2;
            if strcmpi(label,'range')
                smin = value(1);
                smax = value(2);
            elseif strcmpi(label,'showgrid')
                showgrid = value;
            elseif strcmpi(label,'colorbar')
                showcolorbar = value;
            elseif strcmpi(label,'lighting')
                scnlighting = value;
            end
        end
    end
    
    
%     if (nargin > 2)
%         range = sort(range);
%         smin = range(1);
%         smax = range(2);
%     else
%         smin = min(surfdata);
%         smax = max(surfdata);
%     end
else
    surfdata = ones(nnd,1);
end

% Remove non-vertex nodes for display
if length(eltps) == 1   % can only fix uniform meshes
    ieltp = eltp(1);
    if ieltp == 6  % 6-noded triangle
        idx = idx(:,1:3);
    elseif ieltp == 9  % 10-noded triangle
        idx = idx(:,1:3);
    elseif ieltp == 7  % 10-noded tetrahedron
        idx = idx(:,1:4);
    elseif ieltp == 17 % 6-noded surface triangle
        idx = idx(:,1:3);
    end
end

if showgrid
    edgecol = 'black';
else
    edgecol = 'none';
end

if dim==3

    if showsurfdata == true
        h = patch('Vertices', vtx, 'Faces', idx, 'CData',surfdata);
        view([1 1 1]);
        %set(gcf,'Renderer','zbuffer')
        if scnlighting
            shading interp
            lightangle(45,30)
            set(h, ...
                'FaceLighting','phong', ...
                'AmbientStrength',.2,'DiffuseStrength',.8,...
                'SpecularStrength',.0,'SpecularExponent',1,...
                'BackFaceLighting','unlit');
        else
            set(h, 'FaceColor','interp');
            set(findobj(gca,'type','surface'), ...
                'FaceLighting','flat');
        end
    else
        h = patch('Vertices', vtx, 'Faces', idx, 'FaceColor', 'white', 'EdgeColor', 'black');
        view([1 1 1]);
        set(gcf,'Renderer','zbuffer')
    end

else

    if showsurfdata && size(surfdata,1) == size(vtx,1)
        h = patch('Vertices', vtx, 'Faces', idx, 'CData',surfdata, 'FaceColor','interp');
    else
        h = patch('Vertices', vtx, 'Faces', idx, 'FaceVertexCData', surfdata, 'FaceColor', 'flat');
    end
    view([0 0 1]);

end
if showsurfdata && smax > smin
    set(gca,'CLim',[smin,smax]);
end
%if showgrid
    set(h,'EdgeColor',edgecol);
%end
if showcolorbar
    colorbar
end

axis equal
axis tight