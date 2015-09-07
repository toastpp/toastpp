% Convert a Gmsh mesh file to a toastMesh object
% Currently this works only for tetrahedral meshes

function mesh = toastReadGmshMesh (fname)

% Read the file, using the Gmsh-provided Matlab converter
gmsh = load_gmsh(fname);
eltp = 0;

% Extract vertices
vtx = gmsh.POS;

% Look for tetrahedral elements
%tet_idx = find(gmsh.TETS(:,5));
els = gmsh.TETS(1:gmsh.nbTets,:);
if (~isempty(els))
    idx = els(:,1:4);
    eltp = 3;
end

% Look for triangular elements
if eltp == 0
    els = gmsh.TRIANGLES(1:gmsh.nbTriangles,:);
    if (~isempty(els))
        idx = els(:,1:3);
        eltp = 15;
    end
end

if eltp == 0
    error('Gmsh element type not supported by toast importer');
end

% Save region labels
elreg = els(:,end);


% Convert to toast mesh
nel = size(idx,1);
eltp = ones(nel,1)*eltp;
mesh = toastMesh(vtx, idx, eltp);

% Set region labels
mesh.Region(elreg);

end