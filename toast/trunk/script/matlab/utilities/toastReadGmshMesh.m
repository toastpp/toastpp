% Convert a Gmsh mesh file to a toastMesh object
% Currently this works only for tetrahedral meshes

function mesh = toastReadGmshMesh (fname)

% Read the file, using the Gmsh-provided Matlab converter
gmsh = load_gmsh(fname);

% Extract vertices
vtx = gmsh.POS;

% Look only for tetrahedral elements
%tet_idx = find(gmsh.TETS(:,5));
tets = gmsh.TETS(1:gmsh.nbTets,:);

% Save region labels
elreg = tets(:,5);
tets = tets(:,1:4);

% Convert to toast mesh
nel = size(tets,1);
eltp = ones(nel,1)*3; % "4-noded tets"
mesh = toastMesh(vtx, tets, eltp);

% Set region labels
mesh.Region(elreg);

end