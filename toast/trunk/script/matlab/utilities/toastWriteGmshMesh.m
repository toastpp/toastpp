% Write a Toast Mesh in Gmsh format

function toastWriteGmshMesh (mesh, fname, prm)

if nargin > 2 && strcmpi(prm,'surf')
    [vtx,idx,eltp] = mesh.SurfData;
    gmsh.d = 2;
else
    [vtx,idx,eltp] = mesh.Data;
    gmsh.d = size(vtx,2);
end
gmsh.n = size(vtx,2);
gmsh.nNodes = size(vtx,1);
gmsh.nodes = vtx;
gmsh.nElems = size(idx,1);
gmsh.elems = idx;
gmshMeshWriter(gmsh,fname);

end