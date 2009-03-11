function p = toastLabelRegion
%toastLabelRegion     - Assign node labels within a volume defined by a surface.
%
% Synopsis: L = toastLabelRegion(hMesh,vtx,idx,val)
%    hMesh: mesh handle
%    vtx:   matrix of n surface vertices (real n x 3)
%    idx:   index matrix of m surface triangles (int m x 3)
%    val:   label value for internal nodes (int)
%    L:     array of new mesh labels
%
% Given a mesh and a closed surface composed of triangles, this function
% assigns a region label to each mesh node inside the surface. Labels of
% external nodes remain unchanged.
%
% The surface is defined by vertex coordinates (vtx) and a list of triangles,
% each of which is represented by a triplet of node indices (1-based).
% The surface must be closed.
%
% When calling toastLabelRegion repeatedly for concentric shells, start with
% the outermost one and work your way inside, since the labels in the external
% region remain untouched.
%
% The list of modified node labels is returned by the function, but the mesh
% associated with hMesh is also directly updated, so writing out the mesh
% with toastWriteMesh will contain the new label settings.
