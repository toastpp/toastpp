function mesh = gmshMeshReader(gmshFile)
% GMSHMESHREADER Reads mesh definition from the gmshFile file and saves is as a mesh structre
%   mesh = gmshMeshReader(gmshFile)
%   Reads mesh definition from the 'gmshFile' and saves is as a 'mesh' structre. The mesh from Gmsh/toast etc will have the refinement level p =0.
%   Supported Gmsh elements types are: 1, 2, 4 which are (n+1) -simplicies for dimantions n=1,2,3, respectively.
%
%   Fields of the mesh structure:
%    single level mesh fields:
%     n: space dimenstion (1,2,3)
%     p: refinement level
%     nNodes: number of nodes in the mesh up to level p (if mesh if refined this becomes an array)
%     nodes: coordinates of the nodes, node's number is its row number in the array
%     pNodes: list of the ids of the parents of the node (0 for the nodes of the original, p=0, mesh)
%     nElems: number of elements in the mesh up to level p (if mesh if refined this becomes an array)
%     elems: list of nodes in the element, counterclockwise numering, the element's number is its row number in the array
%     oSimps: list of the refined simplicies - the elements of the mesh at the refinement level 0, 
%             including all the nodes (including repetitions) that have been added at refinement levels.
%     nOSimps: no of the original simplicies (elements of the mesh at refinement level=0)
%    additional fields added after refinement:
%     nTotNodesSimp: total number of nodes in one simplex at level p 
%     nNewNodesSimp: no of new nodes in one simplex at level p
%     nTotElemSimp: total number of elements in one simplex at level p 
%     nNewElemSimp: no of new elements in one simplex at level p
%
%   Elements types are (n+1)-simplicies i.e. for n: 
%   n=1: line segment
%   n=2: triangle
%   n=3: tetraheder 

%Open Gmsh file
fid = fopen(gmshFile, 'r');

%Set the refinement level p = 0
p = 0;

%Scrap all tags before $Nodes
tag = '';
while ~strcmp(tag, '$Nodes');
  tag = fgetl(fid);
end

%Read number of nodes
nNodes = str2num(fgetl(fid));

%Read nodes, nodes in Gmsh are always 3D, if the mesh is lower dimenstion one or more coordinates will be equal
nodes = fscanf(fid, '%u %f %f %f\n', [4, nNodes]); nodes = nodes';
nodes = nodes(:, 2:end);
tag = fgetl(fid);

if ~strcmp(tag, '$EndNodes')
  error('gmshMeshReader:Problem with reading in the nodes. Inspect the gmsh file.');
end

%Scrap all tags before $Elements
tag = '';
while ~strcmp(tag, '$Elements');
  tag = fgetl(fid);
end
%Read number of elements. This is temporary number because Gmsh includes all 'sub'elements as well. Those will have to be discarded.
nElems = str2num(fgetl(fid));
elems = zeros(nElems, 4); %Assume maximal number of nodes, 4 in case of tetrahedra
tags = zeros(nElems, 4);  % Assume maximal number of tags 4

maxElType = 0; skip = 0;
%Read elements. Gmsh element definition: el no, el type, no of tags, <tags>, nodes of the element 
for ie = 1:nElems
  %Read the element no, type and no of tags
  nEl = fscanf(fid, '%u', 1);
  gmshElType = fscanf(fid, '%u', 1);
  nElTags = fscanf(fid, '%u', 1);
  %Read tags
  tag = fscanf(fid, '%u', nElTags);
  %read the nodes of the element (no of nodes depends on the element type)
  switch gmshElType
    case 1 %2-node line
      element = fscanf(fid, '%u', 2); scrap = fgetl(fid);
      elType = 2;
    case 2 %3-node triangle
      element = fscanf(fid, '%u', 3); scrap = fgetl(fid);
      elType = 3;
    case 3 %4-node quadrangle
      element = fscanf(fid, '%u', 4); scrap = fgetl(fid);
      elType = 5;
    case 4 %4-node tetratheder
      element = fscanf(fid, '%u', 4); scrap = fgetl(fid);
      elType = 4;
    case 15 %1-node point
      element = fscanf(fid, '%u', 1); scrap = fgetl(fid);
      elType = 1;
    otherwise %other element types are not recognized
      scrap = fgetl(fid); skip = true;
      elType = 0;
  end
  %Find the maximal element type
  if elType > maxElType
    maxElType = elType;
  end
  %Skip unrecongized elements
  if ~skip
    elems(ie, 1:length(element)) = element;
    tags(ie, 1:length(tag)) = tag;
    skip = false;
  end
end
if ~strcmp(fgetl(fid), '$EndElements')
  error('gmshMeshReader:Problem with reading in the elements. Inspect the gmsh file.');
end

%Close Gmsh file
fclose(fid);

%Remove the elements at lower dimensional levels and all not-recongized elements - rows of zeros
perm = sum(elems,2)>0;
elems = elems(perm,:);
tags = tags(perm,:);
%elems = elems(sum(elems, 2)>0, :); %remove unrecongnized

if maxElType <= 4 %Do only for n-simplex meshes
  perm = elems(:,maxElType)>0;
  elems = elems(perm,1:maxElType);
  tags = tags(perm,:);
  %elems = elems(elems(:,maxElType)>0, 1:maxElType); %remove lower dimensions
end
nElems = size(elems,1);

%Remove nodes that do not belong to any element
looseNodeIds = setdiff(1:nNodes, unique(elems(:))); %find loose nodes numbers
nodes = nodes(setdiff(1:nNodes,looseNodeIds),:);
nNodes = size(nodes,1);
%adjust node numbers in elems
for jn = 1:length(looseNodeIds)
  elems(elems > looseNodeIds(jn)) = elems(elems > looseNodeIds(jn))-1;
end

%From the length of elements node field deduce the space dimension n
n = size(elems,2)-1;

%Remove the non-relevant dimensions from the node list (columns which are equal)
firstNode = nodes(1,:);
relDims = (sum(abs(nodes - ones(nNodes,1)*firstNode), 1) > 0); 
%Assign the nodes 
nodes = nodes(:,relDims);

%Set the fields of the structure
mesh.n = n;
mesh.p = p;
mesh.nElems = nElems;
mesh.nNodes = nNodes;
mesh.nodes = nodes;
%Set node's parents to 0 for all nodes of the original mesh
mesh.pNodes = zeros(mesh.nNodes, 2);
mesh.elems = elems;

mesh.tags = tags;

%Set most grand-parent simplex to 0 for all the elements in the original mesh
mesh.oSimps = mesh.elems;
mesh.nOSimps = mesh.nElems;