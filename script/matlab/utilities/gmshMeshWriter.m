function flag = gmshMeshWriter(mesh, gmshFile, values)
% GMSHMESHWRITER Writes mesh in Gmsh ASCII format. Use ONLY WITH SIGNLE LEVEL meshes!
%   flag = gmshMeshWriter(mesh, gmshFile)
%   To extract single level mesh form a refined mesh use getMeshAtLevelP.m

%Predefined fprint formats
MESH_FORMAT = '2.1 0 8';
NODE_FORMAT = '%u %g %g %g';
ELEM_FORMAT = '%u %u %u %u %u';
NODE_DATA_FORMAT1 = '%u %g';
NODE_DATA_FORMAT3 = '%u %g %g %g';
NODE_DATA_FORMAT9 = '%u %g %g %g %g %g %g %g %g %g';

%Values used for output of values to the mesh
NODE_DATA_FORMAT = NODE_DATA_FORMAT1;
viewTitle = 'Solution';
viewTime = 0;
viewTimeStep = 0;
viewNComps = 1; 
viewNNodes = mesh.nNodes(end);
viewPartition = 0;

%Open file
fid = fopen(gmshFile, 'w');
flag = (fid == -1);

%If the dimension of the mesh is not specified (e.g. surface mesh in 3D, d=2), assume mesh.d = mesh.n
if ~isfield(mesh, 'd')
  mesh.d = mesh.n;
end
if ~flag
  if length(mesh.nNodes) > 1
    error('gmshMeshWriter:Multilevel mesh. This writer only works with single level meshes');
  end
  
  switch mesh.d
    case 1
      gmshElType = 1;
    case 2
      gmshElType = 2;
    case 3
      if size(mesh.elems,2) == 3
          gmshElType = 2;
      else
          gmshElType = 4;
      end
  end
  
  
  for id = 1:size(mesh.elems,2) %  1:mesh.d+1
    ELEM_FORMAT = [ELEM_FORMAT, ' %u'];
  end
  ELEM_FORMAT =  [ELEM_FORMAT, '\n'];
  
  fprintf(fid, '$MeshFormat\n');
  fprintf(fid, [MESH_FORMAT '\n']);
  fprintf(fid, '$EndMeshFormat\n');
  
  %Write nodes
  fprintf(fid, '$Nodes\n');
  fprintf(fid, '%u\n', mesh.nNodes);
  for in = 1:mesh.nNodes
    node = zeros(1,3);
    node(1:mesh.n) = mesh.nodes(in,:);
    fprintf(fid, [NODE_FORMAT '\n'], in, node);
  end
  fprintf(fid, '$EndNodes\n');
  
  %write elements
  fprintf(fid, '$Elements\n');
  fprintf(fid, '%u\n', mesh.nElems);
  for ie = 1:mesh.nElems
    element = [gmshElType, 2, 0, 14, mesh.elems(ie, :)];
    fprintf(fid, ELEM_FORMAT, ie, element);
  end
  fprintf(fid, '$EndElements\n');

  %write nodal values
  if nargin >= 3 && ~isempty(values)
    fprintf(fid, '$NodeData\n');
    %write string tags
    fprintf(fid, '%u\n', 1);
    fprintf(fid, [viewTitle '\n']);
    %write real tags
    fprintf(fid, '%u\n', 1);
    fprintf(fid, '%g\n', viewTime); %time value
    %write intereg tags
    fprintf(fid, '%u\n', 4);
    fprintf(fid, '%u\n', viewTimeStep); %timestep
    fprintf(fid, '%u\n', viewNComps);   %number of componenets per node {1,3,9}
    fprintf(fid, '%u\n', viewNNodes);
    fprintf(fid, '%u\n', viewPartition);
    for iv = 1:length(values)
      fprintf(fid, [NODE_DATA_FORMAT '\n'], iv, values(iv,:));
    end
    fprintf(fid, '$EndNodeData\n');
  end
  
  fclose(fid);
end

