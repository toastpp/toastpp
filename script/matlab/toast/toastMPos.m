function p = toastMPos
%toastMPos            - Returns the detector positions defined for a mesh.
%
% Synopsis: MP = toastMPos (hMesh)
%    hMesh: mesh handle
%
% Returns the detector positions of a mesh in an nm x dim matrix, where nm
% is the number of detector positions, and dim is the dimension of the mesh
% (dim = 2 or 3)
% The mesh must previously have loaded a source-detector description file
% (see toastReadQM).
