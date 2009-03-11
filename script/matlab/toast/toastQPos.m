function p = toastQPos
%toastQPos            - Returns the source positions defined for a mesh.
%
% Synopsis: QP = toastQPos (hMesh)
%    hMesh: mesh handle
%
% Returns the source positions of a mesh in an nq x dim matrix, where nq
% is the number of source positions, and dim is the dimension of the mesh
% (dim = 2 or 3)
% The mesh must previously have loaded a source-detector description file
% (see toastReadQM).
%
% Usage example: finding the distances between sources and detectors
%
% qp = toastQPos(hMesh);
% mp = toastMPos(hMesh);
% for i=1:size(qp,1)
%   for j=1:size(mp,1)
%     dst(j,i) = norm(qp(i,:)-mp(j,:))
%   end
% end
