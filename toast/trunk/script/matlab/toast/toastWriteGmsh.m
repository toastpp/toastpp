function toastWriteGmsh (hmesh, fname)
%toastWriteGmsh   - Write a surface mesh in Gmsh format
%
% Synopsis: toastWriteGmsh(hmesh,'fname')
%   hmesh: handle for a toast surface mesh
%   fname: file name (string)
%
% Note:
% - This generates a mesh file that can be read by gmsh (geuz.org/gmsh)
% - Currently, only surface meshes consisting of linear triangular
%   meshes (element type ELTP_TRI3D3) are supported

[vtx idx] = toastMeshData (hmesh);

nnd = size(vtx,1);
nel = size(idx,1);

% for i=1:nel
%     if eltp(i) ~= 16
%         disp('Unexpected element type. Required: 16');
%         return;
%     end
% end

fid = fopen(fname,'w');
fprintf (fid, '$MeshFormat\n');
fprintf (fid, '2.2 0 8\n');
fprintf (fid, '$EndMeshFormat\n');
fprintf (fid, '$Nodes\n');
fprintf (fid, '%d\n', nnd);
for i=1:nnd
    fprintf (fid, '%d %f %f %f\n', i, vtx(i,1), vtx(i,2), vtx(i,3));
end
fprintf (fid, '$EndNodes\n');
fprintf (fid, '$Elements\n');
fprintf (fid, '%d\n', nel);
for i=1:nel
    fprintf (fid, '%d 2 0 %d %d %d\n', i, idx(i,1), idx(i,2), idx(i,3));
end
fprintf (fid, '$EndElements\n');
fclose (fid);
