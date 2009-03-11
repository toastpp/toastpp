function p = toastDeleteMesh
%toastDeleteMesh      - De-allocate a mesh object from dynamic heap.
%
%   Synopsis: toastDeleteMesh(hMesh)
%     hMesh: mesh handle
%
%   Deletes a mesh from dynamic heap. The mesh handle is no longer valid
%   after this call.
%
%   See also toastReadMesh.
