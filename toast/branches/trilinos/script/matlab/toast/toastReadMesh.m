function p = toastReadMesh
%toastReadMesh        - Read a TOAST mesh from file.
%
%   Syntax: hMesh = toastReadMesh('meshfile')
%     meshfile: mesh file name
%     hMesh:    mesh handle
%
%   Reads a TOAST mesh from a file and returns a handle. The mesh handle
%   cannot be used by MATLAB directly but can be passed to other toastXXX
%   functions.
%   After use, the mesh should be de-allocated with toastDeleteMesh.
%
%   See also toastDeleteMesh.
