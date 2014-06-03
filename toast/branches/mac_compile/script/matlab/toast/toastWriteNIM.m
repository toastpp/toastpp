function p = toastWriteNIM
%toastWriteNIM        - Write a nodal parameter array to an image file.
%
% Synposis: toastWriteNIM ('fname', 'meshname', vec)
%     fname:    NIM file name
%     meshname: name of the associated mesh file
%     vec:      real vector of dimension n x 1 containing nodal parameter
%               values, where n is the number of nodes in the mesh.
