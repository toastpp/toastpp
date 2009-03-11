function p = toastReadNIM
%toastReadNIM         - Read a nodal image from file.
%
% Synopsis: nim = toastReadNIM('nimfile', idx)
%    nimfile: NIM file name
%    idx:     image index (>=1, or 0 for last image in file; default: 1)
%
% Reads a nodal image from a file, and returns it in a column vector.
