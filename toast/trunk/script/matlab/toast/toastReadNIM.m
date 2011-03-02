function p = toastReadNIM
%toastReadNIM         - Read a nodal image from file.
%
% Synopsis: nim = toastReadNIM(nimfile, idx)
%           nim = toastReadNIM(nimfile, 'all')
%    nimfile [string]:  NIM file name
%    idx     [integer]: image index (>=1, or 0 for last image in file;
%                       default: 1)
%
% Reads a nodal image from a file, and returns it in a column vector.
% If 'all' is specified as second parameter, all images contained in the
% the file are returned in a matrix, where each column represents one image.
  
