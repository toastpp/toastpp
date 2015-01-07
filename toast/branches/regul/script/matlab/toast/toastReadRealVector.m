function p = toastReadRealVector
%toastReadRealVector  - Read a real-valued TOAST vector from a file.
%
% Synopsis: v = toastReadRealVector('fname', idx)
%    fname: file name
%    idx:   vector index (optional)
%
% Reads a real-valued vector from a file and returns it as a column vector.
% In the file, the vector must be enclosed by square brackets, and the vector
% elements must be separated by white space.
%
% Example:
% [0.2 -3 1.4e2]
%
% If the file contains more than one vector, the optional 'idx' parameter
% can be used to specify which vector to read. Default is idx=1, i.e. read
% the first vector in the file.
%
% If the provided index is out of range, a vector of length 0 is returned.
%
% To read the last vector in a sequence of real vectors stored in 'fname',
% set idx=0.
