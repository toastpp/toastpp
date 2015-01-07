function p = toastReadComplexVector
%toastReadComplexVector - Read a complex TOAST vector from a file.
%
% Synopsis: v = toastReadComplexVector('fname', idx)
%    fname: file name
%    idx:   vector index (optional)
%
% Reads a complex-valued vector from a file and returns it as a column vector.
% In the file, the vector must be enclosed by square brackets, and each
% complex number is defined by a pair of real and imaginary value, enclosed in
% < > and separated by white space.
%
% Example:
% [<0 1> <-0.5 1.3>]
%
% If the file contains more than one vector, the optional 'idx' parameter
% can be used to specify which vector to read. Default is idx=1, i.e. read
% the first vector in the file.
%
% If the provided index is out of range, a vector of length 0 is returned.
%
% To read the last vector in a sequence of complex vectors stored in 'fname',
% set idx=0.
