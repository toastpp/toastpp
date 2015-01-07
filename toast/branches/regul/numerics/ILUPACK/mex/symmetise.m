function p = symmetise(A)
% p = symmetise(A)
% 
% reorder a given nxn matrix A using METIS multilevel nested dissection
% by edges
% 
% input
% -----
% A         nxn matrix
%
% output
% ------
% p         permutation vector. On exit A(p,p) refers to the reordered
%           system
%

p=symilupackmetise(A);