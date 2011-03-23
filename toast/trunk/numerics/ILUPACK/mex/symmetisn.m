function p = symmetisn(A)
% p = symmetisn(A)
% 
% reorder a given nxn matrix A using METIS multilevel nested dissection
% by nodes
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

p=symilupackmetisn(A);
