function [p,D] = symmwmmmd(A)
% [p,D] = symmwmmmd(A)
% 
% reorder and rescale a given nxn SYMMETRIC/HERMITIAN matrix A using symmetric
% maximum weight matching followed by minimum degree
% 
% input
% -----
% A         nxn matrix
%
% output
% ------
% p         permutation vector. On exit D*A(p,p)*D refers to the reordered
%           and rescaled system
%

[p,D]=symmwmilupackmmd(0.5*(abs(A)+abs(A)'));
n=size(A,1);
D=spdiags(D(p),0,n,n);
