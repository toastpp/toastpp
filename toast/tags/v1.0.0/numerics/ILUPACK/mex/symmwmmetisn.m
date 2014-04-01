function [p,D] = symmwmmetisn(A)
% [p,D] = symmwmmetisn(A)
% 
% reorder and rescale a given nxn SYMMETRIC/HERMITIAN matrix A using symmetric
% maximum weight matching followed by METIS multilevel nested dissection by nodes
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

[p,D]=symmwmilupackmetisn(0.5*(abs(A)+abs(A)'));
n=size(A,1);
D=spdiags(D(p),0,n,n);
