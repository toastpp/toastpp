function [pl,pr,Dl,Dr] = mwmrcm(A)
% [pl,pr,Dl,Dr] = mwmrcm(A)
% 
% reorder and rescale a given nxn matrix A using maximum weight
% matching followed by reverse Cuthill-McKee
% 
% input
% -----
% A         n x n  matrix
%
% output
% ------
% pl,pr     left and right permutation vectors
% Dl, Dr    left and right scaling matrices
%
%           On exit B=Dl*A(pl,pr)*Dr would refer to the reordered
%           and rescaled system such that in theory |B(i,i)|=1
%           and |B(i,j)|<=1. In practice, powers of 2 are used for
%           scaling matrices Dl,Dr to avoid rounding errors. For 
%           this reason, the scaled entries fulfill the constraint
%           only within the range that can be achieved by the nearest
%           power of 2 for Dl,Dr
 

[pr,pl,Dr,Dl]=mwmilupackrcm(abs(A));
n=size(A,1);
Dl=spdiags(Dl(pl),0,n,n);
Dr=spdiags(Dr(pr),0,n,n);
