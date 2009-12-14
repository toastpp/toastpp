function p = FDOTAdjOp
%FDOTAdjOp            - Fluorescence adjoint operator
%
%   Syntax: x = FDOTAdjOp (hSolver, y [,q])
%     hSolver: solver handle
%     y:       data vector
%     x:       parameter vector
%     q:       source index (>= 1)
%
%  Calculates the adjoint operation x = Adj(y).
%  If source index q is not provided, the result is accumulated over
%  all sources.
