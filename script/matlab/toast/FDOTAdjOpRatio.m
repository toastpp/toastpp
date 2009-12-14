function p = FDOTAdjOpRatio
%FDOTAdjOpRatio       - Fluorescence adjoint operator, with Born normalisation
%
%   Syntax: x = FDOTAdjOpRatio (hSolver, y [,q])
%     hSolver: solver handle
%     y:       data vector
%     x:       parameter vector
%     q:       source index (>= 1)
%     eps:     Epsilon value for ratio = fluo / (excit + epsilon)
%
%  Calculates the adjoint operation x = Adj(y).
%  If source index q is not provided, the result is accumulated over
%  all sources.
