function p = FDOTFwdOpRatio
%FDOTFwdOpRatio       - Fluorescence forward operator, returns Born ratio
%
%   Syntax: y = FDOTFwdOpRatio (hSolver, x)
%     hSolver: solver handle
%     x:       parameter vector
%     eps:     Epsilon value for ratio = fluo / (excit + epsilon)
%     y:       data vector
%
%  Calculates the forward operation y = Fwd(x)
