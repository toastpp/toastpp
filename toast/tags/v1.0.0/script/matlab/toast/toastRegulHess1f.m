function p = toastRegulHess1f
%toastRegulHess1f     - Applies 1st order regularisation Hessian to vector.
%
% Synopsis: y = toastRegulHess1f (hReg, x, f)
%    hReg: regulariser (handle)
%    x:    coefficient vector (real vector)
%    f:    operand vector (real vector)
%    y:    result of operation (real vector)
%
% Given a regulariser hReg (see toastRegul) and a parameter distribution x,
% this function applies the 1st order Hessian L(x) to a vector f and
% returns the resulting vector.
