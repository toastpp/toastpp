function p = toastRegulHess
%toastRegulHess       - Returns regularisation Hessian for an image.
%
% Synopsis: H = toastRegulHess (hReg, x)
%    hReg: regulariser (handle)
%    x:    coefficient vector (real vector)
%    H:    regularisation Hessian (real sparse matrix)
%
% Given a regulariser hReg (see toastRegul) and a parameter distribution x,
% this function returns the corresponding Hessian matrix in sparse format.
