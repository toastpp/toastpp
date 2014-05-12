function p = toastRegulGradient
%toastRegulGradient   - Returns regularisation gradient for an image
%
% Synopsis: g = toastRegulGradient (hReg, x)
%    hReg: regulariser (handle)
%    x:    coefficient vector (real vector)
%    g:    regularisation gradient vector (real vector)
%
% Given a regulariser hReg (see toastRegul) and a parameter distribution x,
% this function returns the vector representing the corresponding
% regularisation gradient.
