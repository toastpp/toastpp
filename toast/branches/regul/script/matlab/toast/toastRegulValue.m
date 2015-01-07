function p = toastRegulValue
%toastRegulValue      - Returns regularisation value of an image.
%
% Synopsis: v = toastRegulValue (hReg,x)
%    hReg: regulariser (handle)
%    x:    parameter distribution (vector)
%    v:    regularisation value (scalar)
%
% Given a regulariser hReg (see toastRegul) and a parameter distribution x,
% this function returns the corresponding regularisation value.
