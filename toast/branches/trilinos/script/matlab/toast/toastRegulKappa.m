function p = toastRegulKappa
%toastRegulKappa      - Returns diffusivity reference for prior.
%
% Synopsis: ref = toastRegulKappa (hReg,x)
%    hReg: regulariser (handle)
%    x:    parameter distribution (vector)
%    ref:  diffusivity reference images (real vector)
%
% The returned vector contains the concatenated reference images for both
% absorption and diffusion, in the solution basis.
%
% The reference images are constructed from the images passed to
% toastRegul via the KapRefImage parameter.
%
% The reference images allow to define areas where the regularisation
% should permit the existence of edges in the reconstructed images.