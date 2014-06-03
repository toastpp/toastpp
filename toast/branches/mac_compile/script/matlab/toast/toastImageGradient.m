function p = toastImageGradient
%toastImageGradient   - Generate the gradient of an image in g-basis with mask
%
% Synopsis: grad = toastImageGradient (hBasis, img)
%    hBasis: basis mapper handle (see toastSetBasis)
%    img:    image in g-basis format (real or complex vector)
%
%    grad:   image gradient (real or complex n x g matrix, with n=2,3)
%
% Calculates the gradient of the image given in fine resolution
% (g-basis). Any pixel masked out in the basis are removed from the
% calculation, to avoid gradients at the domain boundary. The gradients
% are returned as a dense matrix of dimension n x g, where n is the
% domain dimension (2 or 3). Real and complex images are supported.
  
