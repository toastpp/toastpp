function p = toastMapGridToBasis
%toastMapGridToBasis  - Map an image from fine to coarse basis grid.
%
% Synopsis: img = toastMapGridToBasis (hBasis, sol)
%    hbasis: basis mapper handle (see toastSetBasis)
%    sol:    image in fine raster grid
%    img:    image mapped to coarse raster grid.
%
% Maps an image from the fine ("intermediate") grid basis to the coarse
% basis representation.
% The image can be real or complex-valued.
