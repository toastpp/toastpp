function p = toastMapGridToSol
%toastMapGridToSol    - Map an image from fine grid to inverse solution basis.
%
% Synopsis: sol = toastMapGridToSol (hBasis, img)
%    hbasis: basis mapper handle (see toastSetBasis)
%    img:    image in fine raster grid
%    sol:    image mapped to inverse solution basis.
%
% Maps an image from the fine ("intermediate") grid basis to the solution
% basis representation of the inverse solver (excluding pixels outside
% the support of the mesh).
% The image can be real or complex-valued.
