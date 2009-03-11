function p = toastMapMeshToGrid
%toastMapMeshToGrid   - Map a nodal image to a regular grid
%
% Synopsis: img = toastMapMeshToGrid (hBasis, nim)
%    hbasis: basis mapper handle (see toastSetBasis)
%    nim:    nodal image (real vector)
%    img:    grid image
%
% Maps a nodal image into a regular grid (e.g. for displaying on screen).
% The dimension of the grid raster is defined by the basis mapper.
% The image can be real or complex-valued.
