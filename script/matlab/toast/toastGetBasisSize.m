function p = toastGetBasisSize
%toastGetBasisSize    - Returns the grid dimensions of reconstruction basis.
%
% Synopsis: [bdim gdim] = toastGetBasisSize(hBasis)
%    hBasis: basis mapper handle
%    bdim: b-basis dimensions [1x2 or 1x3 integer]
%    gdim: (optional) g-basis dimensions [1x2 or 1x3 integer]
%
% Returns the grid size of the regular reconstruction basis (b-basis), and
% optionally the grid size of the high-resolution intermediate basis
% (g-basis). If the mapper instance was created without a third argument
% in the call to toastSetBasis, then bdim and gdim are identical.
%
% The returned vectors are integer arrays of dimension 1x2 or 1x3,
% depending on the dimension of the problem.