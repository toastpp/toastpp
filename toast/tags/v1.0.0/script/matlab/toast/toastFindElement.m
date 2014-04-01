function p = toastFindElement
%toastFindElement     - Find the mesh element containing a point.
%
% Synopsis: el = toastFindElement(hMesh, pt)
%    hMesh: mesh handle
%    vtx:   point coordinates (2-D or 3-D coordinate array)
%    el:    index of the element containing pt (>=1, or 0 if pt is outside
%           the mesh).
%
% This function searches through the list of mesh element and checks pt
% against each. It can therefore be slow if applied to a large number of
% points.
%
% For points where the element association is not unique (e.g. on element
% edges) this function returns the first matching element, but it is
% possible that no elements are found due to rounding errors.
