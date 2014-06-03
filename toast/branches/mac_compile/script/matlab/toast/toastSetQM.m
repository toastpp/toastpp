function p = toastSetQM
%toastSetQM           - Assign source and detector positions to a mesh.
%
% Synopsis: toastSetQM(hMesh,Q,M)
%    hMesh:  mesh handle
%    Q:      source array (real matrix nq x dim)
%    M:      detector array (real matrix nm x dim)
%
% Assigns a list of source and detector locations to a mesh.
% Q is a dense nq x dim matrix containing nq source positions, and M is a
% dense nm x dim matrix containing nm detector positions. dim is the
% dimension of the problem (2 or 3) and must correspond to the mesh
% dimension.
%
% This function assumes a dense link list, i.e. all sources are connected
% to all detectors.
%
% Any sources and detectors previously defined by the mesh are overwritten.
%
% The mesh stores only the pointwise positions of sources and detectors,
% not their type and shape. Shape parameters are only specified when the
% source and measurement vectors are constructed with toastQvec and
% toastMvec.
%
% See also: toastReadQM, toastWriteQM, toastGetQM, toastQvec, toastMvec,
%           toastQPos, toastMPos
