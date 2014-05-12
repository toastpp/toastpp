function p = toastSysmatComponent
%toastSysmatComponent - Generate an FEM system matrix term.
%
% Synopsis: S = toastSysmatComponent (hMesh, int_tp, prm)
%           S = toastSysmatComponent (hMesh, int_tp, prm, 'EL')
%
%    hMesh:  mesh handle
%    int_tp: integration type (string; see notes)
%    prm:    nodal parameter (real column vector)
%    'EL'    flag to indicate element basis
%    S:      system matrix (complex, sparse)
%
% Returns a sparse complex n x n matrix (n: number of nodes in the mesh)
% containing system matrix component S. S is constructed by integrating
% the parameter (prm), given by its nodal values and a combination of
% shape functions or shape function derivatives over all mesh elements.
%
% The type of integration is selected via the int_tp string. The
% following are currently supported:
%
% int_tp      Integral
% --------------------------------------------
% 'FF'        \int Fi(r) Fj(r) dr
% 'DD'        \int Di(r) Dj(r) dr
% 'PFF'       \int prm(r) Fi(r) Fj(r) dr
% 'PDD'       \int prm(r) Di(r) Dj(r) dr
% 'BNDPFF'    \int prm(r) Fi(r) Fj(r) dr (over boundary)
%
% where Fi, Di: shape function/shape function derivative for node i
%
% If the 'EL' flag is present, the system matrix is calculated on an
% element basis, rather than a node basis. Length 'n' of all parameter
% vectors in that case must be equal to the number of elements.
% Parameters are considered constant over each element, so parameter
% distribution has discontinuities at element boundaries.
%
% S is a building block for creating the FEM system matrix.

% See also: toastSysmat.
