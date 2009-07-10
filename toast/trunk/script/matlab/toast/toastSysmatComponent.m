function p = toastSysmatComponent
%toastSysmatComponent - Generate an FEM system matrix term.
%
% Synopsis: S = toastSysmatComponent (hMesh, prm, int_tp)
%           S = toastSysmatComponent (hMesh, prm, int_tp, 'EL')
%
%    hMesh:  mesh handle
%    prm:    nodal parameter (real column vector)
%    int_tp: integration type (string; see notes)
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
% 'PDD'       \int prm(r) Di(r) Dj(r) dr
%
% where Di: shape function derivative for node i
%
% If the 'EL' flag is present, the system matrix is calculated on an
% element basis, rather than a node basis. Length 'n' of all parameter
% vectors in that case must be equal to the number of elements.
% Parameters are considered constant over each element, so parameter
% distribution has discontinuities at element boundaries.
%
% S is a building block for creating the FEM system matrix.

% See also: toastSysmat.
