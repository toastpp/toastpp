function p = toastMassmat
%toastMassmat         - Generate an FEM mass matrix for time-domain OT.
%
% Synopsis: B = toastMassmat (hMesh)
%    hMesh: mesh handle
%    B:     mass matrix (sparse, real)
%
% Returns a sparse real n x n matrix (n: number of nodes in the mesh)
% containing mass matrix B. B is composed from element-wise integrals of
% products of shape functions. Unlike the system matrix S, B does not
% depend on optical parameters.
%
% B is required for the evaluation of the time-dependent FEM diffusion
% equation, where it appears in the time derivative term.
%
% See also: toastSysmat.  
