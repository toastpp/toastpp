function p = toastIntGradFGradG
%toastIntGradFGradG   - Product of nodal function gradients over mesh.
%
% Synopsis: r = toastIntFG(hMesh,f,g)
%    hMesh: mesh handle
%    f      first function (length h)
%    g      second function (length h)
%
% Returns the product of the gradients of the two functions, integrated
% over the mesh elements. This function evaluates
%
%   r_i = sum_j sum_k g_j f_k int u_i(r) grad(u_j(r)) grad(u_k(r)) dr
%
% over each element, where u_i, u_j and u_k are the shape functions
% associated with nodes i, j and k.
% Vectors f and g can be both real or complex-valued (but not mixed).
%
% See also: toastIntFG
