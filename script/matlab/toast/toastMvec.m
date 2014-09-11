function p = toastMvec
%toastMvec            - Generate a sparse matrix of measurement column vectors.
%
% Synopsis: M = toastMvec (hMesh, mprof, mwidth)
%    hMesh:  mesh handle
%    mprof:  measurement profile (string: 'Gaussian', 'Cosine', 'TrigBasis')
%    mwidth: measurement radius [mm]
%
% Generates the boundary measurement operators (Direchlet-to-Neumann map)
% for the OT FEM forward problem from measurement specifications stored in a
% mesh. The vectors are returned as column vectors in a sparse m x n matrix,
% where m is the mesh node count, and n is the number of sources. The mesh
% (given by mesh handle hMesh) must contain measurement definitions (see
% toastReadQM).
%
% The boundary operator is normally scaled with a factor of c/2A, where c is
% the speed of light, and A is a term incorporating the refractive index
% mismatch at the medium-air interface. This scaling factor is computed from
% the 'ref' term passed to Mvec (either as a vector of length n of nodal
% values, or as a scalar if the refractive index is homogeneous.
%
% To avoid scaling with the c/2A term altogether, set the 'ref' parameter to 0.
%
% Warning: if 'ref' is omitted, Mvec uses the refractive index values stored
% with the mesh. This behaviour may change in future versions.

