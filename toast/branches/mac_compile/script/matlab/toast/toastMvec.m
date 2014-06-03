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
