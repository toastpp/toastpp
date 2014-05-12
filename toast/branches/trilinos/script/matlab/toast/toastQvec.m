function p = toastQvec
%toastQvec            - Generate a sparse matrix of source column vectors.
%
% Synopsis: Q = toastQvec (hMesh, qtype, qprof, qwidth)
%           Q = toastQvec (hMesh, prm)
%
%    hMesh:  mesh handle
%    qtype:  source type (string: 'Neumann', 'Isotropic')
%    qprof:  source profile (string: 'Point', 'Gaussian', 'Cosine', 'TrigBasis')
%    qwidth: source radius [mm] (real)
%    prm:    source parameter structure (see notes)
%
% Generates the RHS vectors for the FEM forward problem from source
% positions stored in a mesh, and source specifications passed as
% parameters.
%
% Source specifications can be passed as a list of qtype, qprof and
% qwidth, or as a structure prm, which must contain the following fields:
% prm.type, prm.prof and prm.width.
%  
% Neumann sources are defined as a diffuse boundary current (Q is nonzero
% only for boundary nodes), while isotropic sources are defined in the
% mesh volume.
%
% The width parameter defines the size of the sources (1/e radius for
% Gaussian sources, 1/2 radius for trigonometric sources). It is ignored
% for point sources.
%
% The source vectors are returned as column vectors in a sparse m x n
% matrix Q, where m is the mesh node count, and n is the number of
% sources. The mesh (given by mesh handle hMesh) must contain
% source definitions (see toastReadQM).
