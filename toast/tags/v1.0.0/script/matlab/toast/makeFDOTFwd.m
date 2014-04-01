function p = makeFDOTFwd
%makeFDOTFwd          - Create a fluorescence forward solver instance
%
%  Syntax: hSolver = makeFMTFwd (hMesh, hBasis, qvec, hProjList, solverType, solverTol, Mua, Mus, Ref, Freq)
%    hMesh:   mesh handle
%    hBasis:  basis handle
%    qvec:    Matrix of source vectors as given by toastQvec 
%    hProjList: Array of projector handles, as returned by toastMakeProjectorList
%    solverType: FEM solver type = {DIRECT, BICGSTAB, ...}  
%    solverTol: FEM (iterative) solver tolerance 
%    Mua (nodes x 1)
%    Mus (nodes x 1)
%    Ref (nodes x 1)
%    Freq (must be zero for now!)
%    hSolver: solver handle

