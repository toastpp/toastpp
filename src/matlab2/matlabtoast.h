// -*-C++-*-
// ========================================================================
// Declaration of class MatlabToast
// MEX interface for TOAST functions
// ========================================================================

#ifndef __MATLABTOAST_H
#define __MATLABTOAST_H

#include "mex.h"
#include "stoastlib.h"
#include "../common/objmgr.h"
#include "FDOTFwd.h"

#define ASSERTARG(cond,argno,errmsg) AssertArg(cond,__func__,argno,errmsg)
#define ASSERTARG_CHAR(argno) AssertArg_Char(prhs[argno],__func__,(argno+1))
#define GETSTRING_SAFE(idx,pc,len) GetString_Safe(prhs[idx],pc,len,__func__,(idx+1))
#define ASSERTMESH(meshptr,argno) ASSERTARG(meshptr,argno,"Invalid mesh index")
#define GETMESH_SAFE(idx) GetMesh_Safe(prhs[idx],__func__,idx+1)
#define GETBASIS_SAFE(idx) GetBasis_Safe(prhs[idx],__func__,idx+1)

uint64_T Ptr2Handle (void *ptr);
void *Handle2Ptr (uint64_T handle);
void AssertArg (bool cond, const char *func, int argno, const char *errmsg);
void AssertArg_Char(const mxArray *arg, const char *func, int argno);
bool fileExists(const std::string& fileName);
void GetString_Safe (const mxArray *arg, char *pc, int len,
    const char *func, int argno);

#define MESH_INDEXOUTOFRANGE 1
#define MESH_INDEXCLEARED    2
#define MESH_INDEXFORMAT     3

class MatlabToast {
public:
    MatlabToast();
    ~MatlabToast();

	// Methods defined in matlabtoast.cc
    void SetVerbosity (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ThreadCount (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshLin2Quad (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SetQM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void FindElement (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ReadNIM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void WriteNIM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Qvec (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Mvec (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void QPos (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void MPos (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Sysmat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Sysmat_basis (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void Volmat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Bndmat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void BndReflectionTerm (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SetBasis (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ClearBasis (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisSize (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisNLen (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisBLen (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisSLen (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisBuu (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisBvv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisBuv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisBvw (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisDuu (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisDvv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisDuv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisGuv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GetBasisGvv (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void BasisValue (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void GetBasisSupportArea (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void BasisRefine (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void MapBasis (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void MapMeshToBasis (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapMeshToGrid (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapMeshToSol (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapBasisToMesh (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapSolToMesh (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapSolToBasis (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapSolToGrid (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapGridToMesh (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapGridToBasis (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MapGridToSol (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void BasisToMeshMatrix (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshToBasisMatrix (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SampleBasis (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void GridElref (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SampleField (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ImageGradient (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ReadVector (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void WriteVector(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SolutionMask (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void Regul (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ClearRegul (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulValue (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulGradient (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulHDiag (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void RegulHess (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulHess1f (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulKappa (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void RegulSetLocalScaling (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void Fields (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Gradient (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Jacobian (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void JacobianCW (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void Krylov (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void LBFGS (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

    void FDOTAdjOp (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);

    Mesh *GetMesh (const mxArray *idx, int *errid = 0);
    Mesh *GetMesh_Safe (const mxArray *arr, const char *func, int argno);
    Raster *GetBasis (const mxArray *idx, int *errid = 0,
        bool allownull = false);
    Raster *GetBasis_Safe (const mxArray *arr, const char *func, int argno);
    Regularisation *GetRegul (const mxArray *idx);

	// Methods defined in mtMesh.cc
    void MakeMesh (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ReadMesh (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void WriteMesh (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void WriteMeshVtk (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshOpt (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void MeshReorder (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void MeshNodeCount (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void MeshElementCount (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshDimension (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshBB (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void MeshSize (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ClearMesh (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void MeshData (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void SurfData (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ElementNeighbours (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void MarkMeshBoundary (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SysmatSparsityStructure (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SysmatComponent (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void Massmat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ReadQM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void WriteQM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void GetQM (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void DataLinkList (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void MeshRefine (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void SplitElement (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void NodeNeighbour (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void UnwrapPhase (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void IntFG (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);
    void IntGradFGradG (int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]);

    // Methods defined in mtElement.cc
    void ElDof (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ElSize (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ElData (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ElRegion (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ElMat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void ShapeFunc (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void ShapeGrad (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);

protected:
    void GradientReal (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void GradientCplx (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);

private:
    static void ErrorHandler (char *msg);

    Raster **basislist;       // list of basis instantiations
    unsigned int nbasis;      // number of bases

    unsigned int verbosity;   // verbosity level
};

#endif // !__MATLABTOAST_H
