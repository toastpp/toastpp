// TOAST mex driver.
// This is a single entry point for Matlab toast functions that are
// implemented in C++.

#include "mex.h"
#include "toastmex.h"

MatlabToast *mtoast = 0;
MatlabFDOT *mfdot = 0;

void mexClear();

void ProvideToast ()
{
    if (!mtoast) {
	mtoast = new MatlabToast;
	mexAtExit (mexClear);

#ifdef TOAST_THREAD
	Task_Init (0);
#endif
    }
}

void ProvideFDOT ()
{
    ProvideToast();
    if (!mfdot) mfdot = new MatlabFDOT;
}

// =========================================================================
// MEX entry point

void mexClear ()
{
  //    if (mtoast) {
  //	delete mtoast;
  //	mtoast = 0;
  //    }
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned int funcid;

    if (nrhs> 0 && mxIsUint32 (prhs[0])) {
	funcid = *(unsigned int*)mxGetData (prhs[0]);
	nrhs--;
	prhs++;
    }
    else {
        mexErrMsgTxt ("toast: Invalid toast driver call: integer function id expected in first argument.");
    }

    ProvideToast();

    switch (funcid) {
    case TOAST_CLEARMESH:
	mtoast->ClearMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MARKMESHBOUNDARY:
	mtoast->MarkMeshBoundary (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHLIN2QUAD:
	mtoast->MeshLin2Quad (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SETQM:
	mtoast->SetQM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_FINDELEMENT:
	mtoast->FindElement (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_READNIM:
	mtoast->ReadNIM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_WRITENIM:
	mtoast->WriteNIM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_QVEC:
	mtoast->Qvec (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MVEC:
	mtoast->Mvec (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_QPOS:
	mtoast->QPos (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MPOS:
	mtoast->MPos (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SYSMAT:
	mtoast->Sysmat (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SYSMATBASIS:
	mtoast->Sysmat_basis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_VOLMAT:
	mtoast->Volmat (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BNDMAT:
	mtoast->Bndmat (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BNDREFLECTIONTERM:
	mtoast->BndReflectionTerm (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SETBASIS:
	mtoast->SetBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_CLEARBASIS:
	mtoast->ClearBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_GETBASISSIZE:
	mtoast->GetBasisSize (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_NLEN:
        mtoast->GetBasisNLen (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_BLEN:
        mtoast->GetBasisBLen (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_SLEN:
        mtoast->GetBasisSLen (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_VALUE:
	mtoast->BasisValue (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_BUU:
	mtoast->GetBasisBuu (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_BVV:
	mtoast->GetBasisBvv (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_BUV:
	mtoast->GetBasisBuv (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_BVW:
        mtoast->GetBasisBvw (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_DUU:
	mtoast->GetBasisDuu (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_DVV:
	mtoast->GetBasisDvv (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_DUV:
	mtoast->GetBasisDuv (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_SUPPORTAREA:
	mtoast->GetBasisSupportArea (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASIS_REFINE:
	mtoast->BasisRefine (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPBASIS:
	mtoast->MapBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPMESHTOBASIS:
	mtoast->MapMeshToBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPMESHTOGRID:
	mtoast->MapMeshToGrid (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPMESHTOSOL:
	mtoast->MapMeshToSol (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPBASISTOMESH:
	mtoast->MapBasisToMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPSOLTOMESH:
	mtoast->MapSolToMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPSOLTOBASIS:
	mtoast->MapSolToBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPSOLTOGRID:
	mtoast->MapSolToGrid (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPGRIDTOMESH:
	mtoast->MapGridToMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPGRIDTOBASIS:
	mtoast->MapGridToBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MAPGRIDTOSOL:
	mtoast->MapGridToSol (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASISTOMESHMATRIX:
	mtoast->BasisToMeshMatrix (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_BASISSAMPLE:
	mtoast->SampleBasis (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_GRIDELREF:
	mtoast->GridElref (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SAMPLEFIELD:
	mtoast->SampleField (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_IMAGEGRADIENT:
	mtoast->ImageGradient (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_READVECTOR:
	mtoast->ReadVector (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_WRITEVECTOR:
	mtoast->WriteVector (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SOLUTIONMASK:
	mtoast->SolutionMask (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SETVERBOSITY:
	mtoast->SetVerbosity (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_THREADCOUNT:
	mtoast->ThreadCount (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGUL:
	mtoast->Regul (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_CLEARREGUL:
	mtoast->ClearRegul (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULVALUE:
	mtoast->RegulValue (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULGRADIENT:
	mtoast->RegulGradient (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULHDIAG:
	mtoast->RegulHDiag (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULHESS:
	mtoast->RegulHess (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULHESS1F:
	mtoast->RegulHess1f (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULKAPPA:
	mtoast->RegulKappa (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_REGULSETLOCALSCALING:
        mtoast->RegulSetLocalScaling (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_FIELDS:
	mtoast->Fields (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_GRADIENT:
	mtoast->Gradient (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_JACOBIAN:
	mtoast->Jacobian (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_JACOBIANCW:
	mtoast->JacobianCW (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_KRYLOV:
	mtoast->Krylov (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_LBFGS:
        mtoast->LBFGS (nlhs, plhs, nrhs, prhs);
	break;

    // methods defined in mtMesh.cc
    case TOAST_MAKEMESH:
	mtoast->MakeMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_READMESH:
	mtoast->ReadMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_WRITEMESH:
	mtoast->WriteMesh (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_WRITEMESHVTK:
	mtoast->WriteMeshVtk (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHOPT:
	mtoast->MeshOpt (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHREORDER:
	mtoast->MeshReorder (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHNODECOUNT:
	mtoast->MeshNodeCount (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHELEMENTCOUNT:
	mtoast->MeshElementCount (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHDIMENSION:
	mtoast->MeshDimension (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHBB:
	mtoast->MeshBB (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHSIZE:
	mtoast->MeshSize (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHDATA:
	mtoast->MeshData (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SURFDATA:
	mtoast->SurfData (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_ELEMENTNEIGHBOURS:
	mtoast->ElementNeighbours (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SYSMATCOMPONENT:
    mtoast->SysmatComponent (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SPARSITYSTRUCTURE:
    mtoast->SysmatSparsityStructure(nlhs, plhs, nrhs, prhs);
    break;
    case TOAST_MASSMAT:
	mtoast->Massmat (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_READQM:
	mtoast->ReadQM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_WRITEQM:
	mtoast->WriteQM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_GETQM:
	mtoast->GetQM (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_DATALINKLIST:
	mtoast->DataLinkList (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_MESHREFINE:
        mtoast->MeshRefine (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SPLITELEMENT:
	mtoast->SplitElement (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_NODENEIGHBOUR:
	mtoast->NodeNeighbour (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_UNWRAPPHASE:
	mtoast->UnwrapPhase (nlhs, plhs, nrhs, prhs);
	break;

    // methods defined in mtElement.cc
    case TOAST_ELDOF:
        mtoast->ElDof (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_ELEMENTSIZE:
	mtoast->ElSize (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_ELEMENTDATA:
	mtoast->ElData (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_ELEMENTREGION:
	mtoast->ElRegion (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_ELMAT:
	mtoast->ElMat (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SHAPEFUNC:
	mtoast->ShapeFunc (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_SHAPEGRAD:
	mtoast->ShapeGrad (nlhs, plhs, nrhs, prhs);
	break;
// Integral methods added by SP
    case TOAST_INTFG:
	mtoast->IntFG (nlhs, plhs, nrhs, prhs);
	break;
    case TOAST_INTGRADFGRADG:
	mtoast->IntGradFGradG (nlhs, plhs, nrhs, prhs);
	break;

    // FDOT interface
    case FDOT_MAKEFWD:
	ProvideFDOT();
	mfdot->MakeFwd (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_CLEARFWD:
	ProvideFDOT();
	mfdot->DelFwd (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_FWDOP:
	ProvideFDOT();
	mfdot->FwdOp (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_FWDOPRATIO:
	ProvideFDOT();
	mfdot->FwdOpRatio (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_ADJOP:
	ProvideFDOT();
	mfdot->AdjOp (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_ADJOPRATIO:
	ProvideFDOT();
	mfdot->AdjOpRatio (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_EXCIT:
	ProvideFDOT();
	mfdot->Excit (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_SYSMAT:
	ProvideFDOT();
	mfdot->Sysmat (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_MAKEPROJECTORLIST:
	ProvideFDOT();
	mfdot->MakeProjectorList (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_PROJECTTOIMAGE:
	mfdot->ProjectToImage (nlhs, plhs, nrhs, prhs);
	break;
    case FDOT_PROJECTTOFIELD:
	mfdot->ProjectToField (nlhs, plhs, nrhs, prhs);
	break;

    default:
	mexPrintf ("toast: requested function id=%d\n", funcid);
	mexErrMsgTxt ("toast: Invalid toast driver call: function id not recognised.");
	break;
    }
}
