#include "mex.h"
#include "stoastlib.h"
#include "util.h"
#include "mexutil.h"
#include "matlabfdot.h"

// =========================================================================
// =========================================================================
// MatlabFDOT class implementation

MatlabFDOT::MatlabFDOT ()
{
    nfdotfwd = 0;
}

// =========================================================================

MatlabFDOT::~MatlabFDOT ()
{
    if (nfdotfwd) delete []fdotfwdlist;
}

// =========================================================================

FDOTFwd *MatlabFDOT::GetFDOTFwd (const mxArray *idx)
{
    unsigned int fwdid;

    if (mxIsUint32(idx)) {
	fwdid = *(unsigned int*)mxGetData(idx) - 1;
	if (fwdid >= nfdotfwd) {
	    if (nfdotfwd)
		mexPrintf ("GetFDOTFwd: index out of range (expected 1..%d).\n", nfdotfwd);
	    else
		mexPrintf ("GetFDOTFwd: index out of range (no meshes registered).\n");
	    return 0;
	}
	if (!fdotfwdlist[fwdid]) {
	    mexPrintf ("GetFDOTFwd: index points to cleared mesh.\n");
	}
	return fdotfwdlist[fwdid];
    } else {
	mexPrintf ("GetFDOTFwd: Invalid index format (expected uint32).\n");
	return 0;
    }
}

// =========================================================================

void MatlabFDOT::MakeFwd (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFMTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */

    // get mesh pointer from handle
    extern MatlabToast *mtoast;
    QMMesh *mesh = (QMMesh*)mtoast->GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    mesh->MarkBoundary ();
    int nQ, nM;
    int n = mesh->nlen();
    nQ = mesh->nQ;
    nM = mesh->nM;

    // Basis
    Raster *rast = mtoast->GetBasis(prhs[0]);
    ASSERTARG(rast, 2, "Basis not found");

    // QVec
    size_t nr = mxGetM(prhs[2]);
    size_t nc = mxGetN(prhs[2]);
    if ((nc != nQ) || (nr != n))
    {
	cerr << "qVec is wrong size!" << endl;
    }
    RCompRowMatrix qvec(n, nQ);
    RDenseMatrix qvdense(n, nQ);
    CopyTMatrix(qvec, prhs[2]);

    // FWD solver 
    char solverType[256];
    mxGetString(prhs[4], solverType, 256);
    double lin_tol = mxGetScalar(prhs[5]);    
    cout << "Creating forward solver: " << solverType << ", tolerance = " << lin_tol << endl;
    RFwdSolver *FWS = new RFwdSolver (mesh, solverType, lin_tol);
    cout << endl << "Allocating system matrix" << endl;
    FWS->Allocate ();

    // Set optical parameters of FEM solver
    Solution sol(OT_NPARAM, n);
    RVector mua, mus, ref, otc2a(n);
    // Get values from args
    CopyVector (mua, prhs[6]);
    CopyVector (mus, prhs[7]);
    CopyVector (ref, prhs[8]);
    double freq = mxGetScalar (prhs[9]);
    ASSERTARG (n == mua.Dim(), 7, "Mua: wrong size");
    ASSERTARG (n == mus.Dim(), 8, "Mus: wrong size");
    ASSERTARG (n == ref.Dim(), 9, "Ref: wrong size");
    // Set up Solution and calculate system matrix
    sol.SetParam (OT_CMUA,   mua*c0/ref);
    sol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    sol.SetParam (OT_N, ref);
    for (int i = 0; i < otc2a.Dim(); i++)
	otc2a[i] = c0/(2*ref[i]*A_Keijzer(ref[i]));
    sol.SetParam (OT_C2A, otc2a);
    FWS->Reset (sol, freq);

    // Projectors
    Projector ** projPList = new Projector*[nM];
    nr = mxGetM(prhs[3]);
    nc = mxGetN(prhs[3]);
    if ( (nc>1) || (nr<nM) )
    {
	cerr << "Not enough projectors in list" << endl;
    }
    RVector projHandleList(nM);
    CopyVector(projHandleList, prhs[3]);
    for (int i=0; i<nM; ++i) projPList[i] = (Projector*) Handle2Ptr(projHandleList[i]);

    // Create the FDOT forward model
    FDOTFwd * fdot;
    fdot = new FDOTFwd(FWS, *mesh, rast, nQ, qvec, projPList);

    plhs[0] = mxCreateDoubleScalar (Ptr2Handle (fdot));


    // ...
//    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
}

// =========================================================================

void MatlabFDOT::DelFwd (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    unsigned int fwdid;

    if (mxIsUint32(prhs[0])) {
	fwdid = *(unsigned int*)mxGetData(prhs[0]) - 1;
	if (fwdid >= nfdotfwd) {
	    if (nfdotfwd)
		mexPrintf ("DelFDOTFwd: index out of range (expected 1..%d).\n", nfdotfwd);
	    else
		mexPrintf ("DelFDOTFwd: index out of range (no meshes registered).\n");
	    return;
	}
    } else {
	mexPrintf ("DelFDOTFwd: Invalid index format (expected uint32).\n");
	return;
    }

    if (fdotfwdlist[fwdid]) {
	delete fdotfwdlist[fwdid];
	fdotfwdlist[fwdid] = 0;
    } else {
	mexPrintf ("DelFDOTFwd: solver already cleared.\n");
    }
}

// =========================================================================

void MatlabFDOT::FwdOp (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");
    
    // get x vector
    RVector x;
    CopyVector (x, prhs[1]);

    // calculate data vector
    RVector y;
    fsolver->fwdOperator (x, y, false);
    CopyVector (&plhs[0], y);
}

// =========================================================================

void MatlabFDOT::FwdOpRatio (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    if (nrhs < 3)
    {
        mexErrMsgTxt ("FDOTFwdOpRatio: Invalid function parameters");
	return;
    }

    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");
    
    // get x vector
    RVector x;
    CopyVector (x, prhs[1]);

    double eps = mxGetScalar(prhs[2]);

    // calculate data vector
    RVector y;
    fsolver->fwdOperator (x, y, true, eps);
    CopyVector (&plhs[0], y);
}

// =========================================================================

void MatlabFDOT::AdjOp (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");

    // get data vector
    RVector y;
    CopyVector (y, prhs[1]);

    // single-source index
    int q = -1;
    if (nrhs >= 3) {
        q = (int)(mxGetScalar(prhs[2]));
    }

    // calculate parameter vector
    RVector x;
    if (q >= 0) fsolver->adjOperator (y, x, q, false);
    else        fsolver->adjOperator (y, x, false);
    CopyVector (&plhs[0], x);
    /*
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
    */
}

// =========================================================================

void MatlabFDOT::AdjOpRatio (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */
    // get FDOTFwd pointer from handle

    if (nrhs < 3)
    {
        mexErrMsgTxt ("FDOTAdjOpRatio: Invalid function parameters");
	return;
    }

    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");
    
    // get data vector
    RVector y;
    CopyVector (y, prhs[1]);

    // single-source index
    int q = -1;
    if (nrhs >= 4) {
        q = (int)(mxGetScalar(prhs[3]));
    }

    double eps = mxGetScalar(prhs[2]);

    // calculate parameter vector
    RVector x;
    if (q >= 0) fsolver->adjOperator (y, x, q, true, eps);
    else        fsolver->adjOperator (y, x, true, eps);
    CopyVector (&plhs[0], x);
    /*
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
    */
}

// =========================================================================

void MatlabFDOT::Excit (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");
    
    // calculate data vector
    RVector y;
    fsolver->getExcitationImages(y);
    CopyVector (&plhs[0], y);
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
}

// =========================================================================

void MatlabFDOT::Sysmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = GetFDOTFwd (prhs[0]);
    ASSERTARG(fsolver, 1, "DOT fwdsolver not found");
    CopyMatrix( &plhs[0], fsolver->getSysmat() );    

    /*
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
    */
}
