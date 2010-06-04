// Interface:
// RH-1: mesh handle
// RH-2: basis handle
// RH-3: source vector matrix (1 column per source)
// RH-4: Projector handle list (from toastMakeProjectorList)
// RH-5: FEM Solver type
// RH-6: FEM Solver tolerance
// RH-7: Mua
// RH-8: Mus
// RH-9: Ref
// RH-10: Freq
// LH-1: FDOT solver handle

#include "mex.h"
#include "fdotlib.h"
//#include "felib.h"
//#include "stoastlib.h"
#include "util.h"
//#include "projector.h"
//#include "muaSolver.h"
//#include "FDOTFwd.h"
//#include "MLEMSolver.h"
#include "toastmex.h"

#define MAXREGION 100

using namespace std;

void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "makeFDOTFwd: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}


// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFMTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */

    // get mesh pointer from handle
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar(prhs[0]));
    mesh->MarkBoundary ();
    int nQ, nM;
    int n = mesh->nlen();
    nQ = mesh->nQ;
    nM = mesh->nM;

    // Basis
    Raster *rast = (Raster*)Handle2Ptr (mxGetScalar(prhs[1]));

    // QVec
    size_t nr = mxGetM(prhs[2]);
    size_t nc = mxGetN(prhs[2]);
    if ((nc != nQ) || (nr != n)) {
	cerr << "nc=" << nc << ", nQ=" << nQ
	     << ", nr=" << nr << ", n=" << n << endl;
	mexErrMsgTxt("qVec is wrong size!");
    }

    RCompRowMatrix qvec(n, nQ);
    RDenseMatrix qvdense(n, nQ);
    CopyTMatrix(qvec, prhs[2]);

    // FWD solver 
    char solverType[256];
    mxGetString(prhs[4], solverType, 256);
    double lin_tol = mxGetScalar(prhs[5]);    
    cout << "Creating forward solver: " << solverType << ", tolerance = " << lin_tol << endl;

    RFwdSolver *FWS = new RFwdSolver (solverType, lin_tol);
    cout << endl << "Allocating system matrix" << endl;
    FWS->Allocate (*mesh);

    // Set optical parameters of FEM solver
    Solution sol(OT_NPARAM, n);
    RVector mua, mus, ref, otc2a(n);
    // Get values from args
    CopyVector (mua, prhs[6]);
    CopyVector (mus, prhs[7]);
    CopyVector (ref, prhs[8]);
    double freq = mxGetScalar (prhs[9]);
    Assert (n == mua.Dim(), "Mua: wrong size");
    Assert (n == mus.Dim(), "Mus: wrong size");
    Assert (n == ref.Dim(), "Ref: wrong size");
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

    plhs[0] = mxCreateScalarDouble (Ptr2Handle (fdot));


    // ...
//    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
}

// =========================================================================
