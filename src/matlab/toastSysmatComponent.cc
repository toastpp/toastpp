// =========================================================================
// toastSysmatComponent
// Returns a single term of a system matrix
//
// RH parameters:
//     1: mesh handle (double)
//     2: prm NIM (double array
//     6: (optional) 'EL' flag to indicate element basis
// LH parameters:
//     1: system matrix (sparse double matrix)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

// =========================================================================
// Implementation

void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "toastSysmat: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}

void CalcSysmatComponent (QMMesh *mesh, RVector &prm, char *integral_type,
    bool elbasis, mxArray **res)
{
    int n = (elbasis ? mesh->elen() : mesh->nlen());

    // Create forward solver to initialise system matrix
    CFwdSolver FWS (LSOLVER_DIRECT, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    
    FWS.Allocate (*mesh);
    if (!strcasecmp (integral_type, "PDD")) {
	cerr << "Assembling PDD system matrix component" << endl;
	FWS.AssembleSystemMatrixComponent (prm,
	    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    }

    // Return system matrix to MATLAB
    CopyMatrix (res, *FWS.F);
}

void CalcBndSysmat (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    CFwdSolver FWS (LSOLVER_DIRECT, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    FWS.Allocate (*mesh);
    RVector prm(n);
    for (int i = 0; i < n; i++) prm[i] = c0/ref[i];
    AddToSysMatrix (*mesh, *FWS.F, &prm, ASSEMBLE_BNDPFF);
    CopyMatrix (res, *FWS.F);
}

void CalcBndFactors (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    RVector prm(n);
    for (int i = 0; i < n; i++)
	prm[i] = A_Keijzer (ref[i]);
    CopyVector (res, prm);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    RVector prm;
    int n;
    bool elbasis = false;
    char integral_type[32] = "";

    CopyVector (prm, prhs[1]);

    if (mxIsChar(prhs[2]))
	mxGetString (prhs[2], integral_type, 32);


    if (nrhs >= 4 && mxIsChar(prhs[3])) {
	char cbuf[32];
	mxGetString (prhs[3], cbuf, 32);
	elbasis = (strcasecmp (cbuf, "EL") == 0);
    }
    if (elbasis) n = mesh->elen();
    else         n = mesh->nlen();
    Assert (n == prm.Dim(), "Argument 2: wrong size");

    CalcSysmatComponent (mesh, prm, integral_type, elbasis, &plhs[0]);
}
