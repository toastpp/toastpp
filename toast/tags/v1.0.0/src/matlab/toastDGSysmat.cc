// =========================================================================
// toastSysmat
// Generate a system matrix from a qm-mesh and nodal images for
// mua, mus and n.
//
// RH parameters:
//     1: mesh handle (double)
//     2: mua NIM (double array)
//     3: mus NIM (double array)
//     4: n NIM (double array)
//     5: modulation frequency [MHz] (double)
//     6: (optional) 'EL' flag to indicate element basis
// LH parameters:
//     1: system matrix (sparse double matrix)
//     2: (optional) boundary part of system matrix (sparse)
//     3: (optional) boundary matrix prefactors (alpha)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "dgfwdsolver.h"
#include "util.h"

using namespace std;

// =========================================================================
// Implementation

void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "toastDGSysmat: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}

void CalcSysmat (NonconformingMesh *mesh, RVector &mua, RVector &mus, RVector &ref,
		 double freq, bool elbasis, mxArray **res)
{
    int n = (elbasis ? mesh->elen() : mesh->nlen());

    // Set optical coefficients
    Solution sol (OT_NPARAM, n);
    sol.SetParam (OT_CMUA, mua*c0/ref);
    sol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (int i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    mexPrintf("alpha is set to : %f\n", c2a[0]);
    sol.SetParam (OT_C2A, c2a);
    sol.SetParam (OT_N, ref);
    // Create forward solver to initialise system matrix
    CDGFwdSolver DGFWS(LSOLVER_DIRECT, 1e-10);
    //DGFWS.SetDataScaling (DATA_LOG);
    double omega = freq * 2.0*Pi*1e-6;
    cout << "calling allocate function ..." <<endl;
    DGFWS.Allocate (*mesh);
    cout<<"calling assemble system matrix ..."<<endl;
    DGFWS.AssembleSystemMatrix (sol, omega, elbasis);

    // Return system matrix to MATLAB
    CopyMatrix (res, *DGFWS.F);
}


// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    NonconformingMesh *mesh = (NonconformingMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    RVector mua, mus, ref;
    double freq;
    int n;
    bool elbasis = false;

    CopyVector (mua, prhs[1]);
    CopyVector (mus, prhs[2]);
    CopyVector (ref, prhs[3]);
    freq = mxGetScalar (prhs[4]);
    mexPrintf("recieved inputs ...\n");
    if (nrhs >= 6 && mxIsChar(prhs[5])) {
	char cbuf[32];
	mxGetString (prhs[5], cbuf, 32);
	elbasis = (strcasecmp (cbuf, "EL") == 0);
    }
    mexPrintf("elbasis: %d \n ", elbasis);
    if (elbasis) n = mesh->elen();
    else         n = mesh->nlen();
    mexPrintf("n: %d\n", n);
    Assert (n == mua.Dim(), "Argument 2: wrong size");
    Assert (n == mus.Dim(), "Argument 3: wrong size");
    Assert (n == ref.Dim(), "Argument 4: wrong size");
    mexPrintf("calling CalcSysmat function ...\n");
    CalcSysmat (mesh, mua, mus, ref, freq, elbasis, &plhs[0]);

}
