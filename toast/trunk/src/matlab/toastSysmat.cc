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

void CalcSysmat (QMMesh *mesh, RVector &mua, RVector &mus, RVector &ref,
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
    sol.SetParam (OT_C2A, c2a);

    // Create forward solver to initialise system matrix
    CFwdSolver FWS (mesh, LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    double omega = freq * 2.0*Pi*1e-6;
    
    FWS.Allocate ();
    FWS.AssembleSystemMatrix (sol, omega, elbasis);

    // Return system matrix to MATLAB
    CopyMatrix (res, *FWS.F);
}

void CalcBndSysmat (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    CFwdSolver FWS (mesh, LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    FWS.Allocate ();
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

    RVector mua, mus, ref;
    double freq;
    int n;
    bool elbasis = false;

    CopyVector (mua, prhs[1]);
    CopyVector (mus, prhs[2]);
    CopyVector (ref, prhs[3]);
    freq = mxGetScalar (prhs[4]);

    if (nrhs >= 6 && mxIsChar(prhs[5])) {
	char cbuf[32];
	mxGetString (prhs[5], cbuf, 32);
	elbasis = (strcasecmp (cbuf, "EL") == 0);
    }

    if (elbasis) n = mesh->elen();
    else         n = mesh->nlen();
    Assert (n == mua.Dim(), "Argument 2: wrong size");
    Assert (n == mus.Dim(), "Argument 3: wrong size");
    Assert (n == ref.Dim(), "Argument 4: wrong size");

    CalcSysmat (mesh, mua, mus, ref, freq, elbasis, &plhs[0]);


    if (nlhs > 1) { // separately provide boundary integral matrix
	CalcBndSysmat (mesh, ref, &plhs[1]);

	if (nlhs > 2) { // boundary matrix prefactors
	    CalcBndFactors (mesh, ref, &plhs[2]);
	}
    }
}
