// =========================================================================
// toastRegulHess
// Returns the Hessian matrix of the regularisation for a given parameter
// distribution.
//
// RH parameters:
//     1: hReg (handle) regularisation
//     2: x (double vector) solution data
//
// LH parameters:
//     1: H (sparse double matrix) regularisation Hessian
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "toastmex.h"
#include "regul.h"
#include "util.h"

// ==========================================================================
// Build the Hessian of the prior from individual parameter contributions

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x)
{
    int nprm = reg->GetNParam();
    int n = x.Dim();
    int n0 = n/nprm;
    int i, j;

    RCompRowMatrix H, Hi, Hij;
    for (i = 0; i < nprm; i++) {
	for (j = 0; j < nprm; j++) {
	    if (!j) {
		Hi.New(n0,n0);
		if (j==i) reg->SetHess1 (Hi, x, j);
	    } else {
		Hij.New(n0,n0);
		if (j==i) reg->SetHess1 (Hij, x, j);
		Hi = cath (Hi, Hij);
	    }
	}
	if (!i) H = Hi;
	else    H = catv (H, Hi);
    }
    return H;
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = (Regularisation*)Handle2Ptr (mxGetScalar (prhs[0]));

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);
    RCompRowMatrix H (BuildRHessian (reg, x));
    CopyMatrix (&plhs[0], H);

#ifdef UNDEF
    int slen = x.Dim();
    RCompRowMatrix H (slen, slen);
    reg->SetHess1 (H, x, 0);

    CopyMatrix (&plhs[0], H);
#endif
}
