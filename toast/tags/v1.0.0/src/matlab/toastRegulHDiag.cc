// =========================================================================
// toastRegulHDiag
// Obtain diagonal of Hessian of regularisation for given solution x
//
// RH parameters:
//     1: hReg (handle) regularisation
//     2: x (double vector) solution data
//
// LH parameters:
//     1: diag (double vector) diagonal of 2nd derivative of regularisation
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "toastmex.h"
#include "regul.h"
#include "util.h"

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // regularisation
    //int hReg = (int)mxGetScalar (prhs[0]);
    //Regularisation *reg = (Regularisation*)hReg;
    Regularisation *reg = (Regularisation*)Handle2Ptr (mxGetScalar (prhs[0]));

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    RVector diag(x.Dim());
    if (reg) diag = reg->GetHessianDiag (x);
    CopyVector (&plhs[0], diag);
}
