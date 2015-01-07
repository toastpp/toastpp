// =========================================================================
// toastRegulHess1f
// Applies 1st order Hessian L on a vector f at given coefficient vector x
// and returns resulting vector.
//
// RH parameters:
//     1: hReg (handle) regularisation instance handle
//     2: x (double vector) coefficient vector
//     3: f (double vector) operand
//
// LH parameters:
//     1: y (double vector) result of applying L(x) to f
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
    Regularisation *reg = (Regularisation*)Handle2Ptr (mxGetScalar (prhs[0]));

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    // operand vector
    RVector f;
    CopyVector (f, prhs[2]);

    RVector y = reg->GetHess1f (x, f);
    CopyVector (&plhs[0], y);
}
