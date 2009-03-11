// =========================================================================
// toastRegulGradient
// Obtain gradient of regularisation w.r.t. image parameters
//
// RH parameters:
//     1: hReg (handle) regularisation
//     2: x (double vector) solution data
//
// LH parameters:
//     1: grad (double vector) regularisation gradient
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

    RVector grad(x.Dim());
    if (reg) grad = reg->GetGradient (x);

    CopyVector (&plhs[0], grad);
}
