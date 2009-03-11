// =========================================================================
// toastRegulValue
// Obtain prior value for given image
//
// RH parameters:
//     1: hReg (handle) regularisation
//     2: x (double vector) solution data
//
// LH parameters:
//     1: prior (double) prior value
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

    double prior = (reg ? reg->GetValue (x) : 0.0);
    plhs[0] = mxCreateScalarDouble (prior);
}
