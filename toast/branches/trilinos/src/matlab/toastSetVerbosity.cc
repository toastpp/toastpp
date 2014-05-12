// =========================================================================
// toastSetVerbosity
// Set global verbosity level
//
// RH-1: verbosity level (int >= 0)
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1)
        mexErrMsgTxt ("toastSetVerbosity: not enough parameters");
    int verbosity = (int)mxGetScalar (prhs[0]);
    toastVerbosity = verbosity;
}
