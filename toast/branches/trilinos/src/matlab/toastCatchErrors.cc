// =========================================================================
// toastErrorFunc
// Reset TOAST error handler for Matlab
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

// =========================================================================
// Implementation

void MatlabErrorHandler (char *msg)
{
    mexErrMsgTxt (msg);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    SetErrorhandler (MatlabErrorHandler);
}
