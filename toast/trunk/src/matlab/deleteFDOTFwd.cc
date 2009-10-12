// Interface:
// RH-1: hFDOTFwd handle 

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "FDOTFwd.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get FluoSolver pointer from handle
    FDOTFwd *fdot = (FDOTFwd*) Handle2Ptr (mxGetScalar(prhs[0]));
    delete fdot;
}
