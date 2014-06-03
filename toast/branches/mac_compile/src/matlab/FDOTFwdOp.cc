// Interface:
// RH-1: FDOTFwd handle
// RH-2: parameter vector
// LH-1: data vector

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "FDOTFwd.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = (FDOTFwd*)Handle2Ptr (mxGetScalar(prhs[0]));
    
    // get x vector
    RVector x;
    CopyVector (x, prhs[1]);

    // calculate data vector
    RVector y;
    fsolver->fwdOperator (x, y, false);
    CopyVector (&plhs[0], y);
}
