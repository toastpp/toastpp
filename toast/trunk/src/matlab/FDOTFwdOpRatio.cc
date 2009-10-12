// Interface:
// RH-1: FDOTFwd handle
// RH-2: parameter vector
// RH-3: Epsilon value for ratio = fluo / (excit + epsilon)
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
    if (nrhs < 3)
    {
        mexErrMsgTxt ("FDOTFwdOpRatio: Invalid function parameters");
	return;
    }

    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = (FDOTFwd*)Handle2Ptr (mxGetScalar(prhs[0]));
    
    // get x vector
    RVector x;
    CopyVector (x, prhs[1]);

    double eps = mxGetScalar(prhs[2]);

    // calculate data vector
    RVector y;
    fsolver->fwdOperator (x, y, true, eps);
    CopyVector (&plhs[0], y);
}
