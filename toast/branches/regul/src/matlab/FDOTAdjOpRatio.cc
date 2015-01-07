// Interface:
// RH-1: FDOTFwd handle
// RH-2: data vector
// RH-3: Epsilon value for ratio = fluo / (excit + epsilon)
// RH-4: adjoint source index (>= 1) (optional)
// LH-1: paramter vector

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "FDOTFwd.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */
    // get FDOTFwd pointer from handle

    if (nrhs < 3)
    {
        mexErrMsgTxt ("FDOTAdjOpRatio: Invalid function parameters");
	return;
    }

    FDOTFwd *fsolver = (FDOTFwd*)Handle2Ptr (mxGetScalar(prhs[0]));
    
    // get data vector
    RVector y;
    CopyVector (y, prhs[1]);

    // single-source index
    int q = -1;
    if (nrhs >= 4) {
        q = (int)(mxGetScalar(prhs[3]));
    }

    double eps = mxGetScalar(prhs[2]);

    // calculate parameter vector
    RVector x;
    if (q >= 0) fsolver->adjOperator (y, x, q, true, eps);
    else        fsolver->adjOperator (y, x, true, eps);
    CopyVector (&plhs[0], x);
    /*
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
    */
}
