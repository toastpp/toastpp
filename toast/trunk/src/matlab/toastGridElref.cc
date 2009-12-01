// =========================================================================
// toastGridElref
// Returns list of element indices for each grid point in the basis.
//
// RH parameters:
//     1: basis mapper handle
// LH parameters:
//     1: element reference vector
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

using namespace toast;

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int i, n = raster->GLen();
    int *elref = raster->Elref();

    plhs[0] = mxCreateDoubleMatrix (n,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < n; i++) {
      pr[i] = elref[i] + 1;
        // make 1-based. Gridpoints outside mesh support are set to zero
    }
}
