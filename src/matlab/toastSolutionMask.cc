// =========================================================================
// toastSolutionMask
// Returns a permutation vector which extracts the coefficients from a
// basis vector which have support over the domain.
//
// RH parameters:
//     1: basis mapper handle
//
// LH parameters:
//     1: permutation vector
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int slen = raster->SLen();

    plhs[0] = mxCreateDoubleMatrix (1, slen, mxREAL);
    double *pr = mxGetPr (plhs[0]);

    for (int i = 0; i < slen; i++)
	pr[i] = raster->Sol2Basis (i) + 1;
        // "+1" to adjust to matlab's 1-based array indices
}
