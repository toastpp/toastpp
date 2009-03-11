// =========================================================================
// toastDataLinkList
// Returns a permutation vector which allows extraction of the active
// subset of measurements from a full measurement vector.
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: permutation list (integer vector)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar(prhs[0]));

    int q, m, idx;
    plhs[0] = mxCreateDoubleMatrix (1, mesh->nQM, mxREAL);
    double *pr = mxGetPr (plhs[0]);
    for (q = idx = 0; q < mesh->nQ; q++)
	for (m = 0; m < mesh->nM; m++)
	    if (mesh->Connected (q, m))
		pr[idx++] = q*mesh->nM + m + 1;
}
