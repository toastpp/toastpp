// =========================================================================
// toastElementSize
// Returns the sizes (area or volume) of each mesh element
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: element sizes (real array)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //QMMesh *mesh = (QMMesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int n = mesh->elen();
    plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    double *pr = mxGetPr (plhs[0]);
    for (int i = 0; i < n; i++)
	pr[i] = mesh->ElSize(i);
}
