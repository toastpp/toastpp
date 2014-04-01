// =========================================================================
// toastQPos
// Returns the source locations defined in a mesh
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: matrix of source locations (nq x dim)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j;

    //int hMesh = (int)mxGetScalar (prhs[0]);
    //QMMesh *mesh = (QMMesh*)hMesh;
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nq = mesh->nQ;
    int dim = mesh->Dimension();
    
    mxArray *qpos = mxCreateDoubleMatrix (nq, dim, mxREAL);
    double *pr = mxGetPr (qpos);

    for (j = 0; j < dim; j++) {
	for (i = 0; i < nq; i++) {
	    *pr++ = mesh->Q[i][j];
	}
    }

    plhs[0] = qpos;
}
