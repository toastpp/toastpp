// =========================================================================
// toastMPos
// Returns the detectors locations defined in a mesh
//
// RH parameters:
//     1: QM mesh handle
// LH parameters:
//     1: matrix of detector locations (nm x dim)
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

    int nm = mesh->nM;
    int dim = mesh->Dimension();
    
    mxArray *mpos = mxCreateDoubleMatrix (nm, dim, mxREAL);
    double *pr = mxGetPr (mpos);

    for (j = 0; j < dim; j++) {
	for (i = 0; i < nm; i++) {
	    *pr++ = mesh->M[i][j];
	}
    }

    plhs[0] = mpos;
}
