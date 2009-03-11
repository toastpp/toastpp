// =========================================================================
// toastMapSolToMesh
// Map a solution from the inverse basis to an FEM nodal solution
//
// RH parameters:
//     1: basis mapper handle
//     2: raster image (real or complex vector)
// LH parameters:
//     1: nodal image (real or complex vector)
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
    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int nn = raster->mesh().nlen();
    int sn = raster->SLen();

    // nodal image
    int m = mxGetM (prhs[1]);
    int n = mxGetN (prhs[1]);
    if (m*n != sn)
	mexErrMsgTxt ("Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (sn);
	CVector nim (nn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	complex *vptr = img.data_buffer();
	for (int i = 0; i < sn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_SolToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    } else {

	RVector img (sn, mxGetPr(prhs[1]));
	RVector nim (nn);
	raster->Map_SolToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    }
}
