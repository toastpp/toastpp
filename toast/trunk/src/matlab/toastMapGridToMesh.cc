// =========================================================================
// toastMapGridToMesh
// Map a solution from the fine grid basis to an FEM nodal solution
//
// RH parameters:
//     1: basis mapper handle
//     2: raster image in fine basis (real or complex vector)
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
    //int hRaster = (int)mxGetScalar (prhs[0]);
    //Raster *raster = (Raster*)hRaster;
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int nn = raster->mesh().nlen();
    int gn = raster->GLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    if ((int)(m*n) != gn)
	mexErrMsgTxt ("Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector nim (nn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	complex *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_GridToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector nim (nn);
	raster->Map_BasisToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    }
}
