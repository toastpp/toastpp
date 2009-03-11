// =========================================================================
// toastMapGridToBasis
// Map a solution from the fine grid to the regular solution
// grid of the full bounding box
//
// RH parameters:
//     1: basis mapper handle
//     2: image in fine basis grid (real or complex vector)
// LH parameters:
//     1: image in coarse basis grid (real or complex vector)
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
    int gn = raster->GLen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    if ((int)(m*n) != gn)
	mexErrMsgTxt ("Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector bimg(bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	complex *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_GridToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector bimg(bn);
	raster->Map_GridToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    }
}
