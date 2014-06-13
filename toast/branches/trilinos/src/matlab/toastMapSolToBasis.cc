// =========================================================================
// toastMapSolToBasis
// Map a solution from the inverse solution basis to the regular solution
// grid of the full bounding box (including voxels outside mesh support)
//
// RH parameters:
//     1: basis mapper handle
//     2: sparse raster image (real or complex vector)
// LH parameters:
//     1: full raster image (real or complex vector)
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
    int sn = raster->SLen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    if ((int)(m*n) != sn)
	mexErrMsgTxt ("Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (sn);
	CVector bimg(bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	std::complex<double> *vptr = img.data_buffer();
	for (int i = 0; i < sn; i++)
	    vptr[i] = std::complex<double> (pr[i], pi[i]);
	raster->Map_SolToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    } else {

	RVector img (sn, mxGetPr(prhs[1]));
	RVector bimg(bn);
	raster->Map_SolToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    }
}
