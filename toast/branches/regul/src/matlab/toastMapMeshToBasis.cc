// =========================================================================
// toastMapMeshToBasis
// Map a nodal image to the regular solution basis defined by the mapper
//
// RH parameters:
//     1: basis mapper handle
//     2: nodal image (real or complex vector)
// LH parameters:
//     1: raster image (real or complex vector)
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
    //int hRaster = (int)mxGetScalar (prhs[0]);
    //Raster *raster = (Raster*)hRaster;
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int nn = raster->mesh().nlen();
    int bn = raster->BLen();

    // nodal image
    int m = mxGetM (prhs[1]);
    int n = mxGetN (prhs[1]);
    if (m*n != nn)
	mexErrMsgTxt ("Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector nim (nn);
	CVector img (bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	std::complex<double> *vptr = nim.data_buffer();
	for (int i = 0; i < nn; i++)
	    vptr[i] = std::complex<double> (pr[i], pi[i]);
	raster->Map_MeshToBasis (nim, img);
	CopyVector (&plhs[0], img);

    } else {

	RVector nim (nn, mxGetPr (prhs[1]));
	RVector img (bn);
	raster->Map_MeshToBasis (nim, img);
	CopyVector (&plhs[0], img);

    }
}
