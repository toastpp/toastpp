// =========================================================================
// toastMapGridToSol
// Map a solution from the fine grid to the inverse solution basis
// (excluding pixels outside the support of the mesh)
//
// RH parameters:
//     1: basis mapper handle
//     2: image in fine basis grid (real or complex vector)
// LH parameters:
//     1: image in inverse solution basis (real or complex vector)
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
    int gn = raster->GLen();
    int sn = raster->SLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    if ((int)(m*n) != gn)
	mexErrMsgTxt ("Invalid image dimensions");


    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector simg(sn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	std::complex<double> *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++)
	    vptr[i] = std::complex<double> (pr[i], pi[i]);
	raster->Map_GridToSol (img, simg);
	CopyVector (&plhs[0], simg);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector simg(sn);
	raster->Map_GridToSol (img, simg);
	CopyVector (&plhs[0], simg);

    }
}
