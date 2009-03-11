#include "mex.h"
#include "stoastlib.h"
#include "util.h"
#include "toastmex.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Returns the grid dimensions of the basis mapper instance
    // RH parameters:
    //     1: basis mapper handle (pointer)
    // LH parameters:
    //     1: bdim (integer array of size 2 or 3)
    //     2: gdim - optional (integer array of size 2 or 3)

    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    int i, dim = raster->Dim();
    if (nlhs >= 1) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->BDim()[i];
	CopyVector (&plhs[0], tmp);
    }
    if (nlhs >= 2) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->GDim()[i];
	CopyVector (&plhs[1], tmp);
    }
}
