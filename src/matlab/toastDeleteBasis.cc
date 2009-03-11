#include "mex.h"
#include "stoastlib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // De-allocates a basis mapper from dynamic heap.
    // RH parameters:
    //     1: basis mapper handle (pointer)

    //int hRaster = (int)mxGetScalar (prhs[0]);
    //Raster *raster = (Raster*)hRaster;
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    delete raster;
}
