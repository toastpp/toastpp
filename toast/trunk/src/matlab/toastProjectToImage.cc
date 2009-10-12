// Interface:
// RH-1: Projector handle
// RH-2: field vector (nim image)
// LH-1: dense real matrix (imsize x nQ) of projections

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "projector.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get FluoSolver pointer from handle
    Projector * proj = (Projector*)Handle2Ptr (mxGetScalar(prhs[0]));
    
    // get field image (nim format)
    RVector field;
    CopyVector (field, prhs[1]);

    double *pr;

    // calculate projections
    RVector image;
    proj->projectFieldToImage (field, image);

    int n = image.Dim();
    plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    pr = mxGetPr(plhs[0]);

    memcpy (pr, image.data_buffer(), n*sizeof(double));
}
