// Interface:
// RH-1: FDOTFwd handle
// RH-2: field vector (nim image)
// LH-1: dense real matrix (imsize x nQ) of projections

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "FDOTFwd.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = (FDOTFwd*)Handle2Ptr (mxGetScalar(prhs[0]));
    
    // get field image (nim format)
    RVector field;
    CopyVector (field, prhs[1]);

    int nQ = fsolver->getNumSources();
    int n;
    double *pr;

    // calculate projections
    for (int i = 0; i < nQ; i++) {
        Projector *proj = fsolver->getProjector(i);
	RVector image;
	proj->projectFieldToImage (field, image);

	if (!i) {
  	    n = image.Dim();
  	    plhs[0] = mxCreateDoubleMatrix (n, nQ, mxREAL);
	    pr = mxGetPr(plhs[0]);
	}
	
	memcpy (pr, image.data_buffer(), n*sizeof(double));
	pr += n;
    }
}
