#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char cbuf[256];
    double refind = mxGetScalar (prhs[0]);
    mxGetString (prhs[1], cbuf, 256);

    if (!strcasecmp (cbuf, "Keijzer")) {
	plhs[0] = mxCreateDoubleScalar(A_Keijzer(refind));
    } else if (!strcasecmp (cbuf, "Contini")) {
	plhs[0] = mxCreateDoubleScalar(A_Contini(refind));
    } else {
	plhs[0] = mxCreateDoubleScalar(1.0);
    }
}
