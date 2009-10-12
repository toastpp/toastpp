// Interface:
// RH-1: hProjList list of projector handles

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "projector.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get list pf projector handles
    RVector dblList;
    CopyVector(dblList, prhs[0]);
    for (int i=0; i<dblList.Dim(); ++i)
    {
	Projector * proj = (Projector*) Handle2Ptr (dblList[i]);
	delete proj;
    }
}
