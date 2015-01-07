#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    double *pr = mxGetPr(prhs[1]);

    int nlen = mesh->nlen();
    IVector perm(nlen);
    for (int i = 0; i < nlen; i++)
	perm[i] = (int)(pr[i]-0.5); // switch to 0-based
    
    mesh->Reorder(perm);
}
