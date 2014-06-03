#include "mex.h"
#include "felib.h"
#include "util.h"
#include <iostream>
#include <math.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, res, len, *perm, optmode = 0;

    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    optmode = (int)mxGetScalar (prhs[1]);

    len = mesh->nlen();
    perm = new int[len];
    for (i = 0; i < len; i++) perm[i] = i;

    switch (optmode) {
    case 0:
	res = Optimise_MMD (*mesh, perm, 0, len);
	break;
    case 1:
	res = Optimise_MinBandwidth (*mesh, perm, 0, len);
	break;
    case 2:
	res = Optimise_Tinney2 (*mesh, perm, 0, len);
	break;
    }

    if (res)
	mexErrMsgTxt("Optimisation failed");

    plhs[0] = mxCreateDoubleMatrix (len,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < len; i++)
	*pr++ = perm[i]+1; // switch to 1-based

    delete []perm;
}
