// =========================================================================
// toastGetQM
// Extracts the QM link list from a QMMesh and returns it in the form of
// a sparse permutation matrix
//
// RH parameters:
//     1: QM mesh handle (double)
// LH parameters:
//     1: link list permutation matrix (sparse real matrix)
// =========================================================================

#include "mex.h"
#include "toastmex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int q, m, n, idx;
    int nq = mesh->nQ;
    int nm = mesh->nM;
    int nqm = mesh->nQM;

    mxArray *lnk = mxCreateSparse (nm, nq, nqm, mxREAL);
    double  *pr = mxGetPr (lnk);
    mwIndex *ir = mxGetIr (lnk);
    mwIndex *jc = mxGetJc (lnk);

    *jc++ = 0;
    for (q = idx = 0; q < nq; q++) {
	for (m = n = 0; m < nm; m++) {
	    if (!mesh->Connected (q,m)) continue;
	    *pr++ = ++idx;
	    *ir++ = m;
	    n++;
	}
	*jc = *(jc-1) + n;
	jc++;
    }

    plhs[0] = lnk;
}
