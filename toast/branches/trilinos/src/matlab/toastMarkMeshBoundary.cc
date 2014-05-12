// =========================================================================
// toastMarkMeshBoundary
// Flags boundary nodes either automatically from connectivity data, or
// using a supplied boundary type vector.
//
// RH parameters:
//     1: mesh handle
//     2: [optional]: boundary type vector
// LH parameters:
//     1: boundary type vector
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const char id[3] = {0,3,2};
    const double rid[4] = {0,-1,2,1};
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int i, n = mesh->nlen();
    RVector bnd(n);
    if (nrhs > 1) {
	CopyVector (bnd,prhs[1]);
	for (i = 0; i < n; i++)
	    mesh->nlist[i].SetBndTp (id[(int)(bnd[i]+0.5)]);
    } else {
	mesh->MarkBoundary();
    }
    if (nlhs > 0) {
	for (i = 0; i < n; i++)
	    bnd[i] = rid[mesh->nlist[i].BndTp()];
	CopyVector (&plhs[0], bnd);
    }
}
