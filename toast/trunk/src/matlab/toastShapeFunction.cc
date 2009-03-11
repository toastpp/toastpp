// =========================================================================
// toastShapeFunction
// Returns an array of shape function values at a given global coordinate
// for a given mesh element
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: element index (1-based, integer)
//     3: global point coordinates (array 2D or 3D)
// LH parameters:
//     1: array of shape functions of length n (number of nodes in the element)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j;
    double *pr;

    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // element index
    int idx = (int)(floor(mxGetScalar (prhs[1])-0.5));
    if (idx < 0 || idx >= elen) {
	mexPrintf ("Invalid element index\n");
	return;
    }

    // global point
    Point glob(dim);
    const mwSize *gdim = mxGetDimensions (prhs[2]);
    if (gdim[0]*gdim[1] != dim) {
	mexPrintf ("Invalid point dimension\n");
	return;
    }
    pr = mxGetPr (prhs[2]);
    for (i = 0; i < dim; i++) glob[i] = pr[i];

    // calculate shape functions
    RVector fun = mesh->elist[idx]->GlobalShapeF (mesh->nlist, glob);
    int nn = fun.Dim();

    // vertex coordinate list
    mxArray *sf = mxCreateDoubleMatrix (1, nn, mxREAL);
    pr = mxGetPr (sf);
    for (i = 0; i < nn; i++) pr[i] = fun[i];

    plhs[0] = sf;
}
