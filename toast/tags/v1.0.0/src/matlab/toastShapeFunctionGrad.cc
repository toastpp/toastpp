// =========================================================================
// toastShapeFunctionGrad
// Returns a matrix of shape function gradient values at a given global
// coordinate for a given mesh element
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: element index (1-based, integer)
//     3: global point coordinates (array 2D or 3D)
// LH parameters:
//     1: matrix of shape functions of dimension d x n (number of nodes in
//        the element, d: dimension of the problem)
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

    // calculate shape function derivatives
    RDenseMatrix fgrad = mesh->elist[idx]->GlobalShapeD (mesh->nlist, glob);
    int nn = fgrad.nCols();

    // vertex coordinate list
    mxArray *sf = mxCreateDoubleMatrix (dim, nn, mxREAL);
    CopyMatrix (&sf, fgrad);

    plhs[0] = sf;
}
