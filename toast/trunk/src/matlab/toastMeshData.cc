// =========================================================================
// toastMeshData
// Returns the node vertex and element lists of a TOAST mesh in MATLAB
// arrays.
//
// RH parameters:
//     1: mesh handle (pointer)
// LH parameters:
//     1: vertex coordinate list
//     2: element node index list
//     3: element type list
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j;
    double *pr;

    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // vertex coordinate list
    mxArray *vtx = mxCreateDoubleMatrix (nlen, dim, mxREAL);
    pr = mxGetPr (vtx);

    for (i = 0; i < dim; i++)
	for (j = 0; j < nlen; j++)
	    *pr++ = mesh->nlist[j][i];

    // max number of nodes per element
    int nnd = mesh->elist[0]->nNode();
    for (i = 0; i < elen; i++)
	nnd = std::max (nnd, mesh->elist[i]->nNode());

    // element index list
    // (1-based; value 0 indicates unused matrix entry)
    mxArray *idx = mxCreateDoubleMatrix (elen, nnd, mxREAL);
    pr = mxGetPr (idx);
    
    for (i = 0; i < nnd; i++)
	for (j = 0; j < elen; j++)
	    if (i < mesh->elist[j]->nNode())
		*pr++ = mesh->elist[j]->Node[i]+1;
	    else
		*pr++ = 0;

    plhs[0] = vtx;
    plhs[1] = idx;

    if (nrhs >= 3) {
	mxArray *eltp = mxCreateDoubleMatrix (elen, 1, mxREAL);
	pr = mxGetPr (eltp);
	for (j = 0; j < elen; j++) {
	    *pr++ = mesh->elist[j]->Type();
	}
	plhs[2] = eltp;
    }
}
