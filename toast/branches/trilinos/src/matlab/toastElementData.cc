// =========================================================================
// toastElementData
// Returns the node coordinates and element vertex index lists for each
// face of a mesh element
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: element index (>= 1)
// LH parameters:
//     1: element vertices: dense n x d matrix,
//     2: element index list (one row for each face)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j;
    double *pr;

    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    // element index
    int el = (int)mxGetScalar (prhs[1]) - 1;
    if (el < 0 || el >= mesh->elen()) { // invalid index
	plhs[0] = mxCreateDoubleMatrix (0, 0, mxREAL);
	return;
    }

    Element *pel = mesh->elist[el];
    int dim = pel->Dimension();
    int nnd = pel->nNode();
    int nsd = pel->nSide();
    int nsn = pel->nSideNode(0);
    for (i = 1; i < nsd; i++)
	nsn = std::max (nsn, pel->nSideNode(i));

    mxArray *vtx = mxCreateDoubleMatrix (nnd, dim, mxREAL);
    pr = mxGetPr (vtx);
    for (i = 0; i < dim; i++)
	for (j = 0; j < nnd; j++)
	    *pr++ = mesh->nlist[pel->Node[j]][i];

    mxArray *idx = mxCreateDoubleMatrix (nsd, nsn, mxREAL);
    pr = mxGetPr (idx);
    for (i = 0; i < nsn; i++)
	for (j = 0; j < nsd; j++)
	    if (i < pel->nSideNode(j))
		*pr++ = pel->Node[pel->SideNode(j,i)]+1;
	    else
		*pr++ = 0;
    
    plhs[0] = vtx;
    plhs[1] = idx;
}
