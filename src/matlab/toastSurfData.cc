// =========================================================================
// toastSurfData
// Returns the node vertex and element lists of a TOAST mesh in MATLAB
// arrays.
//
// RH parameters:
//     1: mesh handle (pointer)
// LH parameters:
//     1: vertex coordinate list
//     2: element node index list
//     3: permutation index list to extract boundary data from full
//        nodal mesh vector
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k;
    double *pr;

    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    int nbnd = mesh->nlist.NumberOf (BND_ANY);
    int *bndidx = new int[nlen];

    // vertex coordinate list
    mxArray *vtx = mxCreateDoubleMatrix (nbnd, dim, mxREAL);
    pr = mxGetPr (vtx);

    for (i = 0; i < dim; i++)
	for (j = 0; j < nlen; j++)
	    if (mesh->nlist[j].isBnd())
		*pr++ = mesh->nlist[j][i];

    for (j = k = 0; j < nlen; j++)
	bndidx[j] = (mesh->nlist[j].isBnd() ? k++ : -1);
    
    // boundary element index list
    // note: this currently assumes that all elements contain the
    // same number of vertices!

    int nnd = 0, nface, sd, nn, nd, bn, *bndellist, *bndsdlist;
    nface = mesh->BoundaryList (&bndellist, &bndsdlist);
    for (j = 0; j < nface; j++)
	nnd = ::max (nnd, mesh->elist[bndellist[j]]->nSideNode(bndsdlist[j]));
    mxArray *idx = mxCreateDoubleMatrix (nface, nnd, mxREAL);
    pr = mxGetPr (idx);
    
    for (i = 0; i < nnd; i++)
	for (j = 0; j < nface; j++) {
	    Element *pel = mesh->elist[bndellist[j]];
	    sd = bndsdlist[j];
	    nn = pel->nSideNode (sd);
	    if (i < nn) {
		nd = pel->Node[pel->SideNode (sd, i)];
		bn = bndidx[nd]+1;
	    } else bn = 0;
	    *pr++ = bn;
	}

    // generate nodal permutation index list
    mxArray *perm = mxCreateDoubleMatrix (nbnd, 1, mxREAL);
    pr = mxGetPr (perm);
    for (i = 0; i < nlen; i++) {
	if (bndidx[i] >= 0)
	    *pr++ = i+1;
    }

    // cleanup
    delete []bndidx;

    plhs[0] = vtx;
    plhs[1] = idx;
    plhs[2] = perm;
}
