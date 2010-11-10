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
// =========================================================================

#include "mex.h"
#include "nonconformingMesh.h"
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
    NonconformingMesh *mesh = (NonconformingMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    int iedgelen = mesh->iedgelen();
    int bedgelen = mesh->bedgelen();

    // vertex coordinate list
    mxArray *vtx = mxCreateDoubleMatrix (nlen, dim, mxREAL);
    pr = mxGetPr (vtx);

    for (i = 0; i < dim; i++)
	for (j = 0; j < nlen; j++)
	    *pr++ = mesh->nlist[j][i];

    // max number of nodes per element
    int nnd = mesh->elist[0]->nNode();
    for (i = 0; i < elen; i++)
	nnd = ::max (nnd, mesh->elist[i]->nNode());

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

    mxArray *iedge_elist = mxCreateDoubleMatrix(iedgelen, 2, mxREAL);
    pr = mxGetPr(iedge_elist);
    for(i = 0;  i < 2; i++)
	for(j=0; j < iedgelen; j++)
		*pr++ = mesh->iedge_elist[j][i]+1;

    mxArray *iedge_nlist = mxCreateDoubleMatrix(iedgelen, 3, mxREAL);
    pr = mxGetPr(iedge_nlist);
    for(i = 0; i < 3; i++)
	for(j=0; j < iedgelen; j++)
		*pr++ = mesh->iedge_nlist[j][i]+1;

   mxArray *iedge_state = mxCreateDoubleMatrix(iedgelen, 1, mxREAL);
   pr = mxGetPr(iedge_state);
   for(i = 0; i < iedgelen; i++)
		*pr++ = mesh->iedge_state[i];

    mxArray *bedge_elist = mxCreateDoubleMatrix(bedgelen, 1, mxREAL);
    pr = mxGetPr(bedge_elist);
    for(i = 0; i < bedgelen; i++)
		*pr++ = mesh->bedge_elist[i]+1;

    mxArray *bedge_nlist = mxCreateDoubleMatrix(bedgelen, 3, mxREAL);
    pr = mxGetPr(bedge_nlist);
    for(i = 0; i < 3; i++)
	for(j=0; j < bedgelen; j++)
		*pr++ = mesh->bedge_nlist[j][i]+1;

   mxArray *bedge_state = mxCreateDoubleMatrix(bedgelen, 1, mxREAL);
   pr = mxGetPr(bedge_state);
   for(i = 0; i < bedgelen; i++)
		*pr++ = mesh->bedge_state[i];


    plhs[0] = vtx;
    plhs[1] = idx;
    plhs[2] = iedge_elist;
    plhs[3] = iedge_nlist;
    plhs[4] = iedge_state;
    plhs[5] = bedge_elist;
    plhs[6] = bedge_nlist;
    plhs[7] = bedge_state;
}
