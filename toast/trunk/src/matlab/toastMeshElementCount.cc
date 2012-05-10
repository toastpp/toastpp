// =========================================================================
// toastMeshElementCount
// Returns the element count of a TOAST mesh
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: element count (integer)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    plhs[0] = mxCreateDoubleScalar (mesh->elen());
}
