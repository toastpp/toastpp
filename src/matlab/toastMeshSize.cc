// =========================================================================
// toastMeshSize
// Returns the size (area or volume) of a TOAST mesh.
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: mesh size (real)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //QMMesh *mesh = (QMMesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    plhs[0] = mxCreateScalarDouble (mesh->FullSize());
}
