// =========================================================================
// toastMeshData
// Mesh relaxation: moves mesh nodes towards the barycentre of their
// neighbours
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: shift magnitude (real, 0..1)
//     3: number of iterations (integer, >=1)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    double scale = mxGetScalar (prhs[1]);
    int iterations = (int)mxGetScalar (prhs[2]);

    JiggleMesh (mesh, scale, iterations);
}
