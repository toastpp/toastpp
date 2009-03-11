#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // De-allocates a mesh from dynamic heap.
    // RH parameters:
    //     1: mesh handle (pointer)

    //int hMesh = (int)mxGetScalar (prhs[0]);
    //QMMesh *mesh = (QMMesh*)hMesh;
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    delete mesh;
}
