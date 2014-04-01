// =========================================================================
// toastWriteMesh
// Writes a TOAST mesh to a file
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: output mesh file name (string)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    char meshname[256];
    mxGetString (prhs[1], meshname, 256);

    ofstream ofs (meshname);
    ofs << *mesh;

    if (toastVerbosity > 0)
        mexPrintf("Mesh: written to %s\n", meshname);
}
