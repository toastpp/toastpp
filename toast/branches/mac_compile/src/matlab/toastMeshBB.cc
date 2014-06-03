// =========================================================================
// toastMeshBB
// Returns the bounding box of a mesh
//
// RH parameters:
//     1: mesh handle (pointer)
// LH parameters:
//     1: pmin (min corner of the mesh bounding box)
//     2: pmax (max corner of the mesh bounding box)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // mesh handle
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int dim = mesh->Dimension();
    Point pmin(dim), pmax(dim);
    mesh->BoundingBox (pmin, pmax);

    CopyVector (plhs+0, pmin);
    CopyVector (plhs+1, pmax);
}
