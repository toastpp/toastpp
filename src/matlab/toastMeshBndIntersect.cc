// =========================================================================
// toastMeshBndIntersect
// Calculates the intersection of a line with the mesh boundary
//
// RH parameters:
//     1: mesh handle
//     2: internal point (2D or 3D)
//     3: external point (2D or 3D)
// LH parameters:
//     1: intersection point
// =========================================================================

#include "mex.h"
#include "toastmex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int dim = mesh->Dimension();
    Point p1(dim), p2(dim);
    CopyVector (p1, prhs[1]);
    CopyVector (p2, prhs[2]);

    Point s = mesh->BndIntersect (p1, p2);
    CopyVector (&plhs[0], s);
}
