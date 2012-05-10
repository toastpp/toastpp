// =========================================================================
// toastFindElement
// Returns the index of the element containing a point
//
// RH parameters:
//     1: mesh handle
//     2: point (2D or 3D array)
// LH parameters:
//     1: element index (>=1, or 0 if point is outside mesh)
// =========================================================================

#include "mex.h"
#include "toastmex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //QMMesh *mesh = (QMMesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int dim = mesh->Dimension();
    Point pt(dim);
    CopyVector (pt, prhs[1]);
    //cerr << pt[0] << ' ' << pt[1] << ' ' << pt[2] << endl;
    plhs[0] = mxCreateDoubleScalar (mesh->ElFind(pt)+1);
}
