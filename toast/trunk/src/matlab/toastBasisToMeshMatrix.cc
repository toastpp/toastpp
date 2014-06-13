// =========================================================================
// toastBasisToMeshMatrix
// Returns the transformation matrix for mapping from image (b) to mesh (h)
// basis.
//
// RH parameters:
//     1: basis mapper handle
// LH parameters:
//     1: mapping matrix (sparse real)
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

using namespace std;

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));

    // copy mapping matrix
    RCompRowMatrix &map = 
	(RCompRowMatrix&)raster->Basis2MeshMatrix(); // dodgy cast

    CopyMatrix (&plhs[0], map);
}
