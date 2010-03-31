// =========================================================================
// toastMassmat
// Generate a mass matrix from a mesh. The mass matrix is defined as the
// element-wise integral of the product of two shape functions.
// This is independent of optical parameters. 
//
// RH parameters:
//     1: mesh handle
// LH parameters:
//     1: mass matrix (sparse double matrix)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

// =========================================================================
// Implementation

void CalcMassmat (QMMesh *mesh, mxArray **res)
{
    int n = mesh->nlen();

    // Create forward solver to initialise system matrix
    RFwdSolver FWS (LSOLVER_ITERATIVE, 1e-10);
    FWS.AssembleMassMatrix (mesh);

    // Return system matrix to MATLAB
    CopyMatrix (res, *FWS.B);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    CalcMassmat (mesh, &plhs[0]);
}
