// =========================================================================
// toastSetThreadCount
// Set the number of threads to be used by various toast functions (e.g.
// forward solver)
//
// RH parameters:
//     1: number of threads
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "task.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nth = (int)mxGetScalar (prhs[0]);
    Task_Init (nth);
}
