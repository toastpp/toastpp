#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "regul.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // De-allocates a regularisation instance from dynamic heap.
    // RH parameters:
    //     1: regularisation instance handle (pointer)

    Regularisation *regul = (Regularisation*)Handle2Ptr(mxGetScalar (prhs[0]));
    delete regul;
}
