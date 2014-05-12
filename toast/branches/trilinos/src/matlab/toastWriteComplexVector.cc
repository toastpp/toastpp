#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include <iomanip>

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Writes a vector (complex) in TOAST format to a file
    // RH parameters:
    //    1: Vector file name (string)
    //    2: vector (complex double array)

    xAssert(nrhs >= 2, "Expected 2 input arguments");

    char dname[256];
    xAssert (mxGetString (prhs[0], dname, 256) == 0,
	     "Argument 1: string expected");

    xAssert (mxIsNumeric (prhs[1]),
	     "Argument 2: Expected numeric array");

    CVector data;
    CopyVector (data, prhs[1]);
    ofstream ofs(dname);
    ofs.precision(12);
    ofs.setf (ios::scientific);
    ofs << data;
}
