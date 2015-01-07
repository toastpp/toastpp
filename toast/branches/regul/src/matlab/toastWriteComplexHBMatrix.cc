#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include <iomanip>

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Writes a vector (real) in TOAST format to a file
    // RH parameters:
    //     1: Vector file name (string)
    //     2: vector (double array)

    char dname[256];
    mxGetString (prhs[0], dname, 256);
    
    CCompRowMatrix data;
    CopyTMatrix (data, prhs[1]); // WARNING: this is transpose!
    ofstream ofs(dname);
    ofs.precision(12);
    ofs.setf (ios::scientific);
    data.Export (ofs);
}
