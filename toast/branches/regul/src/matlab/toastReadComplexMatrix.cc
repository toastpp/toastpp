#include "mex.h"
#include "toastmex.h"
#include "mathlib.h"
#include "felib.h"

using namespace std;

void ReadMatrix (char *name, CCompRowMatrix &data)
{
    char c;

    ifstream ifs(name);
    ifs >> data;
}    

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Reads a TOAST sparse (complex) from a file
    // RH parameters:
    //     1: Vector file name (string)
    // LH parameters:
    //     1: vector (complex array)

    char dname[256];
    mxGetString (prhs[0], dname, 256);

    CCompRowMatrix data;
    ReadMatrix (dname, data);

    CopyMatrix (&plhs[0], data);
    return;
}
