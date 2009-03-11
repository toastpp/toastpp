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
    //     1: NIM file name (string)
    //     2: mesh name (string)
    //     3: vector (double array)

    char dname[256];
    mxGetString (prhs[0], dname, 256);
    ofstream ofs (dname);
    ofs << "NIM" << endl;

    mxGetString (prhs[1], dname, 256);
    ofs << "Mesh = " << dname << endl;

    ofs << "SolutionType = N/A" << endl;
    
    int n = mxGetN(prhs[2]) * mxGetM(prhs[2]);
    ofs << "ImageSize = " << n << endl;

    ofs << "EndHeader" << endl;

    ofs << "Image 0" << endl;

    double *val = mxGetPr (prhs[2]);
    ofs.precision(12);
    ofs.setf (ios::scientific);
    for (int i = 0; i < n; i++)
	ofs << val[i] << ' ';
    ofs << endl;
}
    
