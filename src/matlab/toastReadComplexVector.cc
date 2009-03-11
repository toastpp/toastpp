#include "mex.h"
#include "mathlib.h"
#include "felib.h"

using namespace std;
using namespace toast;

bool ReadVector (char *name, CVector &data, int idx)
{
    char c;

    ifstream ifs(name);

    if (!idx) { // read last vector
	CVector tmp;
	ifs >> tmp;
	while (ifs.good ()) {
	    data.Copy (tmp);
	    ifs >> tmp;
	}
    } else {
	for (int i = 0; i < idx; i++) {
	    ifs >> data;
	    if (!ifs.good()) {
		cerr << "toastReadComplexVector: index out of range" << endl;
		data.New(0);
		return false;
	    }
	}
    }
    return true;
}    

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Reads a TOAST vector (complex) from a file
    // RH parameters:
    //     1: Vector file name (string)
    //     2: vector index (int, optional, default=1)
    // LH parameters:
    //     1: vector (complex array)

    int idx = 1;
    char dname[256];
    mxGetString (prhs[0], dname, 256);

    if (nrhs > 1) {
	double v = mxGetScalar (prhs[1]);
	idx = (int)(v+0.5);
    }

    CVector data;
    if (idx >= 0)
	ReadVector (dname, data, idx);

    plhs[0] = mxCreateDoubleMatrix (data.Dim(), 1, mxCOMPLEX);
    double *pr = mxGetPr (plhs[0]);
    double *pi = mxGetPi (plhs[0]);

    complex *buf = data.data_buffer();
    for (int i = 0; i < data.Dim(); i++) {
	pr[i] = buf[i].re;
	pi[i] = buf[i].im;
    }
}
