#include "mex.h"
#include "mathlib.h"
#include "felib.h"

using namespace std;

bool ReadVector (char *name, RVector &data, int idx)
{
    char c;

    ifstream ifs(name);

    if (!idx) { // read last vector
	RVector tmp;
	ifs >> tmp;
	while (ifs.good ()) {
	    data.Copy (tmp);
	    ifs >> tmp;
	}
    } else {
	for (int i = 0; i < idx; i++) {
	    ifs >> data;
	    if (!ifs.good()) {
		cerr << "toastReadRealVector: index out of range" << endl;
		data.New(0);
		return false;
	    }
	}
    }
    return true;
}    

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Reads a TOAST vector (real) from a file
    // RH parameters:
    //     1: Vector file name (string)
    //     2: vector index (int, optional, default=1)
    // LH parameters:
    //     1: vector (real array)

    int idx = 1;
    char dname[256];
    mxGetString (prhs[0], dname, 256);

    if (nrhs > 1) {
	double v = mxGetScalar (prhs[1]);
	idx = (int)(v+0.5);
    }

    RVector data;
    if (idx >= 0)
	ReadVector (dname, data, idx);

    plhs[0] = mxCreateDoubleMatrix (data.Dim(), 1, mxREAL);
    double *pr = mxGetPr (plhs[0]);

    double *buf = data.data_buffer();
    for (int i = 0; i < data.Dim(); i++) {
	pr[i] = buf[i];
    }
}
