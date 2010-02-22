#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    RDenseMatrix qvec, mvec;
    RCompRowMatrix lnk;
    char qmname[256];
    int i, j;

    CopyMatrix (qvec, prhs[0]);
    CopyMatrix (mvec, prhs[1]);

    if (!mxIsSparse(prhs[2])) {
	cerr << "Error: toastWriteQM: argument 3: expected sparse matrix\n";
	cerr << "Aborting" << endl;
	return;
    }
    CopyTMatrix (lnk, prhs[2]);

    int dim = qvec.nCols();
    int nq = qvec.nRows();
    int nm = mvec.nRows();
    
    if (nq != lnk.nRows() || nm != lnk.nCols()) {
	cerr << "Error: toastWriteQM: invalid dimension for argument 3:\n";
	cerr << "Was: " << lnk.nCols() << " x " << lnk.nRows() << endl;
	cerr << "Expected: " << nm << " x " << nq << endl;
	cerr << "Aborting" << endl;
	return;
    }

    mxGetString (prhs[3], qmname, 256);

    ofstream ofs(qmname);
    ofs << "QM file" << endl;
    ofs << "Dimension " << dim << endl << endl;
    ofs << "SourceList " << nq << endl;
    for (i = 0; i < nq; i++)
	for (j = 0; j < dim; j++)
	    ofs << qvec(i,j) << (j == dim-1 ? '\n':' ');
    ofs << endl;
    ofs << "MeasurementList " << nm << endl;
    for (i = 0; i < nm; i++)
	for (j = 0; j < dim; j++)
	    ofs << mvec(i,j) << (j == dim-1 ? '\n':' ');
    ofs << endl;

    ofs << "LinkList" << endl;
    int *ci = new int[nm];
    double *idx = new double[nm];

    for (i = 0; i < nq; i++) {
	int nz = lnk.SparseRow (i, ci, idx);
	ofs << nz << ':';
	for (j = 0; j < nz; j++) {
	    if (ci[j])
		ofs << ' ' << ci[j];
	}
	ofs << endl;
    }
}
