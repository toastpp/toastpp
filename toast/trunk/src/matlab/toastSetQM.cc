// =========================================================================
// toastSetQM
// Assigns a set of source and detector locations to a TOAST QMMesh.
//
// RH parameters:
//     1: mesh handle
//     2: source position array (double matrix nq x dim)
//     3: detector position array (double matrix nm x dim)
// =========================================================================


#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"

using namespace std;
using namespace toast;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t i, j, d;
    QMMesh *mesh = (QMMesh*)Handle2Ptr(mxGetScalar(prhs[0]));
    size_t dim = (size_t)mesh->Dimension();

    // Copy source positions
    size_t nq = mxGetM(prhs[1]);
    size_t dimq = mxGetN(prhs[1]);
    d = std::min(dim,dimq);
    if (dim != dimq) {
	cerr << "Warning: toastSetQM: param 2:" << endl;
	cerr << "source array size(" << nq << "," << dimq
	     << ") does not correspond\nto mesh dimension ("
	     << dim << ")" << endl;
    }
    RDenseMatrix Q(nq,dimq);
    CopyMatrix (Q, prhs[1]);
    Point *pQ = new Point[nq];
    for (i = 0; i < nq; i++) {
	pQ[i].New(dim);
	for (j = 0; j < d; j++) pQ[i][j] = Q(i,j);
    }

    // Copy detector positions
    size_t nm = mxGetM(prhs[2]);
    size_t dimm = mxGetN(prhs[2]);
    d = std::min(dim,dimm);
    if (dim != dimm) {
	cerr << "Warning: toastSetQM: param 3:" << endl;
	cerr << "detector array size(" << nm << "," << dimm
	     << ") does not correspond\nto mesh dimension ("
	     << dim << ")" << endl;
    }
    RDenseMatrix M(nm,dimm);
    CopyMatrix (M, prhs[2]);
    Point *pM = new Point[nm];
    for (i = 0; i < nm; i++) {
	pM[i].New(dim);
	for (j = 0; j < d; j++) pM[i][j] = M(i,j);
    }

    // Assign Q and M arrays to mesh
    mesh->SetupQM (pQ, nq, pM, nm);

    mexPrintf ("QM: %d sources, %d detectors, %d measurements\n",
	       mesh->nQ, mesh->nM, mesh->nQM);
}
