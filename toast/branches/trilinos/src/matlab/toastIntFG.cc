// =========================================================================
// toastIntFG: Returns the product of two functions over a mesh basis
// Arguments can be real or complex.
// =========================================================================

#include "mex.h"
#include "felib.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

template<class T>
TVector<T> IntFG (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), Wrong vector size);
    dASSERT(g.Dim() == mesh.nlen(), Wrong vector size);

    int el, nnode, *node, i, j, k, nj, nk, bs;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
        pel   = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
	        nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFFF(i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFFF(i,j,k);
		}
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}


// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int n = mesh->nlen();

    bool isCplx = mxIsComplex (prhs[1]);
    if (isCplx != mxIsComplex (prhs[2]))
	mexErrMsgTxt
	   ("toastIntGradFGradG: Arguments must be both real or both complex");

    if (isCplx) {
	CVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec = IntFG (*mesh, fvec, gvec);
	CopyVector (&plhs[0], rvec);
    } else {
	RVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec= IntFG (*mesh, fvec, gvec);
	CopyVector (&plhs[0],rvec);
    }
}  
