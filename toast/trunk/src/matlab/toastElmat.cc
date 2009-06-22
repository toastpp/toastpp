// =========================================================================
// toastElmat
// Return an element matrix for a single element of a mesh
// The returned matrix contains the integral of products of shape functions
// or shape function derivatives over the element volume.
//
// RH parameters:
//     1: mesh handle
//     2: element index (>= 1)
//     3: integral type string (see below)
//     4+: any additional required parameters
// LH parameters:
//     1: element matrix (dimension and size depend on element type and
//        intgral type)
// =========================================================================

#include "mex.h"
#include "toastmex.h"
#include "util.h"

using namespace std;

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char cbuf[256];
    mxArray *elmat;
    double *pr;

    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int idx = (int)mxGetScalar(prhs[1]) - 1; // convert to zero-based
    Element *pel = mesh->elist[idx];
    int i, j, k, l, nnd = pel->nNode();
    int dim = mesh->Dimension();

    mxGetString (prhs[2], cbuf, 256);

    if (!strcmp(cbuf, "F")) {
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    pr[i] = pel->IntF(i);
    } else if (!strcmp(cbuf, "FF")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntFF(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntFF(i,j);
	}
    } else if (!strcmp(cbuf, "FFF")) {
	mwSize dims[3] = {nnd,nnd,nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    for (j = 0; j < nnd; j++) {
		for (k = 0; k < nnd; k++) {
		    pr[(i*nnd+j)*nnd+k] = pel->IntFFF(i,j,k);
		}
	    }
	}
    } else if (!strcmp(cbuf, "DD")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntDD(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntDD(i,j);
	}
    } else if (!strcmp(cbuf, "FD")) {
	mwSize dims[3] = {nnd, nnd, dim};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++) {
		RVector fd = pel->IntFD(i,j);
		if (fd.Dim() == 0) {
		    plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		    return;
		}
		for (k = 0; k < dim; k++)
		    pr[i + nnd*(j + k*nnd)] = fd[k];
	    }
    } else if (!strcmp(cbuf, "FDD")) {
	mwSize dims[3] = {nnd, nnd, nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++)
		for (k = 0; k < nnd; k++)
		    pr[i+nnd*(j+nnd*k)] = pel->IntFDD(i,j,k);
    } else if (!strcmp(cbuf, "dd")) {
	mwSize dims[4] = {nnd, dim, nnd, dim};
	elmat = mxCreateNumericArray (4, dims, mxDOUBLE_CLASS, mxREAL);
	RSymMatrix intdd = pel->Intdd();
	pr = mxGetPr(elmat);
	for (l = 0; l < dim; l++)
	    for (k = 0; k < nnd; k++)
		for (j = 0; j < dim; j++)
		    for (i = 0; i < nnd; i++)
			*pr++ = intdd(i*dim+j,k*dim+l);

    } else if (!strcmp(cbuf, "BndF")) {

	// for now, only integrals over all boundary sides are supported
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	RVector bndintf = pel->BndIntF();
	for (i = 0; i < pel->nNode(); i++)
	    pr[i] = bndintf[i];

    } else if (!strcmp(cbuf, "BndFF")) {
	int ii, jj, sd;
	if (nrhs > 3) sd = (int)mxGetScalar(prhs[3]) - 1; // side index
	else sd = -1;
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);

	if (sd >= 0) { // integral over a single side
	    for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		i = pel->SideNode(sd,ii);
		pr[i*nnd+i] = pel->BndIntFFSide(i,i,sd);
		for (jj = 0; jj < ii; jj++) {
		    j = pel->SideNode(sd,jj);
		    pr[i*nnd+j] = pr[j*nnd+i] = pel->BndIntFFSide(i,j,sd);
		}
	    }
	} else { // integral over all boundary sides
	    for (sd = 0; sd < pel->nSide(); sd++) {
		if (!pel->IsBoundarySide (sd)) continue;
		for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		    i = pel->SideNode(sd,ii);
		    pr[i*nnd+i] = pel->BndIntFFSide(i,i,sd);
		    for (jj = 0; jj < ii; jj++) {
			j = pel->SideNode(sd,jj);
			pr[i*nnd+j] = pr[j*nnd+i] = pel->BndIntFFSide(i,j,sd);
		    }
		}
	    }
	}
    }
    plhs[0] = elmat;
}
