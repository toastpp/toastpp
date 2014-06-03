// =========================================================================
// toastSampleField
// Given a mesh, a nodal array of basis coefficients, and an array of
// sampling points, this returns the interpolated values of the field at
// the sampling points.
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: field: nodal basis coefficients (real array)
//     3: array of sampling points
// LH parameters:
//     1: array of interpolated field values
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, nsample, el;
    double *pr;

    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // field coefficients
    RVector phi;
    CopyVector (phi, prhs[1]);
    if (phi.Dim() != nlen)
	mexPrintf ("toastSampleField: invalid array dimension\n");

    // list of sampling points
    RDenseMatrix pt;
    CopyMatrix (pt, prhs[2]);
    nsample = pt.nRows();
    if (pt.nCols() != dim)
	mexPrintf ("toastSampleField: invalid array dimension\n");

    // create the output array
    mxArray *sample = mxCreateDoubleMatrix (nsample, 1, mxREAL);
    pr = mxGetPr (sample);

    for (i = 0; i < nsample; i++) {
	Point p(dim);
	for (j = 0; j < dim; j++) p[j] = pt(i,j);
	el = mesh->ElFind (p);
	if (el < 0) { // didn't find an element
	    double dst, dstmin = 1e10;
	    for (j = 0; j < elen; j++) {
		dst = p.Dist (mesh->ElCentre(j));
		if (dst < dstmin) dstmin = dst, el = j;
	    }
	}
	Element *pel = mesh->elist[el];
	RVector fun = pel->GlobalShapeF (mesh->nlist, p);

	double val = 0.0;
	for (j = 0; j < fun.Dim(); j++)
	    val += fun[j] * phi[pel->Node[j]];
	*pr++ = val;
    }

    plhs[0] = sample;
}
