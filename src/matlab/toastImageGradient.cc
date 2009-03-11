// =========================================================================
// toastImageGradient
// Calculate the gradient of an image represented in the fine grid basis
// of the basis mapper. Masked pixels are eliminated from the calculation.
//
// RH parameters:
//     1: basis mapper handle
//     2: image vector in fine grid ('g') format (real or complex)
// LH parameters:
//     1: image gradient (dense g x n matrix, with n=2,3
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;
using namespace toast;

// ============================================================================

template<class T>
void ImageGradient (const IVector &gdim, const RVector &gsize,
    const TVector<T> &im, TVector<T> *grad, const int *mask)
{
    // this should be done by Fourier transform
    int x, y, z, idx, dim = gdim.Dim();
    int nx = gdim[0], ny = gdim[1], nz = (dim >= 3 ? gdim[2]:1);
    int n = nx*ny*nz;

    // x gradient
    double dx = gsize[0]/nx, ix = 1.0/dx, i2x = 0.5*ix;
    TVector<T> &gradx = grad[0];
    gradx.New (n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt   = (!mask || (mask[idx] >= 0));
		bool bleft  = (x > 0 && (!mask || (mask[idx-1] >= 0)));
		bool bright = (x < nx-1 && (!mask || (mask[idx+1] >= 0)));
		if (bleft && bright) {
		    gradx[idx] = (im[idx+1]-im[idx-1]) * i2x;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradx[idx] = 0.0;
		} else if (bleft) {
		    gradx[idx] = (im[idx]-im[idx-1]) * ix;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bright) {
		    gradx[idx] = (im[idx+1]-im[idx]) * ix;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    gradx[idx] = 0.0;
		}
		idx++;
	    }
	}
    }

    // y gradient
    double dy = gsize[1]/ny, iy = 1.0/dy, i2y = 0.5*iy;
    TVector<T> &grady = grad[1];
    grady.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt  = (!mask || (mask[idx] >= 0));
		bool bup   = (y > 0 && (!mask || (mask[idx-nx] >= 0)));
		bool bdown = (y < ny-1 && (!mask || (mask[idx+nx] >= 0)));
		if (bup && bdown) {
		    grady[idx] = (im[idx+nx]-im[idx-nx]) * i2y;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    grady[idx] = 0.0;
		} else if (bup) {
		    grady[idx] = (im[idx]-im[idx-nx]) * iy;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bdown) {
		    grady[idx] = (im[idx+nx]-im[idx]) * iy;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    grady[idx] = 0.0;
		}
		idx++;
	    }
	}
    }
    if (dim < 3) return;

    // z gradient
    double dz = gsize[2]/nz, iz = 1.0/dz, i2z = 0.5*iz;
    int stridez = nx*ny;
    TVector<T> &gradz = grad[2];
    gradz.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt = (!mask || (mask[idx] >= 0));
		bool bfront = (z > 0 && !(mask || (mask[idx-stridez] >= 0)));
		bool bback  = (z < nz-1 && (!mask || (mask[idx+stridez] >= 0)));
	        if (bfront && bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx-stridez]) * i2z;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradz[idx] = 0.0;
		} else if (bfront) {
		    gradz[idx] = (im[idx]-im[idx-stridez]) * iz;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx]) * iz;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    gradz[idx] = 0.0;
		}
		idx++;
	    }
	}
    }
}

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    int glen = raster->GLen();
    int d, dim = gdim.Dim();

    if (mxIsComplex (prhs[1])) {
	CVector img;
	CVector *grad = new CVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	ImageGradient (gdim, gsize, img, grad, elref);
	CDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    } else {
	RVector img;
	RVector *grad = new RVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	ImageGradient (gdim, gsize, img, grad, elref);
	RDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    }
}
