// =========================================================================
// toastSetBasis
// Sets up a mapper between FEM mesh, sampling grid and reconstruction basis.
//
// RH parameters:
//     1: mesh handle (scalar)
//     2: reconstruction basis dimensions (integer vector)
//     3: sampling grid dimensions (integer vector)
// LH parameters:
//     1: mapper handle (pointer)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, arg;
    RDenseMatrix *bb = 0;
    bool bIntermediate = false;

    // mesh
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    // solution basis dimensions
    int dim = mxGetNumberOfElements (prhs[1]);
    if (dim != mesh->Dimension())
	mexErrMsgTxt ("Invalid basis dimensions");
    IVector bdim(dim);
    for (i = 0; i < dim; i++)
	bdim[i] = (int)mxGetPr (prhs[1])[i];

    // sampling grid dimensions
    IVector gdim(dim);
    if (nrhs > 2) {
	arg = 2;
	if (mxGetNumberOfElements (prhs[arg]) == dim) {
	    // argument is intermediate grid dimension
	    for (i = 0; i < dim; i++)
		gdim[i] = (int)mxGetPr (prhs[arg])[i];
	    bIntermediate = true;
	    arg++;
	}
	if (nrhs > arg) {
	    if (mxGetM (prhs[arg]) == dim && mxGetN (prhs[arg]) == 2) {
		// argument is grid bounding box
		bb = new RDenseMatrix (2,dim);
		CopyTMatrix (*bb, prhs[arg]);
		arg++;
	    }
	}
    }

    Raster *raster;
    if (bIntermediate) raster = new Raster (gdim, bdim, mesh, bb);
    else               raster = new Raster (bdim, mesh, bb);
	
    plhs[0] = mxCreateScalarDouble (Ptr2Handle (raster));

    char cbuf[256];
    sprintf (cbuf, "Grid size:  %d [", raster->BLen());
    for (i = 0; i < dim; i++)
	sprintf (cbuf+strlen(cbuf), "%d%c", bdim[i], i==dim-1 ? ']':' ');
    mexPrintf ("%s\n", cbuf);

    if (nrhs > 2) {
	sprintf (cbuf, "Grid2 size: %d [", raster->GLen());
	for (i = 0; i < dim; i++)
	    sprintf (cbuf+strlen(cbuf), "%d%c", gdim[i], i==dim-1 ? ']':' ');
	mexPrintf ("%s\n", cbuf);
    }
    
    mexPrintf ("Sol. size:  %d\n", raster->SLen());
    mexPrintf ("Mesh size:  %d\n", mesh->nlen());
    if (bb) {
	mexPrintf ("Grid bounding box:\n");
	for (i = 0; i < 2; i++) {
	    for (j = 0; j < dim; j++)
		mexPrintf ("  %12g", bb->Get(i,j));
	    mexPrintf ("\n");
	}
	delete bb;
    }
}
