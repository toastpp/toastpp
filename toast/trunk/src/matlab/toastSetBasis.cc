// =========================================================================
// toastSetBasis
// Sets up a mapper between FEM mesh, sampling grid and reconstruction basis.
//
// Interface:
// -----------------
// RH-1: basis type (string):
//       'LINEAR': bilinear pixel or trilinear voxel basis
//       'CUBIC': bicubic pixel or tricubic voxel spline basis
// RH-2: mesh handle
// RH-3: basis dimensions (2D or 3D integer vector)
// RH-4: high-res grid dimensions (2D or 3D integer vector [optional]
// RH-5: bounding box for basis grid [optional, default is mesh BB]
// -----------------
// LH-1: basis handle
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

using namespace std;

// ============================================================================

void QuitOnError (char *msg)
{
    char cbuf[512] = "toastSetBasis: ";
    strcat (cbuf, msg);
    strcat (cbuf, "\n");
    mexErrMsgTxt (cbuf);
}

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char basistype[256] = "LINEAR";
    int i, j, arg = 0, basistp;
    RDenseMatrix *bb = 0;

    // basis type
    if (nrhs <= arg)
	QuitOnError ("Too few parameters provided.");

    if (mxIsChar (prhs[arg])) {
	mxGetString (prhs[arg], basistype, 256);
	arg++;
    }
    if (!strcasecmp (basistype, "LINEAR")) {
	basistp = 0;
    } else if (!strcasecmp (basistype, "CUBIC")) {
	basistp = 1;
    } else {
	basistp = -1;
    }
    if (basistp < 0)
	QuitOnError ("Argument 1: unexpected value.");

    // mesh
    if (nrhs <= arg)
	QuitOnError ("Too few parameters provided.");
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[arg++]));

    // solution basis dimensions
    if (nrhs <= arg)
	QuitOnError ("Too few parameters provided.");
    int dim = mxGetNumberOfElements (prhs[arg]);
    if (dim != mesh->Dimension())
	QuitOnError ("Argument 3: Invalid basis dimensions");
    IVector bdim(dim);
    for (i = 0; i < dim; i++)
	bdim[i] = (int)mxGetPr (prhs[arg])[i];
    arg++;
    IVector gdim(bdim);

    if (nrhs > arg) {
	// sampling grid dimensions
	if (mxGetNumberOfElements (prhs[arg]) == dim) {
	    // argument is intermediate grid dimension
	    for (i = 0; i < dim; i++)
		gdim[i] = (int)mxGetPr (prhs[arg])[i];
	    arg++;
	}
    }
    if (nrhs > arg) {
	if (mxGetM (prhs[arg]) == dim && mxGetN (prhs[arg]) == 2) {
	    // argument is grid bounding box
	    bb = new RDenseMatrix (2,dim);
	    CopyTMatrix (*bb, prhs[arg]);
	    arg++;
	}
    }

    Raster *raster;
    switch (basistp) {
    case 0:
	raster = new Raster_Pixel (bdim, gdim, mesh, bb);
	break;
    case 1:
	raster = new Raster_CubicPixel (bdim, gdim, mesh, bb);
	break;
    }
	
    plhs[0] = mxCreateDoubleScalar (Ptr2Handle (raster));

    char cbuf[256];
    mexPrintf ("basis: Type:       %s-%s\n",
	       mesh->Dimension()==2 ? "BI":"TRI", basistype);
    sprintf (cbuf, "basis: Grid size:  %d [", raster->BLen());
    for (i = 0; i < dim; i++)
	sprintf (cbuf+strlen(cbuf), "%d%c", bdim[i], i==dim-1 ? ']':'x');
    mexPrintf ("%s\n", cbuf);

    mexPrintf ("basis: Sol. size:  %d\n", raster->SLen());
    mexPrintf ("basis: Mesh size:  %d\n", mesh->nlen());
    if (bb) {
	mexPrintf ("basis: Grid bounding box:\n");
	for (i = 0; i < 2; i++) {
	    for (j = 0; j < dim; j++)
		mexPrintf ("  %12g", bb->Get(i,j));
	    mexPrintf ("\n");
	}
	delete bb;
    }
}
