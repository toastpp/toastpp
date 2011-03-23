// =========================================================================
// toastMapBasis2
// Maps a field between two different basis representations
//
// Interface:
// -----------------
// RH-1: basis handle
// RH-2: source and target basis representations (string):
//       'M->G' (mesh to grid)
//       'M->B' (mesh to basis)
//       'M->S' (mesh to solution)
//       'G->M' (grid to mesh)
//       'G->B' (grid to basis)
//       'G->S' (grid to solution)
//       'B->M' (basis to mesh)
//       'B->G' (basis to grid)
//       'B->S' (basis to solution)
//       'S->M' (solution to mesh)
//       'S->G' (solution to grid)
//       'S->B' (solution to basis)
// RH-3: field in source basis representation (real or complex vector)
// -----------------
// LH-1: field in target basis representation
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "util.h"

using namespace std;
using namespace toast;

// ============================================================================

void QuitOnError (char *msg)
{
    char cbuf[512] = "toastMapBasis: ";
    strcat (cbuf, msg);
    strcat (cbuf, "\n");
    mexErrMsgTxt (cbuf);
}

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3)
	QuitOnError ("Not enough parameters provided.");

    // basis handle
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[0]));

    // source and target basis representations
    if (!mxIsChar (prhs[1]))
	QuitOnError ("Argument 2: expected string.");

    char basisstr[256], srcid, tgtid;
    mxGetString (prhs[1], basisstr, 256);
    if (strlen(basisstr) != 4)
	QuitOnError ("Argument 2: format not recognised.");
    if (!strncmp (basisstr+1, "->", 2)) {
	srcid = toupper(basisstr[0]);
	tgtid = toupper(basisstr[3]);
    } else if (!strncmp (basisstr+1, "<-", 2)) {
	srcid = toupper(basisstr[3]);
	tgtid = toupper(basisstr[0]);
    } else {
	QuitOnError ("Argument 2: format not recognised.");
    }
    
    // source field
    mwIndex m = mxGetM (prhs[2]);
    mwIndex n = mxGetN (prhs[2]);
    int nsrc = (int)(m*n);

    if (mxIsComplex (prhs[2])) {

	double *pr = mxGetPr(prhs[2]);
	double *pi = mxGetPi(prhs[2]);
	CVector src(nsrc), tgt;
	toast::complex *vptr = src.data_buffer();
	for (int i = 0; i < nsrc; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	switch (srcid) {
	case 'M':
	    if (raster->mesh().nlen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'G':
		    raster->Map_MeshToGrid (src, tgt);
		    break;
		case 'B':
		    raster->Map_MeshToBasis (src, tgt);
		    break;
		case 'S':
		    raster->Map_MeshToSol (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'G':
	    if (raster->GLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'M':
		    raster->Map_GridToMesh (src, tgt);
		    break;
		case 'B':
		    raster->Map_GridToBasis (src, tgt);
		    break;
		case 'S':
		    raster->Map_GridToSol (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'B':
	    if (raster->BLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
	    } else {
		switch (tgtid) {
		case 'G':
		    raster->Map_BasisToGrid (src, tgt);
		    break;
		case 'S':
		    raster->Map_BasisToSol (src, tgt);
		    break;
		case 'M':
		    raster->Map_BasisToMesh (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'S':
	    if (raster->SLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'B':
		    raster->Map_SolToBasis (src, tgt);
		    break;
		case 'G':
		    raster->Map_SolToGrid (src, tgt);
		    break;
		case 'M':
		    raster->Map_SolToMesh (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	default:
	    QuitOnError ("Argument 2: source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);

    } else {

	RVector src (nsrc, mxGetPr (prhs[2]));
	RVector tgt;

	switch (srcid) {
	case 'M':
	    if (raster->mesh().nlen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'G':
		    raster->Map_MeshToGrid (src, tgt);
		    break;
		case 'B':
		    raster->Map_MeshToBasis (src, tgt);
		    break;
		case 'S':
		    raster->Map_MeshToSol (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'G':
	    if (raster->GLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'M':
		    raster->Map_GridToMesh (src, tgt);
		    break;
		case 'B':
		    raster->Map_GridToBasis (src, tgt);
		    break;
		case 'S':
		    raster->Map_GridToSol (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'B':
	    if (raster->BLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
	    } else {
		switch (tgtid) {
		case 'G':
		    raster->Map_BasisToGrid (src, tgt);
		    break;
		case 'S':
		    raster->Map_BasisToSol (src, tgt);
		    break;
		case 'M':
		    raster->Map_BasisToMesh (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	case 'S':
	    if (raster->SLen() != nsrc) {
		QuitOnError ("Argument 3: source vector unexpected length.");
		return;
	    } else {
		switch (tgtid) {
		case 'B':
		    raster->Map_SolToBasis (src, tgt);
		    break;
		case 'G':
		    raster->Map_SolToGrid (src, tgt);
		    break;
		case 'M':
		    raster->Map_SolToMesh (src, tgt);
		    break;
		default:
		    QuitOnError ("Argument 2: target id not recognised.");
		    return;
		}
	    }
	    break;
	default:
	    QuitOnError ("Argument 2: source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);
    }
}
