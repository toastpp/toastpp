// =========================================================================
// toastVolmat
// Generate a system matrix from a mesh by integrating over elements
//
// RH parameters:
//     1: mesh handle (double)
//     2: integration type (string)
//     3: nodal parameter array, if required
//     6: (optional) 'EL' flag to indicate element basis
// LH parameters:
//     1: system matrix (sparse double matrix)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

// =========================================================================
// Implementation

void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "toastSysmat: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}

void CalcSysmat (Mesh *mesh, const char *intstr, RVector &prm, bool elbasis,
    mxArray **res)
{
    // sysmatrix structure
    int n = mesh->nlen();
    int *rowptr, *colidx, nzero;
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (intstr, "FF")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_FF);
    } else if (!strcasecmp (intstr, "DD")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_DD);
    } else if (!strcasecmp (intstr, "PFF")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_PFF_EL : ASSEMBLE_PFF);
    } else if (!strcasecmp (intstr, "PDD")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_PDD_EL : ASSEMBLE_PDD);
    } else if (!strcasecmp (intstr, "BndPFF")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_BNDPFF_EL : ASSEMBLE_BNDPFF);
    } else if (!strcasecmp (intstr, "BndFF")) {
	RVector dummy(n);
	dummy = 1.0;
	AddToSysMatrix (*mesh, F, &dummy, ASSEMBLE_BNDPFF);
    }
    CopyMatrix (res, F);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    char intstr[256];
    RVector prm;
    bool elbasis = false;

    mxGetString (prhs[1], intstr, 256);
    
    if (nrhs > 2) {
	CopyVector (prm, prhs[2]);
	if (nrhs > 3 && mxIsChar (prhs[3])) {
	    char cbuf[32];
	    mxGetString (prhs[3], cbuf, 32);
	    elbasis = (strcasecmp (cbuf, "EL") == 0);
	}
    }

    CalcSysmat (mesh, intstr, prm, elbasis, &plhs[0]);
}
