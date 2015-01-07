// =========================================================================
// toastSysmatComponent
// Returns a single term of a system matrix
//
// RH parameters:
//     1: mesh handle (double)
//     2: prm NIM (double array)
//     3: integral type (string). Currently supported: FF,DD,PFF,PDD,BNDPFF
//     4: (optional) 'EL' flag to indicate element basis
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

void CalcSysmatComponent (QMMesh *mesh, RVector &prm, char *integral_type,
    bool elbasis, mxArray **res)
{
    int n = mesh->nlen();
    int *rowptr, *colidx, nzero;

    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F(n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (integral_type, "FF")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_FF);
    } else if (!strcasecmp (integral_type, "DD")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_DD);
    } else if (!strcasecmp (integral_type, "PFF")) {
	AddToSysMatrix (*mesh, F, &prm,
            elbasis ? ASSEMBLE_PFF_EL : ASSEMBLE_PFF);
    } else if (!strcasecmp (integral_type, "PDD")) {
	AddToSysMatrix (*mesh, F, &prm,
            elbasis ? ASSEMBLE_PDD_EL : ASSEMBLE_PDD);
    } else if (!strcasecmp (integral_type, "BNDPFF")) {
	AddToSysMatrix (*mesh, F, &prm,
            elbasis ? ASSEMBLE_BNDPFF_EL : ASSEMBLE_BNDPFF);
    }

    // Return system matrix to MATLAB
    CopyMatrix (res, F);
}

void CalcBndSysmat (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    CFwdSolver FWS (mesh, LSOLVER_DIRECT, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    FWS.Allocate ();
    RVector prm(n);
    for (int i = 0; i < n; i++) prm[i] = c0/ref[i];
    AddToSysMatrix (*mesh, *FWS.F, &prm, ASSEMBLE_BNDPFF);
    CopyMatrix (res, *FWS.F);
}

void CalcBndFactors (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    RVector prm(n);
    for (int i = 0; i < n; i++)
	prm[i] = A_Keijzer (ref[i]);
    CopyVector (res, prm);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    RVector prm;
    bool elbasis = false;
    char integral_type[32] = "";

    if (nrhs >= 2 && mxIsChar(prhs[1])) {
	mxGetString (prhs[1], integral_type, 32);
    } else {
	mexErrMsgTxt ("Parameter 2: string expected");
    }


    if (nrhs >= 4 && mxIsChar(prhs[3])) {
	char cbuf[32];
	mxGetString (prhs[3], cbuf, 32);
	elbasis = (strcasecmp (cbuf, "EL") == 0);
    }

    
    if (nrhs >= 3) {
        CopyVector (prm, prhs[2]);
        Assert ( (elbasis ? mesh->elen() : mesh->nlen())
                 == prm.Dim(), (char*)"Argument 3: wrong size");
        
    }

    CalcSysmatComponent (mesh, prm, integral_type, elbasis, &plhs[0]);
}
