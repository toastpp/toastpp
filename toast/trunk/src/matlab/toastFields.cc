// =========================================================================
// toastJacobian
// Calculate direct and adjoint fields for a given distribution of optical
// parameters.
//
// RH parameters:
//     1: mesh handle (double)
//     2: basis mapper handle or 0 to return fields in mesh basis
//     3: source vectors (columns in complex sparse matrix)
//     4: measurement vectors (columns in complex sparse matrix)
//     5: mua (real vector, mesh basis)
//     6: mus (real vector, mesh basis)
//     7: refractive index (real vector, mesh basis)
//     8: modulation frequency [MHz]
//     9: linear solver (string)
//    10: iterative solver tolerance (double) (optional)
// LH parameters:
//     1: matrix of direct fields in solution basis
//     2: matrix of ajdoint fields in solution basis (optional)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace toast;

// =========================================================================
// Implementation

void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, mxArray **dfield, mxArray **afield)
{
    const double c0 = 0.3;
    int i, j, idx, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    slen = (raster ? raster->SLen() : n);

    CVector *dphi, *aphi;
    CFwdSolver FWS (solver, tol);

//  if (FWS.LinSolver() == LSOLVER_DIRECT)
//	mexPrintf ("Using direct solver\n");
//  else
//	mexPrintf ("Using iterative solver, tol=%f\n", FWS.ItSolverTol());

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // Calculate direct and adjoint fields
    FWS.Allocate (*mesh);
    FWS.Reset (msol, omega);
    CVector sphi(slen);

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    FWS.CalcFields (qvec, dphi);
    mxArray *dmx = mxCreateDoubleMatrix (slen, nQ, mxCOMPLEX);
    double *pr = mxGetPr (dmx);
    double *pi = mxGetPi (dmx);

    for (i = idx = 0; i < nQ; i++) {
	if (raster) raster->Map_MeshToSol (dphi[i], sphi);
	else        sphi = dphi[i];
	for (j = 0; j < slen; j++) {
	    pr[idx] = sphi[j].re;
	    pi[idx] = sphi[j].im;
	    idx++;
	}
    }
    delete []dphi;
    *dfield = dmx;

    // build adjoint field vectors
    if (afield) {
	aphi = new CVector[nM];
	for (i = 0; i < nM; i++) aphi[i].New (n);
	FWS.CalcFields (mvec, aphi);
	mxArray *amx = mxCreateDoubleMatrix (slen, nM, mxCOMPLEX);
	pr = mxGetPr (amx);
	pi = mxGetPi (amx);

	for (i = idx = 0; i < nM; i++) {
	    if (raster) raster->Map_MeshToSol (aphi[i], sphi);
	    else        sphi = aphi[i];
	    for (j = 0; j < slen; j++) {
		pr[idx] = sphi[j].re;
		pi[idx] = sphi[j].im;
		idx++;
	    }
	}
	delete []aphi;
	*afield = amx;
    }
}                                                    


// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // mesh
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int n = mesh->nlen();

    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[1]));

    // source vectors
    CCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    CCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // modulation frequency
    double freq = mxGetScalar (prhs[7]);

    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[8], solver, 128);
    if (nrhs >= 10) tol = mxGetScalar (prhs[9]);

    CalcFields (mesh, raster, qvec, mvec, mua, mus, ref, freq,
		solver, tol, &plhs[0], nlhs > 1 ? &plhs[1] : NULL);
}
