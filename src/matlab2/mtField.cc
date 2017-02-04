// ========================================================================
// Implementation of class MatlabToast
// field-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

using namespace std;

// =========================================================================
// Prototypes
// =========================================================================

void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    mxArray **dfield);

// =========================================================================

void MatlabToast::Fields (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    int n = mesh->nlen();

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);

    // source vectors
    CCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[3]));
    RVector mus (n, mxGetPr (prhs[4]));
    RVector ref (n, mxGetPr (prhs[5]));

    // modulation frequency
    double freq = mxGetScalar (prhs[6]);

    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[7], solver, 128);
    if (nrhs >= 9) tol = mxGetScalar (prhs[8]);

    CalcFields (mesh, raster, qvec, mua, mus, ref, freq,
		solver, tol, &plhs[0]);
}

// =========================================================================
// Implementation

void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    mxArray **dfield)
{
    const double c0 = 0.3;
    int i, j, idx, n, dim, nQ, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = qvec.nRows();
    slen = (raster ? raster->SLen() : n);

    CVector *dphi;
    CFwdSolver FWS (mesh, solver, tol);

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
    FWS.SetPrecon (PRECON_ICH);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // Calculate direct and adjoint fields
    FWS.Allocate ();
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
	    pr[idx] = sphi[j].real();
	    pi[idx] = sphi[j].imag();
	    idx++;
	}
    }
    delete []dphi;
    *dfield = dmx;
}                                                    
