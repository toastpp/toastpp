// =========================================================================
// toastJacobian
// Generate a Jacobian (unscaled) for log amplitude and phase data and
// absorption,scattering parameters.
//
// RH parameters:
//     1: mesh handle (double)
//     2: basis mapper handle (0 for none)
//     3: source vectors (columns in complex sparse matrix)
//     4: measurement vectors (columns in complex sparse matrix)
//     5: mua (real vector, mesh basis)
//     6: mus (real vector, mesh basis)
//     7: refractive index (real vector, mesh basis)
//     8: modulation frequency [MHz]
//     9: linear solver (string)
//    10: iterative solver tolerance (double) (optional)
// LH parameters:
//     1: Jacobian matrix (dense double matrix, grid basis)
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace std;

// =========================================================================
// Implementation
// =========================================================================

// Calculate Jacobian from given direct and adjoint fields and boundary
// projection data

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, mxArray **res)
{
    int nQM, slen, ndat, nprm;
    nQM  = mesh->nQM;
    if (raster == RASTER_ELBASIS)
        slen = mesh->elen();
    else if (raster == RASTER_NDBASIS)
        slen = mesh->nlen();
    else
        slen = raster->SLen();
    ndat = nQM * 2;
    nprm = slen * 2;

    RDenseMatrix J(ndat,nprm);

    GenerateJacobian (raster, mesh, dphi, aphi, proj, dscale, J);

    CopyMatrix (res, J);
}

// Calculate Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, mxArray **res)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, nQM, nlen, elen, blen, slen;

    bool elbasis = (raster == RASTER_ELBASIS);
    nlen = mesh->nlen();
    elen = mesh->elen();
    n    = (elbasis ? elen : nlen);
    slen = (raster && raster != RASTER_ELBASIS ? raster->SLen() : n);
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;

    CVector *dphi, *aphi;
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

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (nlen);
    aphi = new CVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (nlen);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega, elbasis);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate projections if required
    CVector *proj = 0;
    DataScale dscale = FWS.GetDataScaling();
    if (dscale == DATA_LOG) {
    	proj = new CVector(nQM);
	*proj = FWS.ProjectAll (mvec, dphi, DATA_LIN);
    	//ProjectAll (*mesh, FWS, mvec, dphi, *proj);
    }

    // Calculate Jacobian
    CalcJacobian (mesh, raster, dphi, aphi, proj, dscale, res);

    delete []dphi;
    delete []aphi;
    if (proj) delete proj;
}                                                                              

void Assert (bool cond, int arg, const char *msg)
{
    if (!cond) {
	char cbuf[256];
	sprintf (cbuf, "toastJacobian: argument %d: %s", arg, msg);
	mexErrMsgTxt (cbuf);
    }
}

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// mesh
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int n   = mesh->nlen();
    int nq  = mesh->nQ;
    int nm  = mesh->nM;
    int nqm = mesh->nQM;
	int nprm = n;

    // raster
    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[1]));


	if (nrhs >= 11 && mxIsChar(prhs[10])) {
	    char cbuf[32];
	    mxGetString (prhs[10], cbuf, 32);
	    if (strcasecmp (cbuf, "EL") == 0 && !raster)
		raster = RASTER_ELBASIS;
		nprm = mesh->elen();
	}

    if (nrhs == 5) {

	// this is the version that provides fields and projections directly
	int i, j;
	double *pr, *pi;
	
	// copy fields
	const mxArray *mx_dphi = prhs[2];
	Assert (mxGetM(mx_dphi) == n, 3, "Unexpected number of rows");
	Assert (mxGetN(mx_dphi) == nq, 3, "Unexpected number of columns");
	Assert (mxIsComplex (mx_dphi), 3, "Must be complex");
	pr  = mxGetPr (mx_dphi);
	pi  = mxGetPi (mx_dphi);
	CVector *dphi = new CVector[nq];
	for (i = 0; i < nq; i++) {
	    dphi[i].New (n);
	    std::complex<double> *v = dphi[i].data_buffer();
	    for (j = 0; j < n; j++)
	        *v++ = std::complex<double>(*pr++, *pi++);
	}

	// copy adjoint fields
	const mxArray *mx_aphi = prhs[3];
	Assert (mxGetM(mx_aphi) == n, 4, "Unexpected number of rows");
	Assert (mxGetN(mx_aphi) == nm, 4, "Unexpected number of columns");
	Assert (mxIsComplex (mx_aphi), 4, "Must be complex");
	pr = mxGetPr (mx_aphi);
	pi = mxGetPi (mx_aphi);
	CVector *aphi = new CVector[nm];
	for (i = 0; i < nm; i++) {
	    aphi[i].New (n);
	    std::complex<double> *v = aphi[i].data_buffer();
	    for (j = 0; j < n; j++)
	        *v++ = std::complex<double>(*pr++, *pi++);
	}

	// copy projections
	const mxArray *mx_proj = prhs[4];
	Assert (mxGetM(mx_proj)*mxGetN(mx_proj) == nqm, 5, "Unexpected size");
	Assert (mxIsComplex(mx_proj), 5, "Must be complex");
	CVector proj(nqm);
	pr = mxGetPr (mx_proj);
	pi = mxGetPi (mx_proj);
	std::complex<double> *v = proj.data_buffer();
	for (i = 0; i < nqm; i++)
	    *v++ = std::complex<double>(*pr++, *pi++);

	CalcJacobian (mesh, raster, dphi, aphi, &proj, DATA_LOG,
		      &plhs[0]);

    } else {

	// this is the version that calculates fields on the fly

	// source vectors
	CCompRowMatrix qvec;
	CopyTMatrix (qvec, prhs[2]);

	// measurement vectors
	CCompRowMatrix mvec;
	CopyTMatrix (mvec, prhs[3]);

	// nodal optical parameters
	RVector mua (nprm, mxGetPr (prhs[4]));
	RVector mus (nprm, mxGetPr (prhs[5]));
	RVector ref (nprm, mxGetPr (prhs[6]));

	// modulation frequency
	double freq = mxGetScalar (prhs[7]);

	// linear solver parameters
	char solver[128];
	double tol = 1e-10;
	mxGetString (prhs[8], solver, 128);
	if (nrhs >= 10) tol = mxGetScalar (prhs[9]);

#ifdef WIN64
	// for now, direct solver doesn't work in WIN64
	if (!stricmp(solver,"direct")) {
		mexWarnMsgTxt("toastJacobian: direct solver not supported. Switching to BiCGSTAB(tol=1e-14)");
		strcpy (solver,"bicgstab");
		tol = 1e-14;
	}
#endif

	CalcJacobian (mesh, raster, qvec, mvec, mua, mus, ref, freq,
		      solver, tol, &plhs[0]);
    }
}
