// ========================================================================
// Implementation of class MatlabToast
// Jacobian-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"
#include "fwdsolver.h"

using namespace std;

unsigned int verb;

// =========================================================================
// Prototypes
// =========================================================================

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, mxArray **res);

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, mxArray **res);

void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    char *solver, double tol, mxArray **res);

// =========================================================================
// Matlab interface
// =========================================================================

// =========================================================================
// Frequency-domain Jacobian

void MatlabToast::Jacobian (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    verb = verbosity;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n   = mesh->nlen();
    int nq  = mesh->nQ;
    int nm  = mesh->nM;
    int nqm = mesh->nQM;

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);

    if (nrhs == 5) {

	// this is the version that provides fields and projections directly
	int i, j;
	double *pr, *pi;
	
	// copy fields
	const mxArray *mx_dphi = prhs[2];
	ASSERTARG(mxGetM(mx_dphi) == n, 3, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_dphi) == nq, 3, "Unexpected number of columns");
        ASSERTARG(mxIsComplex (mx_dphi), 3, "Must be complex");
	pr  = mxGetPr (mx_dphi);
	pi  = mxGetPi (mx_dphi);
	CVector *dphi = new CVector[nq];
	for (i = 0; i < nq; i++) {
	    dphi[i].New (n);
	    std::complex<double> *v = dphi[i].data_buffer();
	    for (j = 0; j < n; j++)
	        *v++ = std::complex<double> (*pr++, *pi++);
	}
	// copy adjoint fields
	const mxArray *mx_aphi = prhs[3];
	ASSERTARG(mxGetM(mx_aphi) == n, 4, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_aphi) == nm, 4, "Unexpected number of columns");
	ASSERTARG(mxIsComplex (mx_aphi), 4, "Must be complex");
	pr = mxGetPr (mx_aphi);
	pi = mxGetPi (mx_aphi);
	CVector *aphi = new CVector[nm];
	for (i = 0; i < nm; i++) {
	    aphi[i].New (n);
	    std::complex<double> *v = aphi[i].data_buffer();
	    for (j = 0; j < n; j++)
	        *v++ = std::complex<double> (*pr++, *pi++);
	}
	// copy projections
	const mxArray *mx_proj = prhs[4];
	ASSERTARG(mxGetM(mx_proj)*mxGetN(mx_proj) == nqm, 5,"Unexpected size");
	ASSERTARG(mxIsComplex(mx_proj), 5, "Must be complex");
	CVector proj(nqm);
	pr = mxGetPr (mx_proj);
	pi = mxGetPi (mx_proj);
	std::complex<double> *v = proj.data_buffer();
	for (i = 0; i < nqm; i++)
	    *v++ = std::complex<double> (*pr++, *pi++);

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

    // bounds check on optical parameters
    ASSERTARG(mxGetM(prhs[4]) == n, 5, "Unexpected size");
    ASSERTARG(mxGetM(prhs[5]) == n, 6, "Unexpected size");
    ASSERTARG(mxGetM(prhs[6]) == n, 7, "Unexpected size");
        
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
	
	CalcJacobian (mesh, raster, qvec, mvec, mua, mus, ref, freq,
		      solver, tol, &plhs[0]);
    }    
}

// ==========================================================================
// CW Jacobian

void MatlabToast::JacobianCW (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    verb = verbosity;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    int n = mesh->nlen();

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);

    // source vectors
    RCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    RCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);
    
    // bounds check on optical parameters
    ASSERTARG(mxGetM(prhs[4]) == n, 5, "Unexpected size");
    ASSERTARG(mxGetM(prhs[5]) == n, 6, "Unexpected size");
    ASSERTARG(mxGetM(prhs[6]) == n, 7, "Unexpected size");

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));
    
    
    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[7], solver, 128);
    if (nrhs >= 9) tol = mxGetScalar (prhs[8]);
	
    CalcJacobianCW (mesh, raster, qvec, mvec, mua, mus, ref,
		  solver, tol, &plhs[0]);
}


// ==========================================================================
// Implementation
// ==========================================================================

// ==========================================================================
// Calculate Jacobian from given direct and adjoint fields and boundary
// projection data

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, mxArray **res)
{
    int nQM, slen, ndat, nprm;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : mesh->nlen());
    ndat = nQM * 2;
    nprm = slen * 2;

    RDenseMatrix J(ndat,nprm);

    ofstream ofs("dbg_m.dat");
    ofs << dphi[0] << endl;
    ofs.close();

    GenerateJacobian (raster, mesh, dphi, aphi, proj, dscale, J);

    CopyMatrix (res, J);
}

// ==========================================================================
// Calculate Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, mxArray **res)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : n);

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
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new CVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);
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

    ofstream ofs("dbg_m.dat");
    ofs << *proj << endl;
    ofs.close();

    // Calculate Jacobian
    CalcJacobian (mesh, raster, dphi, aphi, proj, dscale, res);

    delete []dphi;
    delete []aphi;
    if (proj) delete proj;
}                                                                              

// ==========================================================================
// Calculate CW Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    char *solver, double tol, mxArray **res)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    slen = (raster ? raster->SLen() : n);

    RVector *dphi, *aphi;
    RFwdSolver FWS (mesh, solver, tol);

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

    // build the field vectors
    dphi = new RVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new RVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, 0);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate Jacobian
    int ndat = mesh->nQM;
    int nprm = slen;
    RDenseMatrix J(ndat,nprm);
    GenerateJacobian_cw (raster, mesh, mvec, dphi, aphi, DATA_LOG, &J);

    delete []dphi;
    delete []aphi;

    CopyMatrix (res, J);
}                                                                              
