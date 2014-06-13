// ========================================================================
// Implementation of class MatlabToast
// Jacobian-related methods
// ========================================================================

#include "matlabtoast.h"
#include "toastmex.h"

using namespace std;

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

void CalcJacobianCW (const Raster &raster, const QMMesh &mesh,
    const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J);

void CalcJacobianCW (const QMMesh &mesh,
    const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J);

void Project (const QMMesh &mesh, int q, const RVector &phi,
    RVector &proj);

RVector IntFG (const Mesh &mesh, const RVector &f, const RVector &g);

void PMDF_mua (const RVector &dphi, const RVector &aphi, RVector &pmdf);
RVector PMDF_mua (const RVector &dphi, const RVector &aphi, double proj);
CVector PMDF_mua (const CVector &dphi, const CVector &aphi);
void PMDF_mua (const CVector &pmdf, std::complex<double> proj,
    RVector &pmdf_mod, RVector &pmdf_arg);

// =========================================================================
// Matlab interface
// =========================================================================

// =========================================================================
// Frequency-domain Jacobian

void MatlabToast::Jacobian (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
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
    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    int n = mesh->nlen();

    // raster
    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    // source vectors
    RCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    RCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[7], solver, 128);
    if (nrhs >= 10) tol = mxGetScalar (prhs[8]);
	
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
    int i, n, dim, nQ, nM, nQM, blen, slen;

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
    int i, n, dim, nQ, nM, nQM, blen, slen;

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

    bool logparam = true; // make user-definable

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
    if (raster)
	CalcJacobianCW (*raster, *mesh, dphi, aphi, logparam, J);
    else
	CalcJacobianCW (*mesh, dphi, aphi, logparam, J);

    delete []dphi;
    delete []aphi;

    CopyMatrix (res, J);
}                                                                              

// ============================================================================

void CalcJacobianCW (const Raster &raster, const QMMesh &mesh,
    const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J)
{
    int i, j, jj, k, idx, dim, nQ, nM, nQM, slen, glen;
    const RGenericSparseMatrix &B = raster.Mesh2GridMatrix();
    const IVector &gdim = raster.GDim();
    const IVector &bdim = raster.BDim();
    const RVector &gsize = raster.GSize();
    dim  = raster.Dim();
    glen = raster.GLen();
    slen = raster.SLen();
    nQ   = mesh.nQ;
    nM   = mesh.nM;
    nQM  = mesh.nQM;

    RVector pmdf_mua(glen);
    RVector pmdf_basis(slen);
    RVector proj;
    RVector cdfield(glen);
    RVector *cafield = new RVector[nM];
    RVector *cdfield_grad = new RVector[dim];
    RVector *cafield_grad = new RVector[dim];

    // resample all measurement fields to fine grid
    //cout << "Allocating " << glen*nM*8*2 << " bytes for cafield" << endl;
    for (i = 0; i < nM; i++) {
        cafield[i].New (glen);
	raster.Map_MeshToGrid (aphi[i], cafield[i]);
    }

    for (i = idx = 0; i < nQ; i++) {

        // resample field and field gradient for source i to fine grid
	raster.Map_MeshToGrid (dphi[i], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad,
		       raster.Elref());

	if (logparam) {
	    proj.New (mesh.nQMref[i]);
	    Project (mesh, i, dphi[i], proj);
	}

	for (j = jj = 0; j < nM; j++) {
	    if (!mesh.Connected (i,j)) continue;

	    // measurement field gradient
	    ImageGradient (gdim, gsize, cafield[j], cafield_grad,
			   raster.Elref());

	    PMDF_mua (cdfield, cafield[j], pmdf_mua);
	    if (logparam) pmdf_mua /= proj[jj];

	    // map into solution basis
	    raster.Map_GridToSol (pmdf_mua, pmdf_basis);
	    for (k = 0; k < slen; k++)
	        J(idx,k) = pmdf_basis[k];

	    idx++;
	    jj++;
	}
    }
    delete []cafield;
    delete []cdfield_grad;
    delete []cafield_grad;
}

// ============================================================================
// This version doesn't use base mapping and works directly on the mesh basis

void CalcJacobianCW (const QMMesh &mesh,
    const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J)
{
    cerr << "Jacobian: using mesh basis" << endl;
    cerr << "Dim: " << J.nRows() << " x " << J.nCols() << endl;

    int i, j, jj, k, n, idx, dim, nQ, nM, nQM, slen, glen;
    n    = mesh.nlen();
    nQ   = mesh.nQ;
    nM   = mesh.nM;
    nQM  = mesh.nQM;

    RVector pmdf_mua (n);
    RVector cdfield(n);
    RVector proj;

    for (i = idx = 0; i < nQ; i++) {

	proj.New (mesh.nQMref[i]);
	Project (mesh, i, dphi[i], proj);

	for (j = jj = 0; j < nM; j++) {

	    if (!mesh.Connected (i,j)) continue;
	    pmdf_mua = IntFG (mesh, dphi[i], aphi[j]);
	    if (logparam) pmdf_mua /= proj[jj];

	    // map into solution basis
	    for (k = 0; k < n; k++)
	        J(idx,k) = -pmdf_mua[k];

	    idx++;
	    jj++;
	}
    }
}

// ============================================================================
// generate a projection from a field - real case

void Project (const QMMesh &mesh, int q, const RVector &phi,
    RVector &proj)
{
    int i, m, el, nv, dim, nd, in;
    double dphi;
    double c2a;

    for (i = 0; i < mesh.nQMref[q]; i++) {
	m = mesh.QMref[q][i];
	el = mesh.Mel[m];
	nv = mesh.elist[el]->nNode();
	dim = mesh.elist[el]->Dimension();
	dphi = 0.0;
	
	RVector fun = mesh.elist[el]->GlobalShapeF (mesh.nlist, mesh.M[m]);
	for (in = 0; in < nv; in++) {
	    nd = mesh.elist[el]->Node[in];
	    c2a = mesh.plist[nd].C2A();   // reflection parameter
	    dphi += phi[nd]*fun[in]*c2a;
	}
	proj[i] = dphi;
    }
}

// ============================================================================

RVector IntFG (const Mesh &mesh, const RVector &f, const RVector &g)
{
    dASSERT(f.Dim() == mesh.nlen(), Wrong vector size);
    dASSERT(g.Dim() == mesh.nlen(), Wrong vector size);

    int el, nnode, *node, i, j, k, nj, nk, bs;
    double sum;
    Element *pel;
    RVector tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
        pel   = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
	        nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFFF(i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFFF(i,j,k);
		}
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}

// ============================================================================
// PMDF calculations

// absorption PMDF (real)
void PMDF_mua (const RVector &dphi, const RVector &aphi, RVector &pmdf)
{
    int i, len = pmdf.Dim();
    for (i = 0; i < len; i++)
        pmdf[i] = -dphi[i]*aphi[i];
}

// Intensity PMDF for absorption, given the direct and adjoint fields
// and measurement proj for the given source-detector pair
RVector PMDF_mua (const RVector &dphi, const RVector &aphi, double proj)
{
    return (dphi*aphi) / proj;
}

// Complex PMDF for absorption, given the direct and adjoint fields
CVector PMDF_mua (const CVector &dphi, const CVector &aphi)
{
    return -(dphi*aphi);
}

// Extract modulation amplitude PMDF and phase PMDF from complex PMDF
void PMDF_mua (const CVector &pmdf, std::complex<double> proj,
    RVector &pmdf_mod, RVector &pmdf_arg)
{
    double idenom = 1.0/norm(proj);
    for (int i = 0; i < pmdf.Dim(); i++) {
        pmdf_mod[i] = (pmdf[i].real()*proj.real() + pmdf[i].imag()*proj.imag())
	    * idenom;
	pmdf_arg[i] = (pmdf[i].imag()*proj.real() - pmdf[i].real()*proj.imag())
	    * idenom;
    }
}

