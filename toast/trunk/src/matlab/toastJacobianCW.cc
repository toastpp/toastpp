// =========================================================================
// toastJacobianCW
// Generate a Jacobian (unscaled) for log intensity (CW) data and
// absorption coefficient.
//
// RH parameters:
//     1: mesh handle (double)
//     2: basis mapper handle (double)
//     3: source vectors (columns in complex sparse matrix)
//     4: measurement vectors (columns in complex sparse matrix)
//     5: mua (real vector, mesh basis)
//     6: mus (real vector, mesh basis)
//     7: refractive index (real vector, mesh basis)
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
using namespace toast;

// =========================================================================
// local prototypes

//void Project (const QMMesh &mesh, int q, const RVector &phi,
//    RVector &proj);
void PMDF_mua (const RVector &dphi, const RVector &aphi, RVector &pmdf);
void GenerateJacobian (const Raster &raster, const QMMesh &mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J);
void GenerateJacobian (const QMMesh &mesh, const RCompRowMatrix &mvec,
    const RVector *dphi, const RVector *aphi,
    bool logparam, RDenseMatrix &J);

// =========================================================================
// Implementation

void CalcJacobian (QMMesh *mesh, Raster *raster,
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
    RFwdSolver FWS (solver, tol);

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
    FWS.Allocate (*mesh);
    FWS.Reset (msol, 0);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate Jacobian
    int ndat = mesh->nQM;
    int nprm = slen;
    RDenseMatrix J(ndat,nprm);
    if (raster)
	GenerateJacobian (*raster, *mesh, mvec, dphi, aphi, logparam, J);
    else
	GenerateJacobian (*mesh, mvec, dphi, aphi, logparam, J);

    delete []dphi;
    delete []aphi;

    CopyMatrix (res, J);
}                                                                              

// ============================================================================
// generate a projection from a field - complex case

#ifdef UNDEF
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
#endif

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
void PMDF_mua (const CVector &pmdf, complex proj,
    RVector &pmdf_mod, RVector &pmdf_arg)
{
    double idenom = 1.0/(proj.re*proj.re + proj.im*proj.im);
    for (int i = 0; i < pmdf.Dim(); i++) {
        pmdf_mod[i] = (pmdf[i].re*proj.re + pmdf[i].im*proj.im) * idenom;
	pmdf_arg[i] = (pmdf[i].im*proj.re - pmdf[i].re*proj.im) * idenom;
    }
}

// ============================================================================

void GenerateJacobian (const Raster &raster, const QMMesh &mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
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

	if (logparam)
	    proj = ProjectSingle (&mesh, i, mvec, dphi[i], DATA_LIN);

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
// This version doesn't use base mapping and works directly on the mesh basis

void GenerateJacobian (const QMMesh &mesh, const RCompRowMatrix &mvec,
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

	if (logparam)
	    proj = ProjectSingle (&mesh, i, mvec, dphi[i], DATA_LIN);

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
	
    CalcJacobian (mesh, raster, qvec, mvec, mua, mus, ref,
		  solver, tol, &plhs[0]);
}
