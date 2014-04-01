// ==========================================================================
// SolverLM: Levenberg-Marquardt
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "regul.h"
#include "solverlm.h"
#include "dnsmatrix_mpi.h"
#include "jacobian_mpi.h"
#include "timing.h"
//#include <unistd.h>
#include <time.h>
//#include <sys/times.h>
//#include <sys/param.h>
#include "timing.h"

#define DEBUG_FRECHET
#define BYHAND

using namespace std;
using namespace toast;

// ==========================================================================
// external references

extern int bWriteJ;
extern int g_imgfmt;
extern int bOutputUpdate;
extern int bOutputGradient;
extern double clock0;
extern char g_meshname[256];
extern int g_nimsize;
extern double g_refind;

extern void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh);

// ==========================================================================
// local prototypes

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
    int nprm);

RVector RescaleHessian (const RDenseMatrixMPI &Jmod,
    const RDenseMatrixMPI &Jarg, const RVector &x, RCompRowMatrix *RHess);

bool LineSearchWithPrior (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, const RVector &p0, const RVector &p,
    const Regularisation *reg, const RMatrix *cov = 0);

bool CheckRange (const Solution &sol);

RDenseMatrix mul (const RMatrix &A, const RMatrix &B)
{
    int i, j, k;
    double sum;
    int ar = A.nRows();
    int ac = A.nCols();
    int br = B.nRows();
    int bc = B.nCols();
    xASSERT (ac == br, Incompatible dimensions);
    RDenseMatrix tmp (ar, bc);
    for (i = 0; i < ar; i++) {
	for (j = 0; j < bc; j++) {
	    sum = 0.0;
	    for (k = 0; k < ac; k++)
		sum += A.Get(i,k) * B.Get(k,j);
	    tmp(i,j) = sum;
	}
    }
    return tmp;
}

// ==========================================================================
// parameter block for definition of distributed Jacobian

#ifdef TOAST_MPI
struct JAC_DATA_MPI {
    int np, rank;        // number of processes, current process
    int rtot;            // total number of rows in complete Jacobian
    int r0, r1;          // first and last+1 rows of Jacobian to calculate
};
#endif // TOAST_MPI

// ==========================================================================
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const RDenseMatrixMPI *Jmod;   // distributed Jacobian (log amp part)
    const RDenseMatrixMPI *Jarg;   // distributed Jacobian (phase part)
    const RVector *M;              // normalisation diagonal matrix
    const double *lambda;          // diagonal scaling factor
    const Regularisation *reg;     // regularisation
    const RCompRowMatrix *RHess;   // Hessian of prior (1 matrix per parameter)
};

// ==========================================================================
// This data structure defines the Jacobian implicitly

struct JAC_DATA {
    const QMMesh *mesh;            // mesh
    const Raster *raster;          // basis mapper
    const CFwdSolver *FWS;         // forward solver
          CVector **dphi;          // fields for each source in mesh basis
          CVector **aphi;          // adjoint fields
    const CCompRowMatrix *mvec;    // measurement vectors
    const RVector *dat;            // measurements
    const RVector *sd;             // sd
    const RVector *prmscl;         // parameter scaling vector
    const Regularisation *reg;     // regularisation
    const RCompRowMatrix *RHess;   // Hessian of prior (1 matrix per parameter)
    const RVector *M;              // normalisation diagonal matrix
    const double *lambda;          // diagonal scaling factor
          double omega;            // modulation frequency
};

// ==========================================================================

#ifdef UNDEF
static RSymMatrix Hess_full (const HESS_DATA &data, const RVector &x)
{
    // Computes the explicit Hessian from its components, defined as
    // (J^T J + M P'' M + lambda I)
    // where J is the normalised Jacobian,
    // P is the second derivative of the prior (un-normalised),
    // M is the diagonal normalisation matrix,
    // lambda is a scalar, and I is identity

    // unpack the data
    const RMatrix *J          =  data.J;
    const Regularisation *reg =  data.reg;
    const RVector &M          = *data.M;
    const double lambda       = *data.lambda;
    int m = J->nRows(), n = J->nCols();

    RSymMatrix hess (ATA (*J));

    // add prior to Hessian
    int i, j, k, nz, *colidx = new int[n];
    double *val = new double[n];
    for (i = 0; i < n; i++) {
	nz = reg->GetHessianRow (x, i, colidx, val);
	for (k = 0; k < nz && (j=colidx[k]) <= i; k++)
	    hess(i,j) += val[k] * M[i]*M[j];
    }
    delete []colidx;
    delete []val;

    hess.AddDiag (lambda);
    return hess;
}
#endif

// ==========================================================================

#ifdef UNDEF
static RVector Hess_diag (const HESS_DATA &data, const RVector &x)
{
    // Computes the diagonal of the Hessian given implicitly as
    // (J^T J + P + lambda D)
    // where J and P are full, D is the diagonal of J^T J + P, and
    // lambda is a scalar

    // unpack the data
    const RMatrix *J            =  data.J;
    const RCompRowMatrix *RHess =  data.RHess;
    const RVector &M            = *data.M;
    const double lambda         = *data.lambda;
    int m = J->nRows(), n = J->nCols();

    RVector diag = ATA_diag (*J);

    // add prior to Hessian
    diag += RHess->Diag();

    diag += lambda;
    return diag;
}
#endif

// ==========================================================================

static RVector JTJx_clbk (const RVector &x, void *context)
{
    // Computes (J^T J + M P'' M + lambda I) x
    // where J is the column-normalised Jacobian,
    // P'' is the second derivative of the prior term,
    // M is the diagonal scaling matrix, and
    // lambda is a scalar

    // unpack the data
    HESS_DATA *data             = (HESS_DATA*)context;
    const RDenseMatrixMPI *Jmod =  data->Jmod;
    const RDenseMatrixMPI *Jarg =  data->Jarg;
    const RCompRowMatrix *RHess =  data->RHess;
    const RVector &M            = *data->M;
    const double lambda         = *data->lambda;
    int m2 = Jmod->nRows();
    int m = m2*2;
    int n = Jmod->nCols();

    RVector Px(n);

    // add prior to Hessian
    if (RHess) {
	int i, j, k, nz, *colidx = new int[n];
	double *val = new double[n];
	for (i = 0; i < n; i++) {
	    nz = RHess->SparseRow (i, colidx, val);
	    for (k = 0; k < nz; k++) {
		j = colidx[k];
		Px[i] += val[k] * x[j] * M[i]*M[j];
	    }
	}
	delete []colidx;
	delete []val;
    }

    RVector y(m);
    RVector ymod(y,0,m2), yarg(y,m2,m2);
    RVector z(n), zmod(n), zarg(n);
    Jmod->Ax(x,ymod);
    Jarg->Ax(x,yarg);
    Jmod->ATx(ymod,zmod);
    Jarg->ATx(yarg,zarg);
    z = zmod+zarg;

    return z + Px + lambda*x;
}

// ==========================================================================
// Generate the Frechet derivatives: direct method
// This version only requires the availability of the direct fields
// (dphi_h) passed in mesh basis
// Note: parameter and data scaling not yet implemented

RVector FrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &mvec, const CVector *dphi_h,
    const RVector &sd, const RVector &x)
{
    int i, j, q;
    int dim = mesh.Dimension();
    int nq = mesh.nQ;
    int nqm = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int ofs_lnmod = 0;
    int ofs_phase = nqm;
    CVector qdelta_g(glen);
    CVector qdelta_h(nlen);
    CVector phi_h_delta(nlen);
    CVector dphi_p(slen);
    CVector dphi_g(glen);
    CVector *dphi_g_grad = new CVector[dim];
    CVector *dphi_g_gradgrad = new CVector[dim];
    CVector tmp(slen);
    RVector y_h_delta (nqm*2);

    for (i = 0; i < dim; i++) {
	dphi_g_grad[i].New (glen);
	dphi_g_gradgrad[i].New (glen);
    }

    RVector alpha_h(x, 0, slen);
    RVector beta_h(x, slen, slen);
    RVector alpha_g(glen), beta_g(glen);
    CVector calpha_g(glen), cbeta_g(glen);
    raster.Map_SolToGrid (alpha_h, alpha_g);
    raster.Map_SolToGrid (beta_h, beta_g);
    SetReal (calpha_g, alpha_g);
    SetReal (cbeta_g, beta_g);

    for (q = 0; q < nq; q++) {
	// build up perturbation rhs (qdelta_g)
	raster.Map_MeshToGrid (dphi_h[q], dphi_g);
	qdelta_g = -calpha_g * dphi_g;
	ImageGradient (raster.GDim(), raster.GSize(), dphi_g, dphi_g_grad, 
		       raster.Elref());
	for (i = 0; i < dim; i++) {
	    dphi_g_grad[i] *= cbeta_g;
	    // note: improve divergence calculation
	    ImageGradient (raster.GDim(), raster.GSize(), dphi_g_grad[i],
			   dphi_g_gradgrad, raster.Elref());
	    qdelta_g += dphi_g_gradgrad[i];
	}

	// to make the results compatible with the definition of the
	// explicit Jacobian, we need to scale with voxel size here
	double scl = 1.0;
	RVector gsize = raster.GSize();
	IVector bdim = raster.BDim();
	for (i = 0; i < gsize.Dim(); i++) scl *= gsize[i]/bdim[i];
	qdelta_g *= toast::complex(1.0/scl,0);

	// remap to mesh basis
	raster.Map_GridToMesh (qdelta_g, qdelta_h);

	// scale with element support size
	CVector escl(nlen);
	for (i = 0; i < mesh.elen(); i++) {
	    Element *pel = mesh.elist[i];
	    for (j = 0; j < pel->nNode(); j++) {
		escl[pel->Node[j]].re += pel->IntF(j);
	    }
	}
	qdelta_h *= escl;

	// solve for phi_delta
	FWS.CalcField (qdelta_h, phi_h_delta);

	// boundary projections
	int nm = mesh.nQMref[q];
	CVector projdelta = ProjectSingle (&mesh, q, mvec, phi_h_delta);

	if (FWS.GetDataScaling() == DATA_LOG) {
	    CVector proj (nm);
	    Project_cplx (mesh, q, dphi_h[q], proj);
	    for (i = 0; i < nm; i++) {
		double scl = norm2 (proj[i]);
		toast::complex c = projdelta[i]*conj(proj[i]);
		projdelta[i] = c/scl;
	    }
	}

	// map data types into real vector
	for (i = 0; i < nm; i++) {
	    toast::complex proji = projdelta[i];
	    y_h_delta[ofs_lnmod++] = proji.re;
	    y_h_delta[ofs_phase++] = proji.im;
	}
    }
    // scale with sd
    y_h_delta /= sd;

    delete []dphi_g_grad;
    delete []dphi_g_gradgrad;

    return y_h_delta;
}

// ==========================================================================
// Generate the Frechet derivatives:
// This version assumes the availability of all direct and adjoint fields
// (dphi_h and aphi_h) passed in mesh basis

RVector FrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &mvec, const CVector *dphi_h,
    const CVector *aphi_h, const RVector &sd, const RVector &x)
{
    if (!aphi_h)
	return FrechetDerivative (mesh,raster,FWS,mvec, dphi_h, sd, x);
    // no adjoint fields - use direct method

    tic();
    int i, q, m, mm, ofs;
    int n    = x.Dim();
    int nq   = mesh.nQ;
    int nqm  = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    CVector res(nqm);
    RVector rres(nqm*2);

    // multiply the direct fields with the rhs before
    // multiplying in the adjoint fields

    CVector pd(glen), pa(glen);
    RVector xa(x, 0, n/2), xb(x, n/2, n/2);
    RVector xag(glen), xbg(glen);
    raster.Map_SolToGrid (xa, xag);
    raster.Map_SolToGrid (xb, xbg);
    CVector xacg(glen), xbcg(glen);
    SetReal (xacg, xag);
    SetReal (xbcg, xbg);

    CVector *pd_grad = new CVector[dim];
    CVector *pa_grad = new CVector[dim];
    CVector pas(slen), *pas_grad = new CVector[dim];
    CVector *pds_alpha = new CVector[nq];
    CVector **pds_alpha_grad = new CVector*[nq];
    for (i = 0; i < nq; i++) pds_alpha_grad[i] = new CVector[dim];
    CVector *proj = new CVector[nq];
    int *qofs = new int[nq];
    CVector pd_alpha, *pd_alpha_grad = new CVector[dim];

    for (i = 0; i < dim; i++)
	pas_grad[i].New(slen);

    // loop over sources
    for (q = 0; q < nq; q++) {
	qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	raster.Map_MeshToGrid (dphi_h[q], pd);
	ImageGradient (gdim, gsize, pd, pd_grad, raster.Elref());

	if (FWS.GetDataScaling() == DATA_LOG) {
	    proj[q].New (mesh.nQMref[q]);
	    Project_cplx (mesh, q, dphi_h[q], proj[q]);
	}
	pd_alpha = pd * xacg;
	pds_alpha[q].New (slen);
	raster.Map_GridToSol (pd_alpha, pds_alpha[q]);
	for (i = 0; i < dim; i++) {
	    pd_alpha_grad[i] = pd_grad[i] * xbcg;
	    pds_alpha_grad[q][i].New(slen);
	    raster.Map_GridToSol (pd_alpha_grad[i], pds_alpha_grad[q][i]);
	}
    }
    // loop over detectors
    for (m = 0; m < mesh.nM; m++) {
	raster.Map_MeshToGrid (aphi_h[m], pa);
	raster.Map_GridToSol (pa, pas);
	ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());
	for (i = 0; i < dim; i++)
	    raster.Map_GridToSol (pa_grad[i], pas_grad[i]);
	
	for (q = 0; q < nq; q++) {
	    if (!mesh.Connected (q,m)) continue;
	    mm = mesh.QMref[q][m];
	    ofs = qofs[q] + mm;
	    res[ofs] = -pds_alpha[q] & pas;
	    for (i = 0; i < dim; i++)
		res[ofs] -= pds_alpha_grad[q][i] & pas_grad[i];
	    if (FWS.GetDataScaling() == DATA_LOG) {
		double scl = norm2 (proj[q][mm]);
		toast::complex c = res[ofs]*conj(proj[q][mm]);
		res[ofs] = c/scl;
	    }
	}
    }
    // copy result to real vector and scale with sd
    for (i = 0; i < nqm; i++) {
	rres[i]     = re(res[i]) / sd[i];
	rres[i+nqm] = im(res[i]) / sd[i+nqm];
    }

    delete []pd_grad;
    delete []pa_grad;
    delete []pas_grad;
    delete []pds_alpha;
    for (i = 0; i < nq; i++) delete []pds_alpha_grad[i];
    delete []pds_alpha_grad;
    delete []proj;
    delete []qofs;
    delete []pd_alpha_grad;

    cerr << "Frechet timing: " << toc() << endl;
    return rres;
}

// ==========================================================================
// Adjoint Frechet derivative:
// Direct method (requires only direct fields)

RVector AdjointFrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &mvec, const CVector *dphi_h,
    const RVector &sd, const RVector &y)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    double term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    CVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    CVector dgrad (slen);
    ofs_mod = 0;         // data offset for Mod data
    ofs_arg = mesh.nQM;  // data offset for Arg data
    RVector rres (slen*2);
    RVector grad_cmua(rres, 0, slen);       // mua part of grad
    RVector grad_ckappa (rres, slen, slen); // kappa part of grad
    
    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster.Map_MeshToGrid (dphi_h[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh.nQMref[q];

	RVector y_mod (y, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector b_mod(n);
	b_mod = y_mod/s_mod;

	RVector y_arg (y, ofs_arg, n);
	RVector s_arg (sd, ofs_arg, n);
	RVector b_arg(n);
	b_arg = y_arg/s_arg;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	Project_cplx (mesh, q, dphi_h[q], cproj);
	wqa = toast::complex(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].re;
	    double ip = cproj[idx].im;
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = /* -2.0 * */ b_mod[idx] /* / (ype[idx]*s_mod[idx]) */;
	    wqa += qs * toast::complex (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = /* -2.0 * */ b_arg[idx] /* / (ype[idx]*s_arg[idx]) */;
	    wqa += qs * toast::complex (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster.Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster.Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster.Map_GridToSol (gk, dgrad);
	grad_ckappa -= Re(dgrad);

	ofs_mod += n; // step to next source
	ofs_arg += n;
	delete []cdfield_grad;
	delete []cafield_grad;
    }
    return rres;
}

// ==========================================================================
// Adjoint Frechet derivative:
// This method uses direct and adjoint fields. If adjoint fields are not
// available, falls back to direct method above

RVector AdjointFrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &mvec, const CVector *dphi_h,
    const CVector *aphi_h, const RVector &sd, const RVector &y)
{
    if (!aphi_h) 
	return AdjointFrechetDerivative (mesh, raster, FWS, mvec, dphi_h,
					 sd, y);
    // no adjoint fields - use direct method

    tic();
    int i, q, m, mm, ofs;
    int nq   = mesh.nQ;
    int nqm  = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector rres(slen*2);
    RVector ysd = y/sd; // data scaling

    CVector pd(glen), *pds = new CVector[nq], pa(glen), pas(slen);
    CVector das(slen), dbs(slen);
    CVector db(glen);

    CVector *pd_grad = new CVector[dim];
    CVector **pds_grad = new CVector*[nq];
    for (i = 0; i < nq; i++) pds_grad[i] = new CVector[dim];
    CVector *pa_grad = new CVector[dim];
    CVector *pas_grad = new CVector[dim];
    CVector *proj = new CVector[nq];
    int *qofs = new int[nq];

    for (i = 0; i < dim; i++)
	pas_grad[i].New(slen);

    // loop over sources
    for (q = ofs = 0; q < nq; q++) {
	qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	raster.Map_MeshToGrid (dphi_h[q], pd);
	ImageGradient (gdim, gsize, pd, pd_grad, raster.Elref());

	if (FWS.GetDataScaling() == DATA_LOG) {
	    proj[q].New (mesh.nQMref[q]);
	    Project_cplx (mesh, q, dphi_h[q], proj[q]);
	}
	pds[q].New(slen);
	raster.Map_GridToSol (pd, pds[q]);
	for (i = 0; i < dim; i++) {
	    pds_grad[q][i].New(slen);
	    raster.Map_GridToSol (pd_grad[i], pds_grad[q][i]);
	}
    }

    // loop over detectors
    for (m = 0; m < mesh.nM; m++) {
	raster.Map_MeshToGrid (aphi_h[m], pa);
	ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());
	raster.Map_GridToSol (pa, pas);
	for (i = 0; i < dim; i++)
	    raster.Map_GridToSol (pa_grad[i], pas_grad[i]);

	for (q = 0; q < nq; q++) {
	    if (!mesh.Connected (q,m)) continue;
	    mm = mesh.QMref[q][m];
	    ofs = qofs[q] + mm;
	    das = pds[q] * pas;
	    dbs.Clear();
	    for (i = 0; i < dim; i++)
		dbs += pds_grad[q][i] * pas_grad[i];
	    double yre = ysd[ofs], yim = ysd[ofs+nqm];
	    if (FWS.GetDataScaling() == DATA_LOG) {
		// rescale for log data (lnamp + phase)
		toast::complex logscl = proj[q][mm]/norm2 (proj[q][mm]);
		toast::complex logdat = toast::complex(yre,yim) * logscl;
		yre = logdat.re;
		yim = logdat.im;
	    }
	    for (i = 0; i < slen; i++) {
		// apply data to PMDF; write to real vector
		rres[i]      -= das[i].re*yre + das[i].im*yim;
		rres[i+slen] -= dbs[i].re*yre + dbs[i].im*yim;
	    }
	}
    }
    delete []pds;
    delete []pd_grad;
    for (i = 0; i < nq; i++) delete []pds_grad[i];
    delete []pds_grad;
    delete []pa_grad;
    delete []pas_grad;
    delete []proj;
    delete []qofs;

    cerr << "Adjoint Frechet timing: " << toc() << endl;
    return rres;
}

// ==========================================================================
// Generate diagonal scaling matrix M that scales the diagonal of the
// Hessian to 1. This version is for implicit Jacobian given in terms of
// direct and adjoint fields

RVector ImplicitHessianScaling (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &mvec, const CVector *dphi_h,
    const CVector *aphi_h, const RVector &sd, const RVector &prmscl)
{
    int i, q, m, mm, ofs;
    int nq   = mesh.nQ;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    int nqm  = mesh.nQM;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector M(slen*2);
    tic();

    CVector *pd = new CVector[nq];
    CVector pa(glen);
    CVector das(slen), dbs(slen);
    CVector db(glen);

    CVector **pd_grad = new CVector*[nq];
    for (i = 0; i < nq; i++) pd_grad[i] = new CVector[dim];
    CVector *pa_grad = new CVector[dim];
    CVector *proj = new CVector[nq];
    int *qofs = new int[nq];

    // loop over sources
    for (q = ofs = 0; q < nq; q++) {
	qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	pd[q].New(glen);
	for (i = 0; i < dim; i++) pd_grad[q][i].New(glen);
	raster.Map_MeshToGrid (dphi_h[q], pd[q]);
	ImageGradient (gdim, gsize, pd[q], pd_grad[q], raster.Elref());

	if (FWS.GetDataScaling() == DATA_LOG) {
	    proj[q].New (mesh.nQMref[q]);
	    Project_cplx (mesh, q, dphi_h[q], proj[q]);
	    
	}
    }

    // loop over detectors
    for (m = 0; m < mesh.nM; m++) {
	raster.Map_MeshToGrid (aphi_h[m], pa);
	ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());

	for (q = 0; q < nq; q++) {
	    if (!mesh.Connected (q,m)) continue;
	    mm = mesh.QMref[q][m];
	    ofs = qofs[q] + mm;
	    CVector da = pd[q] * pa;
	    raster.Map_GridToSol (da, das);
	    db.Clear();
	    for (i = 0; i < dim; i++)
		db += pd_grad[q][i] * pa_grad[i];
	    raster.Map_GridToSol (db, dbs);
	    if (FWS.GetDataScaling() == DATA_LOG) {
		double scl = norm2 (proj[q][mm]);
		das *= (conj(proj[q][mm])/scl);
		dbs *= (conj(proj[q][mm])/scl);
	    }
	    for (i = 0; i < slen; i++) { // data scaling
		das[i].re /= sd[ofs];
		dbs[i].re /= sd[ofs];
		das[i].im /= sd[ofs+nqm];
		dbs[i].im /= sd[ofs+nqm];
	    }
	    for (i = 0; i < slen; i++) {
		M[i]      += norm2 (das[i]);
		M[i+slen] += norm2 (dbs[i]);
	    }
	}
    }
    delete []pd;
    for (i = 0; i < nq; i++) delete []pd_grad[i];
    delete []pd_grad;
    delete []pa_grad;
    delete []proj;
    delete []qofs;

    M *= prmscl*prmscl; // parameter scaling
    return M;
}

// ==========================================================================
// Generate diagonal scaling matrix M that scales the diagonal of the
// Hessian to 1. This version solves one direct and adjoint field

RVector ExplicitHessianScaling1 (const QMMesh &mesh, const Raster &raster,
    const CFwdSolver &FWS, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec,  const RVector &sd, const RVector &prmscl)
{
    int i, q, m;
    int nq   = mesh.nQ;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    int nqm  = mesh.nQM;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector M(slen*2);
    tic();

    CVector pd, pa(glen);
    CVector das(slen), dbs(slen);
    CVector db(glen);

    CVector *pd_grad = new CVector[dim];
    CVector *pa_grad = new CVector[dim];
    CVector *proj = new CVector[nq];
    int *qofs = new int[nq];

    CVector qsum(nlen), msum(nlen), dphi(nlen), aphi(nlen);

    cout << "Using Explicit Hessian Scaling\n";
    // loop over sources
    for (q = 0; q < nq; q++) {
      qsum += qvec.Row(q);
    }
    // loop over detectors
    for (m = 0; m < mesh.nM; m++) {
      msum += mvec.Row(m);
    }
    FWS.CalcField (qsum, dphi);
    FWS.CalcField (msum, aphi);

    //	qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
    pd.New(glen);
    for (i = 0; i < dim; i++) pd_grad[i].New(glen);
    raster.Map_MeshToGrid (dphi, pd);
    ImageGradient (gdim, gsize, pd, pd_grad, raster.Elref());

    raster.Map_MeshToGrid (aphi, pa);
    ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());
    CVector da = pd * pa;
    raster.Map_GridToSol (da, das);
    db.Clear();
    for (i = 0; i < dim; i++)
	    db += pd_grad[i] * pa_grad[i];
    raster.Map_GridToSol (db, dbs);
    for (i = 0; i < slen; i++) {
		M[i]      += norm2 (das[i]);
		M[i+slen] += norm2 (dbs[i]);
    }
    delete []pd_grad;
    delete []pa_grad;
    delete []proj;
    delete []qofs;

    M *= prmscl*prmscl; // parameter scaling
    return M;
}

// ==========================================================================

#ifdef UNDEF // think about this!

RVector ApproxHessianScaling (const QMMesh &mesh, const Raster &raster,
    CFwdSolver &FWS, const CCompRowMatrix &qvec, const CCompRowMatrix &mvec)
{
    int i, j, idx;
    int dim  = raster.Dim();
    int glen = raster.GLen();
    int slen = raster.SLen();
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    CVector proj(mesh.nQM);
    CVector q(qvec.nCols());
    CVector m(mvec.nCols());
    CVector projsum(mvec.nCols());
    CVector dphisum (mesh.nlen());
    CVector aphisum (mesh.nlen());
    CVector *dphi = new CVector[mesh.nQ];
    for (i = 0; i < mesh.nQ; i++) {
	dphi[i].New (mesh.nlen());
	q += qvec.Row(i);
    }
    FWS.CalcField (mesh, q, dphisum);
    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
    ProjectAll (mesh, FWS, mvec, dphi, proj);
    for (i = idx = 0; i < mesh.nQ; i++) {
	for (j = 0; j < mesh.nM; j++) {
	    if (!mesh.Connected (i,j)) continue;
	    projsum[j] += proj[idx++];
	}
    }
    for (j = 0; j < mesh.nM; j++) {
	m += mvec.Row(j) / projsum[j];
    }
    FWS.CalcField (mesh, m, aphisum);
    
    CVector dgrid(glen);
    CVector agrid(glen);
    raster.Map_MeshToGrid (dphisum, dgrid);
    raster.Map_MeshToGrid (aphisum, agrid);

    RVector dgrid_re_grad[dim], dgrid_im_grad[dim];
    RVector agrid_re_grad[dim], agrid_im_grad[dim];
    ImageGradient (gdim, gsize, Re(dgrid), dgrid_re_grad, raster.Elref());
    ImageGradient (gdim, gsize, Im(dgrid), dgrid_im_grad, raster.Elref());
    ImageGradient (gdim, gsize, Re(agrid), agrid_re_grad, raster.Elref());
    ImageGradient (gdim, gsize, Im(agrid), agrid_re_grad, raster.Elref());

    RVector pmdf1 (glen*2);
    RVector pmdf1_mua (pmdf1, 0, glen);
    RVector pmdf1_kap (pmdf1, glen, glen);
    RVector pmdf2 (glen*2);
    RVector pmdf2_mua (pmdf2, 0, glen);
    RVector pmdf2_kap (pmdf2, glen, glen);

    PMDF_mua_Re (dphisum, aphisum, pmdf1_mua);
    PMDF_kappa_Re (dgrid_re_grad, dgrid_im_grad, agrid_re_grad, agrid_im_grad,
		   pmdf1_kap, dim);
    PMDF_mua_Im (dphisum, aphisum, pmdf2_mua);
    PMDF_kappa_Im (dgrid_re_grad, dgrid_im_grad, agrid_re_grad, agrid_im_grad,
		   pmdf2_kap, dim);

    if (FWS.datatype = MEAS_FMOD_FARG) {
	toast::complex scal = toast::complex(1,0); // don't know yet
	RVector pmdf_mod(glen*2), pmdf_arg(glen*2);
	PMDF_Mod (pmdf1, pmdf2, scal, pmdf_mod);
	PMDF_Arg (pmdf1, pmdf2, scal, pmdf_arg);
	pmdf1 = pmdf_mod, pmdf2 = pmdf_arg;
    }

    RVector M(slen*2);
    RVector pmdf_basis(slen);
    raster.Map_GridToSol (pmdf1_mua, pmdf_basis);
    for (i = 0; i < slen; i++) M[i] = pmdf_basis[i];
    raster.Map_GridToSol (pmdf2_mua, pmdf_basis);
    for (i = 0; i < slen; i++) M[i] += pmdf_basis[i];
    raster.Map_GridToSol (pmdf1_kap, pmdf_basis);
    for (i = 0; i < slen; i++) M[slen+i] = pmdf_basis[i];
    raster.Map_GridToSol (pmdf2_kap, pmdf_basis);
    for (i = 0; i < slen; i++) M[slen+i] += pmdf_basis[i];

    delete []dphi;
}

#endif

// ==========================================================================

RVector DistanceField (const Raster &raster)
{
    int i, j, k, n, idx;

    const Mesh *mesh = &raster.mesh();
    IVector gdim = raster.GDim();
    int dim  = gdim.Dim();
    int slen = raster.SLen();
    int glen = raster.GLen();
    int nlen = mesh->nlen();
    int nbnd = mesh->nbnd();

    // first mark all boundary nodes
    int *bndidx = new int[nbnd];
    for (i = j = 0; i < nlen; i++)
	if (mesh->nlist[i].isBnd()) bndidx[j++] = i;

    RVector distfield(glen);
    Point pt(dim);
    Point bbmin(dim), bbmax(dim);
    mesh->BoundingBox (bbmin, bbmax);
    int ni = gdim[0];
    int nj = gdim[1];
    int nk = (dim > 2 ? gdim[2]:1);

    double di = (bbmax[0]-bbmin[0])/(double)(ni-1);
    double dj = (bbmax[1]-bbmin[1])/(double)(nj-1);
    double dk = (dim > 2 ? (bbmax[2]-bbmin[2])/(double)(nk-1) : 1.0);

    double dstmin, dst;

    for (k = idx = 0; k < nk; k++) {
	if (dim > 2) pt[2] = k*dk+bbmin[2];
	for (j = 0; j < nj; j++) {
	    pt[1] = j*dj+bbmin[1];
	    for (i = 0; i < ni; i++) {
		pt[0] = i*di+bbmin[0];

		dstmin = 1e10;
		for (n = 0; n < nbnd; n++) {
		    dst = pt.Dist (mesh->nlist[bndidx[n]]);
		    if (dst < dstmin) dstmin = dst;
		}
		distfield[idx++] = dstmin;
	    }
	}
    }

    delete []bndidx;
    return distfield;
}

// ==========================================================================

RVector DistScaling (const Raster &raster, const RVector &logx)
{
    RVector x = exp(logx);

    double c = 0.3/1.4; // for now - GENERALISE!
    double scal_mua = 1e4;
    double scal_mus = 1e2;

    int slen = raster.SLen();
    RVector gdist = DistanceField (raster);
    gdist = -exp(-gdist*0.2)+1.1;

    RVector sdist (slen);
    raster.Map_GridToSol (gdist, sdist);

    RVector xmua (x, 0, slen);
    RVector xkap (x, slen, slen);
    RVector M (slen*2);
    RVector Mmua (M, 0, slen);
    RVector Mkap (M, slen, slen);

    RVector smua(slen), smus(slen);
    smua = xmua/c*scal_mua;
    smus = inv(xkap*3.0)*c*scal_mus;
    
    Mmua = sdist*smua;
    Mkap = sdist*smus;

    return inv(sqr(M));
}

// ==========================================================================

static RVector Frechet_clbk (const RVector &x, void *context)
{
    // Jacobian-free version of the Hx callback function

    // unpack the data
    JAC_DATA *data = (JAC_DATA*)context;
    const QMMesh *mesh = data->mesh;         // mesh
    const Raster *raster = data->raster;     // basis mapper
    const CFwdSolver *FWS = data->FWS;       // forward solver
          CVector **dphi_h = data->dphi;     // fields in mesh basis
	  CVector **aphi_h = data->aphi;     // adjoint fields
    const CCompRowMatrix *mvec = data->mvec; // measurement vectors
    const RVector *dat = data->dat;         // data
    const RVector *sd = data->sd;            // sd
    const RVector *prmscl = data->prmscl;    // parameter scaling vector
    const RCompRowMatrix *RHess = data->RHess;
    const RVector &M  = *data->M;
    const double lambda = *data->lambda;

    RVector xscl(x);
    xscl *= M;        // apply scaling Jx -> JMx
    xscl *= *prmscl;  // apply param scaling dy/dx -> dy/d(log x)

    // Frechet derivative z = F' x
    RVector y = FrechetDerivative (*mesh, *raster, *FWS, *mvec,
        *dphi_h, *aphi_h, *sd, xscl);

    // Adjoint Frechet derivative F'* z
    RVector z = AdjointFrechetDerivative (*mesh, *raster, *FWS, *mvec,
	*dphi_h, *aphi_h, *sd, y);

    z *= M;        // apply scaling J^T y -> M J^T y
    z *= *prmscl;  // apply param scaling dy/dx -> dy/d(log x)


#ifdef UNDEF
    // temporary!!!
    RVector z2 = AdjointFrechetDerivative (*mesh, *raster, *FWS, *mvec,
	*dphi_h, 0, *sd, *prmscl * M, y);
    ofstream ofs("dbg1.dat");
    ofs << z << endl;
    ofstream ofs2("dbg2.dat");
    ofs2 << z2 << endl;
#endif

    int n = raster->SLen()*2;
    RVector Px(n);

    // add prior to Hessian
    if (RHess) {
	int i, j, k, nz, *colidx = new int[n];
	double *val = new double[n];
	for (i = 0; i < n; i++) {
	    nz = RHess->SparseRow (i, colidx, val);
	    for (k = 0; k < nz; k++) {
		j = colidx[k];
		Px[i] += val[k] * x[j] * M[i]*M[j];
	    }
	}
	delete []colidx;
	delete []val;
    }

    return z + Px + lambda*x;
}

// ==========================================================================
SolverLM::SolverLM (ParamParser *_pp): Solver (_pp)
{
    nrmax = 50;
    itmax = 1000;
    gmres_tol = 1e-8;
    gn_tol = 1e-8;
    lambda0 = 10.0;
    lambda_scale = 4.0;
    do_linesearch = true;
    hscale = LM_HSCALE_IMPLICIT;
    alpha_min = 0.0;
    alpha0 = 1.0;
    Jfill = 1.0;
    reg = NULL;
    modelerr = false;
}

SolverLM::~SolverLM ()
{
    if (reg) delete reg;
}

// ==========================================================================

void SolverLM::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    bool Use_precon = true;

    // initialisations
    int i, inr;
    const QMMesh *mesh = FWS.MeshPtr();
    int nqm  = mesh->nQM;
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int nprm = bsol.nActive();
    int n    = slen*2;
    double of, fmin, errstart, err0, err00, err1;
    double lambda = lambda0, alpha = -1.0;

    RVector r(n), s(n), d(n);
    RVector rmod(n), rarg(n);
    RVector h(n), dold(n), M(n);

    CVector *dphi = 0, *aphi = 0, *aphi_hscale = 0;
    RVector proj(ndat);
    RVector prmscl;
    RVector sdmod(sd,0,nqm), sdarg(sd,nqm,nqm);
    RDenseMatrixMPI Jmod(nqm,n), Jarg(nqm,n);
    
    RMatrix *cov;

    Regularisation *reg;
    RCompRowMatrix RHess;

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    // initial solution (scaled)

    RVector x(x0), x1(x0);
    // current solution, trial solution (scaled);

    // reset forward solver
    LOGOUT("Resetting forward solver");
    FWS.Reset (msol, omega);

    // calculate direct fields and projections
    dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);
    FWS.CalcFields (qvec, dphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

    // calculate adjoint fields
    if ((precon != LM_PRECON_GMRES_DIRECT) || hscale == LM_HSCALE_IMPLICIT) {
	aphi = new CVector[mesh->nM];
	for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);
	FWS.CalcFields (mvec, aphi);
	if (hscale == LM_HSCALE_IMPLICIT)
	    aphi_hscale = aphi;
	if (precon == LM_PRECON_GMRES_DIRECT) aphi = 0;
    }

    int i_count = 0; // iteration counter

    // create the regularisation instance
    reg = Regularisation::Create (pp, &x0, &raster);

    // add model error means to data
    if (modelerr) {
	(RVector)data -= merr;
	cov = &mcov;
    } else {
	cov = 0;
    }

    HESS_DATA hdata = {&Jmod, &Jarg, &M, &lambda, reg, &RHess};
    // set the components for the implicit definition of the Hessian

    bool logparam = !strcmp (pscaler->ScalingMethod(), "LOG");
    JAC_DATA jdata = {mesh, &raster, &FWS, &dphi, &aphi, &mvec, &data, &sd,
		      &prmscl, reg, &RHess, &M, &lambda, omega};
    // set the components for the implicit definition of the Jacobian

    of  = ObjectiveFunction::get_value (data, proj, sd, cov);
    of += reg->GetValue (x1);
    // actually this will be 0 since x = x0 = x1 here
    err0 = errstart = of;

    err00 = 0; // this is previous error

    //    test_biscale();
    LOGOUT("Starting error: %f", errstart);

    LOGOUT("Iteration 0  CPU %f  OF %f", toc(clock0), err0);

    // LM iteration loop
    for (inr = 0; (!nrmax || inr < nrmax) &&
	   err0 > gn_tol*errstart && fabs(err00-err0) > gn_tol; inr++) {

	RVector kap = bsol.GetParam (OT_CKAPPA);
	prmscl = pscaler->JScale (bsol.GetActiveParams());
	RVector b = (data-proj)/sd;     // negative residual
	RVector bmod(b,0,nqm);
	RVector barg(b,nqm,nqm);

	// update Hessian of prior
	RHess = BuildRHessian (reg, x, nprm);

	// Generate Jacobian
	LOGOUT("Calculating Jacobian ...");
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
			  FWS.GetDataScaling(), Jmod, Jarg);

	// apply data space rescaling T
	Jmod.RowScale (inv(sdmod));
	Jarg.RowScale (inv(sdarg));

	// apply solution space rescaling S
	Jmod.ColScale (bsol.GetActiveParams());
	Jarg.ColScale (bsol.GetActiveParams());

	// apply Hessian normalisation M
	M = RescaleHessian (Jmod, Jarg, x, &RHess);
	Jmod.ColScale (M);
	Jarg.ColScale (M);
	    
	// calculate Gradient
	Jmod.ATx (bmod,rmod);
	Jarg.ATx (barg,rarg);
	r = rmod+rarg;

	LOGOUT("Gradient norm: %f", l2norm(r));
	
	// add prior to gradient
	RVector gpsi = reg->GetGradient (x);
	gpsi *= M; // apply Hessian normalisation
	r -= gpsi;
	LOGOUT("Gradient norm with penalty: %f", l2norm(r));

	bool KeepGoing = true;
	i_count = 0;
	err00 = err0;
	while (i_count < itmax && KeepGoing) {
	    LOGOUT("LM iteration %d", i_count);

	    // solve Hh = r

	    switch (precon) {
	    case LM_PRECON_NONE:    // assume H = I
	        h = r;
		break;
#ifdef UNDEF
	    case LM_PRECON_HDIAG:   // assume H = diag(H)
		LOGOUT ("Solving Hessian: LM-DIAG ...");
	        h = r / Hess_diag (hdata, x);
		break;
	    case LM_PRECON_CH: {    // solve with Cholesky factorisation
		LOGOUT ("Solving Hessian: LM-CH ...");
	        RSymMatrix hess = Hess_full (hdata, x);// need explicit Hessian
		CHdecomp (hess);
		h = CHsubst(hess, r);
	        } break;
	    case LM_PRECON_ICH:     // Incomplete Cholesky factorisation
	        cerr << "LM ICH preconditioner not yet implemented" << endl;
		exit(1);
#endif
	    case LM_PRECON_GMRES: {  // GMRES solver with implicit Hessian
		LOGOUT ("Solving Hessian: LM-GMRES ...");
	        double tol = gmres_tol;
		static RPrecon_Diag precon;
		precon.ResetFromDiagonal (M);
		GMRES (JTJx_clbk, &hdata, r, h, tol, &precon, 20);
		//GMRES (JTJx_clbk, &hdata, r, h, tol);
		//GMRES (Hess_full (hdata, x), r, h, tol); // debug only
	        } break;
	    case LM_PRECON_BICGSTAB: { // BiCGSTAB with implicit Hessian
		LOGOUT ("Solving Hessian: LM-BICGSTAB ...");
		double tol = gmres_tol;
		static RPrecon_Diag precon;
		precon.ResetFromDiagonal (M);
		BiCGSTAB (JTJx_clbk, &hdata, r, h, tol, &precon);
	        } break;
	    case LM_PRECON_GMRES_JACOBIANFREE:
	    case LM_PRECON_GMRES_DIRECT: {
		double tol = gmres_tol;
		static RPrecon_Diag prec;
		prec.ResetFromDiagonal (M);
		GMRES (Frechet_clbk, &jdata, r, h, tol, &prec, 20);
		//GMRES (Frechet_clbk, &jdata, r, h, tol);
		} break;

	    }
	    d = h;
	    double stepsize = alpha0;
	    LOGOUT("Step size: %f", stepsize);

#ifdef NEWTON_CG
	    double beta;
	    if(inr >  0) {
	        beta = (h&r)/hrold;
		d += dold*beta;
	    }
	    dold = d;
	    hrold = h&r;
#endif

	    bool status_valid = false;
	    if (do_linesearch) {
		static double alpha_ls = stepsize;
		//if (alpha < 0.0) alpha = stepsize; // initialise step length
		if (!LineSearchWithPrior (FWS, raster, pscaler, qvec, mvec,
                  data, sd, omega, d, err0, fmin, alpha_ls, msol, x0, x,
		  reg, cov)) {
		    lambda *= lambda_scale;
		    if (lambda_scale == 1.0) KeepGoing = false;
		    LOGOUT ("No decrease in line search");
		    continue;
		}
		if (alpha_ls >= alpha_min) {
		    alpha = alpha_ls;
		    err1 = fmin;
		    status_valid = true;
		} else {
		    alpha = alpha_min;
		}
	    } else {
		alpha = stepsize; // fixed step length
	    }

	    x1 = x + d*alpha;  // trial solution	

	    if (bOutputUpdate) {
		static int count = 0;
		static double muamin, muamax, kappamin, kappamax;
		char cbuf[256];
		RVector tmp = pscaler->Unscale(d*alpha);
		RVector tmp1 (tmp, 0, slen);
		RVector tmp2 (tmp, slen, slen);
		RVector vg(glen);
		raster.Map_SolToGrid (tmp1, vg);
		if (!count) ImageScale (vg, muamin, muamax);
		sprintf (cbuf, "images/update_mua_%03d.ppm", count);
		WritePPM (vg, raster.GDim(), &muamin, &muamax, cbuf);
		raster.Map_SolToGrid (tmp2, vg);
		if (!count) ImageScale (vg, kappamin, kappamax);
		sprintf (cbuf, "images/update_kappa_%03d.ppm", count);
		WritePPM (vg, raster.GDim(), &kappamin, &kappamax, cbuf);
		count++;
	    }

	    raster.Map_ActiveSolToMesh (pscaler->Unscale(x1), msol);
	    if (!status_valid) {
		FWS.Reset (msol, omega);
		FWS.CalcFields (qvec, dphi);
		if (aphi) FWS.CalcFields (mvec, aphi);
		else if (aphi_hscale)
		    FWS.CalcFields (mvec, aphi_hscale);
		proj = FWS.ProjectAll_real (mvec, dphi);
		err1 = ObjectiveFunction::get_value (data, proj, sd, cov);
	    }

	    if(err1 < err0) {  // accept update
	        //err00 = err0;
		err0 = err1;
		x = x1;
		bsol.SetActiveParams (pscaler->Unscale(x));
		raster.Map_SolToMesh (bsol, msol);
		FWS.Reset (msol, omega);
		FWS.CalcFields (qvec, dphi);
		if (aphi) FWS.CalcFields (mvec, aphi);
		else if (aphi_hscale)
		    FWS.CalcFields (mvec, aphi_hscale);
		proj = FWS.ProjectAll_real (mvec, dphi);

		switch (g_imgfmt) {
		case IMGFMT_NIM:
		    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
		    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
		    break;
		case IMGFMT_PGM:
		case IMGFMT_PPM: {
		    char cbuf[256];
		    bool col = (g_imgfmt == IMGFMT_PPM);
		    double smin, smax, c = 0.3/g_refind;
		    RVector mua(slen), mus(slen), vb(blen);
		    mua = bsol.GetParam (OT_CMUA) / c;
		    smin = vmin(mua), smax = vmax(mua);
		    raster.Map_SolToBasis (mua, vb);
		    sprintf (cbuf, "recon_mua_%03d.ppm", inr);
		    WritePixmap (vb, raster.BDim(), &smin, &smax, cbuf, col);
		    mus = inv(bsol.GetParam (OT_CKAPPA))*(c/3.0) - mua;
		    smin = vmin(mus), smax = vmax(mus);
		    raster.Map_SolToBasis (mus, vb);
		    sprintf (cbuf, "recon_mus_%03d.ppm", inr);
		    WritePixmap (vb, raster.BDim(), &smin, &smax, cbuf, col);
		    } break;
		case IMGFMT_RAW:
		    Solution rsol(OT_NPARAM, blen);
		    raster.Map_SolToBasis (bsol, rsol, true);
		    // need to unscale here
		    rsol.WriteImg_mua (inr+1, "recon_mua.raw");
		    rsol.WriteImg_mus (inr+1, "recon_mus.raw");
		    break;
		}
		lambda /= lambda_scale;
		KeepGoing = false;
	    } else {          // reject update
		if (lambda_scale == 1.0) KeepGoing = false;
		// if we don't allow lambda adjustment
		// then this is the end of it
	        lambda *= lambda_scale;
	    }
	    if (lambda) {
		LOGOUT("LM param lambda: %g", lambda);
	    }
	    i_count++;
	}
	LOGOUT("Iteration %d  CPU %f  OF %f", inr+1, toc(clock0), err0);
    } // end of NR loop;

    // final residuals
    double rd = ObjectiveFunction::get_value (data, proj, sd, cov);
    double rp = reg->GetValue (x);
    double tau = reg->GetTau();
    LOGOUT("Residuals  TAU %g  DATA %g  PRIOR/TAU %g  IT %d", tau, rd,
		rp/tau, inr);
}

void SolverLM::ReadParams (ParamParser &pp)
{
    char cbuf[256], c;
    bool def = false;

    // 1. === MAX NUMBER OF ITERATIONS ===

    if (!pp.GetInt ("ITMAX", nrmax) || nrmax < 0) do {
	cout << "\nMax number of GN iterations (0 for unlimited):\n>> ";
	cin >> nrmax;
    } while (nrmax < 0);

    // 2. === NONLINEAR SOLVER CONVERGENCE CRITERION ===

    if (!pp.GetReal ("NONLIN_TOL", gn_tol) ||
	gn_tol <= 0.0) do {
	    cout << "\nSelect convergence criterion for "
		 << "GN nonlinear solver (>0):\n>> ";
	    cin >> gn_tol;
	} while (gn_tol <= 0.0);

    // 3. === PRECONDITIONER ===

    if (pp.GetString ("LM_PRECON", cbuf)) {
        if (!strcasecmp (cbuf, "NONE"))
	    precon = LM_PRECON_NONE,  def = true;
	else if (!strcasecmp (cbuf, "HDIAG"))
	    precon = LM_PRECON_HDIAG, def = true;
	else if (!strcasecmp (cbuf, "CH"))
	    precon = LM_PRECON_CH,    def = true;
	else if (!strcasecmp (cbuf, "ICH"))
	    precon = LM_PRECON_ICH,   def = true;
	else if (!strcasecmp (cbuf, "GMRES"))
	    precon = LM_PRECON_GMRES, def = true;
	else if (!strcasecmp (cbuf, "BICGSTAB"))
	    precon = LM_PRECON_BICGSTAB, def = true;
	else if (!strcasecmp (cbuf, "GMRES_JFREE"))
	    precon = LM_PRECON_GMRES_JACOBIANFREE, def = true;
	else if (!strcasecmp (cbuf, "GMRES_DIRECT"))
	    precon = LM_PRECON_GMRES_DIRECT, def = true;
    }
    while (!def) {
        int cmd;
	cout << "\nSelect LM preconditioner:\n";
	cout << "(0) None\n";
	cout << "(1) Diagonal of Hessian\n";
	cout << "(2) Cholesky factorisation of full Hessian\n";
	cout << "(3) Incomplete CH factorisation of sparse Hessian\n";
	cout << "(4) BiCGSTAB solver (implicit Hessian)\n";
	cout << "(5) GMRES solver (implicit Hessian)\n";
	cout << "(6) GMRES solver (implicit Jacobian)\n";
	cout << "(7) GMRES solver (direct method)\n";
	cout << "[0|1|2|3|4|5|6|7] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: precon = LM_PRECON_NONE,  def = true; break;
	case 1: precon = LM_PRECON_HDIAG, def = true; break;
	case 2: precon = LM_PRECON_CH,    def = true; break;
	case 3: precon = LM_PRECON_ICH,   def = true; break;
	case 4: precon = LM_PRECON_BICGSTAB, def = true; break;
	case 5: precon = LM_PRECON_GMRES, def = true; break;
	case 6: precon = LM_PRECON_GMRES_JACOBIANFREE, def = true; break;
	case 7: precon = LM_PRECON_GMRES_DIRECT, def = true; break;
	}
    }

    // === HESSIAN SCALING STRATEGY
    def = false;
    if (precon != LM_PRECON_GMRES_DIRECT) {
	hscale = LM_HSCALE_IMPLICIT, def = true;
    } else if (pp.GetString ("LM_HESS_SCALING", cbuf)) {
	if (!strcasecmp (cbuf, "LM_HSCALE_NONE"))
	    hscale = LM_HSCALE_NONE, def = true;
	else if (!strcasecmp (cbuf, "LM_HSCALE_IMPLICIT"))
	    hscale = LM_HSCALE_IMPLICIT, def = true;
	else if (!strcasecmp (cbuf, "LM_HSCALE_EXPLICIT"))
	    hscale = LM_HSCALE_EXPLICIT, def = true;
	else if (!strcasecmp (cbuf, "LM_HSCALE_DISTANCE"))
	    hscale = LM_HSCALE_DISTANCE, def = true;;
    }
    while (!def) {
	int cmd;
	cout << "\nSelect Hessian scaling strategy:\n";
	cout << "(0) None\n";
	cout << "(1) Implicit (will generate adjoint fields\n";
	cout << "(2) Explicit\n";
	cout << "(3) Boundary distance scaling\n";
	cout << "[0|1|2|3] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: hscale = LM_HSCALE_NONE,     def = true; break;
	case 1: hscale = LM_HSCALE_IMPLICIT, def = true; break;
	case 2: hscale = LM_HSCALE_EXPLICIT, def = true; break;
	case 3: hscale = LM_HSCALE_DISTANCE, def = true; break;
	}
    }

    // 4. === GMRES CONVERGENCE CRITERION ===

    if (precon == LM_PRECON_GMRES || 
        precon == LM_PRECON_BICGSTAB ||
        precon == LM_PRECON_GMRES_JACOBIANFREE ||
	precon == LM_PRECON_GMRES_DIRECT) {
        if (!pp.GetReal ("LM_GMRES_TOL", gmres_tol) ||
	  gmres_tol <= 0.0) do {
	    cout << "\nSelect convergence criterion for Hessian solver (>0):\n";
	    cout << ">> ";
	    cin >> gmres_tol;
	} while (gmres_tol <= 0.0);
    }

    // 5. === LINE SEARCH ===

    if (!pp.GetBool ("LM_LINESEARCH", do_linesearch)) {
	char cmd;
	cout << "\nPerform line search in GN solver?\n";
	cout << "(y|n) >> ";
	cin >> cmd;
	do_linesearch = (toupper (cmd) == 'Y');
    }
    if (do_linesearch) {
	if (!pp.GetReal ("LM_MIN_STEPLENGTH", alpha_min) || alpha_min < 0.0)
	    do {
		cout << "\nMin. step length for line search\n";
		cout << "(>=0, 0:accept all step lengths):\n>> ";
		cin >> alpha_min;
	    } while (alpha_min < 0.0);
    }
    if (!pp.GetReal ("LM_INIT_STEPLENGTH", alpha0) || alpha0 <= 0.0)
	do {
	    cout << "\nInitial/fixed step length (>0)\n";
	    cin  >> alpha0;
	} while (alpha0 <= 0.0);

    // 6. === CONTROL PARAMETER LAMBDA ===

    if (!pp.GetReal ("LM_LAMBDA_0", lambda0) || lambda0 < 0.0) do {
	cout << "\nInitial value of control parameter lambda (>=0):\n";
	cout << ">> ";
	cin >> lambda0;
    } while (lambda0 < 0.0);
    if (!pp.GetReal ("LM_LAMBDA_SCALE", lambda_scale) || lambda_scale < 1.0) {
	cout << "\nScaling factor for control parameter lambda adjustment\n";
	cout << "(>=1, 1:keep lambda fixed):\n>> ";
	cin >> lambda_scale;
    } while (lambda_scale < 1.0);

    // 9. === JACOBIAN FILL FRACTION ===

    if (!pp.GetReal ("JACOBIAN_FILL", Jfill)) {
	cout << "\nJacobian fill fraction (1 for dense):\n";
	cout << ">> ";
	cin >> Jfill;
    }

    // 10. === COVARIANCE DATA ===

    if (!pp.GetBool ("MODEL_CORRECTION", modelerr)) {
	cout << "\nUse covariance data for model correction?\n";
	cout << "[y|n]: >> ";
	do {
	    cin >> c;
	} while (c != 'y' && c != 'n');
	modelerr = (c == 'y');
    }

    if (modelerr) {

	if (!pp.GetReal ("MODEL_LAMBDA", mcov_lambda)) {
	    cout << "Covariance diagonal scaling factor lambda: ";
	    cin >> mcov_lambda;
	}

	RDenseMatrix tmp;
	int i, j, n;
	if (!pp.GetString ("FMOD_COVARIANCE", mcov_fmod_fname)) {
	    cout << "Covariance file for modulation: \n>> ";
	    cin >> mcov_fmod_fname;
	}
	ifstream ifs1 (mcov_fmod_fname);
	ifs1 >> tmp;
	n = tmp.nRows();
	for (i = 0; i < n; i++) tmp(i,i) += mcov_lambda*mcov_lambda;
	tmp = inverse (tmp);
	mcov.New (2*n,2*n);
	for (i = 0; i < n; i++)
	    for (j = 0; j <= i; j++)
		mcov(i,j) = tmp(i,j);

	if (!pp.GetString ("FARG_COVARIANCE", mcov_farg_fname)) {
	    cout << "Covariance file for phase: \n>> ";
	    cin >> mcov_farg_fname;
	}
	ifstream ifs2 (mcov_farg_fname);
	ifs2 >> tmp;
	if (n != tmp.nRows()) {
	    cerr << "inconsistent dimension" << endl;
	    exit (0);
	}
	for (i = 0; i < n; i++) tmp(i,i) += mcov_lambda*mcov_lambda;
	tmp = inverse (tmp);
	for (i = 0; i < n; i++)
	    for (j = 0; j <= i; j++)
		mcov(n+i,n+j) = tmp(i,j);
	    
	CHdecomp (mcov);

	if (!pp.GetString ("FMOD_MODEL_ERROR", merr_fmod_fname)) {
	    cout << "Model error data for modulation: \n>> ";
	    cin >> merr_fmod_fname;
	}
	merr.New(2*n);
	RVector merr_fmod(merr,0,n);
	ifstream ifs3 (merr_fmod_fname);
	ifs3 >> merr_fmod;
	
	if (!pp.GetString ("FARG_MODEL_ERROR", merr_farg_fname)) {
	    cout << "Model error data for phase: \n>> ";
	    cin >> merr_farg_fname;
	}
	RVector merr_farg(merr,n,n);
	ifstream ifs4 (merr_farg_fname);
	ifs4 >> merr_farg;
    }
}

void SolverLM::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "LM");
    pp.PutReal ("NONLIN_TOL", gn_tol);
    pp.PutInt ("ITMAX", nrmax);

    switch (precon) {
    case LM_PRECON_NONE:
	pp.PutString ("LM_PRECON", "NONE");
	break;
    case LM_PRECON_HDIAG:
	pp.PutString ("LM_PRECON", "HDIAG");
	break;
    case LM_PRECON_CH:
	pp.PutString ("LM_PRECON", "CH");
	break;
    case LM_PRECON_ICH:
	pp.PutString ("LM_PRECON", "ICH");
	break;
    case LM_PRECON_BICGSTAB:
	pp.PutString ("LM_PRECON", "BICGSTAB");
	pp.PutReal ("LM_GMRES_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES:
	pp.PutString ("LM_PRECON", "GMRES");
	pp.PutReal ("LM_GMRES_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES_JACOBIANFREE:
	pp.PutString ("LM_PRECON", "GMRES_JFREE");
	pp.PutReal ("LM_GMRES_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES_DIRECT:
	pp.PutString ("LM_PRECON", "GMRES_DIRECT");
	pp.PutReal ("LM_GMRES_TOL", gmres_tol);
	break;
    }

    switch (hscale) {
    case LM_HSCALE_IMPLICIT:
	pp.PutString ("LM_HESS_SCALING", "LM_HSCALE_IMPLICIT");
	break;
    case LM_HSCALE_EXPLICIT:
	pp.PutString ("LM_HESS_SCALING", "LM_HSCALE_EXPLICIT");
	break;
    case LM_HSCALE_DISTANCE:
	pp.PutString ("LM_HESS_SCALING", "LM_HSCALE_DISTANCE");
	break;
    }

    pp.PutBool ("LM_LINESEARCH", do_linesearch);
    if (do_linesearch)
	pp.PutReal ("LM_MIN_STEPLENGTH", alpha_min);
    pp.PutReal ("LM_INIT_STEPLENGTH", alpha0);

    pp.PutReal ("LM_LAMBDA_0", lambda0);
    pp.PutReal ("LM_LAMBDA_SCALE", lambda_scale);
    pp.PutReal ("JACOBIAN_FILL", Jfill);

    pp.PutBool ("MODEL_CORRECTION", modelerr);
    if (modelerr) {
	pp.PutReal ("MODEL_LAMBDA", mcov_lambda);
	pp.PutString ("FMOD_COVARIANCE", mcov_fmod_fname);
	pp.PutString ("FARG_COVARIANCE", mcov_farg_fname);
	pp.PutString ("FMOD_MODEL_ERROR", merr_fmod_fname);
	pp.PutString ("FARG_MODEL_ERROR", merr_farg_fname);
    }
}

// ==========================================================================

bool LineSearchWithPrior (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, const RVector &p0, const RVector &p,
    const Regularisation *reg, const RMatrix *cov)
{
    const int MAXIT = 16;
    double x0 = 0.0, xm, fm, fmd, fmp, fb, fbd, fbp, fmind, fminp,xb = lambda;
    RVector proj(data.Dim());
    
    RVector p1 = p+grad*xb;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    if (CheckRange (meshsol)) {
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
	fbp = reg->GetValue (p1);
	fb  = fbd + fbp;
	LOGOUT("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g", xb,fbd,fbp,fb);
    } else {
	LOGOUT ("Parameters out of range in trial step");
	fb = f0*4.0; // force reduction in step size
    }

    if (fb < f0) { // increase interval
        xm = xb; fm = fb;
	xb *= 2.0;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	if (CheckRange (meshsol)) {
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
	    fbp = reg->GetValue (p1);
	    fb  = fbd + fbp;
	    LOGOUT("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
		     xb, fbd, fbp, fb);
	} else {
	    LOGOUT ("Parameters out of range in trial step");
	    fb = fm*4.0; // stop growing the interval
	}

	while (fb < fm) {
	    x0 = xm; f0 = fm;
	    xm = xb; fm = fb;
	    xb *= 2.0;
	    p1 = p+grad*xb;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    if (CheckRange (meshsol)) {
		proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
		fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
		fbp = reg->GetValue (p1);
		fb  = fbd + fbp;
		LOGOUT("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
			 xb, fbd, fbp, fb);
	    } else {
		LOGOUT ("Parameters out of range in trial step");
		fb = fm*4.0; // stop growing the interval
	    }
	}
    } else { // decrease interval
        xm = 0.5*xb;
	p1 = p+grad*xm;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	if (CheckRange (meshsol)) {
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fmd = ObjectiveFunction::get_value (data, proj, sd, cov);
	    fmp = reg->GetValue (p1);
	    fm  = fmd + fmp;
	    LOGOUT("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
		     xm, fmd, fmp, fm);
	} else {
	    LOGOUT ("Parameters out of range in trial step");
	    fm = f0*4.0; // force reduction of interval
	}
	int itcount = 0;
	while (fm > f0) {
  	    if (++itcount > MAXIT) return false;
	    xb = xm; fb = fm;
	    xm = 0.5*xb;
	    p1 = p+grad*xm;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    if (CheckRange (meshsol)) {
		proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
		fmd = ObjectiveFunction::get_value (data, proj, sd, cov);
		fmp = reg->GetValue (p1);
		fm  = fmd + fmp;
		LOGOUT("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
			 xm, fmd, fmp, fm);
	    } else {
		LOGOUT ("Parameters out of range in trial step");
		fm = f0*4.0; // force reduction of interval
	    }
	}
    }
    // quadratic interpolation
    double a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
    double b = (f0-fb)/(x0-xb) - a*(x0+xb);
    lambda = -b/(2.0*a);
    p1 = p+grad*lambda;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fmind = ObjectiveFunction::get_value (data, proj, sd, cov);
    fminp = reg->GetValue (p1);
    fmin  = fmind + fminp;
    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
    }
    LOGOUT("Lsearch final: STEP %g OF %g PRIOR %g TOTAL %g",
		lambda, fmind, fminp, fmin);
    // restimate tau 
    //    tau = fmind/fminp;
 
    return true;
}

bool CheckRange (const Solution &sol)
{
    bool inrange = true;

    const double MIN_CMUA = 0;
    const double MAX_CMUA = 0.1;
    const double MIN_CKAPPA = 0;
    const double MAX_CKAPPA = 10;

    double vmin, vmax;
    sol.Extents (OT_CMUA, vmin, vmax);
    if (vmin < MIN_CMUA || vmax > MAX_CMUA) {
	cerr << "WARNING: " << vmin << " < CMUA < " << vmax
	     << " in trial solution" << endl;
	inrange = false;
    }
    sol.Extents (OT_CKAPPA, vmin, vmax);
    if (vmin < MIN_CKAPPA || vmax > MAX_CKAPPA) {
	cerr << "WARNING: " << vmin << " < CKAPPA < " << vmax
	     << " in trial solution" << endl;
	inrange = false;
    }
    return inrange;
}

// ==========================================================================
// Build the Hessian of the prior from individual parameter contributions

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
    int nprm)
{
    int n = x.Dim();
    int n0 = n/nprm;
    int i, j;

    RCompRowMatrix H, Hi, Hij;
    for (i = 0; i < nprm; i++) {
	for (j = 0; j < nprm; j++) {
	    if (!j) {
		Hi.New(n0,n0);
		if (j==i) reg->SetHess1 (Hi, x, j);
	    } else {
		Hij.New(n0,n0);
		if (j==i) reg->SetHess1 (Hij, x, j);
		Hi = cath (Hi, Hij);
	    }
	}
	if (!i) H = Hi;
	else    H = catv (H, Hi);
    }
    return H;
}

// ==========================================================================

RVector RescaleHessian (const RDenseMatrixMPI &Jmod,
    const RDenseMatrixMPI &Jarg, const RVector &x, RCompRowMatrix *RHess)
{
    int i, n = Jmod.nCols();
    RVector M(n);

    // J^T J contribution
    M += Jmod.ColSumSq();
    M += Jarg.ColSumSq();

    // Psi'' contribution
    if (RHess) M += RHess->Diag();

    for (i = 0; i < n; i++) {
	M[i] = 1.0/sqrt (M[i]);
    }

    return M;
}
