// ==========================================================================
// SolverLM: Levenberg-Marquardt
// CW, multi-wavelength version
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "mwsolution.h"
#include "solverlm2_cw_mw.h"
#include "supertoast_util.h"
#include "supertoast_cw_mw.h"
#include "timing.h"
#include <time.h>
#include "timing.h"

#define MAX_NOFWLENGTH 50

using namespace std;

// ==========================================================================
// external references

extern int g_imgfmt;
extern char g_prefix[256];
extern double clock0;
extern double g_refind;

// Implicit structure of multi-wavelength Jacobian
// only the mua and kappa versions of J for all wavelengths are stored,
// together with the extinction coefficients
struct MWJacobian {
    const MWsolution *msol;        // MW solution in mesh basis
    const Solution *bsol;          // solution in reconstruction basis
    RDenseMatrix Jmua[MAX_NOFWLENGTH]; // mua part of Jacobian
    RDenseMatrix Jkap[MAX_NOFWLENGTH]; // kappa part of Jacobian
    RVector scale_A[MAX_NOFWLENGTH];
    RVector scale_b[MAX_NOFWLENGTH];
    RVector prmscale;
    RVector hscale;
    RVector datascale;

    RVector Ax (const RVector &x) const;
    RVector ATx (const RVector &y) const;
    RVector RescaleHessian (const RVector &x, RCompRowMatrix *RHess) const;
};

// ==========================================================================
// local prototypes

void GenerateJacobian (RFwdSolverMW &FWS, const Raster &raster,
    const QMMesh *mesh,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    DataScale dscale, MWJacobian *J);

static RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
    int nprm);

static RVector RescaleHessian (const RMatrix *J, const RVector &x,
    RCompRowMatrix *RHess);

static RVector Mergewavels (const RVector &b, int r0, int len, int nofwavel);

static RDenseMatrix mul (const RMatrix &A, const RMatrix &B)
{
    int i, j, k;
    double sum;
    int ar = A.nRows();
    int ac = A.nCols();
    int br = B.nRows();
    int bc = B.nCols();
    xASSERT (ac == br, "Incompatible dimensions");
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
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const MWJacobian *J;         // Jacobian
    //const RVector *M;            // normalisation diagonal matrix
    const double *lambda;        // diagonal scaling factor
    const Regularisation *reg;   // regularisation
    const RCompRowMatrix *RHess; // Hessian of prior (1 matrix per parameter)
    const RVector *Acoeff;       // scattering prefactor
    const RVector *bcoeff;       // scattering power
    const RDenseMatrix *extcoef; // extinction coefficients
    int activeWl;                // active wavelength index
};

// ==========================================================================
// This data structure defines the Jacobian implicitly

struct JAC_DATA {
    const QMMesh *mesh;            // mesh
    const Raster *raster;          // basis mapper
    const RFwdSolverMW *FWS;       // forward solver
          RVector **dphi;          // fields for each source in mesh basis
          RVector **aphi;          // adjoint fields
    const RCompRowMatrix *mvec;    // measurement vectors
    const RVector *dat;            // measurements
    const RVector *sd;             // sd
    const RVector *prmscl;         // parameter scaling vector
    const Regularisation *reg;     // regularisation
    const RCompRowMatrix *RHess;   // Hessian of prior (1 matrix per parameter)
    const RVector *M;              // normalisation diagonal matrix
    const double *lambda;          // diagonal scaling factor
    const MWsolution *sol;         // solution instance
    const RDenseMatrix *excoef;    // extinction coefficients
};

// ==========================================================================

static RSymMatrix Hess_full (const HESS_DATA &data, const RVector &x)
{
    ERROR_UNDEF;
    return RSymMatrix();

#ifdef UNDEF
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
#endif
}

// ==========================================================================

static RVector Hess_diag (const HESS_DATA &data, const RVector &x)
{
    ERROR_UNDEF;
    return RVector();

#ifdef UNDEF
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
#endif
}

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
    const MWJacobian *J         =  data->J;
    const RCompRowMatrix *RHess =  data->RHess;
    const RVector &M            = J->hscale;
    const double lambda         = *data->lambda;
    int n = x.Dim();

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

    RVector JTJx = J->ATx (J->Ax (x));
   
    return JTJx + Px + lambda*x;
}

// ==========================================================================
// Generate the Frechet derivatives: direct method
// This version only requires the availability of the direct fields
// (dphi_h) passed in mesh basis
// Note: parameter and data scaling not yet implemented

static RVector FrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const RCompRowMatrix &mvec, const RVector *dphi_h,
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
    RVector qdelta_g(glen);
    RVector qdelta_h(nlen);
    RVector phi_h_delta(nlen);
    RVector dphi_p(slen);
    RVector dphi_g(glen);
    RVector *dphi_g_grad = new RVector[dim];
    RVector *dphi_g_gradgrad = new RVector[dim];
    RVector tmp(slen);
    RVector y_h_delta (nqm*2);

    for (i = 0; i < dim; i++) {
	dphi_g_grad[i].New (glen);
	dphi_g_gradgrad[i].New (glen);
    }

    RVector alpha_h(x, 0, slen);
    RVector beta_h(x, slen, slen);
    RVector alpha_g(glen), beta_g(glen);
    //RVector calpha_g(glen), cbeta_g(glen);
    raster.Map_SolToGrid (alpha_h, alpha_g);
    raster.Map_SolToGrid (beta_h, beta_g);
    //SetReal (calpha_g, alpha_g);
    //SetReal (cbeta_g, beta_g);

    for (q = 0; q < nq; q++) {
	// build up perturbation rhs (qdelta_g)
	raster.Map_MeshToGrid (dphi_h[q], dphi_g);
	qdelta_g = -alpha_g * dphi_g;
	ImageGradient (raster.GDim(), raster.GSize(), dphi_g, dphi_g_grad, 
		       raster.Elref());
	for (i = 0; i < dim; i++) {
	    dphi_g_grad[i] *= beta_g;
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
	qdelta_g *= 1.0/scl;

	// remap to mesh basis
	raster.Map_GridToMesh (qdelta_g, qdelta_h);

	// scale with element support size
	RVector escl(nlen);
	for (i = 0; i < mesh.elen(); i++) {
	    Element *pel = mesh.elist[i];
	    for (j = 0; j < pel->nNode(); j++) {
		escl[pel->Node[j]] += pel->IntF(j);
	    }
	}
	qdelta_h *= escl;

	// solve for phi_delta
	FWS.CalcField (qdelta_h, phi_h_delta);

	// boundary projections
	int nm = mesh.nQMref[q];
	RVector projdelta = ProjectSingle (&mesh, q, mvec, phi_h_delta);

	if (FWS.GetDataScaling() == DATA_LOG) {
	    RVector proj = FWS.ProjectSingle (q, mvec, dphi_h[q], DATA_LIN);
	    for (i = 0; i < nm; i++) {
		double scl = norm2 (proj[i]);
		double c = projdelta[i]*proj[i];
		projdelta[i] = c/scl;
	    }
	}

	// map data types into real vector
	for (i = 0; i < nm; i++) {
	    double proji = projdelta[i];
	    y_h_delta[ofs_lnmod++] = proji;
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

static RVector FrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const MWsolution &sol, const RCompRowMatrix &mvec,
    const RVector *dphi_h, const RVector *aphi_h, const RVector &sd,
    const RVector &x, const RDenseMatrix &excoef)
{
    if (!aphi_h)
	return FrechetDerivative (mesh,raster,FWS,mvec, dphi_h, sd, x);
    // no adjoint fields - use direct method

    tic();
    int i, j, w, q, m, ofs, mofs;
    int n    = x.Dim();
    int nq   = mesh.nQ;
    int nm   = mesh.nM;
    int nqm  = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    int nofwavel = sol.nofwavel;
    int nchromo = excoef.nCols();
    int nprm    = sol.nParam();
    bool bFactorA = sol.IsActive(nchromo);
    bool bPowerb  = sol.IsActive(nchromo+1);
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector res(nqm*nofwavel);

    // multiply the direct fields with the rhs before
    // multiplying in the adjoint fields

    RVector pd(glen), pa(glen);

    RVector *pd_grad = new RVector[dim];
    RVector *pa_grad = new RVector[dim];
    RVector pas(slen), *pas_grad = new RVector[dim];
    RVector *pds = new RVector[nq];
    RVector **pds_grad = new RVector*[nq];
    for (i = 0; i < nq; i++) pds_grad[i] = new RVector[dim];
    RVector *proj = new RVector[nq];
    //int *qofs = new int[nq];
    RVector sA, sb;

    for (i = 0; i < dim; i++)
	pas_grad[i].New(slen);

    // loop over wavelengths
    for (w = 0; w < nofwavel; w++) {
	const RVector *wdphi_h = dphi_h + nq*w;
	const RVector *waphi_h = aphi_h + nm*w;
	if (bFactorA)
	    raster.Map_MeshToSol (sol.GetJacobianCoeff_A(w), sA);
	if (bPowerb)
	    raster.Map_MeshToSol (sol.GetJacobianCoeff_b(w), sb);

	// loop over sources
	for (q = 0; q < nq; q++) {
	    //qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	    raster.Map_MeshToGrid (wdphi_h[q], pd);
	    ImageGradient (gdim, gsize, pd, pd_grad, raster.Elref());
	    
	    if (FWS.GetDataScaling() == DATA_LOG) {
		proj[q].New (mesh.nQMref[q]);
		proj[q] = FWS.ProjectSingle (q, mvec, wdphi_h[q], DATA_LIN);
	    }
	    pds[q].New (slen);
	    raster.Map_GridToSol (pd, pds[q]);
	    for (i = 0; i < dim; i++) {
		pds_grad[q][i].New(slen);
		raster.Map_GridToSol (pd_grad[i], pds_grad[q][i]);
	    }
	}
	// loop over detectors
	for (m = 0; m < mesh.nM; m++) {
	    raster.Map_MeshToGrid (waphi_h[m], pa);
	    raster.Map_GridToSol (pa, pas);
	    ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());
	    for (i = 0; i < dim; i++)
		raster.Map_GridToSol (pa_grad[i], pas_grad[i]);
	    
	    for (q = 0; q < nq; q++) {
		if (!mesh.Connected (q,m)) continue;
		ofs = mesh.QMofs[q][m];
		mofs = ofs-mesh.Qofs[q];
		//mm = mesh.QMref[q][m];
		//ofs = w*nqm + qofs[q] + mm;
		res[ofs+w*nqm] = 0;
		int prmofs = 0;
		// sum over chromophores
		for (j = 0; j < nchromo; j++) {
		    RVector xj(x,prmofs,slen);
		    res[ofs+w*nqm] -= ((pds[q]*xj) & pas) * excoef(w,j);
		    prmofs += slen;
		}
		// scattering prefactor
		if (bFactorA) {
		    RVector xj(x,prmofs,slen);
		    xj *= sA;
		    for (i = 0; i < dim; i++)
			res[ofs+w*nqm] -= (pds_grad[q][i]*xj) & pas_grad[i];
		    prmofs += slen;
		}
		// scattering power
		if (bPowerb) {
		    RVector xj(x,prmofs,slen);
		    xj *= sb;
		    for (i = 0; i < dim; i++)
			res[ofs+w*nqm] -= (pds_grad[q][i]*xj) & pas_grad[i];
		    prmofs += slen;
		}
		if (FWS.GetDataScaling() == DATA_LOG) {
		    double scl = norm2 (proj[q][mofs]);
		    double c = res[ofs+w*nqm]*proj[q][mofs];
		    res[ofs+w*nqm] = c/scl;
		}
	    }
	}
    }

    // scale with sd
    res /= sd;

    delete []pd_grad;
    delete []pa_grad;
    delete []pas_grad;
    delete []pds;
    for (i = 0; i < nq; i++) delete []pds_grad[i];
    delete []pds_grad;
    delete []proj;
    //delete []qofs;

    cerr << "Frechet timing: " << toc() << endl;
    return res;
}

// ==========================================================================
// Adjoint Frechet derivative:
// Direct method (requires only direct fields)

static RVector AdjointFrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const MWsolution &sol, const RCompRowMatrix &mvec,
    const RVector *dphi_h, const RVector &sd, const RVector &y,
    const RDenseMatrix &excoef)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    double term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    RVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    RVector dgrad (slen);
    ofs_mod = 0;         // data offset for Mod data
    ofs_arg = mesh.nQM;  // data offset for Arg data
    RVector rres (slen*2);
    RVector grad_cmua(rres, 0, slen);       // mua part of grad
    RVector grad_ckappa (rres, slen, slen); // kappa part of grad
    
    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        RVector cdfield (glen);
        RVector *cdfield_grad = new RVector[dim];
	raster.Map_MeshToGrid (dphi_h[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh.nQMref[q];

	RVector y_mod (y, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector b_mod(n);
	b_mod = y_mod/s_mod;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	RVector cproj(n);
	cproj = FWS.ProjectSingle (q, mvec, dphi_h[q], DATA_LIN);
	wqa = 0.0;
	wqb = 0.0;

	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const RVector qs = mvec.Row(m);
	    double rp = cproj[idx];
	    double dn = 1.0/(rp*rp);

	    // amplitude term
	    term = /* -2.0 * */ b_mod[idx] /* / (ype[idx]*s_mod[idx]) */;
	    wqa += qs * (term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	RVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);

	RVector cafield(glen);
	RVector *cafield_grad = new RVector[dim];
	raster.Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster.Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua -= dgrad;

	// diffusion contribution
	// multiply complex field gradients
	RVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster.Map_GridToSol (gk, dgrad);
	grad_ckappa -= dgrad;

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

static RVector AdjointFrechetDerivative (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const MWsolution &sol, const RCompRowMatrix &mvec,
    const RVector *dphi_h, const RVector *aphi_h, const RVector &sd,
    const RVector &y, const RDenseMatrix &excoef)
{
    if (!aphi_h) 
	return AdjointFrechetDerivative (mesh, raster, FWS, sol, mvec, dphi_h,
					 sd, y, excoef);
    // no adjoint fields - use direct method

    tic();
    int i, j, w, q, m, ofs, mofs;
    int nq   = mesh.nQ;
    int nm   = mesh.nM;
    int nqm  = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    int nofwavel = sol.nofwavel;
    int nchromo = excoef.nCols();
    int nprm    = sol.nParam();
    int nactive = sol.nActive();
    bool bFactorA = sol.IsActive(nchromo);
    bool bPowerb  = sol.IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;

    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector rres(slen*nactive);
    RVector ysd = y/sd; // data scaling

    RVector pd(glen), *pds = new RVector[nq], pa(glen), pas(slen);
    RVector das(slen), dbs(slen);
    RVector db(glen);
    RVector sA, sb;

    RVector *pd_grad = new RVector[dim];
    RVector **pds_grad = new RVector*[nq];
    for (i = 0; i < nq; i++) pds_grad[i] = new RVector[dim];
    RVector *pa_grad = new RVector[dim];
    RVector *pas_grad = new RVector[dim];
    RVector *proj = new RVector[nq];
    //int *qofs = new int[nq];

    for (i = 0; i < dim; i++)
	pas_grad[i].New(slen);

    // loop over wavelengths
    for (w = 0; w < nofwavel; w++) {
	const RVector *wdphi_h = dphi_h + nq*w;
	const RVector *waphi_h = aphi_h + nm*w;

	// loop over sources
	for (q = ofs = 0; q < nq; q++) {
	    //qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	    raster.Map_MeshToGrid (wdphi_h[q], pd);
	    ImageGradient (gdim, gsize, pd, pd_grad, raster.Elref());
	    
	    if (FWS.GetDataScaling() == DATA_LOG) {
		proj[q].New (mesh.nQMref[q]);
		proj[q] = FWS.ProjectSingle (q, mvec, wdphi_h[q], DATA_LIN);
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
	    raster.Map_MeshToGrid (waphi_h[m], pa);
	    ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());
	    raster.Map_GridToSol (pa, pas);
	    for (i = 0; i < dim; i++)
		raster.Map_GridToSol (pa_grad[i], pas_grad[i]);
	    
	    for (q = 0; q < nq; q++) {
		if (!mesh.Connected (q,m)) continue;
		ofs = mesh.QMofs[q][m];
		mofs = ofs-mesh.Qofs[q];
		//mm = mesh.QMref[q][m];
		//ofs = nqm*w + qofs[q] + mm;
		das = pds[q] * pas;
		dbs.Clear();
		for (i = 0; i < dim; i++)
		    dbs += pds_grad[q][i] * pas_grad[i];
		double y = ysd[ofs+w*nqm];
		if (FWS.GetDataScaling() == DATA_LOG) {
		    // rescale for log data (lnamp + phase)
		    double logscl = proj[q][mofs]/norm2 (proj[q][mofs]);
		    double logdat = y * logscl;
		    y = logdat;
		}

		// chromophore blocks
		int prmofs = 0;
		for (j = 0; j < nchromo; j++) {
		    double ex = excoef (w,j);
		    for (i = 0; i < slen; i++)
			rres[prmofs+i] -= das[i]*ex*y;
		    prmofs += slen;
		}
		// scattering prefactor
		if (bFactorA) {
		    raster.Map_MeshToSol (sol.GetJacobianCoeff_A(w), sA);
		    for (i = 0; i < slen; i++)
			rres[prmofs+i] -= dbs[i]*sA[i]*y;
		    prmofs += slen;
		}
		// scattering power
		if (bPowerb) {
		    raster.Map_MeshToSol (sol.GetJacobianCoeff_b(w), sb);
		    for (i = 0; i < slen; i++)
			rres[prmofs+i] -= dbs[i]*sb[i]*y;
		    prmofs += slen;
		}
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
    //delete []qofs;

    cerr << "Adjoint Frechet timing: " << toc() << endl;
    return rres;
}

// ==========================================================================
// Generate diagonal scaling matrix M that scales the diagonal of the
// Hessian to 1. This version is for implicit Jacobian given in terms of
// direct and adjoint fields

static RVector ImplicitHessianScaling (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const MWsolution &sol, const RCompRowMatrix &mvec,
    const RVector *dphi_h, const RVector *aphi_h, const RVector &sd,
    const RVector &prmscl, const RDenseMatrix &excoef)
{
    int i, j, w, q, m, ofs, mofs;
    int nq   = mesh.nQ;
    int nm   = mesh.nM;
    int nqm  = mesh.nQM;
    int nlen = mesh.nlen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = mesh.Dimension();
    int nprm = sol.nParam();
    int nofwavel = sol.nofwavel;
    int nchromo  = excoef.nCols();
    bool bFactorA = sol.IsActive(nchromo);
    bool bPowerb  = sol.IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;

    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    RVector M(slen*sol.nActive());

    RVector *pd = new RVector[nq];
    RVector pa(glen);
    RVector das(slen), dbs(slen);
    RVector db(glen);
    RVector sA, sb;

    RVector **pd_grad = new RVector*[nq];
    for (i = 0; i < nq; i++) pd_grad[i] = new RVector[dim];
    RVector *pa_grad = new RVector[dim];
    RVector *proj = new RVector[nq];
    //int *qofs = new int[nq];
    
    // loop over wavelengths
    for (w = 0; w < nofwavel; w++) {
	const RVector *wdphi_h = dphi_h + w*nq;
	const RVector *waphi_h = aphi_h + w*nm;
	if (bFactorA)
	    raster.Map_MeshToSol (sol.GetJacobianCoeff_A(w), sA);
	if (bPowerb)
	    raster.Map_MeshToSol (sol.GetJacobianCoeff_b(w), sb);

	// loop over sources
	for (q = 0; q < nq; q++) {
	    //qofs[q] = (q ? qofs[q-1]+mesh.nQMref[q-1] : 0);
	    pd[q].New(glen);
	    for (i = 0; i < dim; i++) pd_grad[q][i].New(glen);
	    raster.Map_MeshToGrid (wdphi_h[q], pd[q]);
	    ImageGradient (gdim, gsize, pd[q], pd_grad[q], raster.Elref());

	    if (FWS.GetDataScaling() == DATA_LOG) {
		proj[q].New (mesh.nQMref[q]);
		proj[q] = FWS.ProjectSingle (q, mvec, wdphi_h[q], DATA_LIN);
	    }
	}

	// loop over detectors
	for (m = 0; m < nm; m++) {
	    raster.Map_MeshToGrid (waphi_h[m], pa);
	    ImageGradient (gdim, gsize, pa, pa_grad, raster.Elref());

	    for (q = 0; q < nq; q++) {
		if (!mesh.Connected (q,m)) continue;
		ofs = mesh.QMofs[q][m];
		mofs = ofs-mesh.Qofs[q];
		//mm = mesh.QMref[q][m];
		//ofs = qofs[q] + mm;
		RVector da = pd[q] * pa;
		raster.Map_GridToSol (da, das);
		db.Clear();
		for (i = 0; i < dim; i++)
		    db += pd_grad[q][i] * pa_grad[i];
		raster.Map_GridToSol (db, dbs);
		if (FWS.GetDataScaling() == DATA_LOG) {
		    double scl = norm2 (proj[q][mofs]);
		    das *= (proj[q][mofs]/scl);
		    dbs *= (proj[q][mofs]/scl);
		}
		for (i = 0; i < slen; i++) { // data scaling
		    das[i] /= sd[ofs+w*nqm];
		    dbs[i] /= sd[ofs+w*nqm];
		}
		// chromophore blocks
		int prmofs = 0;
		for (j = 0; j < nchromo; j++) {
		    double ex = excoef (w,j);
		    for (i = 0; i < slen; i++)
			M[prmofs+i] += norm2(das[i]*ex);
		    prmofs += slen;
		}
		// scattering prefactor
		if (bFactorA) {
		    for (i = 0; i < slen; i++)
			M[prmofs+i] += norm2(dbs[i]*sA[i]);
		    prmofs += slen;
		}
		// scattering power
		if (bPowerb) {
		    for (i = 0; i < slen; i++)
			M[prmofs+i] += norm2(dbs[i]*sb[i]);
		    prmofs += slen;
		}
	    }
	}
    }
    delete []pd;
    for (i = 0; i < nq; i++) delete []pd_grad[i];
    delete []pd_grad;
    delete []pa_grad;
    delete []proj;
    //delete []qofs;

    M *= prmscl*prmscl; // parameter scaling
    return M;
}

// ==========================================================================
// Generate diagonal scaling matrix M that scales the diagonal of the
// Hessian to 1. This version solves one direct and adjoint field

static RVector ExplicitHessianScaling1 (const QMMesh &mesh, const Raster &raster,
    const RFwdSolverMW &FWS, const RCompRowMatrix &qvec,
    const RCompRowMatrix &mvec,  const RVector &sd, const RVector &prmscl)
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

    RVector pd, pa(glen);
    RVector das(slen), dbs(slen);
    RVector db(glen);

    RVector *pd_grad = new RVector[dim];
    RVector *pa_grad = new RVector[dim];
    RVector *proj = new RVector[nq];
    int *qofs = new int[nq];

    RVector qsum(nlen), msum(nlen), dphi(nlen), aphi(nlen);

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
    RVector da = pd * pa;
    raster.Map_GridToSol (da, das);
    db.Clear();
    for (i = 0; i < dim; i++)
	    db += pd_grad[i] * pa_grad[i];
    raster.Map_GridToSol (db, dbs);
    for (i = 0; i < slen; i++) {
		M[i]      += norm2 (das[i]);
		//M[i+slen] += norm2 (dbs[i]);
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
    CFwdSolverMW &FWS, const CCompRowMatrix &qvec, const CCompRowMatrix &mvec)
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

static RVector DistanceField (const Raster &raster)
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

static RVector DistScaling (const Raster &raster, const RVector &logx)
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
    const RFwdSolverMW *FWS = data->FWS;     // forward solver
          RVector **dphi_h = data->dphi;     // fields in mesh basis
	  RVector **aphi_h = data->aphi;     // adjoint fields
    const RCompRowMatrix *mvec = data->mvec; // measurement vectors
    const RVector *dat = data->dat;         // data
    const RVector *sd = data->sd;            // sd
    const RVector *prmscl = data->prmscl;    // parameter scaling vector
    const RCompRowMatrix *RHess = data->RHess;
    const RVector &M  = *data->M;
    const double lambda = *data->lambda;
    const MWsolution *sol = data->sol;
    const RDenseMatrix *excoef = data->excoef;

    RVector xscl(x);
    xscl *= M;        // apply scaling Jx -> JMx
    xscl *= *prmscl;  // apply param scaling dy/dx -> dy/d(log x)

    // Frechet derivative z = F' x
    RVector y = FrechetDerivative (*mesh, *raster, *FWS, *sol, *mvec,
        *dphi_h, *aphi_h, *sd, xscl, *excoef);

    // Adjoint Frechet derivative F'* z
    RVector z = AdjointFrechetDerivative (*mesh, *raster, *FWS, *sol, *mvec,
        *dphi_h, *aphi_h, *sd, y, *excoef);

    z *= M;        // apply scaling J^T y -> M J^T y
    z *= *prmscl;  // apply param scaling dy/dx -> dy/d(log x)


    int n = x.Dim();
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
SolverLM2_CW_MW::SolverLM2_CW_MW (ParamParser *_pp): Solver_CW (_pp)
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
    reg = NULL;
}

SolverLM2_CW_MW::~SolverLM2_CW_MW ()
{
    if (reg) delete reg;
}

// ==========================================================================

void SolverLM2_CW_MW::Solve (RFwdSolverMW &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, MWsolution &msol,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec)
{
    RDenseMatrix extcoef = msol.extcoef;
    // extinction coefficients for all chromophores at all wavelengths

    int nlambda = msol.nofwavel;   // number of wavelengths
    int nch =     extcoef.nCols(); // number of chromophores

    const QMMesh *mesh = FWS.MeshPtr();

    bool Use_precon = true;
    bool Gradient_descent = false;

    // initialisations
    int i, inr;
    int nq   = mesh->nQ;
    int nm   = mesh->nM;
    int nqm  = mesh->nQM;
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int nprm = bsol.nActive();
    int n    = bsol.ActiveDim();
    double of, of_value, of_prior, fmin, errstart, err0, err00, err1;
    double lambda = lambda0, alpha = -1.0;

    RVector r(n), s(n), d(n);
    RVector h(n), dold(n);

    RVector *dphi = 0, *aphi = 0, *aphi_hscale = 0;
    RVector proj(ndat);
    RVector prmscl;
    
    Regularisation *reg;
    RCompRowMatrix RHess;
    
    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    // initial solution (scaled)
    
    RVector x(x0), x1(x0);
    // current solution, trial solution (scaled);

    // allocate field vectors for all wavelengths
    dphi = new RVector[nq*nlambda];
    for (i = 0; i < nq*nlambda; i++) dphi[i].New (nlen);
    if ((precon != LM_PRECON_GMRES_DIRECT) || hscale == LM_HSCALE_IMPLICIT) {
	aphi = new RVector[nm*nlambda];
	for (i = 0; i < nm*nlambda; i++) aphi[i].New (nlen);
    }

    // reset forward solver
    LOGOUT("Resetting forward solver");
    for (i = 0; i < nlambda; i++) {
	FWS.Reset (*msol.swsol[i], 0);
	FWS.CalcFields (qvec, dphi+i*nq);
	RVector proj_i(proj, i*nqm, nqm);
	proj_i = FWS.ProjectAll (mvec, dphi+i*nq);

	if (aphi) // calculate adjoint fields
	    FWS.CalcFields (mvec, aphi+i*nm);
    }
    if (aphi) {
	if (hscale == LM_HSCALE_IMPLICIT)
	    aphi_hscale = aphi;
	if (precon == LM_PRECON_GMRES_DIRECT) aphi = 0;
    }
    
    int i_count = 0; // iteration counter
    
    // create the regularisation instance
    reg = Regularisation::Create (pp, &x0, &raster);
    
    // Set up the context data for the line search callback
    OF_CLBK_DATA ofdata;
    ofdata.fws = &FWS;
    ofdata.raster = &raster;
    ofdata.pscaler = pscaler;
    ofdata.meshsol = &msol;
    ofdata.reg = reg;
    ofdata.qvec = &qvec;
    ofdata.mvec = &mvec;
    ofdata.data = &data;
    ofdata.sd = &sd;

    HESS_DATA hdata;
    hdata.lambda = &lambda;
    hdata.reg = reg;
    hdata.RHess = &RHess;
    // set the components for the implicit definition of the Hessian
    
    bool logparam = !strcmp (pscaler->ScalingMethod(), "LOG");
    //JAC_DATA jdata = {mesh, &raster, &FWS, &dphi, &aphi, &mvec, &data, &sd,
    //		      &prmscl, reg, &RHess, &M, &lambda, &msol, &extcoef};
    // set the components for the implicit definition of the Jacobian
    
    of_value = ObjectiveFunction::get_value (data, proj, sd);
    of_prior = reg->GetValue (x1);
    of = of_value + of_prior;
    // actually this will be 0 since x = x0 = x1 here
    err0 = errstart = of;
    
    err00 = 0; // this is previous error
    
    //    test_biscale();
    LOGOUT("Starting error: %f", errstart);
    
    LOGOUT("Iteration 0  CPU %f  OF %g [LH %g PR %g]", toc(clock0),
	   of, of_value, of_prior);

    MWJacobian J;
    J.msol = &msol;
    J.bsol = &bsol;

    hdata.J = &J;

    // LM iteration loop
    for (inr = 0; (!nrmax || inr < nrmax) &&
	     err0 > gn_tol*errstart && fabs(err00-err0) > gn_tol; inr++) {

	RVector kap = bsol.GetParam (OT_CKAPPA);
	prmscl = pscaler->JScale (bsol.GetActiveParams());
	RVector b = (data-proj)/sd;     // negative residual
	
	// update Hessian of prior
	RHess = BuildRHessian (reg, x, nprm);
	
	// generate Jacobian for all wavelengths
	LOGOUT("Calculating Jacobian ...");
	GenerateJacobian (FWS, raster, mesh, qvec, mvec,
			  FWS.GetDataScaling(), &J);

	// apply solution space rescaling S
	J.prmscale = pscaler->JScale (bsol.GetActiveParams());
	
	// apply data space rescaling T
	J.datascale = (inv(sd));
	
	// apply Hessian normalisation M
	J.hscale = J.RescaleHessian (x, &RHess);
	
	// calculate Gradient
	r = J.ATx (b);   // gradient

	LOGOUT("Gradient norm: %f", l2norm(r));
	
	// add prior to gradient
	RVector gpsi = reg->GetGradient (x);
	gpsi *= J.hscale; // apply Hessian normalisation
	r -= gpsi;
	LOGOUT("Gradient norm with penalty: %f", l2norm(r));
	cerr << "T(reg_gradient)=" << toc() << endl;
	tic();

	bool KeepGoing = true;
	i_count = 0;
	err00 = err0;
	while (i_count < itmax && KeepGoing) {
	    LOGOUT("LM iteration %d", i_count);

	    if (Gradient_descent) {

		h = r;

	    } else {  // solve Hh = r

		switch (precon) {
		case LM_PRECON_NONE:    // assume H = I
		    h = r;
		    break;
		case LM_PRECON_HDIAG:   // assume H = diag(H)
		    LOGOUT ("Solving Hessian: LM-DIAG ...");
		    h = r / Hess_diag (hdata, x);
		    break;
		case LM_PRECON_CH: {    // solve with Cholesky factorisation
		    LOGOUT ("Solving Hessian: LM-CH ...");
		    RSymMatrix hess = Hess_full (hdata, x);// explicit Hessian
		    CHdecomp (hess);
		    h = CHsubst(hess, r);
	            } break;
		case LM_PRECON_ICH:     // Incomplete Cholesky factorisation
		    cerr << "LM ICH preconditioner not implemented" << endl;
		    exit(1);
		case LM_PRECON_PCG: {
		    LOGOUT ("Solving Hessian: LM-PCG ...");
		    double tol = gmres_tol;
		    static RPrecon_Diag precon;
		    precon.ResetFromDiagonal (J.hscale);
		    PCG (JTJx_clbk, &hdata, r, h, tol, &precon);
		    } break;
		case LM_PRECON_BICGSTAB: {
		    LOGOUT ("Solving Hessian: LM-BICGSTAB ...");
		    double tol = gmres_tol;
		    static RPrecon_Diag precon;
		    precon.ResetFromDiagonal (J.hscale);
		    BiCGSTAB (JTJx_clbk, &hdata, r, h, tol, &precon);
		    } break;
		case LM_PRECON_GMRES: {  // GMRES solver with implicit Hessian
		    LOGOUT ("Solving Hessian: LM-GMRES ...");
		    double tol = gmres_tol;
		    static RPrecon_Diag precon;
		    precon.ResetFromDiagonal (J.hscale);
		    GMRES (JTJx_clbk, &hdata, r, h, tol, &precon, 20);
		    //GMRES (JTJx_clbk, &hdata, r, h, tol);
		    //GMRES (Hess_full (hdata, x), r, h, tol); // debug only
	            } break;
		}
	    
	    }

	    d = h;
	    double stepsize = alpha0;

#ifdef NEWTON_CG
	    double beta;
	    if(inr >  0) {
	        beta = (h&r)/hrold;
		d += dold*beta;
	    }
	    dold = d;
	    hrold = h&r;
#endif

	    cerr << "T(solve_H)=" << toc() << endl;
	    tic();

	    bool status_valid = false;
	    if (do_linesearch) {
		static double alpha_ls = stepsize;
		//if (alpha < 0.0) alpha = stepsize; // initialise step length
		if (LineSearch (x, d, alpha_ls, err0, of_clbk, &alpha_ls,
				&fmin, &ofdata) != 0 ||
		    fabs(err00-fmin) <= gn_tol) {
		    LOGOUT ("No decrease in line search");
		    lambda *= lambda_scale;
		    if (lambda_scale == 1.0) {
			if (Gradient_descent) KeepGoing = false;
			else {
			    Gradient_descent = true, alpha_ls = stepsize;
			    LOGOUT ("Reverting to gradient direction");
			}
		    }
		    continue;
		}
		Gradient_descent = false;
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
	    cerr << "T(linesearch)=" << toc() << endl;
	    tic();

	    x1 = x + d*alpha;  // trial solution	

	    raster.Map_ActiveSolToMesh (pscaler->Unscale(x1), msol);
	    msol.RegisterChange();

	    if (1/*!status_valid*/) {
		for (i = 0; i < nlambda; i++) {
		    FWS.Reset (*msol.swsol[i], 0);
		    FWS.CalcFields (qvec, dphi+i*nq);
		    if (aphi)
			FWS.CalcFields (mvec, aphi+i*nm);
		    else if (aphi_hscale)
			FWS.CalcFields (mvec, aphi_hscale+i*nm);
		    RVector proj_i(proj, i*nqm, nqm);
		    proj_i = FWS.ProjectAll (mvec, dphi+i*nq);
		}
		err1 = ObjectiveFunction::get_value (data, proj, sd);
	    }

	    if(err1 < err0) {  // accept update
	        //err00 = err0;
		err0 = err1;
		x = x1;
		bsol.SetActiveParams (pscaler->Unscale(x));
		raster.Map_SolToMesh (bsol, msol);
		msol.RegisterChange();
		for (i = 0; i < nlambda; i++) {
		    FWS.Reset (*msol.swsol[i], 0);
		    FWS.CalcFields (qvec, dphi+i*nq);
		    if (aphi)
			FWS.CalcFields (mvec, aphi+i*nm);
		    else if (aphi_hscale)
		    	FWS.CalcFields (mvec, aphi_hscale+i*nm);
		    RVector proj_i(proj, i*nqm, nqm);
		    proj_i = FWS.ProjectAll (mvec, dphi+i*nq);
		}

		if (g_imgfmt != IMGFMT_RAW) {
		    for (i = 0; i < msol.nParam(); i++) {
			char fname[256];
			if (msol.IsActive(i)) {
			    if (i < msol.nmuaChromo) 
			        sprintf (fname,"%sreconChromophore_%d.nim",
				    g_prefix,i+1);
			    if (i == msol.nmuaChromo)
			        sprintf (fname,"%sreconScatPrefactor_A.nim",
				    g_prefix);
			    if (i == msol.nmuaChromo + 1) 
			        sprintf (fname,"%sreconScatPower_b.nim",
				    g_prefix);
			    msol.WriteImgGeneric (inr+1, fname, i);
			}
		    }
		}
		if (g_imgfmt != IMGFMT_NIM) {
		    Solution gsol(msol.nParam(), blen);
		    raster.Map_SolToBasis (bsol, gsol, true);
		    for (i = 0; i < msol.nParam(); i++) {
			char fname[256];
			if (msol.IsActive(i)) {
			    if (i < msol.nmuaChromo)
				sprintf (fname, "%sreconChromophore_%d.raw",
				    g_prefix,i+1);
			    else if (i == msol.nmuaChromo)
			        sprintf (fname, "%sreconScatPrefactor_A.raw",
				    g_prefix);
			    else if (i == msol.nmuaChromo+1)
			        sprintf (fname,"%sreconScatPower_b.raw",
				    g_prefix);
			    else
				xERROR("Invalid parameter index during output");
			    gsol.WriteImgGeneric (inr+1, fname, i);
			}
		    }
		}

#ifdef UNDEF
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
		}
#endif
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
	    cerr << "T(update)=" << toc() << endl;
	    tic();

	    i_count++;
	}
        of_prior = reg->GetValue (x);
	LOGOUT("Iteration %d  CPU %f  OF %g [LH %g PR %g]",
	       inr+1, toc(clock0), err0, err0-of_prior, of_prior);
    } // end of NR loop;

    if (err0 <= gn_tol*errstart)
        LOGOUT ("LM solver convergence criterion satisfied");
    else if (fabs(err00-err0) <= gn_tol)
        LOGOUT ("LM solver stopped (insufficient improvement)");
    else if (nrmax && inr == nrmax)
        LOGOUT ("LM solver iteration limit reached");
    else
        LOGOUT ("LM solver terminated");

    // final residuals
    double rd = ObjectiveFunction::get_value (data, proj, sd);
    double rp = reg->GetValue (x);
    double tau = reg->GetTau();
    LOGOUT("Residuals  TAU %g  DATA %g  PRIOR/TAU %g  IT %d", tau, rd,
		rp/tau, inr);
}

void SolverLM2_CW_MW::ReadParams (ParamParser &pp)
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
	else if (!strcasecmp (cbuf, "PCG"))
	    precon = LM_PRECON_PCG, def = true;
	else if (!strcasecmp (cbuf, "BICGSTAB"))
	    precon = LM_PRECON_BICGSTAB, def = true;
	else if (!strcasecmp (cbuf, "GMRES"))
	    precon = LM_PRECON_GMRES, def = true;
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
	cout << "(4) PCG solver (implicit Hessian)\n";
	cout << "(5) BiCGSTAB solver (implicit Hessian)\n";
	cout << "(6) GMRES solver (implicit Hessian)\n";
	cout << "(7) GMRES solver (implicit Jacobian)\n";
	cout << "(8) GMRES solver (direct method)\n";
	cout << "[0|1|2|3|4|5|6|7|8] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: precon = LM_PRECON_NONE,     def = true; break;
	case 1: precon = LM_PRECON_HDIAG,    def = true; break;
	case 2: precon = LM_PRECON_CH,       def = true; break;
	case 3: precon = LM_PRECON_ICH,      def = true; break;
	case 4: precon = LM_PRECON_PCG,      def = true; break;
	case 5: precon = LM_PRECON_BICGSTAB, def = true; break;
	case 6: precon = LM_PRECON_GMRES,    def = true; break;
	case 7: precon = LM_PRECON_GMRES_JACOBIANFREE, def = true; break;
	case 8: precon = LM_PRECON_GMRES_DIRECT, def = true; break;
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

    // 4. === PRECON CONVERGENCE CRITERION ===

    if (precon == LM_PRECON_PCG || precon == LM_PRECON_BICGSTAB ||
        precon == LM_PRECON_GMRES || precon == LM_PRECON_GMRES_JACOBIANFREE ||
	precon == LM_PRECON_GMRES_DIRECT) {
      if (!(pp.GetReal ("LM_PRECON_TOL", gmres_tol) ||
	    pp.GetReal ("LM_GMRES_TOL", gmres_tol)) ||
	  gmres_tol <= 0.0) do {
	    cout << "\nSelect LM GMRES convergence criterion (>0):\n";
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
}

void SolverLM2_CW_MW::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "LM2");
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
    case LM_PRECON_PCG:
        pp.PutString ("LM_PRECON", "PCG");
	pp.PutReal ("LM_PRECON_TOL", gmres_tol);
	break;
    case LM_PRECON_BICGSTAB:
        pp.PutString ("LM_PRECON", "BICGSTAB");
	pp.PutReal ("LM_PRECON_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES:
	pp.PutString ("LM_PRECON", "GMRES");
	pp.PutReal ("LM_PRECON_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES_JACOBIANFREE:
	pp.PutString ("LM_PRECON", "GMRES_JFREE");
	pp.PutReal ("LM_PRECON_TOL", gmres_tol);
	break;
    case LM_PRECON_GMRES_DIRECT:
	pp.PutString ("LM_PRECON", "GMRES_DIRECT");
	pp.PutReal ("LM_PRECON_TOL", gmres_tol);
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
}

// ==========================================================================
// Build the Hessian of the prior from individual parameter contributions

static RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
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

static RVector RescaleHessian (const  RMatrix *J, const RVector &x,
    RCompRowMatrix *RHess)
{
    int i, j, m = J->nRows(), n = J->nCols();
    RVector M(n);

    // J^T J contribution
    for (j = 0; j < m; j++) {
	RVector r = J->Row (j);
	for (i = 0; i < n; i++) M[i] += r[i]*r[i];
    }
    // Psi'' contribution
    if (RHess) M += RHess->Diag();

    for (i = 0; i < n; i++) {
	M[i] = 1.0/sqrt (M[i]);
    }

    return M;
}

#ifdef TOAST_MPI
static RVector RescaleHessianPart (RDenseMatrix &Jpart, const RVector &x,
    RCompRowMatrix *RHess, int rank)
{
    int i, j, k, nz;
    int m = Jpart.nRows(), n = Jpart.nCols();
    RVector M(n), Mpart(n);

    // J^T J contribution
    for (j = 0; j < m; j++) {
	RVector r = Jpart.Row(j);
	for (i = 0; i < n; i++) Mpart[i] += r[i]*r[i];
    }
    MPI_Reduce ((void*)Mpart.data_buffer(), (void*)M.data_buffer(),
	n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // add contributions

    if (!rank) {
	// Psi'' contribution
	if (RHess) M += RHess->Diag();

	for (i = 0; i < n; i++) 
	    M[i] = 1.0/sqrt (M[i]);
    }
    MPI_Bcast ((void*)M.data_buffer(), n, MPI_DOUBLE, 0,
	MPI_COMM_WORLD);
    return M;
    
}
#endif

// =================================================================
// select a part of a vector from all wavelength parts and merge 
static RVector Mergewavels (const RVector &b, int r0, int len, int nofwavel)
{
    int i, j, wlofs, k = 0;
    int swlen = b.Dim()/nofwavel;

    RVector bparts (len*nofwavel);
     
    for (i = 0; i < nofwavel; i++) {
	wlofs = i * swlen;
	for (j = 0; j < len; j++) {
	    bparts[k] = b[wlofs+r0+j];
	    k++;
	}
    }
    return bparts;
}

// ==========================================================================
// This function builds the full multiwavelength Jacobian, containing
// all chromophores

void GenerateJacobian (RFwdSolverMW &FWS, const Raster &raster,
    const QMMesh *mesh,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    DataScale dscale, MWJacobian *J)
{
    int i, j;
    int dim     = raster.Dim();
    int nlen    = mesh->nlen();
    int slen    = raster.SLen();
    int nq      = mesh->nQ;
    int nm      = mesh->nM;
    int nqm     = mesh->nQM;
    const MWsolution *msol = J->msol;
    const Solution *bsol = J->bsol;
    const RDenseMatrix *excoef = &msol->extcoef;
    int nchromo = excoef->nCols();
    int nlambda = excoef->nRows();
    int nprm    = bsol->nParam();
    bool bFactorA = bsol->IsActive(nchromo);
    bool bPowerb  = bsol->IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;

    // create field buffers
    RVector *dphi = new RVector[nq];
    for (i = 0; i < nq; i++) dphi[i].New(nlen);
    RVector *aphi = new RVector[nm];
    for (i = 0; i < nm; i++) aphi[i].New(nlen);

    J->Jmua[i].New(nqm,slen);
    J->Jkap[i].New(nqm,slen);

    // loop over wavelengths
    for (i = 0; i < nlambda; i++) {

	// Generic mua/kappa Jacobian
	FWS.Reset (*msol->swsol[i]);
	FWS.CalcFields (qvec, dphi);
	FWS.CalcFields (mvec, aphi);
	GenerateJacobian_cw (&raster, mesh, mvec, dphi, aphi, dscale,
			     &J->Jmua[i], &J->Jkap[i]);

	if (bFactorA)
	    raster.Map_MeshToSol (msol->GetJacobianCoeff_A(i), J->scale_A[i]);
	if (bPowerb)
  	    raster.Map_MeshToSol (msol->GetJacobianCoeff_b(i), J->scale_b[i]);
    }

    // cleanup
    delete []dphi;
    delete []aphi;
}


RVector MWJacobian::Ax (const RVector &x) const
{
    int i, j, ii, jj, xidx, yidx;
    double e, Jij, prms;

    int nchromo = msol->extcoef.nCols();
    int nlambda = msol->extcoef.nRows();
    int slen    = Jmua[0].nCols();
    int ndat    = Jmua[0].nRows();
    int nqm     = ndat;

    bool bA = bsol->IsActive (nchromo);
    bool bb = bsol->IsActive (nchromo+1);
    int n   = bsol->ActiveDim();
    int m   = ndat*nlambda;

    RVector y(m);

    // loop over chromophores
    for (j = 0; j < nchromo; j++) {
        for (i = 0; i < nlambda; i++) {
	    e = msol->extcoef(i,j);
	    const RDenseMatrix &J = Jmua[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = j*slen+jj;
		prms = e*prmscale[xidx]*hscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    y[yidx] += Jij*x[xidx];
		}
	    }
	}
    }
    // scattering parameters
    if (bA) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = nchromo*slen+jj;
		prms = scale_A[i][jj]*prmscale[xidx]*hscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    y[yidx] += Jij*x[xidx];
		}
	    }
	}
    }
    if (bb) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = (nchromo+1)*slen+jj;
		prms = scale_b[i][jj]*prmscale[xidx]*hscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    y[yidx] += Jij*x[xidx];
		}
	    }
	}
    }
    return y;
}

RVector MWJacobian::ATx (const RVector &y) const
{
    int i, j, ii, jj, xidx, yidx;
    double e, Jij, prms;

    int nchromo = msol->extcoef.nCols();
    int nlambda = msol->extcoef.nRows();
    int slen    = Jmua[0].nCols();
    int ndat    = Jmua[0].nRows();
    int nqm     = ndat;

    bool bA = bsol->IsActive (nchromo);
    bool bb = bsol->IsActive (nchromo+1);
    int n   = bsol->ActiveDim();

    RVector x(n);

    // loop over chromophores
    for (j = 0; j < nchromo; j++) {
        for (i = 0; i < nlambda; i++) {
	    e = msol->extcoef(i,j);
	    const RDenseMatrix &J = Jmua[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = j*slen+jj;
		prms = e*prmscale[xidx]*hscale[xidx];
	        for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    x[xidx] += Jij*y[yidx];   
	        }
	    }
        }
    }
    // scattering parameters
    if (bA) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = nchromo*slen+jj;
		prms = scale_A[i][jj]*prmscale[xidx]*hscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    x[xidx] += Jij*y[yidx];
		}
	    }
	}
    }
    if (bb) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = (nchromo+1)*slen+jj;
		prms = scale_b[i][jj]*prmscale[xidx]*hscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    Jij = J(ii,jj)*prms*datascale[yidx];
		    x[xidx] += Jij*y[yidx];
		}
	    }
	}
    }
    return x;
}


RVector MWJacobian::RescaleHessian (const RVector &x, RCompRowMatrix *RHess)
    const
{
    int i, j, ii, jj, xidx, yidx;
    double e, prms, v;

    int nchromo = msol->extcoef.nCols();
    int nlambda = msol->extcoef.nRows();
    int slen    = Jmua[0].nCols();
    int ndat    = Jmua[0].nRows();
    int nqm     = ndat;
    int n   = bsol->ActiveDim();
    bool bA = bsol->IsActive (nchromo);
    bool bb = bsol->IsActive (nchromo+1);

    RVector M(n);

    // loop over chromophores
    for (j = 0; j < nchromo; j++) {
        for (i = 0; i < nlambda; i++) {
	    e = msol->extcoef(i,j);
	    const RDenseMatrix &J = Jmua[i];
	    for (jj = 0; jj < slen; jj++) {
		xidx = j*slen+jj;
	        prms = e * prmscale[xidx];
	        for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    v = J(ii,jj)*prms*datascale[yidx];
		    M[xidx] += v*v;
		}
	    }
	}
    }
    // scattering parameters
    if (bA) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = nchromo*slen+jj;
		prms = scale_A[i][jj]*prmscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    v = J(ii,jj)*prms*datascale[yidx];
		    M[xidx] += v*v;
		}
	    }
	}
    }
    if (bb) {
        for (i = 0; i < nlambda; i++) {
	    const RDenseMatrix &J = Jkap[i];
	    for (jj = 0; jj < slen; jj++) {
	        xidx = nchromo*slen+jj;
		prms = scale_b[i][jj]*prmscale[xidx];
		for (ii = 0; ii < nqm; ii++) {
		    // log amplitude part
		    yidx = i*nqm+ii;
		    v = J(ii,jj)*prms*datascale[yidx];
		    M[xidx] += v*v;
		}
	    }
	}
    }
    if (RHess) M += RHess->Diag();

    for (i = 0; i < n; i++) {
	M[i] = 1.0/sqrt (M[i]);
    }

    return M;
}
