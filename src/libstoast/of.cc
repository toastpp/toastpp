// NOTE: Most of this is obsolete
// Has been replaced by class Regularisation
// This file is in need of re-work.

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "util.h"
#include "timing.h"
#include <fstream>

//#define OUTPUT_FIELDS

using namespace std;

// fix these extern references
template<class T> extern void ImageGradient (const IVector &dim,
    const RVector &size, const TVector<T> &im, TVector<T> *grad,
    const int *mask);

extern char g_meshname[256];
extern int g_nimsize;

// =========================================================================
// class ObjectiveFunction

ObjectiveFunction::ObjectiveFunction (const Raster *_raster)
{
    data = 0;
    sd = 0;
    scale = 0;
    apply_prior = false;
    hyperp[0] = hyperp[1] = 1.0;
    nbg_rowptr = 0;
    nbg_colidx = 0;
    raster = _raster;
}

void ObjectiveFunction::SetTVParam (double _tv_beta2,
    double _hypermua, double _hyperkappa)
{
    tv_beta2 = _tv_beta2;
    hyperp[0] = _hypermua;
    hyperp[1] = _hyperkappa;
    apply_prior = true;
}

double ObjectiveFunction::get_posterior (const RVector *proj) const
{
    
    RVector fy = (*data - *proj) / *sd;
    if (scale) fy = *scale * fy;
    return dot (fy, fy);

#ifdef UNDEF
    int i, dim = data->Dim();
    double temp, sum = 0.0;
    for (i = 0; i < dim; i++) {
        temp = ((*data)[i] - (*proj)[i]) / (*sd)[i];
	sum += temp*temp;
    }
    return sum;
#endif
}

double ObjectiveFunction::get_prior (PRIOR_OLD prior, const RVector &x,
    const RVector &xp, const RVector &x0, double tau,
    const RCompRowMatrix &LTL)
{
    switch (prior) {
    case PRIOR_LAPLACIAN_OLD:
	return penalty_tikh1 (x, x0, tau, tau, LTL);
    case PRIOR_DIAG_OLD:
	return tau * penalty_tikh0 (x, x0, xp);
    default:
	return 0.0;
    }
}

#ifdef UNDEF
double ObjectiveFunction::get_prior (const Solution *sol) const
{
    double prior = 0.0;
    if (!sol || !apply_prior) return 0.0;

#if REGULARISATION == REGULARISATION_GMRF  // Markov random field
    int i, j, k, c, r1, r2, dim;
    double h, v1, v2, edge;
    const double *hp = hyper;

    for (i = 0; i < 3; i++) {
        Solution::PARAM prm = (Solution::PARAM)i;
	if (sol->IsActive(prm)) {
	    RVector val(sol->GetParam (prm));
	    h = *hp++;
	    dim = val.Dim();
	    for (j = 0; j < dim; j++) {
	        v1 = val[j];
		r1 = nbg_rowptr[j];
		r2 = nbg_rowptr[j+1];
		for (k = r1; k < r2; k++) {
		    c = nbg_colidx[k];
		    if (c == j) continue; 
		    v2 = val[c];
		    edge = v1-v2;
		    prior += h * pow (fabs(edge), gmrf_exp);
		}
	    }
	}
    }
#elif REGULARISATION == REGULARISATION_TV  // total variation

    RVector bval (raster->BLen());
    IVector bdim = raster->BDim();
    RVector gsize = raster->GSize();
    const int *elref = raster->Elref();
    int i, j, k, idx, dim = bdim.Dim();
    double dx, dy, dz, Dxf, Dyf, Dzf, psi, sum;
    double *hp = hyper;
    dx = gsize[0]/(bdim[0]-1);
    dy = gsize[1]/(bdim[1]-1);
    if (dim > 2) dz = gsize[2]/(bdim[2]-1);

    for (i = 0; i < 3; i++) {
        Solution::PARAM prm = (Solution::PARAM)i;
	if (!sol->IsActive (prm)) continue;
	RVector val(sol->GetParam (prm));
	raster->Map_SolToBasis (val, bval);
	sum = 0.0;

	if (dim > 2) {
	    for (k = 1; k < bdim[2]; k++) {
	        for (j = 1; j < bdim[1]; j++) {
		    for (i = 0; i < bdim[0]; i++) {
		        idx = k*bdim[0]*bdim[1] + j*bdim[0] + i;
			if (elref[idx] < 0) continue;
			if (elref[idx-1] >= 0)
			    Dxf = (bval[idx]-bval[idx-1])/dx;
			else
			    Dxf = 0.0;
			if (elref[idx-bdim[0]] >= 0)
			    Dyf = (bval[idx]-bval[idx-bdim[0]])/dy;
			else
			    Dyf = 0.0;
			if (elref[idx-bdim[0]*bdim[1]] >= 0)
			    Dzf = (bval[idx]-bval[idx-bdim[0]*bdim[1]])/dz;
			else
			    Dzf = 0.0;
			psi = 2.0*sqrt(Dxf*Dxf + Dyf*Dyf + Dzf*Dzf + tv_beta2);
			sum += psi * dx*dy*dz;
		    }
		}
	    }
	} else {
	    // not implemented
	}
	prior += sum * *hp++;
    }

#ifdef UNDEF
    RVector bval(raster->BLen());
    RVector bgrad[3];
    int i, j, k, dim = raster->BDim().Dim();
    int blen = raster->BLen();
    double term, sum, g, *hp = hyper;

    for (i = 0; i < dim; i++) bgrad[i].New(blen);

    for (i = 0; i < 3; i++) {
        Solution::PARAM prm = (Solution::PARAM)i;
	if (sol->IsActive (prm)) {
	    RVector val(sol->GetParam (prm));
	    raster->Map_SolToBasis (val, bval);
	    ImageGradient (raster->BDim(), raster->GSize(), bval, bgrad,
			   raster->BasisToSol());
	    for (sum = 0.0, j = 0; j < blen; j++) {
	        term = tv_beta2;
		for (k = 0; k < dim; k++) {
		    g = bgrad[k][j];
		    term += g*g;
		}
		sum += sqrt(term);
	    }
	    prior += sum * *hp++;
	}
    }

#endif


#endif

    return prior;
}
#endif // UNDEF

double ObjectiveFunction::get_value (const RVector &data, const RVector &proj,
    const RVector &sd, const RMatrix *scale)
{
    RVector fy = (data - proj) / sd;
    if (scale) fy = *scale * fy;
    return dot (fy, fy);

#ifdef UNDEF
    int i, dim = data.Dim();
    double temp, sum = 0.0;
    for (i = 0; i < dim; i++) {
        temp = (data[i] - proj[i]) / sd[i];
	sum += temp*temp;
    }
    return sum;
#endif
}

double ObjectiveFunction::get_value (const RVector *proj, const Solution *sol)
    const
{
    double of;
    of = get_posterior (proj);
    //if (apply_prior) of += get_prior (sol);
    return of;
}

void ObjectiveFunction::add_gradient_data (RVector &grad, const Raster &raster,
    const CFwdSolver &FWS, const RVector &proj, CVector *dphi,
    const CCompRowMatrix &mvec) const
{
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs_mod, ofs_arg;
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
    RVector grad_cmua(grad, 0, slen);       // mua part of grad
    RVector grad_ckappa (grad, slen, slen); // kappa part of grad
    
    double tm_mesh2grid = 0.0;
    double tm_grid2sol = 0.0;
    double tm_gradient = 0.0;
    double tm_innerloop = 0.0;

    LOGOUT1_INIT_PROGRESSBAR ("AddGradient", 50, mesh.nQ);
    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	tic();
	raster.Map_MeshToGrid (dphi[q], cdfield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);
	tm_gradient += toc();

        n = mesh.nQMref[q];

	RVector y_mod (*data, ofs_mod, n);
	RVector s_mod (*sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector y_arg (*data, ofs_arg, n);
	RVector s_arg (*sd, ofs_arg, n);
	RVector ypm_arg (proj, ofs_arg, n);
	RVector b_arg(n);
	b_arg = (y_arg-ypm_arg)/s_arg;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	//Project_cplx (mesh, q, dphi[q], cproj);
	cproj = ProjectSingle (&mesh, q, mvec, dphi[q], DATA_LIN);
	wqa = std::complex<double>(0,0);
	wqb = 0.0;

	tic();
	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * std::complex<double> (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa += qs * std::complex<double> (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}
	tm_innerloop += toc();

	// adjoint field and gradient
	CVector wphia (mesh.nlen());

	FWS.CalcField (wqa, wphia);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	tic();
	raster.Map_MeshToGrid (wphia, cafield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);
	tm_gradient += toc();

	// absorption contribution
	tic();
	raster.Map_GridToSol (cdfield * cafield, dgrad);

	tm_grid2sol += toc();
	grad_cmua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	tic();
	raster.Map_GridToSol (gk, dgrad);

	tm_grid2sol += toc();
	grad_ckappa -= Re(dgrad);

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
	LOGOUT1_PROGRESS (q);
    }
}

void ObjectiveFunction::add_gradient_prior (const Solution &sol,
    RVector &grad) const
{
#if REGULARISATION == REGULARISATION_GMRF
    int prm, j, k, c, r1, r2, dim, ofs = 0;
    double h, v1, v2, edge, tmp;
    const double *hp = hyperp;

    for (prm = 0; prm < sol.nParam(); prm++) {
        if (sol.IsActive(prm)) {
 	    RVector val(sol.GetParam (prm));
	    h = *hp++;
	    dim = val.Dim();
	    for (j = 0; j < dim; j++) {
	        v1 = val[j];
	        r1 = nbg_rowptr[j];
		r2 = nbg_rowptr[j+1];
		for (k = r1; k < r2; k++) {
		    c = nbg_colidx[k];
		    if (c == j) continue; 
		    v2 = val[c];
		    edge = v1-v2;
		    tmp = pow (fabs(edge), gmrf_exp-1.0);
		    grad[j+ofs] += h*gmrf_exp * (edge >= 0.0 ? tmp:-tmp);
		}
	    }
	    ofs += dim;
	}
    }
#elif REGULARISATION == REGULARISATION_TV
    IVector bdim = raster->BDim();
    int i, j, k, m, mm, p, idx, ofs = 0, blen = raster->BLen();
    int slen = raster->SLen();
    int dim = bdim.Dim();
    int blenm = 1;
    for (i = 0; i < dim; i++) blenm *= bdim[i]-1;
    RVector bval(blen);
    RVector gsize = raster->GSize();
    RVector x(blenm*dim);
    RVector g(blen), gs(slen);
    RCompRowMatrix Dop(blenm*dim,blen);
    RCompRowMatrix psi(blenm*dim,blenm*dim);
    double sum, term, dx, dy, dz;
    dx = gsize[0]/(bdim[0]-1);
    dy = gsize[1]/(bdim[1]-1);
    if (dim > 2) dz = gsize[2]/(bdim[2]-1);

    // build gradient operator matrix (should be done externally)
    int *rowptr1, *colidx1, *rowptr2, *colidx2;
    double *val1, *val2;
    rowptr1 = new int[blenm*dim+1];
    colidx1 = new int[blenm*dim*2];
    val1    = new double[blenm*dim*2];
    rowptr2 = new int[blenm*dim+1];
    colidx2 = new int[blenm*dim];
    val2    = new double[blenm*dim];

    for (p = 0; p < 3; p++) {
        Solution::PARAM prm = (Solution::PARAM)p;
	if (!sol.IsActive(prm)) continue;
	RVector sval(sol.GetParam (prm));
	raster->Map_SolToBasis (sval, bval);

	if (dim > 2) {
	  
	    for (m = mm = 0, k = 1; k < bdim[2]; k++) {
	        for (j = 1; j < bdim[1]; j++) {
		    for (i = 1; i < bdim[0]; i++) {
		        sum = 0.0;
			idx = k*bdim[0]*bdim[1]+ j*bdim[0] + i;
			// Dx part
			colidx1[m]   = idx-1;
			val1[m]      = -dx;
			colidx1[m+1] = idx;
			val1[m+1]    = dx;
			term = (bval[idx]-bval[idx-1])/dx;
			sum += term*term;
			// Dy part
			colidx1[m+blenm*2]   = idx-bdim[0];
			val1[m+blenm*2]      = -dy;
			colidx1[m+blenm*2+1] = idx;
			val1[m+blenm*2+1]    = dy;
			term = (bval[idx]-bval[idx-bdim[0]])/dy;
			sum += term*term;
			// Dz part
			colidx1[m+blenm*4]   = idx-bdim[0]*bdim[1];
			val1[m+blenm*4]      = -dz;
			colidx1[m+blenm*4+1] = idx;
			val1[m+blenm*4+1]    = dz;
			term = (bval[idx]-bval[idx-bdim[0]*bdim[1]])/dz;
			sum += term*term;
			
			colidx2[mm] = mm;
			val2[mm] = 1.0/sqrt(sum+tv_beta2);
			mm++;
			m += 2;
		    }
		}
	    }

	} else { // dim <= 2

	    // not implemented yet
 
	}

	Dop.Initialise (rowptr1, colidx1, val1);
	psi.Initialise (rowptr2, colidx2, val2);
	Dop.Ax (bval, x);
	psi.Ax (x, x);
	Dop.ATx (x, g);
	raster->Map_BasisToSol (g, gs);
	for (i = 0; i < gs.Dim(); i++)
	    grad[ofs+i] += gs[i];
	ofs += gs.Dim();
    }
    delete []rowptr1;
    delete []colidx1;
    delete []val1;
    delete []rowptr2;
    delete []colidx2;
    delete []val2;

#ifdef UNDEF
    RVector bval(raster->BLen()), gp(raster->BLen());
    RVector bgrad[3];
    int i, j, k;
    double term, num, g, h, *hp = hyper;

    for (i = 0; i < dim; i++) bgrad[i].New (blen);

    for (i = 0; i < 3; i++) {
        Solution::PARAM prm = (Solution::PARAM)i;
	if (sol.IsActive(prm)) {
	    RVector val(sol.GetParam (prm));
	    raster->Map_SolToBasis (val, bval);
	    



	    ImageGradient (raster->BDim(), raster->GSize(), bval, bgrad,
			   raster->BasisToSol());
	    for (d = 0; d < dim; d++) {
	        for (j = 0; j < blen; j++) {
		    for (t = 0, k = 0; k < dim; k++)
		       t += bgrad[j][k] * bgrad[j][k];
		    psid = 1.0/sqrt(t + tv_beta2);
		    gp[j] = psid*bgrad[d];
		}
		
	    }

	    h = *hp++;
	    for (j = 0; j < blen; j++) {
	        term = tv_beta2;
		num  = 0.0;
		for (k = 0; k < dim; k++) {
		    g = bgrad[k][j];
		    term += g*g;
		    num += g;
		}
		gp[j] += h * num / sqrt(term);
	    }
	}
    }
    raster->Map_BasisToSol (gp, grad);
#endif

#endif
}

void ObjectiveFunction::get_gradient (const Raster &raster,
    const CFwdSolver &FWS, const RVector &proj, CVector *dphi,
    const CCompRowMatrix &mvec, const Solution *sol, RVector &grad) const
{
    grad.Clear();
    add_gradient_data (grad, raster, FWS, proj, dphi, mvec);
    if (apply_prior && sol)
	add_gradient_prior (*sol, grad);
}

// =========================================================================
// methods for priors

double penalty_tikh0 (const RVector &x, const RVector &x0, const RVector &xs)
{
  /* simple tikhonov norm of solution */
    int i, dim = x.Dim();
    double temp, sum = 0.0;
    for (i = 0; i < dim; i++) {
        temp = (x[i] - x0[i]) / xs[i];
	sum += temp*temp;
    }
    return sum;
}

void penalty_gradient_tikh0_add (const RVector &x, const RVector &x0,
    const RVector &xs, RVector &grad, const double tau)

{
  /* simple tikhonov norm of solution */
    int i, dim = x.Dim();
    for (i = 0; i < dim; i++) {
	 grad[i] -= tau*(x[i] - x0[i]) / (xs[i]*xs[i]);
    }
}
void penalty_gradient_tikh0_rescale_add (const RVector &x, const RVector &x0,
    RVector &grad, const double tau)

{
  /* simple tikhonov norm of solution */
    int i, dim = x.Dim();
    for (i = 0; i < dim; i++) {
      grad[i] -= tau*(x[i] - x0[i]);
    }
}
double penalty_tikh1 (const RVector &x, const RVector &x0, 
 const double tau1, const double tau2, const RCompRowMatrix& LTL)
{
  /* first tikhonov norm of solution */
    int dim = x.Dim();
    RVector diff(x-x0), x1,x2, xL(dim), x1L, x2L;
    x1.Relink (diff, 0, dim/2);
    x2.Relink (diff, dim/2, dim/2);
    x1L.Relink (xL, 0, dim/2);
    x2L.Relink (xL, dim/2, dim/2);

    x1L = LTL * x1; x1L *= tau1;
    x2L = LTL * x2; x2L *= tau2;
    return xL&diff;
}
void penalty_gradient_tikh1_add (const RVector &x, const RVector &x0, 
RVector &grad, const double tau1, const double tau2, const RCompRowMatrix& LTL)
{
  /* first tikhonov norm of solution */
    int dim = x.Dim();
    RVector diff(x-x0), x1,x2, xL(dim), x1L, x2L;
    x1.Relink (diff, 0, dim/2);
    x2.Relink (diff, dim/2, dim/2);
    x1L.Relink (xL, 0, dim/2);
    x2L.Relink (xL, dim/2, dim/2);

    x1L = LTL * x1; x1L *= tau1;
    x2L = LTL * x2; x2L *= tau2;
    grad -= xL;
}


// =========================================================================
// MS: variation of penalty-tikh1:
// use mean value of current image as reference base line

double penalty_tikh1_avg (const RVector &x, const double tau,
    const RCompRowMatrix &LTL)
{
    RVector xm = x-mean(x);
    return (xm & (LTL * xm)) * tau;
}

void penalty_gradient_tikh1_avg_add (const RVector &x, RVector &grad,
    const double tau, const RCompRowMatrix &LTL)
{
    RVector xm = x-mean(x);
    grad -= (LTL * xm) * tau;
}
