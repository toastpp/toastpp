// ==========================================================================
// SolverLM: Levenberg-Marquardt
// ==========================================================================

#include "util.h"
#include "raster.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "regul.h"
#include "solverlm.h"
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <sys/param.h>

// ==========================================================================
// external references

extern int bWriteJ;
extern int g_imgfmt;
extern int bOutputUpdate;
extern int bOutputGradient;
extern clock_t clock0;
extern char g_meshname[256];
extern int g_nimsize;
extern double g_refind;

void GenerateJacobian (const Raster &raster,
    const CVector *dphi, const CVector *aphi, Measurement meas,
    RMatrix *J, const RVector &sd, const RVector &kap, int fill,
    RVector *colsum = 0);

extern void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh);

// ==========================================================================
// local prototypes

RVector RescaleHessian (const RMatrix *J, const RCompRowMatrix &Lap,
    double tau);

bool LineSearchWithPrior (FwdSolver &FWS, const Raster &raster,
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
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const RMatrix *J;                // Jacobian
    const RVector *M;                // normalisation diagonal matrix
    const double *lambda;            // diagonal scaling factor
    const Regularisation *reg;       // regularisation
};

// ==========================================================================

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

// ==========================================================================

static RVector Hess_diag (const HESS_DATA &data, const RVector &x)
{
    // Computes the diagonal of the Hessian given implicitly as
    // (J^T J + P + lambda D)
    // where J and P are full, D is the diagonal of J^T J + P, and
    // lambda is a scalar

    // unpack the data
    const RMatrix *J           = data.J;
    const Regularisation *reg  =  data.reg;
    const RVector &M           = *data.M;
    const double lambda        = *data.lambda;
    int m = J->nRows(), n = J->nCols();

    RVector diag = ATA_diag (*J);

    // add prior to Hessian
    diag += reg->GetHessianDiag (x);

    diag += lambda;
    return diag;
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
    HESS_DATA *data           = (HESS_DATA*)context;
    const RMatrix *J          =  data->J;
    const Regularisation *reg =  data->reg;
    const RVector &M          = *data->M;
    const double lambda       = *data->lambda;
    int m = J->nRows(), n = J->nCols();

    RVector Px(n);

    // add prior to Hessian
    int i, j, k, nz, *colidx = new int[n];
    double *val = new double[n];
    for (i = 0; i < n; i++) {
	nz = reg->GetHessianRow (x, i, colidx, val);
	for (k = 0; k < nz; k++) {
	    j = colidx[k];
	    Px[i] += val[k] * x[j] * M[i]*M[j];
	}
    }
    delete []colidx;
    delete []val;

    return ATx (*J, Ax (*J, x)) + Px + lambda*x;
}

// ==========================================================================
SolverLM::SolverLM (): Solver()
{
    itmax = 1000;
    gmres_tol = 1e-8;
    lambda0 = 10.0;
    lambda_scale = 4.0;
    do_linesearch = true;
    alpha_min = 0.0;
    Jfill = 1.0;
    modelerr = false;
}

void SolverLM::Solve (FwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega,
    double ftol)
{
    bool Use_precon = true;

    // initialisations
    struct tms tm;
    int i,j, inr, nrmax = 50;
    const QMMesh &mesh = raster.mesh();
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh.nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int n    = slen*2;
    int Jn;
    double of, fmin, errstart, err0, err00, err1;
    double lambda = lambda0, alpha = -1.0;
    double hrold;

    RVector r(n), s(n), d(n);
    RVector h(n), dold(n), M(n);

    RVector proj(ndat);
    RVector colsum(n);
    RMatrix *J;
    RMatrix *cov;

    if (Jfill == 1.0) {
	J = new RDenseMatrix(ndat,n);
    } else {
	Jn = (int)(slen*Jfill);
	J = new RCompRowMatrix(ndat,n);
    }
    //RDenseMatrix J(ndat,n);
    RCompRowMatrix Graph(slen,slen);
    Regularisation *reg;

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    // initial solution (scaled)

    RVector x(x0), x1(x0);
    // current solution, trial solution (scaled);

    // allocate field vectors
    CVector *dphi = new CVector[mesh.nQ];
    for (i = 0; i < mesh.nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh.nM];
    for (i = 0; i < mesh.nM; i++) aphi[i].New (nlen);

    int i_count = 0; // iteration counter
    LOGOUT("Resetting forward solver");
    FWS.Reset (msol, omega);
    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
    FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
    ProjectAll (mesh, FWS, mvec, dphi, proj);

    LOGOUT("Building neighbour graph");
    int *rowptr, *colidx, nzero, row, col, i0, i1, k;

//#define OLD_GRAPH
#ifdef OLD_GRAPH
    raster.NeighbourGraph (rowptr, colidx, nzero);

    double v, valptr[nzero];
    for (row = 0; row < slen; row++) {
        i0 = rowptr[row]; i1 = rowptr[row+1];
	for (i = i0; i < i1; i++) {
	    col = colidx[i];
	    valptr[i] = (col == row ? i1-i0-1 : -1); // Discrete Laplacian
	}
    }
    Graph.Initialise (rowptr, colidx, valptr);
    Graph.PrintFillinGraph ("fillin.pgm", 600, true, false);
    delete []rowptr;
    delete []colidx;
#else
    RVector kappa(blen);
    kappa = 1.0;
    Regularisation::CreateHessian (&raster, kappa, Graph);
#endif

    switch (prior) {
    case PRIOR_DIAG:
	reg = new Tikhonov0 (tau, &x0, &x0);
	break;
    case PRIOR_LAPLACIAN:
	reg = new Tikhonov1 (tau, &x0, &Graph);
	break;
    default:
	reg = new NullRegularisation ();
	break;
    }

    // add model error means to data
    if (modelerr) {
	(RVector)data -= merr;
	cov = &mcov;
    } else {
	cov = 0;
    }

    HESS_DATA hdata = {J, &M, &lambda, reg};
    // set the components for the implicit definition of the Hessian

    of  = ObjectiveFunction::get_value (data, proj, sd, cov);
    of += ObjectiveFunction::get_prior (prior, x1, x, x0, tau, Graph);
    // actually this will be 0 since x = x0 = x1 here
    err0 = errstart = of;

    err00 = 0; // this is previous error

    //    test_biscale();
    LOGOUT_1PRM("Starting error: %f", errstart);

    times (&tm);
    LOGOUT_2PRM("Iteration 0  CPU %f  OF %f",
        (double)(tm.tms_utime-clock0)/(double)HZ, err0);

    // LM iteration loop
    for (inr = 0; inr < nrmax &&
	   err0 > ftol*errstart && fabs(err00-err0) > ftol; inr++) {

	RVector kap = bsol.GetParam (OT_CKAPPA);

	// Generate Jacobian
	LOGOUT("Calculating Jacobian ...");
	GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap, Jn,
			  &colsum);

	//ofstream ofs ("jac.dat");
	//ofs << *J << endl;
	//exit (0);

	// apply model error correction
	if (modelerr) {
	    *J = mul (mcov, *J);
	}

	// apply solution space rescaling S
	pscaler->ScaleJacobian (bsol.GetActiveParams(), *J);

	// apply Hessian normalisation M
	M = RescaleHessian (J, Graph, tau);
	J->ColScale (M);

	if (bWriteJ) WriteJacobian (J, raster, mesh);

	// calculate Gradient
	RVector b = (data-proj)/sd;             // negative residual
	J->ATx (b,r);                            // gradient
	if (bOutputGradient) {
	    RVector r1(r,0,slen);
	    RVector r2(r,slen,slen);
	    RVector rm(raster.mesh().nlen());
	    raster.Map_SolToMesh (r1, rm);
	    WriteImage (rm, inr, "gradient_mua.nim");
	    raster.Map_SolToMesh (r2, rm);
	    WriteImage (rm, inr, "gradient_mus.nim");
	}

//#define TEMP
#ifdef TEMP
	char cbuf[256];
	RVector Jr(Ax(J,r));
	RVector JTJr(ATx(J,Jr));
	RVector JTJr_mua (JTJr, 0, slen);
	RVector JTJr_kappa (JTJr, slen, slen);
	RVector outv (glen);
	RVector bb(data-proj);
	RVector b1(bb, 0, ndat/2);
	RVector b2(bb, ndat/2, ndat/2);
	for (int pw = 0; pw < 10; pw++) {
	    raster.Map_SolToGrid (JTJr_mua, outv);
	    sprintf (cbuf, "images/JTJn_mua_%03d.ppm", pw+1);
	    WritePPM (outv, raster.GDim(), 0, 0, cbuf);
	    sprintf (cbuf, "images/JTJn_kappa_%03d.ppm", pw+1);
	    raster.Map_SolToGrid (JTJr_kappa, outv);
	    WritePPM (outv, raster.GDim(), 0, 0, cbuf);
	    Jr = Ax(J,JTJr);
	    JTJr = ATx(J,Jr);

	    sprintf (cbuf, "images/JJTnb_mod_%03d.ppm", pw);
	    WriteData_pixmap (b1, raster.mesh(), false, false, cbuf);
	    sprintf (cbuf, "images/JJTnb_arg_%03d.ppm", pw);
	    WriteData_pixmap (b2, raster.mesh(), false, false, cbuf);
	    RVector JTJJT = ATx(J,bb);
	    RVector JTJJT_mua (JTJJT, 0, slen);
	    RVector JTJJT_kappa (JTJJT, slen, slen);
	    raster.Map_SolToGrid (JTJJT_mua, outv);
	    sprintf (cbuf, "images/JTJJTnb_mua_%03d.ppm", pw);
	    WritePPM (outv, raster.GDim(), 0, 0, cbuf);
	    raster.Map_SolToGrid (JTJJT_kappa, outv);
	    sprintf (cbuf, "images/JTJJTnb_kappa_%03d.ppm", pw);
	    WritePPM (outv, raster.GDim(), 0, 0, cbuf);
	    bb = Ax(J,ATx(J,bb));
	}
	exit(0);
#endif

	LOGOUT_1PRM("Gradient norm: %f", l2norm(r));
	
	// add prior to gradient
	RVector gpsi = reg->GetGradient (x);
	gpsi *= M; // apply Hessian normalisation
	r -= gpsi;
	LOGOUT_1PRM("Gradient norm with penalty: %f", l2norm(r));

	bool KeepGoing = true;
	i_count = 0;
	err00 = err0;
	while (i_count < itmax && KeepGoing) {
	    LOGOUT_1PRM ("LM iteration %d", i_count);

	    // solve Hh = r

	    switch (precon) {
	    case LM_PRECON_NONE:    // assume H = I
	        h = r;
		break;
	    case LM_PRECON_HDIAG:   // assume H = diag(H)
		LOGOUT ("Calculating preconditioner: LM-DIAG");
	        h = r / Hess_diag (hdata, x);
		break;
	    case LM_PRECON_CH: {    // solve with Cholesky factorisation
		LOGOUT ("Calculating preconditioner: LM-CH");
	        RSymMatrix hess = Hess_full (hdata, x);// need explicit Hessian
		CHdecomp (hess);
		h = CHsubst(hess, r);
	        } break;
	    case LM_PRECON_ICH:     // Incomplete Cholesky factorisation
	        cerr << "LM ICH preconditioner not yet implemented" << endl;
		exit(1);
	    case LM_PRECON_GMRES: {  // GMRES solver with implicit Hessian
		LOGOUT ("Calculating preconditioner: LM-GMRES");
	        double tol = gmres_tol;
		GMRES (JTJx_clbk, &hdata, r, h, tol);
		//GMRES (Hess_full (hdata, x), r, h, tol); // sanity check
	        } break;
	    }
	    d = h;
	    double stepsize = l2norm(h)/l2norm(r); // Newton step/SD step
	    stepsize = 0.1;
	    LOGOUT_1PRM ("Step size: %f", stepsize);

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
		alpha = 1.0; // fixed step length
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
		FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
		FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
		ProjectAll (mesh, FWS, mvec, dphi, proj);
		err1 = ObjectiveFunction::get_value (data, proj, sd, cov);
	    }

	    if(err1 < err0) {  // accept update
	        //err00 = err0;
		err0 = err1;
		x = x1;
		bsol.SetActiveParams (pscaler->Unscale(x));
		raster.Map_SolToMesh (bsol, msol);
		FWS.Reset (msol, omega);
		FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
		FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
		ProjectAll (mesh, FWS, mvec, dphi, proj);

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
		    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
		    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
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
	    LOGOUT_2PRM("error0: %f, error1: %f", err0, err1);
	    LOGOUT_3PRM("lambda: %g, alpha: %f, tau: %f", lambda, alpha, tau);
	    i_count++;
	}
	times (&tm);
	LOGOUT_3PRM("Iteration %d  CPU %f  OF %f",
	    inr+1, (double)(tm.tms_utime-clock0)/(double)HZ, err0);
    } // end of NR loop;

    // final residuals
    double rd = ObjectiveFunction::get_value (data, proj, sd, cov);
    double rp = reg->GetValue (x);
    LOGOUT_4PRM("Residuals  TAU %g  DATA %g  PRIOR/TAU %g  IT %d", tau, rd, rp/tau, inr);
}

void SolverLM::ReadParams (ParamParser &pp)
{
    char cbuf[256], cbuf2[256], c;
    bool def = false;

    // 1. === PRECONDITIONER ===

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
    }
    while (!def) {
        int cmd;
	cout << "\nSelect LM preconditioner:\n";
	cout << "(0) None\n";
	cout << "(1) Diagonal of Hessian\n";
	cout << "(2) Cholesky factorisation of full Hessian\n";
	cout << "(3) Incomplete CH factorisation of sparse Hessian\n";
	cout << "(4) GMRES solver with implicit Hessian\n";
	cout << "[0|1|2|3|4] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: precon = LM_PRECON_NONE,  def = true; break;
	case 1: precon = LM_PRECON_HDIAG, def = true; break;
	case 2: precon = LM_PRECON_CH,    def = true; break;
	case 3: precon = LM_PRECON_ICH,   def = true; break;
	case 4: precon = LM_PRECON_GMRES, def = true; break;
	}
    }

    // 2. === GMRES CONVERGENCE CRITERION ===

    if (precon == LM_PRECON_GMRES) {
        if (!pp.GetReal ("LM_GMRES_TOL", gmres_tol) ||
	  gmres_tol <= 0.0) do {
	    cout << "\nSelect LM GMRES convergence criterion (>0):\n";
	    cout << ">> ";
	    cin >> gmres_tol;
	} while (gmres_tol <= 0.0);
    }

    // 3. === LINE SEARCH ===

    if (!pp.GetBool ("LM_LINESEARCH", do_linesearch)) {
	char cmd;
	cout << "\nPerform line search in LM solver?\n";
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

    // 4. === CONTROL PARAMETER LAMBDA ===

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

    // 5. === PRIOR ===

    def = false;
    if (pp.GetString ("LM_PRIOR", cbuf)) {
	if (!strcasecmp (cbuf, "NONE"))
	    prior = PRIOR_NONE, def = true;
	else if (!strcasecmp (cbuf, "LAPLACIAN"))
	    prior = PRIOR_LAPLACIAN, def = true;
	else if (!strcasecmp (cbuf, "DIAG"))
	    prior = PRIOR_DIAG, def = true;
    }
    while (!def) {
	int cmd;
	cout << "\nSelect a regularisation method:\n";
	cout << "(0) None\n";
	cout << "(1) Laplacian\n";
	cout << "(2) Diag\n";
	cout << "[0|1|2] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: prior = PRIOR_NONE, def = true;      break;
	case 1: prior = PRIOR_LAPLACIAN, def = true; break;
	case 2: prior = PRIOR_DIAG, def = true;      break;
	}
    }

    // 6. === REGULARISATION PARAMETER TAU ===

    if (!pp.GetReal ("LM_TAU", tau)) {
	cout << "\nSelect LM regularisation parameter tau:\n";
	cout << ">> ";
	cin >> tau;
    }

    // 7. === JACOBIAN FILL FRACTION ===

    if (!pp.GetReal ("JACOBIAN_FILL", Jfill)) {
	cout << "\nJacobian fill fraction (1 for dense):\n";
	cout << ">> ";
	cin >> Jfill;
    }

    // 8. === COVARIANCE DATA ===

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
    case LM_PRECON_GMRES:
	pp.PutString ("LM_PRECON", "GMRES");
	pp.PutReal ("LM_GMRES_TOL", gmres_tol);
	break;
    }

    pp.PutBool ("LM_LINESEARCH", do_linesearch);
    if (do_linesearch)
	pp.PutReal ("LM_MIN_STEPLENGTH", alpha_min);

    pp.PutReal ("LM_LAMBDA_0", lambda0);
    pp.PutReal ("LM_LAMBDA_SCALE", lambda_scale);
    pp.PutReal ("JACOBIAN_FILL", Jfill);

    switch (prior) {
    case PRIOR_NONE:
	pp.PutString ("LM_PRIOR", "NONE");
	break;
    case PRIOR_LAPLACIAN:
	pp.PutString ("LM_PRIOR", "LAPLACIAN");
	break;
    case PRIOR_DIAG:
	pp.PutString ("LM_PRIOR", "DIAG");
	break;
    }

    pp.PutReal ("LM_TAU", tau);

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

bool LineSearchWithPrior (FwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, const RVector &p0, const RVector &p,
    const Regularisation *reg, const RMatrix *cov)
{
    const int MAXIT = 6;
    double x0 = 0.0, xm, fm, fmd, fmp, fb, fbd, fbp, fmind, fminp,xb = lambda;
    RVector proj(data.Dim());
    
    RVector p1 = p+grad*xb;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    if (CheckRange (meshsol)) {
	FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
	fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
	fbp = reg->GetValue (p1);
	fb  = fbd + fbp;
	LOGOUT_4PRM("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g", xb,fbd,fbp,fb);
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
	    FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
	    fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
	    fbp = reg->GetValue (p1);
	    fb  = fbd + fbp;
	    LOGOUT_4PRM ("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
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
		FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
		fbd = ObjectiveFunction::get_value (data, proj, sd, cov);
		fbp = reg->GetValue (p1);
		fb  = fbd + fbp;
		LOGOUT_4PRM ("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
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
	    FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
	    fmd = ObjectiveFunction::get_value (data, proj, sd, cov);
	    fmp = reg->GetValue (p1);
	    fm  = fmd + fmp;
	    LOGOUT_4PRM ("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
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
		FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
		fmd = ObjectiveFunction::get_value (data, proj, sd, cov);
		fmp = reg->GetValue (p1);
		fm  = fmd + fmp;
		LOGOUT_4PRM ("Lsearch: STEP %g OF %g PRIOR %g TOTAL %g",
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
    FWS.ProjectAll (qvec, mvec, meshsol, omega, proj);
    fmind = ObjectiveFunction::get_value (data, proj, sd, cov);
    fminp = reg->GetValue (p1);
    fmin  = fmind + fminp;
    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
    }
    LOGOUT_4PRM("Lsearch final: STEP %g OF %g PRIOR %g TOTAL %g",
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
    const double MAX_CKAPPA = 1;

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

RVector RescaleHessian (const RMatrix *J, const RCompRowMatrix &Lap,
    double tau)
{
    int i, j, k, nz;
    int m = J->nRows(), n = J->nCols(), n2 = Lap.nCols();
    int *colidx = new int[n];
    double *val = new double[n];
    RVector M(n);

    // J^T J contribution
    for (j = 0; j < m; j++) {
	RVector r = J->Row (j);
	for (i = 0; i < n; i++) M[i] += r[i]*r[i];
    }
    // Psi'' contribution
    for (i = 0; i < n; i++) {
	nz = Lap.SparseRow (i%n2, colidx, val);
	for (k = 0; k < nz; k++)
	    if (colidx[k] == (i%n2)) M[i] += val[k] * tau;
    }
    for (i = 0; i < n; i++) {
	M[i] = 1.0/sqrt (M[i]);
    }
    delete []colidx;
    delete []val;
    return M;
}
