// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "mwsolution.h"
#include "solverpcg_mw.h"
#include "supertoast_mw.h"
#include "timing.h"
#include <time.h>

using namespace std;

#define PCG_RESET_INTERVAL 10
#define DJTJ_LIMIT 1e-8
//#define RESCALE_HESSIAN

// ==========================================================================
// external references

extern int g_imgfmt;
extern char g_prefix[256];
extern double clock0;

// ==========================================================================
// local prototypes

void ATA_diag (const RMatrix &A, RVector &diag);
// Builds the diagonal of ATA

void ATA_sparse (const Raster &raster, const RDenseMatrix &A,
    RCompRowMatrix &ATA_L, RVector &ATA_d);
// Creates sparse version ATA (using only neighbours obtained from raster)
// and returns its CH factorisation in ATA_L (lower triangle) and ATA_d
// (diagonal)

void ATA_dense (const Raster &raster, const RDenseMatrix &a,
    RSymMatrix &ata);

void MW_get_gradient (const Raster &raster, CFwdSolverMW &FWS,
		      CVector *dphi,
		      const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd, double omega);

// ==========================================================================

SolverPCG_MW::SolverPCG_MW (ParamParser *_pp): Solver_MW (_pp)
{
    itmax = 50;
    cg_tol = 1e-8;
    cg_delta = 1e-5;
    alpha = 0.0; // "auto"
    precon = PCG_PRECON_NONE;
}

// This is an implementation of preconditioned nonlinear CG
// from the Shewchuk paper, B5 (pg.53) but without Secant method

void SolverPCG_MW::Solve (CFwdSolverMW &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, MWsolution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    const int reset_count = PCG_RESET_INTERVAL;

    // initialisations
    int i, res;
    const QMMesh *mesh = FWS.MeshPtr();
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int nprm = bsol.nActive();
    int n    = bsol.ActiveDim();
    int nofwavel = msol.nofwavel;
    bool pvalid;
    double delta_new, delta_old, delta_mid, delta_0, delta_d, beta, alpha0;
    double of, ofp, of_value, of_prior, fmin;
    double gamma = 1.0;
    RVector r(n), r0(n), s(n), d(n), M(n);
    RVector proj(ndat);
    RDenseMatrix J;
    RSymMatrix JTJ;
    RCompRowMatrix JTJ_L;
    RVector JTJ_d;

    Regularisation *reg;

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    // initial solution (scaled)

    RVector x(x0);
    // current solution (scaled)

    RVector mu0(n/2), kap0(n/2);
    //    x0.Relink(mu0,0,n/2); 
    //    x0.Relink(kap0,n/2,n/2); 
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    //cout << "Initial kappa\n" << kap0 << endl;

    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    // Start of Shewchuk implementation

    int i_count = 0; // iteration counter
    int k_count = 0; // reset counter

    LOGOUT("Generating fields and gradient");
    proj = FWS.ProjectAll_wavel_real (qvec, mvec, msol, omega);

    //ofstream ofs ("dbg_rec.dat");
    //ofs << proj << endl;
    //exit (0);

    // initialise regularisation instance
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
    ofdata.omega = omega;
    ofdata.data = &data;
    ofdata.sd = &sd;

    // r = -f'(x)
    walltic();
    MW_get_gradient (raster, FWS, dphi, qvec, mvec, &msol, r, data, sd, omega);
    LOGOUT("gradient time=%g", walltoc());
    pscaler->ScaleGradient (bsol.GetActiveParams(), r);
    r += reg->GetGradient (x0);
    r0 = r;
    r = -r;

    of_value = ObjectiveFunction::get_value (data, proj, sd);
    of_prior = reg->GetValue (x0);
    of = of_value + of_prior;
    ofp = of * 10+cg_delta; // make sure stopping criterion is not satisfied

    LOGOUT("Iteration 0  CPU %f  OF %g [LH %g PR %g]", toc(clock0),
	   of, of_value, of_prior);

    // apply preconditioner
    switch (precon) {
    case PCG_PRECON_NONE:
        s = r;
	break;
    case PCG_PRECON_DIAGJTJ:
        s = r/M; // element-wise division
	break;
    case PCG_PRECON_SPARSEJTJ:
        CholeskySolve (JTJ_L, JTJ_d, r, s);
	break;
    case PCG_PRECON_FULLJTJ:
        s = CHsubst (JTJ, r);
	break;
    }

    d = s;
    delta_new = r & d;                 // r^t M^-1 r
    delta_0 = delta_new;

    while (i_count < itmax                       // iteration limit
	   && delta_new > cg_tol*cg_tol*delta_0  // convergence criterion
	   && (ofp-of)/of > cg_delta) {          // stopping criterion
        delta_d = d & d;

	if (!alpha) { // initialise step length
	    alpha = of / l2norm (d);
	    LOGOUT("Initial step length reset to %f", alpha);
	}

	// line search. this replaces the Secant method of the Shewchuk code
	alpha0 = alpha;
	res = LineSearch (x, d, alpha0, of, of_clbk, &alpha, &fmin, &ofdata);
	if (res != 0) {
	    LOGOUT ("** Line search failed. Resetting.");
	    d = r;
	    res = LineSearch (x, d, alpha0, of, of_clbk, &alpha, &fmin,&ofdata);
	    if (res != 0) {
	        LOGOUT ("** Line search failed after reset. Terminating.");
		break;
	    }
	}

	x += d*alpha; // update scaled solution
	bsol.SetActiveParams (pscaler->Unscale(x));
	raster.Map_SolToMesh (bsol, msol, true);
	msol.RegisterChange();

	// output parameter images
	if (g_imgfmt != IMGFMT_RAW) {
	    for (i = 0; i < msol.nParam(); i++) {
	        char fname[256];
		if (msol.IsActive(i)) {
		    if (i < msol.nmuaChromo) 
		        sprintf (fname,"%sreconChromophore_%d.nim",
				 g_prefix,i+1);
		    else if (i == msol.nmuaChromo)
		        sprintf (fname,"%sreconScatPrefactor_A.nim",
				 g_prefix);
		    else if (i == msol.nmuaChromo + 1) 
		        sprintf (fname,"%sreconScatPower_b.nim",
				 g_prefix);
		    else
		        xERROR("Invalid parameter index during output");
		    msol.WriteImgGeneric (i_count+1, fname, i);
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
		        sprintf (fname,"%sreconChromophore_%d.raw",
				 g_prefix,i+1);
		    else if (i == msol.nmuaChromo)
		        sprintf (fname,"%sreconScatPrefactor_A.raw",
				 g_prefix);
		    else if (i == msol.nmuaChromo + 1) 
		        sprintf (fname,"%sreconScatPower_b.raw",
				 g_prefix);
		    else
		        xERROR("Invalid parameter index during output");
		    gsol.WriteImgGeneric (i_count+1, fname, i);
		}
	    }
	}

	// r = -f'(x)
	LOGOUT ("Generating fields and gradient");
	proj = FWS.ProjectAll_wavel_real (qvec, mvec, msol, omega);

	walltic();
	MW_get_gradient (raster, FWS, dphi, qvec, mvec, &msol, r, data, sd,
			 omega);
	LOGOUT("gradient time=%g", walltoc());
	pscaler->ScaleGradient (bsol.GetActiveParams(), r);
	r += reg->GetGradient(x);

#ifdef RESCALE_HESSIAN
	RVector x1(x);
	RVector S(x1-x0);
	RVector Y(r-r0);
	gamma = (Y&S) / (Y&Y);
	LOGOUT("Hessian scale ", gamma);
	x0 = x1;
	r0 = r;
#endif

	r = -r;
	ofp = of;
	of_value = ObjectiveFunction::get_value (data, proj, sd);
	of_prior = reg->GetValue(x);
	of = of_value + of_prior;
	delta_old = delta_new;
	delta_mid = r & s;

	k_count++;

	// apply preconditioner
	switch (precon) {
	case PCG_PRECON_NONE:
	    s = r;
	    break;
	case PCG_PRECON_DIAGJTJ:
	    s = (r/M)*gamma; // element-wise division
	    break;
	case PCG_PRECON_SPARSEJTJ:
	    CholeskySolve (JTJ_L, JTJ_d, r, s);
	    break;
	case PCG_PRECON_FULLJTJ:
	    s = CHsubst (JTJ, r);
	    break;
	}

	delta_new = r & s;
	beta = (delta_new - delta_mid) / delta_old;
	if (k_count == reset_count || beta <= 0.0) {
	    d = s;
	    k_count = 0;
	} else {
	    d = s + d * beta;
	}
	i_count++;
	LOGOUT ("Iteration %d  CPU %f  OF %g [LH %g PR %g]",
		i_count, toc(clock0), of, of_value, of_prior);

#ifdef DO_PROFILE
	LOGOUT("Solver time: %f", solver_time);
#endif
    }

    if (delta_new <= cg_tol*cg_tol*delta_0)
        LOGOUT ("PCG solver convergence criterion satisfied");
    else if ((ofp-of)/of <= cg_delta)
        LOGOUT ("PCG solver stopped (insufficient improvement)");
    else if (i_count >= itmax)
        LOGOUT ("PCG solver iteration limit reached");
    else
        LOGOUT ("PCG solver terminated");
}

void SolverPCG_MW::ReadParams (ParamParser &pp)
{
    char cbuf[256];
    bool def = false;

    // === SOLVER CONVERGENCE CRITERION ===
    if (!pp.GetReal ("NONLIN_TOL", cg_tol) || cg_tol <= 0.0) {
        cout << "\nPCG nonlinear solver -----------------------------------\n";
	cout << "Enter the convergence criterion epsilon:\n";
	cout << "<r_k,d_k> < epsilon^2 * <r_0,d_0>\n";
	cout << "for gradient r and preconditioned gradient d\n";
        do {
	    cout << "\nNONLIN_TOL (float, >0):\n>> ";
	    cin >> cg_tol;
	} while (cg_tol <= 0.0);
    }

    // === STOPPING CRITERION ===
    if (!pp.GetReal ("NONLIN_DELTA", cg_delta)) {
        cout << "\nCG nonlinear solver ------------------------------------\n";
	cout << "Enter stopping criterion delta. This stops the\n";
	cout << "reconstruction if [f(x_{k-1})-f(x)]/f(x) < delta\n";
	do {
	    cout << "\nNONLIN_DELTA (float, >=0):\n>> ";
	    cin >> cg_delta;
	} while (cg_delta < 0.0);
    }

    // === MAX ITERATION COUNT ===
    if (!(pp.GetInt ("PCG_ITMAX", itmax) ||
	  pp.GetInt ("NONLIN_ITMAX", itmax)) || itmax <= 0) {
        cout << "\nPCG nonlinear solver -----------------------------------\n";
	cout << "Enter the maximum number of iterations to be performed if\n";
	cout << "the convergence criterion is not satisfied:\n";
	do {
  	    cout << "\nNONLIN_ITMAX (int, >0):\n>> ";
	    cin >> itmax;
	} while (itmax <= 0);
    }

    // === PRECONDITIONER ===
    if (pp.GetString ("PCG_PRECON", cbuf)) {
	if (!strcasecmp (cbuf, "NONE"))
	    precon = PCG_PRECON_NONE,      def = true;
	else if (!strcasecmp (cbuf, "DIAGJTJ"))
	    precon = PCG_PRECON_DIAGJTJ,   def = true;
	else if (!strcasecmp (cbuf, "SPARSEJTJ"))
	    precon = PCG_PRECON_SPARSEJTJ, def = true;
	else if (!strcasecmp (cbuf, "FULLJTJ"))
	    precon = PCG_PRECON_FULLJTJ,   def = true;
    }
    if (!def) {
        cout << "\nPCG nonlinear solver -----------------------------------\n";
	cout << "Select the preconditioner for the nonlinear CG solver:\n\n";
	cout << "(0) None\n";
	cout << "(1) Diagonal of Hessian\n";
	cout << "(2) Sparse Hessian\n";
	cout << "(3) Full Hessian\n";
	while (!def) {
	    cout << "\nPCG_PRECON [0|1|2|3] >> ";
	    int cmd;
	    cin >> cmd;
	    switch (cmd) {
	    case 0: precon = PCG_PRECON_NONE,      def = true; break;
	    case 1: precon = PCG_PRECON_DIAGJTJ,   def = true; break;
	    case 2: precon = PCG_PRECON_SPARSEJTJ, def = true; break;
	    case 3: precon = PCG_PRECON_FULLJTJ,   def = true; break;
	    }
	}
    }

    // === INITIAL STEP LENGTH FOR LINE SEARCH ===
    if (!pp.GetReal ("LS_INIT_STEPLENGTH", alpha) || alpha < 0.0) {
        cout << "\nPCG nonlinear solver -----------------------------------\n";
	cout << "Select the initial step length for the line search\n";
	do {
	    cout << "\nLS_INIT_STEPLENGTH (float, >=0, 0=auto):\n>> ";
	    cin >> alpha;
	} while (alpha < 0.0);
    }
}

void SolverPCG_MW::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "PCG");
    pp.PutReal ("NONLIN_TOL", cg_tol);
    pp.PutReal ("NONLIN_DELTA", cg_delta);
    pp.PutInt ("NONLIN_ITMAX", itmax);

    switch (precon) {
    case PCG_PRECON_NONE:
	pp.PutString ("PCG_PRECON", "NONE");
	break;
    case PCG_PRECON_DIAGJTJ:
	pp.PutString ("PCG_PRECON", "DIAGJTJ");
	break;
    case PCG_PRECON_SPARSEJTJ:
	pp.PutString ("PCG_PRECON", "SPARSEJTJ");
	break;
    case PCG_PRECON_FULLJTJ:
	pp.PutString ("PCG_PRECON", "FULLJTJ");
	break;
    }

    pp.PutReal ("LS_INIT_STEPLENGTH", alpha);
}

// ==========================================================================

void ATA_diag (const RMatrix &A, RVector &diag)
{
    int i, j;
    int n = A.nRows();
    int m = A.nCols();
    double sum, Ai;

    for (j = 0; j < m; j++) {
        sum = 0.0;
	for (i = 0; i < n; i++) {
	    Ai = A(i,j);
	    sum += Ai*Ai;
	}
	diag[j] = sum;
    }
}

// ==========================================================================

void ATA_sparse (const Raster &raster, const RDenseMatrix &a,
    RCompRowMatrix &ata_L, RVector &ata_d)
{
    int i, i0, i1, k, row, col;
    idxtype *rowptr, *rowptr2, *colidx, *colidx2;
    int nzero, nzero2;
    int slen = raster.SLen();
    int n    = slen*2; // mua and kappa
    int nr   = a.nRows();

    LOGOUT ("Building neighbour graph");
    raster.NeighbourGraph (rowptr, colidx, nzero);

    // since we need mua and kappa blocks in ata we need to
    // duplicate the nonzero structure in both rows and columns
    LOGOUT ("Building ATA sparsity pattern");
    nzero2 = nzero*4;
    rowptr2 = new idxtype[n+1];
    colidx2 = new idxtype[nzero2];
    for (i = 0; i <= slen; i++)
        rowptr2[i] = rowptr[i]*2;
    for (i = 1; i <= slen; i++)
        rowptr2[i+slen] = rowptr2[i+slen-1] + rowptr2[i] - rowptr2[i-1];
    for (i = 0; i < slen; i++) {
        int j, c, nc = rowptr[i+1] - rowptr[i];
	for (j = 0; j < nc; j++) {
	    c = colidx[rowptr[i]+j];
	    colidx2[rowptr2[i]+j]         = c;
	    colidx2[rowptr2[i]+j+nc]      = c+slen;
	    colidx2[rowptr2[i+slen]+j]    = c;
	    colidx2[rowptr2[i+slen]+j+nc] = c+slen;
	}
    }
    delete []rowptr;
    delete []colidx;

    LOGOUT ("Building sparse ATA");
    double v, *valptr2 = new double[nzero2];
    for (row = 0; row < n; row++) {
        //cerr << row << endl;
        i0 = rowptr2[row]; i1 = rowptr2[row+1];
	for (i = i0; i < i1; i++) {
	    col = colidx2[i];
	    for (v = 0.0, k = 0; k < nr; k++)
	        v += a(k,row)*a(k,col);
	    valptr2[i] = v;
	}
    }

    RCompRowMatrix ata(n, n, rowptr2, colidx2, valptr2);

    // sanity check: test symmetry
    for (i = 0; i < n; i++) {
        idxtype c1 = rowptr2[i];
	idxtype c2 = rowptr2[i+1];
	for (int j = c1; j < c2; j++) {
  	    idxtype c = colidx2[j];
	    double v = valptr2[j];
	    if (v != ata(c,i))
	        LOGOUT ("JTJ not symmetric!");
	}
    }

    ata.PrintFillinGraph ("fillin.pgm", 600, true, false);
    delete []rowptr2;
    delete []colidx2;
    delete []valptr2;

    double ata_min, ata_max;
    for (i = 0; i < n; i++) {
        double ata_ii = ata(i,i);
	if (!i || ata_ii < ata_min) ata_min = ata_ii;
	if (!i || ata_ii > ata_max) ata_max = ata_ii;
    }
    LOGOUT("ATA diagonal range %f to %f", ata_min, ata_max);


#ifdef RESCALE_JTJ
    // rescale from diagonal norm
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        ata_ii = ata(i,i);
	sum += ata_ii * ata_ii;
    }
    sum = sqrt (sum) / n;
    sum *= JTJ_SCALE;
    for (i = 0; i < n; i++)
        ata(i,i) += sum;
    LOGOUT("Added %f to ATA diagonal", sum);
#endif

    ata.CalculateIncompleteCholeskyFillin (rowptr, colidx);

    ata_L.New (n, n);
    ata_L.Initialise (rowptr, colidx);
    ata_d.New (n);

    while (!IncompleteCholeskyFactorize (ata, ata_L, ata_d, true)) {
        for (i = 0; i < n; i++) ata(i,i) *= 2.0;
	LOGOUT ("*** ATA not positive definite. Scaling diagonal up.");
    }

    delete []rowptr;
    delete []colidx;
}

// ==========================================================================

void ATA_dense (const Raster &raster, const RDenseMatrix &a,
		RSymMatrix &ata)
{
    ata = ATA(a);

#ifdef RESCALE_JTJ
    // rescale from diagonal norm
    double ata_ii, sum = 0.0;
    int i, j, n = ata.nRows();
    for (i = 0; i < n; i++) {
        ata_ii = ata(i,i);
	sum += ata_ii * ata_ii;
    }
    sum = sqrt (sum) / n;
    sum *= JTJ_SCALE;
    for (i = 0; i < n; i++)
        ata(i,i) += sum;
    LOGOUT("Added %f to ATA diagonal", sum);
#endif

    LOGOUT ("Calculating CH decomposition of ATA ...");
    if (!CHdecomp (ata, true)) { // not positive definite
        LOGOUT ("*** ATA not positive definite. Aborting.");
	exit (1);
    }
}

// ==========================================================================
// note: threaded version disabled for now
// Doesn't seem to provide any runtime benefit in the current form
// Also may be responsible for occasional (non-reproducible) CTDs

#if THREAD_LEVEL==2

struct SINGLEGRADIENT_THREADDATA {
    const QMMesh *mesh;
    const Raster *raster;
    CFwdSolverMW *fws;
    const CCompRowMatrix *mvec;
    const RVector *proj;
    const RVector *data;
    const RVector *sd;
    CVector *dphi;
    RVector *grad;
};

void single_gradient_data_engine (task_data *td)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg;
    int itask = td->proc;
    int ntask = td->np;
    SINGLEGRADIENT_THREADDATA *thdata = (SINGLEGRADIENT_THREADDATA*)td->data;
    const QMMesh *mesh = thdata->mesh;
    const Raster *raster = thdata->raster;
    CFwdSolverMW *fws = thdata->fws;
    const CCompRowMatrix *mvec = thdata->mvec;
    const RVector &proj = *thdata->proj;
    const RVector &data = *thdata->data;
    const RVector &sd = *thdata->sd;
    CVector *dphi = thdata->dphi;
    int nq = mesh->nQ;
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    std::complex<double> term;
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    RVector *grad = thdata->grad;
    RVector grad_cmua_loc (slen);
    RVector grad_ckap_loc (slen);
    CVector wqa (mesh->nlen());
    RVector wqb (mesh->nlen());
    CVector dgrad (slen);
    ofs_mod = 0;          // data offset for Mod data
    ofs_arg = mesh->nQM;  // data offset for Arg data
    for (i = 0; i < q0; i++) {
	ofs_mod += mesh->nQMref[i];
	ofs_arg += mesh->nQMref[i];
    }
    
    for (q = q0; q < q1; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector b_mod(n);
	RVector b_arg(n);
	RVector s_mod(n);
	RVector s_arg(n);
	for (i = 0; i < n; i++) {
	    s_mod[i] = sd[ofs_mod+i];
	    s_arg[i] = sd[ofs_arg+i];
 	    b_mod[i] = (data[ofs_mod+i]-proj[ofs_mod+i])/s_mod[i];
	    b_arg[i] = (data[ofs_arg+i]-proj[ofs_arg+i])/s_arg[i];
	}

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);

	cproj = fws->ProjectSingle (q, *mvec, dphi[q], DATA_LIN);

	wqa = std::complex<double>(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec->Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    term = std::complex<double>(
	        -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]),
		-2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]));

	    // amplitude term
	    wqa += qs * std::complex<double> (term.real()*rp*dn,
					      -term.real()*ip*dn);

	    // phase term
	    wqa += qs * std::complex<double> (-term.imag()*ip*dn,
					      -term.imag()*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh->nlen());
	fws->CalcField (wqa, wphia);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua_loc -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckap_loc -= Re(dgrad);

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
    }

    // assemble into global gradient vector
    Task::UserMutex_lock();
    {
        RVector grad_cmua(*grad, 0, slen);    // mua part of grad
	RVector grad_ckap(*grad, slen, slen); // kappa part of grad
	grad_cmua += grad_cmua_loc;
	grad_ckap += grad_ckap_loc;
    }
    Task::UserMutex_unlock();
}

RVector single_gradient_data (const Raster &raster,
    CFwdSolverMW &FWS, const RVector &proj, CVector *dphi,
    const CCompRowMatrix &mvec, const RVector &data,
    const RVector &sd)
{
    const QMMesh &mesh = *FWS.meshptr;
    RVector grad(raster.SLen()*2);
    
    struct SINGLEGRADIENT_THREADDATA {
	const QMMesh *mesh;
	const Raster *raster;
	CFwdSolverMW *fws;
	const CCompRowMatrix *mvec;
	const RVector *proj;
	const RVector *data;
	const RVector *sd;
	CVector *dphi;
	RVector *grad;
    } thdata;
    thdata.mesh = &mesh;
    thdata.raster = &raster;
    thdata.fws = &FWS;
    thdata.mvec = &mvec;
    thdata.proj = &proj;
    thdata.data = &data;
    thdata.sd = &sd;
    thdata.dphi = dphi;
    thdata.grad = &grad;

    Task::Multiprocess (single_gradient_data_engine, &thdata);

    // hack: make this general
    //double c = 0.3/1.4;
    return grad/* *c */; // dy/dcmua -> dy/dmua,  dy/dckappa -> dy/dkappa
    // MS120119: removed scaling by c
}

#else

RVector single_gradient_data (const Raster &raster,
    CFwdSolverMW &FWS, const RVector &proj, CVector *dphi,
    const CCompRowMatrix &mvec, const RVector &data,
    const RVector &sd)
{
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    std::complex<double> term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    RVector grad(slen*2);
    CVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    CVector dgrad (slen);
    ofs_mod = 0;         // data offset for Mod data
    ofs_arg = mesh.nQM;  // data offset for Arg data
    RVector grad_cmua(grad, 0, slen);    // mua part of grad
    RVector grad_ckap(grad, slen, slen); // kappa part of grad

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

	RVector y_mod (data, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector y_arg (data, ofs_arg, n);
	RVector s_arg (sd, ofs_arg, n);
	RVector ypm_arg (proj, ofs_arg, n);
	RVector b_arg(n);
	b_arg = (y_arg-ypm_arg)/s_arg;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	cproj = FWS.ProjectSingle (q, mvec, dphi[q], DATA_LIN);
	wqa = std::complex<double>(0,0);
	wqb = 0.0;

	tic();
	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    term = std::complex<double>(
		-2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]),
		-2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]));

	    // amplitude term
	    wqa += qs * std::complex<double> (term.real()*rp*dn,
					      -term.real()*ip*dn);

	    // phase term
	    wqa += qs * std::complex<double> (-term.imag()*ip*dn,
					      -term.imag()*rp*dn);

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
	grad_ckap -= Re(dgrad);

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
	LOGOUT1_PROGRESS (q);
    }

    // hack: make this general
    //double c = 0.3/1.4;
    return grad/* *c */; // dy/dcmua -> dy/dmua,  dy/dckappa -> dy/dkappa
    // MS120119: removed scaling by c
}

#endif

// ==========================================================================

void MW_get_gradient (const Raster &raster, CFwdSolverMW &FWS,
		      CVector *dphi,
		      const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd, double omega)
{
    int i, j;
    const QMMesh *mesh = FWS.MeshPtr();
    RDenseMatrix extcoef = msol->extcoef;
    int nqm     = mesh->nQM;
    int ndat_i  = nqm*2;
    int slen    = raster.SLen();
    int nchromo = extcoef.nCols();
    int nlambda = extcoef.nRows();
    int nprm    = msol->nParam();
    bool bFactorA = msol->IsActive(nchromo);
    bool bPowerb  = msol->IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;
    RVector proj_i, data_i(ndat_i), sd_i(ndat_i);
    RVector sA, sb;

    grad.Clear();
    for (i = 0; i < nlambda; i++) {
	FWS.Reset (*msol->swsol[i], omega);
	FWS.CalcFields (qvec, dphi);
	proj_i = FWS.ProjectAll_real (mvec, dphi);
	for (j = 0; j < nqm; j++) {
	    data_i[j] = data[i*2*nqm+j];
	    data_i[j+nqm] = data[(i*2+1)*nqm+j];
	    sd_i[j] = sd[i*2*nqm+j];
	    sd_i[j+nqm] = sd[(i*2+1)*nqm+j];
	}
	//data_i = RVector(data,i*ndat_i,ndat_i);
	//sd_i = RVector(sd,i*ndat_i,ndat_i);
	RVector sgrad = single_gradient_data (raster, FWS, proj_i, dphi, mvec, data_i, sd_i);
	RVector sgrad_mua(sgrad,0,slen);
	RVector sgrad_kap(sgrad,slen,slen);

	// chromophore contributions
	for (j = 0; j < nchromo; j++) {
	    RVector grad_ch(grad, j*slen, slen);
	    grad_ch += sgrad_mua*extcoef(i,j);
	}

	// scattering prefactor
	if (bFactorA) {
	    raster.Map_MeshToSol (msol->GetJacobianCoeff_A(i), sA);
	    RVector grad_A(grad, nchromo*slen, slen);
	    grad_A += sgrad_kap*sA;
	}

	// scattering power
	if (bPowerb) {
	    raster.Map_MeshToSol (msol->GetJacobianCoeff_b(i), sb);
	    RVector grad_b(grad, (nchromo+1)*slen, slen);
	    grad_b += sgrad_kap*sb;
	}
    }
}
