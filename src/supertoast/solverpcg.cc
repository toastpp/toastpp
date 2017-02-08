// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "solverpcg.h"
#include "supertoast.h"
#include "supertoast_util.h"
#include "timing.h"
#include <time.h>

using namespace std;

#define PCG_RESET_INTERVAL 10
#define DJTJ_LIMIT 1e-8
//#define RESCALE_HESSIAN

// ==========================================================================
// external references

extern int g_imgfmt;
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

// ==========================================================================

SolverPCG::SolverPCG (ParamParser *_pp): Solver (_pp)
{
    itmax = 50;
    cg_tol = 1e-8;
    cg_delta = 1e-5;
    alpha = 0.0; // "auto"
    precon = PCG_PRECON_NONE;
}

// This is an implementation of preconditioned nonlinear CG
// from the Shewchuk paper, B5 (pg.53) but without Secant method

void SolverPCG::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
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
    int n    = slen*2;
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
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);

    // Start of Shewchuk implementation

    int i_count = 0; // iteration counter
    int k_count = 0; // reset counter

    LOGOUT("Generating fields and gradient");
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

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
    OF.get_gradient (raster, FWS, proj, dphi, mvec, 0/* &bsol*/, r);

    pscaler->ScaleGradient (bsol.GetActiveParams(), r);
    r += reg->GetGradient (x0);
    r0 = r;
    r = -r;

    of_value = OF.get_posterior (&proj);
    of_prior = reg->GetValue (x0);
    of = of_value + of_prior;
    ofp = of * 10+cg_delta; // make sure stopping criterion is not satisfied

    if (precon != PCG_PRECON_NONE) {
        // calculate preconditioner M
	LOGOUT("Calculating Jacobian ...");
	J.New (ndat, n);
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
            FWS.GetDataScaling(), J);
	J.RowScale (inv(sd));
	pscaler->ScaleJacobian (bsol.GetActiveParams(), J);
	switch (precon) {
	case PCG_PRECON_FULLJTJ:
	    // Full Hessian preconditioner setup
	    LOGOUT ("Generating dense JTJ");
	    ATA_dense (raster, J, JTJ);
	    LOGOUT ("Using preconditioner FULLJTJ");
	    break;
	case PCG_PRECON_SPARSEJTJ:
	    // Incomplete CH preconditioner setup
	    LOGOUT ("Generating sparse JTJ");
	    ATA_sparse (raster, J, JTJ_L, JTJ_d);
	    LOGOUT ("Using preconditioner SPARSEJTJ");
	    break;
	case PCG_PRECON_DIAGJTJ:
	    // Diagonal Hessian preconditioner setup
	    LOGOUT ("Calculating diagonal of JTJ");
	    ATA_diag (J, M);
	    LOGOUT("Range: %f to %f", vmin(M), vmax(M));
#ifdef DJTJ_LIMIT
	    M.Clip (DJTJ_LIMIT, 1e50);
	    LOGOUT("Cutoff at %f", DJTJ_LIMIT);
#endif // DJTJ_LIMIT
	    LOGOUT ("Using preconditioner DIAGJTJ");
	    break;
	}
	LOGOUT("Precon reset interval: %d", reset_count);
    } else {
	LOGOUT ("Using preconditioner NONE");
    }

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
	if (g_imgfmt != IMGFMT_RAW) {
	    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
	    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	}
	if (g_imgfmt != IMGFMT_NIM) {
	    Solution rsol(OT_NPARAM, blen);
	    raster.Map_SolToBasis (bsol, rsol, true);
	    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
	    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	}

	// r = -f'(x)
	LOGOUT ("Generating fields and gradient");
	FWS.Reset (msol, omega);
	//for (i = 0; i < mesh->nQ; i++) dphi[i].Clear();
	FWS.CalcFields (qvec, dphi);
	FWS.CalcFields (mvec, aphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
	OF.get_gradient (raster, FWS, proj, dphi, mvec, &bsol, r);
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
	of_value = OF.get_posterior (&proj);
	of_prior = reg->GetValue(x);
	of = of_value + of_prior;
	delta_old = delta_new;
	delta_mid = r & s;

	k_count++;

	if (precon != PCG_PRECON_NONE && k_count == reset_count) {
	    // re-calculate preconditioner and reset CG
	    LOGOUT ("Calculating Jacobian ...");
	    GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
                FWS.GetDataScaling(), J);
	    J.RowScale (inv(sd)); // data space rescaling
	    pscaler->ScaleJacobian (bsol.GetActiveParams(), J);
	    switch (precon) {
	    case PCG_PRECON_FULLJTJ:
	        LOGOUT ("Generating dense JTJ ...");
		ATA_dense (raster, J, JTJ);
		break;
	    case PCG_PRECON_SPARSEJTJ:
  	        LOGOUT ("Generating sparse JTJ ...");
		ATA_sparse (raster, J, JTJ_L, JTJ_d);
		break;
	    case PCG_PRECON_DIAGJTJ:
	        LOGOUT ("Calculating diagonal of JTJ ...");
		ATA_diag (J, M);
		LOGOUT("Range %f to %f", vmin(M), vmax(M));
#ifdef DJTJ_LIMIT
		M.Clip (DJTJ_LIMIT, 1e50);
		LOGOUT("Cutoff at %f", DJTJ_LIMIT);
#endif // DJTJ_LIMIT
	    }
	}

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

void SolverPCG::ReadParams (ParamParser &pp)
{
    char cbuf[256];
    bool def = false;

    // === CONVERGENCE CRITERION ===
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

// ==========================================================================

void SolverPCG::WriteParams (ParamParser &pp)
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
