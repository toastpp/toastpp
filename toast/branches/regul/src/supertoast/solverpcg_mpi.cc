// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "dnsmatrix_mpi.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "of.h"
#include "solverpcg.h"
#include "jacobian_mpi.h"
#include "timing.h"
#include <time.h>

using namespace std;

#define PCG_RESET_INTERVAL 10
#define DJTJ_LIMIT 1e-8
//#define RESCALE_HESSIAN

//#define OUTPUT_PMDF

// ==========================================================================
// external references

extern int bWriteGrad;
extern int bWriteJ;
extern int g_imgfmt;
extern double clock0;

extern char g_meshname[256];
extern int g_nimsize;

void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh);
void WriteImage (const RVector &nim, int imgno, char *nimname);

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

bool LineSearch (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    const ObjectiveFunction &OF, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &p, const RVector &grad, double f0,
    double &fmin, double &lambda, Solution &meshsol, RVector &proj,
    bool &proj_valid);

// ==========================================================================

SolverPCG::SolverPCG (ParamParser *_pp): Solver (_pp)
{
    itmax = 50;
    cg_tol = 1e-8;
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
    int i;
    const QMMesh *mesh = FWS.MeshPtr();
    int dim  = raster.Dim();
    int nqm  = mesh->nQM;
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int n    = slen*2;
    bool pvalid;
    double delta_new, delta_old, delta_mid, delta_0, delta_d, beta;
    double of, of_value, of_prior, fmin, alpha = -1.0;
    double gamma = 1.0;
    RVector r(n), r0(n), s(n), d(n), M(n);
    RVector proj(ndat);
    RVector sdmod(sd, 0, nqm), sdarg(sd, nqm, nqm);
    RDenseMatrixMPI Jmod, Jarg;

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
    FWS.CalcFields (*mesh, mesh->nQ, qvec, dphi);
    FWS.CalcFields (*mesh, mesh->nM, mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

#ifdef OUTPUT_INITIAL_PROJECTION
    ofstream ofs1 ("init_fmod.fem");
    RVector proj1(proj, 0, ndat/2);
    ofs1 << proj1 << endl;
    ofstream ofs2 ("init_farg.fem");
    RVector proj2(proj, ndat/2, ndat/2);
    ofs2 << proj2 << endl;
#endif

    // r = -f'(x)
    OF.get_gradient (raster, FWS, proj, dphi, mvec, &bsol, r);

    pscaler->ScaleGradient (bsol.GetActiveParams(), r);
    r0 = r;

    if (bWriteGrad) {
        RVector r_mua (r, 0, slen);
	RVector r_mus (r, slen, slen);
	RVector rim(blen);
	raster.Map_SolToBasis (r_mua, rim);
	WriteImage (rim, 0, "grad_mua.raw");
	raster.Map_SolToBasis (r_mus, rim);
	WriteImage (rim, 0, "grad_mus.raw");
    }

    r = -r;

    of_value = OF.get_posterior (&proj);
    of_prior = 0.0; // OF.get_prior (&bsol); // FIX!
    of = of_value + of_prior;

    if (precon != PCG_PRECON_NONE) {

        // calculate preconditioner M
	LOGOUT("Calculating Jacobian ...");
	Jmod.New (slen, n);
	Jarg.New (slen, n);
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
	    FWS.GetDataScaling(), Jmod, Jarg);

	// Data scaling
	Jmod.RowScale (inv(sdmod));
	Jarg.RowScale (inv(sdarg));
	
	// Parameter scaling
	Jmod.ColScale (bsol.GetActiveParams());
	Jarg.ColScale (bsol.GetActiveParams());

#ifdef UNDEF // NOT IMPLEMENTED
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
#endif
	LOGOUT("Precon reset interval: %d", reset_count);
    } else {
	LOGOUT ("Using preconditioner NONE");
    }

    LOGOUT("Iteration 0  CPU %f  OF %f  (prior %f)",
        toc(clock0), of, of_prior);

#ifdef UNDEF // NOT IMPLEMENTED
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
#else
    s = r;
#endif

    d = s;
    delta_new = r & d;                 // r^t M^-1 r
    delta_0 = delta_new;

    while (i_count < itmax && delta_new > cg_tol*cg_tol*delta_0) {
        delta_d = d & d;

	if (alpha < 0.0) { // initialise step length
	    alpha = of / l2norm (d);
	    LOGOUT("Initial step length reset to %f", alpha);
	}
	// line search. this replaces the Secant method of the Shewchuk code
	if (LineSearch (FWS, raster, pscaler, OF, qvec, mvec, data, sd, omega,
            x, d, of, fmin, alpha, msol, proj, pvalid)) {

	    x += d*alpha; // update scaled solution
	    bsol.SetActiveParams (pscaler->Unscale(x));
	    raster.Map_SolToMesh (bsol, msol, true);
	    switch (g_imgfmt) {
	    case IMGFMT_NIM:
	        msol.WriteImg_mua (i_count+1, "recon_mua.nim");
		msol.WriteImg_mus (i_count+1, "recon_mus.nim");
		break;
	    case IMGFMT_RAW:
	        Solution rsol(OT_NPARAM, blen);
		raster.Map_SolToBasis (bsol, rsol, true);
		rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
		rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
		break;
	    }

	} else {
	    LOGOUT ("** Line search failed. Resetting.");
	    k_count = reset_count-1; // force reset
	    i_count--; // don't increment iteration counter
	}

	// r = -f'(x)
	LOGOUT ("Generating fields and gradient");
	FWS.Reset (msol, omega);
	//for (i = 0; i < mesh->nQ; i++) dphi[i].Clear();
	FWS.CalcFields (*mesh, mesh->nQ, qvec, dphi);
	FWS.CalcFields (*mesh, mesh->nM, mvec, aphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
	OF.get_gradient (raster, FWS, proj, dphi, mvec, &bsol, r);
	pscaler->ScaleGradient (bsol.GetActiveParams(), r);
	if (bWriteGrad) {
	    RVector r_mua (r, 0, slen);
	    RVector r_mus (r, slen, slen);
	    RVector rim(blen);
	    raster.Map_SolToBasis (r_mua, rim);
	    WriteImage (rim, i_count+1, "grad_mua.raw");
	    raster.Map_SolToBasis (r_mus, rim);
	    WriteImage (rim, i_count+1, "grad_mus.raw");
	}

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
	of_value = OF.get_posterior (&proj);
	of_prior = 0.0; // OF.get_prior (&bsol); // FIX!
	of = of_value + of_prior;
	delta_old = delta_new;
	delta_mid = r & s;

	k_count++;

	if (precon != PCG_PRECON_NONE && k_count == reset_count) {
	    // re-calculate preconditioner and reset CG
	    LOGOUT ("Calculating Jacobian ...");
	    GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
	        FWS.GetDataScaling(), Jmod, Jarg);

	    // Data scaling
	    Jmod.RowScale (inv(sdmod));
	    Jarg.RowScale (inv(sdarg));

	    // Parameter scaling
	    Jmod.ColScale (bsol.GetActiveParams());
	    Jarg.ColScale (bsol.GetActiveParams());

#ifdef UNDEF // NOT IMPLEMENTED
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
#endif
	}

#ifdef UNDEF // NOT IMPLEMENTED
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
#else
	s = r;
#endif

	delta_new = r & s;
	beta = (delta_new - delta_mid) / delta_old;
	if (k_count == reset_count || beta <= 0.0) {
	    d = s;
	    k_count = 0;
	} else {
	    d = s + d * beta;
	}
	i_count++;
	LOGOUT("Iteration %d  CPU %f  OF %f  (prior %f)",
	    i_count, toc(clock0), of, of_prior);

#ifdef DO_PROFILE
	LOGOUT("Solver time: %f", solver_time);
#endif
    }
}

void SolverPCG::ReadParams (ParamParser &pp)
{
    char cbuf[256];
    bool def = false;

    // 1. === NONLINEAR SOLVER CONVERGENCE CRITERION ===

    if (!pp.GetReal ("NONLIN_TOL", cg_tol) ||
	cg_tol <= 0.0) do {
	    cout << "\nSelect convergence criterion for "
		 << "PCG nonlinear solver (>0):\n>> ";
	    cin >> cg_tol;
	} while (cg_tol <= 0.0);

    // 2. === PRECONDITIONER ===

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
    while (!def) {
	int cmd;
	cout << "\nSelect PCG preconditioner:\n";
	cout << "(0) None\n";
	cout << "(1) Diagonal of Hessian\n";
	cout << "(2) Sparse Hessian\n";
	cout << "(3) Full Hessian\n";
	cout << "[0|1|2|3] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: precon = PCG_PRECON_NONE,      def = true; break;
	case 1: precon = PCG_PRECON_DIAGJTJ,   def = true; break;
	case 2: precon = PCG_PRECON_SPARSEJTJ, def = true; break;
	case 3: precon = PCG_PRECON_FULLJTJ,   def = true; break;
	}
    }

    // 3. === MAX ITERATION COUNT ===

    if (!pp.GetInt ("PCG_ITMAX", itmax) || itmax <= 0) {
	do {
	    cout << "\nMax number of PCG iterations (>0):\n";
	    cout << ">> ";
	    cin >> itmax;
	} while (itmax <= 0);
    }
}

void SolverPCG::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "PCG");
    pp.PutReal ("NONLIN_TOL", cg_tol);

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

    pp.PutInt ("PCG_ITMAX", itmax);
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
    int *rowptr, *rowptr2, *colidx, *colidx2, nzero, nzero2;
    int slen = raster.SLen();
    int n    = slen*2; // mua and kappa
    int nr   = a.nRows();

    LOGOUT ("Building neighbour graph");
    raster.NeighbourGraph (rowptr, colidx, nzero);

    // since we need mua and kappa blocks in ata we need to
    // duplicate the nonzero structure in both rows and columns
    LOGOUT ("Building ATA sparsity pattern");
    nzero2 = nzero*4;
    rowptr2 = new int[n+1];
    colidx2 = new int[nzero2];
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
        int c1 = rowptr2[i];
	int c2 = rowptr2[i+1];
	for (int j = c1; j < c2; j++) {
  	    int c = colidx2[j];
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

bool LineSearch (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    const ObjectiveFunction &OF, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &p, const RVector &grad, double f0,
    double &fmin, double &lambda, Solution &meshsol, RVector &proj,
    bool &proj_valid)
{
    LOGOUT ("Linesearch Start");
    const int MAXIT = 16;
    double x0 = 0.0, xm, fm, fb, xb = lambda;
    proj_valid = false;
    
    RVector p1 = p+grad*xb;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    while (!meshsol.Valid()) {
	LOGOUT ("** Invalid nodal parameters. Reducing step size");
	xb *= 0.5;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    }
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fb = ObjectiveFunction::get_value (data, proj, sd);
    LOGOUT("Step  %f  OF %f", xb, fb);

    if (fb < f0) { // increase interval
        xm = xb; fm = fb;
	xb *= 2.0;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	if (!meshsol.Valid()) {
	    LOGOUT ("** Invalid nodal parameters. Truncating step size.");
	    fb = fm*10.0; // invalidate this step
	} else {
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xb, fb);
	}
	while (fb < fm) {
	    x0 = xm; f0 = fm;
	    xm = xb; fm = fb;
	    xb *= 2.0;
	    p1 = p+grad*xb;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xb, fb);
	}
    } else { // decrease interval
        xm = 0.5*xb;
	p1 = p+grad*xm;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	fm = ObjectiveFunction::get_value (data, proj, sd);
	LOGOUT("Step  %f  OF %f", xm, fm);
	int itcount = 0;
	while (fm > f0) {
  	    if (++itcount > MAXIT) return false;
	    xb = xm; fb = fm;
	    xm = 0.5*xb;
	    p1 = p+grad*xm;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fm = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xm, fm);
	}
    }
    // quadratic interpolation
    double a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
    double b = (f0-fb)/(x0-xb) - a*(x0+xb);
    lambda = -b/(2.0*a);
    p1 = p+grad*lambda;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fmin = ObjectiveFunction::get_value (data, proj, sd);
    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
    } else proj_valid = true;
    LOGOUT("Final %f  OF %f", lambda, fmin);
    LOGOUT ("Linesearch End");
    return true;
}
