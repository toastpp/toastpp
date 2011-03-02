// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "mwsolution.h"
#include "solverpcg_cw_mw_mpi.h"
//#include "supertoast_cw_mw.h"
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

void MW_get_gradient (const Raster &raster, RFwdSolverMW &FWS,
		      RVector *dphi,
		      const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd);

// ==========================================================================

SolverPCG_CW_MW_MPI::SolverPCG_CW_MW_MPI (ParamParser *_pp)
: Solver_CW_MW_MPI (_pp)
{
    itmax = 50;
    cg_tol = 1e-8;
    alpha = 0.0; // "auto"
    precon = PCG_PRECON_NONE;
}

// This is an implementation of preconditioned nonlinear CG
// from the Shewchuk paper, B5 (pg.53) but without Secant method

void SolverPCG_CW_MW_MPI::Solve (RFwdSolverMW &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, MWsolution &msol,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec)
{
    const int reset_count = PCG_RESET_INTERVAL;

    // initialisations
    int i;
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
    double delta_new, delta_old, delta_mid, delta_0, delta_d, beta;
    double of, of_value, of_prior, fmin;
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

    RVector *dphi = new RVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    // Start of Shewchuk implementation

    int i_count = 0; // iteration counter
    int k_count = 0; // reset counter

    LOGOUT("Generating fields and gradient");
    for (i = 0; i < nofwavel; i++) {
	FWS.Reset (*msol.swsol[i], 0);
	FWS.CalcFields (qvec, dphi);
	RVector proj_i (proj, i*mesh->nQM, mesh->nQM);
	proj_i = FWS.ProjectAll (mvec, dphi);
    }

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

    // r = -f'(x)
    MW_get_gradient (raster, FWS, dphi, qvec, mvec, &msol, r, data, sd);

    pscaler->ScaleGradient (bsol.GetActiveParams(), r);
    r += reg->GetGradient (x0);
    r0 = r;

    r = -r;

    of_value = ObjectiveFunction::get_value (data, proj, sd);
    of_prior = reg->GetValue (x0);
    of = of_value + of_prior;

#ifdef UNDEF
    if (precon != PCG_PRECON_NONE) {
        // calculate preconditioner M
	LOGOUT("Calculating Jacobian ...");
	J.New (ndat, n);
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
            FWS.GetDataScaling(), J);
	J.RowScale (inv(sd));
	pscaler->ScaleJacobian (bsol.GetActiveParams(), J);
	if (bWriteJ) WriteJacobian (&J, raster, *mesh);
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
	    LOGOUT_2PRM ("Range: %f to %f", vmin(M), vmax(M));
#ifdef DJTJ_LIMIT
	    M.Clip (DJTJ_LIMIT, 1e50);
	    LOGOUT_1PRM ("Cutoff at %f", DJTJ_LIMIT);
#endif // DJTJ_LIMIT
	    LOGOUT ("Using preconditioner DIAGJTJ");
	    break;
	}
	LOGOUT_1PRM ("Precon reset interval: %d", reset_count);
    } else {
	LOGOUT ("Using preconditioner NONE");
    }
#endif


    LOGOUT_3PRM ("Iteration 0  CPU %f  OF %g  (prior %g)",
        toc(clock0), of, of_prior);

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

    while (i_count < itmax && delta_new > cg_tol*cg_tol*delta_0) {
        delta_d = d & d;

	if (!alpha) { // initialise step length
	    alpha = of / l2norm (d);
	    LOGOUT_1PRM ("Initial step length reset to %f", alpha);
	}
	// line search. this replaces the Secant method of the Shewchuk code
	if (LineSearch (x, d, alpha, of, of_clbk, &alpha, &fmin, &ofdata)==0) {
	    x += d*alpha; // update scaled solution
	    bsol.SetActiveParams (pscaler->Unscale(x));
	    raster.Map_SolToMesh (bsol, msol, true);
	    msol.RegisterChange();
	    switch (g_imgfmt) {
	    case IMGFMT_NIM:
		for (i = 0; i < msol.nParam(); i++) {
		    char fname[256];
		    if (msol.IsActive(i)) {
			if (i < msol.nmuaChromo) 
			    sprintf (fname,"reconChromophore_%d.nim",i+1);
			else if (i == msol.nmuaChromo)
			    sprintf (fname,"reconScatPrefactor_A.nim");
			else if (i == msol.nmuaChromo + 1) 
			    sprintf (fname,"reconScatPower_b.nim");
			msol.WriteImgGeneric (i_count+1, fname, i);
		    }
		}
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
	    bsol.SetActiveParams (pscaler->Unscale(x));
	    raster.Map_SolToMesh (bsol, msol, true);
	    msol.RegisterChange();
	    k_count = reset_count-1; // force reset
	    i_count--; // don't increment iteration counter
	}

	// r = -f'(x)
	LOGOUT ("Generating fields and gradient");
	for (i = 0; i < nofwavel; i++) {
	    FWS.Reset (*msol.swsol[i], 0);
	    FWS.CalcFields (qvec, dphi);
	    RVector proj_i (proj, i*mesh->nQM, mesh->nQM);
	    proj_i = FWS.ProjectAll (mvec, dphi);
	}

	MW_get_gradient (raster, FWS, dphi, qvec, mvec, &msol, r, data, sd);
	pscaler->ScaleGradient (bsol.GetActiveParams(), r);
	r += reg->GetGradient(x);

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
	LOGOUT_1PRM ("Hessian scale ", gamma);
	x0 = x1;
	r0 = r;
#endif

	r = -r;
	of_value = ObjectiveFunction::get_value (data, proj, sd);
	of_prior = reg->GetValue(x);
	of = of_value + of_prior;
	delta_old = delta_new;
	delta_mid = r & s;

	k_count++;

#ifdef UNDEF
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
		LOGOUT_2PRM ("Range %f to %f", vmin(M), vmax(M));
#ifdef DJTJ_LIMIT
		M.Clip (DJTJ_LIMIT, 1e50);
		LOGOUT_1PRM ("Cutoff at %f", DJTJ_LIMIT);
#endif // DJTJ_LIMIT
	    }
	}
#endif

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
	LOGOUT_4PRM ("Iteration %d  CPU %f  OF %g  (prior %g)",
	    i_count, toc(clock0), of, of_prior);

#ifdef DO_PROFILE
	LOGOUT_1PRM ("Solver time: %f", solver_time);
#endif
    }
}

void SolverPCG_CW_MW_MPI::ReadParams (ParamParser &pp)
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

    // 4. === INITIAL STEP LENGTH FOR LINE SEARCH ===
    
    if (!pp.GetReal ("LS_INIT_STEPLENGTH", alpha) || alpha < 0.0) do {
	cout << "\nSelect initial step length for line search (0=auto):\n>> ";
	cin >> alpha;
    } while (alpha < 0.0);
}

void SolverPCG_CW_MW_MPI::WriteParams (ParamParser &pp)
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
    LOGOUT_2PRM ("ATA diagonal range %f to %f", ata_min, ata_max);


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
    LOGOUT_1PRM ("Added %f to ATA diagonal", sum);
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
    LOGOUT_1PRM ("Added %f to ATA diagonal", sum);
#endif

    LOGOUT ("Calculating CH decomposition of ATA ...");
    if (!CHdecomp (ata, true)) { // not positive definite
        LOGOUT ("*** ATA not positive definite. Aborting.");
	exit (1);
    }
}

// ==========================================================================

RVector single_gradient_data (const Raster &raster,
    RFwdSolverMW &FWS, const RVector &proj, RVector *dphi,
    const RCompRowMatrix &mvec, const RVector &data,
    const RVector &sd)
{
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    double term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    RVector grad(slen*2);
    RVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    RVector dgrad (slen);
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
        RVector cdfield (glen);
        RVector *cdfield_grad = new RVector[dim];
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

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	RVector cproj(n);
	cproj = FWS.ProjectSingle (q, mvec, dphi[q], DATA_LIN);
	wqa = 0.0;
	wqb = 0.0;

	tic();
	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const RVector qs = mvec.Row(m);
	    double rp = cproj[idx];
	    double dn = 1.0/(rp*rp);

	    // amplitude term
	    term = -/* 2.0 * */ b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * (term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}
	tm_innerloop += toc();

	// adjoint field and gradient
	RVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);

	RVector cafield(glen);
	RVector *cafield_grad = new RVector[dim];
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
	grad_cmua -= dgrad;

	// diffusion contribution
	// multiply complex field gradients
	RVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	tic();
	raster.Map_GridToSol (gk, dgrad);

	tm_grid2sol += toc();
	grad_ckappa -= dgrad;

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
	LOGOUT1_PROGRESS (q);
    }
    return grad;
}

// ==========================================================================

void MW_get_gradient (const Raster &raster, RFwdSolverMW &FWS,
		      RVector *dphi,
		      const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd)
{
    int i, j;
    const QMMesh *mesh = FWS.MeshPtr();
    RDenseMatrix extcoef = msol->extcoef;
    int nqm     = mesh->nQM;
    int slen    = raster.SLen();
    int nchromo = extcoef.nCols();
    int nlambda = extcoef.nRows();
    int nprm    = msol->nParam();
    bool bFactorA = msol->IsActive(nchromo);
    bool bPowerb  = msol->IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;
    RVector proj_i(mesh->nQM);
    RVector data_i(mesh->nQM);
    RVector sd_i(mesh->nQM);
    RVector sA, sb;

    grad.Clear();
    for (i = 0; i < nlambda; i++) {
	FWS.Reset (*msol->swsol[i], 0);
	FWS.CalcFields (qvec, dphi);
	proj_i = FWS.ProjectAll (mvec, dphi);
	data_i = RVector(data,i*nqm,nqm);
	sd_i = RVector(sd,i*nqm,nqm);
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
