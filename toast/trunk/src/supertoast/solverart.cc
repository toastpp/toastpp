// ==========================================================================
// SolverART: ART and BlockART nonlinear solvers
// ==========================================================================

#include "stoastlib.h"
#include "pscaler.h"
#include "solverart.h"
#include "fwdsolver.h"
#include "jacobian.h"
#include "util.h"
#include "timing.h"

using namespace std;

// ==========================================================================
// external references

extern int bWriteJ;
extern char g_meshname[256];
extern int g_imgfmt;
extern double clock0;

// ==========================================================================
// local prototypes

static RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
			      int nprm);

struct ART_DATA {
    RDenseMatrix *J;
    double tau;
};

// ==========================================================================
// ART algorithm (inner loop, run once over all measurement to obtain
// an update

// Callback function: return a row of the Jacobian
typedef RVector (*clbkART_Ai)(int i, void *context);

void ART_loop (const RVector &b, RVector &x, double lambda,
    clbkART_Ai Ai, void *context)
{
    int i, ii, m = b.Dim(), n = x.Dim();

    // build a row permutation vector
    IVector perm(m);
    for (i = 0; i < m; i++) perm[i] = i;
    for (i = 0; i < 2*m; i++) {
	int i1 = (int)(m * (rand()/(RAND_MAX+1.0)));
	int i2 = (int)(m * (rand()/(RAND_MAX+1.0)));
	int tmp = perm[i1];
	perm[i1] = perm[i2]; perm[i2] = tmp;
    }

    RVector Rgrad(n);
    RDenseMatrix *J = ((ART_DATA*)context)->J;
    double tau = ((ART_DATA*)context)->tau;

    for (ii = 0; ii < m; ii++) {
	i = perm[ii];
	RVector ai = Ai(i, J);

	x += lambda/((ai & ai) + tau) * (ai * (b[i] - (ai & x)));
    }
}


// ==========================================================================
// ART callback function

RVector ART_Ai (int i, void *context)
{
    RDenseMatrix *J = (RDenseMatrix*)context;
    return J->Row(i);
}


// ==========================================================================

SolverART::SolverART (ParamParser *_pp): Solver (_pp)
{
    itmax = 0;
    itmax_art = 10;
    tol = 1e-8;
    lambda = 1000;
}

// ==========================================================================

void SolverART::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    cout << "Starting ART solver ..." << endl;

    int i, q, m, r, ai, iter;
    const QMMesh *mesh = FWS.MeshPtr();
    int nprm = 2;
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int blen = raster.BLen();
    int slen = raster.SLen();
    int n    = slen*nprm;
    int nq   = mesh->nQ;
    int nm   = mesh->nM;

    double err0; // initial objective value
    double err;  // current objective value
    double errp; // objective value in previous iteration

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    RVector x(x0), x1(x);

    RVector proj(ndat);
    RVector upd(n);

    // create Jacobian instance
    RDenseMatrix *J = new RDenseMatrix (ndat, n);

    // create the regularisation instance
    Regularisation *reg = Regularisation::Create (pp, &x0, &raster);
    double tau = reg->GetTau();

    // allocate field vectors
    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);

    LOGOUT("Resetting forward solver");
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

    err0 = ObjectiveFunction::get_value (data, proj, sd);
    err0 += reg->GetValue(x);
    err  = err0;
    errp = 1e10;

    LOGOUT("Iteration 0  CPU %f  OF %f", toc(clock0), err0);

    ART_DATA art_data = {J, tau};

    // LM iteration loop
    for (iter = 0; !itmax || iter < itmax; iter++) {

	// subtract baseline from data
	RVector b = (data-proj)/sd; // negative residual

	// generate Jacobian
	LOGOUT("Calculating Jacobian ...");
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
			  FWS.GetDataScaling(), *J);

	// apply data space rescaling T
	J->RowScale (inv(sd));

	// apply solution space rescaling S
	pscaler->ScaleJacobian (bsol.GetActiveParams(), *J);

	// ART iteration
	upd.Clear();

	LOGOUT("Starting ART iteration ...");
	for (ai = 0; ai < itmax_art; ai++)
	    ART_loop (b, upd, 0.001, ART_Ai, (void*)&art_data);
	LOGOUT("Finished!");

	// update solution
	x1 = x + upd;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(x1), msol);
	FWS.Reset (msol, omega);
	FWS.CalcFields (qvec, dphi);
	FWS.CalcFields (mvec, aphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
	double err1 = ObjectiveFunction::get_value (data, proj, sd);
	err1 += reg->GetValue(x1);

	if (err1 < errp) {
	    errp = err;
	    err = err1;
	    x = x1;
	    bsol.SetActiveParams (pscaler->Unscale(x));
	    raster.Map_SolToMesh (bsol, msol);
	    lambda *= 2.0;

	    switch (g_imgfmt) {
	    case IMGFMT_NIM:
		msol.WriteImg_mua (iter+1, "recon_mua.nim");
		msol.WriteImg_mus (iter+1, "recon_mus.nim");
		break;
	    case IMGFMT_RAW:
		Solution rsol(OT_NPARAM, blen);
		raster.Map_SolToBasis (bsol, rsol, true);
		rsol.WriteImg_mua (iter+1, "recon_mua.raw");
		rsol.WriteImg_mus (iter+1, "recon_mus.raw");
		break;
	    }
	    LOGOUT("Iteration %d  CPU %f  OF %f", iter+1, toc(clock0), err);

	    // test for convergence
	    if (err < err0*tol) {
		LOGOUT("Solution converged to tolerance");
		break;
	    } else if (fabs (err-errp) < tol) {
		LOGOUT("Improvement below tolerance limit");
		break;
	    }

	} else {

	    lambda *= 0.5;
	    raster.Map_ActiveSolToMesh(pscaler->Unscale(x), msol);
	    FWS.Reset (msol, omega);
	    FWS.CalcFields (qvec, dphi);
	    FWS.CalcFields (mvec, aphi);
	    proj = FWS.ProjectAll_real (mvec, dphi);
	    LOGOUT ("No improvement, reducing lambda");
	    iter--;

	}
    }

    delete J;
    delete reg;
}

// ==========================================================================

void SolverART::ReadParams (ParamParser &pp)
{
    // === MAX NUMBER OF (OUTER) ITERATIONS ===
    if (!pp.GetInt ("ITMAX", itmax) || itmax < 0) do {
	cout << "\nMax number of nonlinear iterations (0 for unlimited):\n>> ";
	cin >> itmax;
    } while (itmax < 0);

    // === MAX NUMBER of (INNER) ART ITERATIONS ===
    if (!pp.GetInt ("ITMAX_ART", itmax_art) || itmax_art <= 0) do {
	cout << "\nMax number of (inner) ART iterations (>0):\n>> ";
	cin >> itmax_art;
    } while (itmax_art <= 0);

    // === NONLINEAR SOLVER CONVERGENCE CRITERION ===
    if (!pp.GetReal ("NONLIN_TOL", tol) ||
	tol <= 0.0) do {
	    cout << "\nConvergence criterion for "
		 << "ART nonlinear solver (>0):\n>> ";
	    cin >> tol;
	} while (tol <= 0.0);
}

// ==========================================================================

void SolverART::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "ART");
    pp.PutInt ("ITMAX", itmax);
    pp.PutInt ("ITMAX_ART", itmax_art);
    pp.PutReal ("NONLIN_TOL", tol);
}

// ==========================================================================

#ifdef UNDEF
void SolverART::BlockSolve (const RVector &b, const RDenseMatrix &A,
    RVector &x, const double mu)
{
    // Solves under-determined problem b = Ax
    // by calculating x = A^T (AA^T + mu I)^-1 b

    RSymMatrix AAt = AAT(A);
    for (int i = 0; i < AAt.nRows(); i++)
	AAt(i,i) += mu;
    CHdecomp (AAt);
    RVector xb = CHsubst (AAt, b);
    x = transpose(A) * xb;
}
#endif

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
