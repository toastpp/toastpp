// ==========================================================================
// SolverLBFGS: LBFGS method (Limited-memory Broyden-Fletcher-Goldfarb-Shanno)
// This version interfaces to the external liblbfgs library
// ==========================================================================

#include "stoastlib.h"
#include "solverlbfgs.h"
#include "lbfgs.h"
#include "util.h"
#include "timing.h"

using namespace std;

// ==========================================================================
// external references

extern int g_imgfmt;
extern double clock0;

// ==========================================================================
// Data structure to provide context for LBFGS library callback functions

struct LBFGS_DATA {
    const ObjectiveFunction *of;
    const Raster *raster;
    CFwdSolver *fws;
    Regularisation *reg;
    const CCompRowMatrix *qvec, *mvec;
    const Scaler *pscaler;
    Solution *msol;
    Solution *bsol;
    double omega;
};

// ==========================================================================

SolverLBFGS::SolverLBFGS (ParamParser *_pp): Solver (_pp)
{
    itmax = 50;
    epsilon = 1e-6;
    delta = 1e-5;
    history = 5;
}

// ==========================================================================
// Callback function for LBFGS library solver:
// Evaluate objective function and gradient

static lbfgsfloatval_t evaluate (void *instance,
    const lbfgsfloatval_t *x_lbfgs,
    lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
    // Calculate objective function and gradient for given solution x
    LBFGS_DATA *ldata = (LBFGS_DATA*)instance;
    const ObjectiveFunction *of = ldata->of;
    const Raster *raster = ldata->raster;
    CFwdSolver *fws = ldata->fws;
    Regularisation *reg = ldata->reg;
    const CCompRowMatrix *qvec = ldata->qvec;
    const CCompRowMatrix *mvec = ldata->mvec;
    const Scaler *pscaler = ldata->pscaler;
    Solution *msol = ldata->msol;
    Solution *bsol = ldata->bsol;
    double omega = ldata->omega;
    const QMMesh *mesh = fws->MeshPtr();
    int i, nlen = mesh->nlen();
    
    RVector x(n, (double*)x_lbfgs);
    RVector r(n);
    RVector proj;
    double f;

    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    raster->Map_ActiveSolToMesh (pscaler->Unscale(x), *msol);
    fws->Reset (*msol, omega);
    fws->CalcFields (*qvec, dphi);
    proj = fws->ProjectAll_real (*mvec, dphi);
    bsol->SetActiveParams (pscaler->Unscale(x));

    if (g) {
        of->get_gradient (*raster, *fws, proj, dphi, *mvec, bsol, r);
	pscaler->ScaleGradient (bsol->GetActiveParams(), r);
	r += reg->GetGradient(x);
	memcpy (g, r.data_buffer(), n*sizeof(double));
    }

    f = of->get_posterior (&proj);
    f += reg->GetValue (x);

    delete []dphi;

    return f;
}

// ==========================================================================
// Callback function for LBFGS library solver:
// iteration progress output (write images and echo objective function)

static int progress (void *instance, const lbfgsfloatval_t *x_lbfgs,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls)
{
    LBFGS_DATA *ldata = (LBFGS_DATA*)instance;
    Solution *msol = ldata->msol;
    Solution *bsol = ldata->bsol;
    const Raster *raster = ldata->raster;
    const Scaler *pscaler = ldata->pscaler;
    Regularisation *reg = ldata->reg;

    RVector x(n, (double*)x_lbfgs);
    bsol->SetActiveParams (pscaler->Unscale(x));
    raster->Map_SolToMesh (*bsol, *msol, true);

    if (g_imgfmt != IMGFMT_RAW) {
	msol->WriteImg_mua (k, "recon_mua.nim");
	msol->WriteImg_mus (k, "recon_mus.nim");
    }
    if (g_imgfmt != IMGFMT_NIM) {
	int blen = raster->BLen();
	Solution rsol(OT_NPARAM, blen);
	raster->Map_SolToBasis (*bsol, rsol, true);
	rsol.WriteImg_mua (k, "recon_mua.raw");
	rsol.WriteImg_mus (k, "recon_mus.raw");
    }
    double pr = reg->GetValue (x);
    LOGOUT ("Iteration %d  CPU %f  OF %g [LH %g PR %g]",
	    k, toc(clock0), fx, fx-pr, pr);
    return 0;
}

// ==========================================================================

void SolverLBFGS::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    // Limited memory BFGS solver
    // Uses identity or diag(JTJ) as initial Hessian

    // Ref: Byrd et. al. "Representations of Quasi-Newton matrices and
    // their use in limited memory methods", Mathematical Programming 63(4),
    // pp. 129-156 (1996)

    LOGOUT ("SOLVER: L-BFGS");

    int n = raster.SLen()*2;
    int ret;
    lbfgs_parameter_t param;

    lbfgs_parameter_init (&param);
    param.m = history;
    param.epsilon = epsilon;
    param.delta = delta;
    param.past = 1;
    param.max_iterations = itmax;

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    Regularisation *reg = Regularisation::Create (pp, &x0, &raster);

    LBFGS_DATA ldata;
    ldata.of = &OF;
    ldata.raster = &raster;
    ldata.fws = &FWS;
    ldata.reg = reg;
    ldata.qvec = &qvec;
    ldata.mvec = &mvec;
    ldata.pscaler = pscaler;
    ldata.msol = &msol;
    ldata.bsol = &bsol;
    ldata.omega = omega;

    // find initial guess for minimizer
    RVector x = pscaler->Scale (bsol.GetActiveParams());

    lbfgsfloatval_t *x_lbfgs = (lbfgsfloatval_t*)x.data_buffer();
    lbfgsfloatval_t fx;

    // initial objective function
    double of = (double)evaluate (&ldata, x_lbfgs, 0, x.Dim(), 0);
    double fp = reg->GetValue (x);
    LOGOUT("Iteration 0  CPU %f  OF %g [LH %g PR %g]", toc(clock0),
	   of, of-fp, fp);

    ret = lbfgs (n, x_lbfgs, &fx, evaluate, progress, &ldata, &param);

    delete reg;
}

// ==========================================================================

void SolverLBFGS::ReadParams (ParamParser &pp)
{
    // === CONVERGENCE LIMIT ===
    if (!(pp.GetReal ("LBFGS_TOL", epsilon) ||
	  pp.GetReal ("NONLIN_TOL", epsilon)) || epsilon <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the convergence criterion epsilon:\n";
	cout << "|g(x)| / max(1,|x|) < epsilon\n";
	cout << "with parameter vector x and gradient g.\n";
	    do {
		cout << "\nNONLIN_TOL (float, >0):\n>> ";
		cin >> epsilon;
	    } while (epsilon <= 0);
    }

    // === STOPPING CRITERION ===
    if (!pp.GetReal ("NONLIN_DELTA", delta)) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the stopping criterion delta. This stops the\n";
	cout << "reconstruction if [f(x_(k-1))-f(x)]/f(x) < delta\n";
        do {
	    cout << "\nNONLIN_DELTA (float, >=0):\n>> ";
	    cin >> delta;
	} while (delta < 0.0);
    }

    // === MAX ITERATION COUNT ===
    if (!(pp.GetInt ("LBFGS_ITMAX", itmax) ||
	  pp.GetInt ("NONLIN_ITMAX", itmax)) || itmax <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the maximum number of iterations to be performed if\n";
	cout << "the convergence and tolerance limits are not satisfied:\n";
	do {
	    cout << "\nNONLIN_ITMAX (int, >0):\n>> ";
	    cin >> itmax;
	} while (itmax <= 0);
    }

    // === NUMBER OF PREVIOUS VECTORS TO STORE ===
    if (!pp.GetInt ("LBFGS_HISTORY", history) || history <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the number of corrections to approximate the inverse\n";
	cout << "Hessian matrix. Larger values require more memory, but can\n";
	cout << "improve convergence rates. Values less than 3 are not\n";
	cout << "recommended.\n";
	do {
	    cout << "\nLBFGS_HISTORY (int, >=3):\n>> ";
	    cin >> history;
	} while (history <= 0);
    }
}

// ==========================================================================

void SolverLBFGS::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "LBFGS");
    pp.PutInt ("NONLIN_ITMAX", itmax);
    pp.PutReal ("NONLIN_TOL", epsilon);
    pp.PutReal ("NONLIN_DELTA", delta);
    pp.PutInt ("LBFGS_HISTORY", history);
}
