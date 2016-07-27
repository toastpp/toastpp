#include "matlabtoast.h"
#include "mexutil.h"
#include "util.h"
#include "timing.h"
#include "lbfgs.h"

// ==========================================================================
// LBFGS solver routines
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
    mxArray *func;
};

// ==========================================================================
// Callback function for LBFGS library solver:
// Evaluate objective function and gradient

static lbfgsfloatval_t evaluate (void *instance,
    const lbfgsfloatval_t *x_lbfgs,
    lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
    int i;
    LBFGS_DATA *ldata = (LBFGS_DATA*)instance;
    mxArray *clbkArg[2];
    clbkArg[0] = ldata->func;
    
    clbkArg[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *pr = mxGetPr(clbkArg[1]);
    for (i = 0; i < n; i++)
	pr[i] = x_lbfgs[i];
    
    mxArray *clbkRes[2];

    mexCallMATLAB (2, clbkRes, 2, clbkArg, "feval");

    double of = mxGetScalar(clbkRes[0]);
    pr = mxGetPr(clbkRes[1]);
    for (i = 0; i < n; i++)
	g[i] = pr[i];
    
    return of;
}

// ==========================================================================
// Callback function for LBFGS library solver:
// iteration progress output (write images and echo objective function)

static int progress (void *instance, const lbfgsfloatval_t *x_lbfgs,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls)
{
    return 0;
}

// ============================================================================
// LBFGS driver function

void Solve_LBFGS (mxArray *func, RVector &x0, Regularisation *reg,
    lbfgs_parameter_t param)
{
    mexPrintf("L-BFGS:\n");
    mexPrintf("  history: %d\n", param.m);
    mexPrintf("  epsilon: %g\n", param.epsilon);
    mexPrintf("  maxit:   %d\n", param.max_iterations);

    int ret;
    int n = x0.Dim();
    RVector x(x0);
    lbfgsfloatval_t *x_lbfgs = (lbfgsfloatval_t*)x.data_buffer();
    lbfgsfloatval_t fx;

    LBFGS_DATA ldata;
    ldata.func = func;
    ldata.reg = reg;

    ret = lbfgs (n, x_lbfgs, &fx, evaluate, progress, &ldata, &param);
}

// ============================================================================
// MATLAB entry point

void MatlabToast::LBFGS (int nlhs, mxArray *plhs[],
			 int nrhs, const mxArray *prhs[])
{
    lbfgs_parameter_t param;
    mxArray *clbkFunc, *field;

    if (nrhs < 3)
	mexErrMsgTxt ("toastLBFGS: insufficient number of parameters");

    // callback function handle
    if (mxGetClassID(prhs[0]) != mxFUNCTION_CLASS)
	mexErrMsgTxt
           ("toastLBFGS: Argument 1: invalid type (function handle expected)");
    clbkFunc = (mxArray*)prhs[0];

    // initial solution vector
    RVector x0;
    CopyVector (x0, prhs[1]);

    // context parameters
    lbfgs_parameter_init (&param);

    if (mxGetClassID(prhs[2]) != mxSTRUCT_CLASS)
	mexErrMsgTxt
	    ("toastLBFGS: Argument 3: invalid type (structure expected)");

    field = mxGetField (prhs[2], 0, "history");
    if (field) {
	if (mxIsNumeric(field))
	    param.m = (int)mxGetScalar(field);
	else
	    mexErrMsgTxt ("toastLBFGS: Argument 3, field 'history': invalid type (numeric expected)");
    } else {
	param.m = 5;
    }

    field = mxGetField (prhs[2], 0, "epsilon");
    if (field) {
	if (mxIsNumeric(field))
	    param.epsilon = mxGetScalar(field);
	else
	    mexErrMsgTxt ("toastLBFGS: Argument 3, field 'epsilon': invalid type (numeric expected)");
    } else {
	param.epsilon = 1e-6;
    }

    field = mxGetField (prhs[2], 0, "itmax");
    if (field) {
	if (mxIsNumeric(field))
	    param.max_iterations = (int)mxGetScalar(field);
	else
	    mexErrMsgTxt ("toastLBFGS: Argument 3, field 'itmax': invalid type (numeric expected)");
    } else {
	param.max_iterations = 50;
    }

    Regularisation *reg = NULL;
    field = mxGetField (prhs[2], 0, "hreg");
    if (field)
	reg = (Regularisation*)Handle2Ptr (mxGetScalar(field));

    Solve_LBFGS (clbkFunc, x0, reg, param);
}
