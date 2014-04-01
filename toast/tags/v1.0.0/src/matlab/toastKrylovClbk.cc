// =========================================================================
// toastKrylovClbk
// Callback function for parametric Krylov solver Hz=g
// This function computes Hx for any x
// Solve linear system Hz = g where H is given by components
// H = J(x)^T J(x) + psi''(x)
//
// RH parameters:
//     1: x (real vector)
//     2: J (dense double matrix)
//     3: g (real vector)
//     4: M (real vector: diagonal of Hessian rescaling matrix)
//     5: lambda (real): LM control parameter
//     6: hReg (regularisation handle)
//     7: tol (real): convergence criterion
// LH parameters:
//     1: z (real vector)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "toastmex.h"
#include "regul.h"
#include "util.h"

// ==========================================================================
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const RMatrix *J;                // Jacobian
    const RVector *M;                // normalisation diagonal matrix
    const double *lambda;            // diagonal scaling factor
    const Regularisation *reg;       // regularisation
    const RCompRowMatrix *RHess;     // Hessian of regularisation operator
};

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
    const RMatrix *J            =  data->J;
    const RCompRowMatrix *RHess =  data->RHess;
    const RVector &M            = *data->M;
    const double lambda         = *data->lambda;
    int m = J->nRows(), n = J->nCols();

    RVector Px(n);

    // add prior to Hessian
    if (RHess) {
	int i, j, k, nz, *colidx = new int[n];
	double *val = new double[n];
	for (i = 0; i < n; i++) {
	    nz = RHess->SparseRow (i, colidx, val);
	    //nz = reg->GetHessianRow (x, i, colidx, val);
	    for (k = 0; k < nz; k++) {
		j = colidx[k];
		Px[i] += val[k] * x[j] * M[i]*M[j];
	    }
	}
	delete []colidx;
	delete []val;
    }

    return ATx (*J, Ax (*J, x)) + Px + lambda*x;
}

// ==========================================================================
// Build the Hessian of the prior from individual parameter contributions

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
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

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m = mxGetM (prhs[1]);
    int n = mxGetN (prhs[1]);

    // copy current solution
    RVector x(n);
    CopyVector (x, prhs[0]);

    // copy Jacobian
    RDenseMatrix J(m,n);
    CopyMatrix (J, prhs[1]);

    // copy Gradient
    RVector g(n);
    CopyVector (g, prhs[2]);

    // copy M
    RVector M(n);
    CopyVector (M, prhs[3]);

    // copy lambda
    double lambda = mxGetScalar (prhs[4]);

    // regularisation handle
    Regularisation *reg = (Regularisation*)Handle2Ptr (mxGetScalar (prhs[5]));

    // copy tolerance
    double tol = mxGetScalar (prhs[6]);

    RVector z(n);
    int nprm = 2; // for now
    RCompRowMatrix RHess = BuildRHessian (reg, x, nprm);
    
    HESS_DATA hdata = {&J, &M, &lambda, reg, &RHess};
    static RPrecon_Diag precon;
    precon.ResetFromDiagonal (M);
    GMRES (JTJx_clbk, &hdata, g, z, tol, &precon, 30);

    // copy result back to MATLAB
    CopyVector (&plhs[0], z);
}
