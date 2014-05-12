// ==========================================================================
// SolverLM: Levenberg-Marquardt
// Frequency domain, multi-wavelength version

#ifndef __SOLVERLM_MW_H
#define __SOLVERLM_MW_H

#include "solver_mw.h"
#include "of.h"

class SolverLM_MW: public Solver_MW {
public:
    SolverLM_MW (ParamParser *_pp = NULL);
    ~SolverLM_MW ();
    SOLVER Type() { return SOLVER_LM; }
    void Solve (CFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

    enum LM_PRECON {       // Hessian inversion method
	LM_PRECON_NONE,    //   no preconditioner
	LM_PRECON_HDIAG,   //   diagonal of Hessian
	LM_PRECON_CH,      //   Cholesky factorisation using explicit Hessian
	LM_PRECON_ICH,     //   Incomplete Cholesky with sparse Hessian
	LM_PRECON_GMRES,   //   "matrix-less" Krylov subspace method (GMRES)
	LM_PRECON_GMRES_JACOBIANFREE,// implicit Jacobian (using d+a fields)
	LM_PRECON_GMRES_DIRECT       // implicit Jacobian (direct method)
    } precon;
       
    enum LM_HESS_SCALING { // Hessian scaling method
	LM_HSCALE_NONE,       // no scaling
	LM_HSCALE_IMPLICIT,   // implicit scaling (requires adjoint fields)
	LM_HSCALE_EXPLICIT,   // explicit scaling
	LM_HSCALE_DISTANCE    // boundary distance scaling
    } hscale;

private:
    Regularisation *reg;   // regularisation instance
    PRIOR_OLD prior;           // regularisation method
    //double tau;            // regularisation parameter
    double lambda0;        // initial value of control parameter
    double lambda_scale;   // scaling factor for lambda adjustment
    double alpha_min;      // minimum step length for line search
    double alpha0;         // initial (DGN) or fixed (LM) step length
    double gn_tol;         // convergence criterion for (outer) GN solver
    double gmres_tol;      // convergence criterion for (inner) GMRES solver
    int nrmax;             // max GN iterations in outer loop
    int itmax;             // max LM iterations in inner loop
    bool do_linesearch;    // perform linesearch for parameter alpha (including
                           // prior)
    double Jfill;          // Jacobian fill fraction (1 for dense)

    // parameters for incorporating model error
    bool modelerr;
    RSymMatrix mcov;
    RVector merr;
    double mcov_lambda;
    char mcov_fmod_fname[128];
    char mcov_farg_fname[128];
    char merr_fmod_fname[128];
    char merr_farg_fname[128];
};

#endif // !__SOLVERLM_MW_H
