// ==========================================================================
// SolverLM: Levenberg-Marquardt

#ifndef __SOLVERLM_H
#define __SOLVERLM_H

#include "solver.h"
#include "of.h"

class SolverLM: public Solver {
public:
    SolverLM ();
    SOLVER Type() { return SOLVER_LM; }
    void Solve (FwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega,
       double ftol);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

    enum LM_PRECON {       // Hessian inversion method
	LM_PRECON_NONE,    //   no preconditioner
	LM_PRECON_HDIAG,   //   diagonal of Hessian
	LM_PRECON_CH,      //   Cholesky factorisation using explicit Hessian
	LM_PRECON_ICH,     //   Incomplete Cholesky with sparse Hessian
	LM_PRECON_GMRES    //   "matrix-less" Krylov subspace method (GMRES)
    } precon;
       
private:
    PRIOR prior;           // regularisation method
    double tau;            // regularisation parameter
    double lambda0;        // initial value of control parameter
    double lambda_scale;   // scaling factor for lambda adjustment
    double alpha_min;      // minimum step length for line search
    double gmres_tol;      // convergence criterion for GMRES solver
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

#endif // !__SOLVERLM_H
