// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients

#ifndef __SOLVERPCG_H
#define __SOLVERPCG_H

#include "solver.h"

class SolverPCG: public Solver {
public:
    SolverPCG (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_PCG; }
    void Solve (CFwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;                  // max number of CG iterations
    double cg_tol;              // CG convergence criterion
    double cg_delta;            // CG stopping criterion
    double alpha;               // step length for line search

    enum PCG_PRECON {           // preconditioners for PCG solver
	PCG_PRECON_NONE,        //   no preconditioner
	PCG_PRECON_DIAGJTJ,     //   diagonal of JTJ
	PCG_PRECON_SPARSEJTJ,   //   sparse JTJ
	PCG_PRECON_FULLJTJ      //   complete JTJ
    } precon;
};

#endif // !__SOLVER_PCG
