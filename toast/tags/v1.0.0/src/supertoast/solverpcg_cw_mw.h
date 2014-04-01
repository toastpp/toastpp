// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients

#ifndef __SOLVERPCG_CW_MW_H
#define __SOLVERPCG_CW_MW_H

#include "solver_cw.h"

class SolverPCG_CW_MW: public Solver_CW {
public:
    SolverPCG_CW_MW (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_PCG; }
    void Solve (RFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const RCompRowMatrix &qvec, const RCompRowMatrix &mvec);
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

#endif // !__SOLVERPCG_CW_MW
