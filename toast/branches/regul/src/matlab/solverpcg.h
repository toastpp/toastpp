// ==========================================================================
// SolverPCG: preconditioned nonlinear conjugate gradients

#ifndef __SOLVERPCG_H
#define __SOLVERPCG_H

#include "solver.h"

class SolverPCG: public Solver {
public:
    SolverPCG () {}
    SOLVER Type() { return SOLVER_PCG; }
    void Solve (FwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega,
       double ftol);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;                  // max number of CG iterations

    enum PCG_PRECON {           // preconditioners for PCG solver
	PCG_PRECON_NONE,        //   no preconditioner
	PCG_PRECON_DIAGJTJ,     //   diagonal of JTJ
	PCG_PRECON_SPARSEJTJ,   //   sparse JTJ
	PCG_PRECON_FULLJTJ      //   complete JTJ
    } precon;
};

#endif // !__SOLVER_PCG
