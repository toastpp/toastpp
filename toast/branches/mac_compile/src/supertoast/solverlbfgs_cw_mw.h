// ==========================================================================
// SolverBFGS: LBFGS method (Limited-memory Broyden-Fletcher-Goldfarb-Shanno)

#ifndef __SOLVERLBFGS_CW_MW_H
#define __SOLVERLBFGS_CW_MW_H

#include "solver_cw.h"

class SolverLBFGS_CW_MW: public Solver_CW {
public:
    SolverLBFGS_CW_MW (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_LBFGS; }
    void Solve (RFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const RCompRowMatrix &qvec, const RCompRowMatrix &mvec);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;   // max. iteration count
    int history; // number of previous steps to be stored
    double epsilon; // convergence limit
    // convergence criterion is: ||g|| < epsilon * max(1,||x||)
    double delta;   // stopping criterion:
    // (f(x_past) - f(x))/f(x) < delta
};

#endif // !__SOLVERBFGS_CW_MW_H

