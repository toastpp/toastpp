// ==========================================================================
// SolverBFGS: LBFGS method (Limited-memory Broyden-Fletcher-Goldfarb-Shanno)

#ifndef __SOLVERLBFGS_MW_H
#define __SOLVERLBFGS_MW_H

#include "solver_mw.h"

class SolverLBFGS_MW: public Solver_MW {
public:
    SolverLBFGS_MW (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_LBFGS; }
    void Solve (CFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;   // max. iteration count
    int history; // number of previous steps to be stored
    double epsilon; // convergence limit
    // stopping criterion is: ||g|| < epsilon * max(1,||x||)
    double delta;   // stopping criterion:
    // (f(x_past) - f(x))/f(x) < delta
};

#endif // !__SOLVERBFGS_CW_MW_H

