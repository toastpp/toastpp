// ==========================================================================
// SolverBFGS: BFGS method (Broyden-Fletcher-Goldfarb-Shanno)

#ifndef __SOLVERBFGS_H
#define __SOLVERBFGS_H

#include "solver.h"

class SolverBFGS: public Solver {
public:
    SolverBFGS (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_BFGS; }
    void Solve (CFwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;
};

#endif // !__SOLVERBFGS_H

