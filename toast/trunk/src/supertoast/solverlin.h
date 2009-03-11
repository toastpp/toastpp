// ==========================================================================
// SolverLIN: linear solver for difference reconstructions

#ifndef __SOLVERLIN_H
#define __SOLVERLIN_H

#include "solver.h"

class SolverLIN: public Solver {
public:
    SolverLIN (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_LINEAR; }
    void Solve (CFwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);
};

#endif // !__SOLVERLIN_H

