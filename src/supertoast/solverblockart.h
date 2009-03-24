// ==========================================================================
// SolverART: ART nonlinear inverse solver

#ifndef __SOLVERBLOCKART_H
#define __SOLVERBLOCKART_H

#include "solver.h"

class SolverBlockART: public Solver {
public:
    SolverBlockART (ParamParser *_pp = NULL);
    SOLVER Type() { return SOLVER_BLOCKART; }
    void Solve (CFwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;     // max number of outer (nonlinear) iterations. 0=unlimited
    int itmax_art; // max number of inner (ART) iterations
    double tol;    // nonlinear solver tolerance criterion
    double lambda; // ART relaxation parameter
};

#endif // !__SOLVERBLOCKART_H

