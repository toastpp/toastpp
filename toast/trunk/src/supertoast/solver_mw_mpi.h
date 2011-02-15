// -*-C++-*-
// ==========================================================================
// Nonlinear solver classes
// ==========================================================================

#ifndef __SOLVER_MW_MPI_H
#define __SOLVER_MW_MPI_H

#include "stoastlib.h"
#include "solver.h"

// ==========================================================================
// Solver virtual base class

class Solver_MW_MPI {
public:
    Solver_MW_MPI (ParamParser *_pp = NULL): pp(_pp) {}
    virtual SOLVER Type() = 0;
    virtual void Solve (CFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
       double omega) = 0;
    virtual void ReadParams (ParamParser &pp) {}
    virtual void WriteParams (ParamParser &pp) {}
    static Solver_MW_MPI *Create (SOLVER solver);
    static Solver_MW_MPI *Create (ParamParser *_pp);

protected:
    ParamParser *pp;
};

#endif // !__SOLVER_MW_MPI_H
