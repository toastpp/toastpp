// -*-C++-*-
// ==========================================================================
// Nonlinear solver classes
// ==========================================================================

#ifndef __SOLVER_CW_MW_MPI_H
#define __SOLVER_CW_MW_MPI_H

#include "stoastlib.h"
#include "solver.h"

// ==========================================================================
// Solver virtual base class

class Solver_CW_MW_MPI {
public:
    Solver_CW_MW_MPI (ParamParser *_pp = NULL): pp(_pp) {}
    virtual SOLVER Type() = 0;
    virtual void Solve (RFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const RCompRowMatrix &qvec, const RCompRowMatrix &mvec) = 0;
    virtual void ReadParams (ParamParser &pp) {}
    virtual void WriteParams (ParamParser &pp) {}
    static Solver_CW_MW_MPI *Create (SOLVER solver);
    static Solver_CW_MW_MPI *Create (ParamParser *_pp);

protected:
    ParamParser *pp;
};

#endif // !__SOLVER_CW_MW_MPI_H
