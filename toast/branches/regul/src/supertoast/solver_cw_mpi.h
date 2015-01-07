// -*-C++-*-
// ==========================================================================
// Nonlinear solver classes
// ==========================================================================

#ifndef __SOLVER_CW_H
#define __SOLVER_CW_H

#include "mathlib.h"
#include "felib.h"
#include "pparse.h"
#include "solver.h"
#include "fwdsolver_mw.h"

class FwdSolverMW;
class Raster;
class ObjectiveFunction;
class Scaler;
class Solution;

// ==========================================================================
// Solver virtual base class

class Solver_CW {
public:
    Solver_CW (ParamParser *_pp = NULL): pp(_pp) {}
    virtual SOLVER Type() = 0;
    virtual void Solve (RFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const RCompRowMatrix &qvec, const RCompRowMatrix &mvec) = 0;
    virtual void ReadParams (ParamParser &pp) {}
    virtual void WriteParams (ParamParser &pp) {}
    static Solver_CW *Create (SOLVER solver);
    static Solver_CW *Create (ParamParser *_pp);

protected:
    ParamParser *pp;
};

#endif // !__SOLVER_CW_H
