// -*-C++-*-
// ==========================================================================
// Nonlinear solver classes
// ==========================================================================

#ifndef __SOLVER_H
#define __SOLVER_H

#include "stoastlib.h"

// ==========================================================================
// List of currently implemented solvers

enum SOLVER {
    SOLVER_PCG,             //   preconditioned conjugate gradient
    SOLVER_LM,              //   Levenberg-Marquardt
    SOLVER_LM_UNDER,        //   underdetermined Levenberg-Marquardt
    SOLVER_BFGS,            //   Broyden-Fletcher-Goldfarb-Shanno
    SOLVER_LBFGS,           //   limited memory version of BFGS
    SOLVER_ART,             //   algebraic reconstruction technique
    SOLVER_BLOCKART,        //   Block-ART solver
    SOLVER_LINEAR           //   linear solver (only for difference recon)
};

// ==========================================================================
// Solver virtual base class

class Solver {
public:
    Solver (ParamParser *_pp = NULL): pp(_pp) {}
    virtual SOLVER Type() = 0;
    virtual void Solve (CFwdSolver &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, Solution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
       double omega) = 0;
    virtual void ReadParams (ParamParser &pp) {}
    virtual void WriteParams (ParamParser &pp) {}
    static Solver *Create (SOLVER solver);
    static Solver *Create (ParamParser *_pp);

protected:
    ParamParser *pp;
};

#endif // !__SOLVER_H
