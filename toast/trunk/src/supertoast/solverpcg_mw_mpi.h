// ==========================================================================
// SolverPCG: Preconditioned nonlinear conjugate gradients
// Frequency domain, multi-wavelength version

#ifndef __SOLVERPCG_MW_MPI_H
#define __SOLVERPCG_MW_MPI_H

#include "solver_mw_mpi.h"
#include "of.h"

class SolverPCG_MW_MPI: public Solver_MW_MPI {
public:
    SolverPCG_MW_MPI (ParamParser *_pp = NULL);
    ~SolverPCG_MW_MPI ();
    SOLVER Type() { return SOLVER_PCG; }
    void Solve (CFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

private:
    int itmax;                  // max number of CG iterations
    double cg_tol;              // CG convergence criterion
    double alpha;               // step length for line search

    enum PCG_PRECON {           // preconditioners for PCG solver
	PCG_PRECON_NONE,        //   no preconditioner
	PCG_PRECON_DIAGJTJ,     //   diagonal of JTJ
	PCG_PRECON_SPARSEJTJ,   //   sparse JTJ
	PCG_PRECON_FULLJTJ      //   complete JTJ
    } precon;
};

#endif // !__SOLVERPCG_MW_MPI_H
