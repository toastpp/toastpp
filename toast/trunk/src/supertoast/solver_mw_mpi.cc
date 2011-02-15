#include "stoastlib.h"
#include "solver_mw_mpi.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "of.h"
#include "solverlm_mw_mpi.h"
#include "solverpcg_mw_mpi.h"
//#include "solverlin.h"

using namespace std;

// ==========================================================================
// Solver virtual base class

Solver_MW_MPI *Solver_MW_MPI::Create (SOLVER solver)
{
    switch (solver) {
    case SOLVER_PCG: return new SolverPCG_MW_MPI; break;
    case SOLVER_LM:  return new SolverLM_MW_MPI; break;
    default:         return 0;
    }
}

Solver_MW_MPI *Solver_MW_MPI::Create (ParamParser *_pp)
{
    char cbuf[256];
    Solver_MW_MPI *s = 0;

    if (_pp->GetString ("SOLVER", cbuf)) {
	if (!strcasecmp (cbuf, "PCG"))
	    s = new SolverPCG_MW_MPI (_pp);
	else if (!strcasecmp (cbuf, "LM"))
	    s = new SolverLM_MW_MPI (_pp);
	//else if (!strcasecmp (cbuf, "LINEAR"))
	//    s = new SolverLIN (_pp);
    }
    while (!s) {
	int cmd;
	cout << "\nSelect main solver:\n";
	cout << "(1) PCG (preconditioned conjugate gradient)\n";
	cout << "(2) LM (Levenberg-Marquardt)\n";
	//cout << "(3) LM_UNDER (underdetermined LM)\n";
	//cout << "(4) BFGS (Broyden-Fletcher-Goldfarb-Shanno)\n";
	//cout << "(5) LBFGS (limited memory version of BFGS)\n";
	//cout << "(6) LINEAR (difference reconstructions only)\n";
	cout << "[1|2|3|4|5|6] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: s = new SolverPCG_MW_MPI (_pp); break;
	case 2: s = new SolverLM_MW_MPI (_pp); break;
	//case 6: s = new SolverLIN_MPI (_pp); break;
	}
    }
    s->ReadParams (*_pp);
    s->WriteParams (*_pp);
    return s;
}
