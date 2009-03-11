#include "solver.h"
#include "raster.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "of.h"
#include "solverlm.h"
#include "solverpcg.h"
#include "solverlin.h"

// ==========================================================================
// Solver virtual base class

Solver *Solver::Create (SOLVER solver)
{
    switch (solver) {
    case SOLVER_PCG: return new SolverPCG; break;
    case SOLVER_LM:  return new SolverLM; break;
    default:         return 0;
    }
}

Solver *Solver::Create (ParamParser &pp)
{
    char cbuf[256];
    Solver *s = 0;

    if (pp.GetString ("SOLVER", cbuf)) {
	if (!strcasecmp (cbuf, "PCG"))
	    s = new SolverPCG;
	else if (!strcasecmp (cbuf, "LM"))
	    s = new SolverLM;
	else if (!strcasecmp (cbuf, "LINEAR"))
	    s = new SolverLIN;
    }
    while (!s) {
	int cmd;
	cout << "\nSelect main solver:\n";
	cout << "(1) PCG (preconditioned conjugate gradient)\n";
	cout << "(2) LM (Levenberg-Marquardt)\n";
	cout << "(3) LM_UNDER (underdetermined LM)\n";
	cout << "(4) BFGS (Broyden-Fletcher-Goldfarb-Shanno)\n";
	cout << "(5) LBFGS (limited memory version of BFGS)\n";
	cout << "(6) LINEAR (difference reconstructions only)\n";
	cout << "[1|2|3|4|5|6] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: s = new SolverPCG; break;
	case 2: s = new SolverLM; break;
	case 6: s = new SolverLIN; break;
	}
    }
    s->ReadParams (pp);
    s->WriteParams (pp);
    return s;
}
