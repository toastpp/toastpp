#include "stoastlib.h"
#include "solver.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "of.h"
#include "solverlm.h"
#include "solverpcg.h"
#include "solverbfgs.h"
#include "solverlbfgs.h"
#include "solverart.h"
#include "solverblockart.h"
#include "solverlin.h"

using namespace std;

// ==========================================================================
// Solver virtual base class

Solver *Solver::Create (SOLVER solver)
{
    switch (solver) {
    case SOLVER_PCG:      return new SolverPCG; break;
    case SOLVER_LM:       return new SolverLM; break;
    case SOLVER_BFGS:     return new SolverBFGS; break;
    case SOLVER_LBFGS:    return new SolverLBFGS; break;
    case SOLVER_ART:      return new SolverART; break;
    case SOLVER_BLOCKART: return new SolverBlockART; break;
    case SOLVER_LINEAR:   return new SolverLIN; break;
    default:              return 0;
    }
}

Solver *Solver::Create (ParamParser *_pp)
{
    char cbuf[256];
    Solver *s = 0;

    if (_pp->GetString ("SOLVER", cbuf)) {
	if (!strcasecmp (cbuf, "PCG"))
	    s = new SolverPCG (_pp);
	else if (!strcasecmp (cbuf, "LM"))
	    s = new SolverLM (_pp);
	else if (!strcasecmp (cbuf, "BFGS"))
	    s = new SolverBFGS (_pp);
	else if (!strcasecmp (cbuf, "LBFGS"))
	    s = new SolverLBFGS (_pp);
	else if (!strcasecmp (cbuf, "ART"))
	    s = new SolverART (_pp);
	else if (!strcasecmp (cbuf, "BLOCKART"))
	    s = new SolverBlockART (_pp);
	else if (!strcasecmp (cbuf, "LINEAR"))
	    s = new SolverLIN (_pp);
    }
    while (!s) {
	int cmd;
	cout << "\nSelect main solver:\n";
	cout << "(1) PCG (preconditioned conjugate gradient)\n";
	cout << "(2) LM (Levenberg-Marquardt)\n";
	cout << "(3) LM_UNDER (underdetermined LM)\n";
	cout << "(4) BFGS (Broyden-Fletcher-Goldfarb-Shanno)\n";
	cout << "(5) LBFGS (limited memory version of BFGS)\n";
	cout << "(6) ART (algebraic reconstruction technique)\n";
	cout << "(7) Block-ART (block version of ART)\n";
	cout << "(8) LINEAR (difference reconstructions only)\n";
	cout << "[1|2|3|4|5|6|7|8] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: s = new SolverPCG (_pp); break;
	case 2: s = new SolverLM (_pp); break;
	case 4: s = new SolverBFGS (_pp); break;
	case 5: s = new SolverLBFGS (_pp); break;
	case 6: s = new SolverART (_pp); break;
	case 7: s = new SolverBlockART (_pp); break;
	case 8: s = new SolverLIN (_pp); break;
	}
    }
    s->ReadParams (*_pp);
    s->WriteParams (*_pp);
    return s;
}
