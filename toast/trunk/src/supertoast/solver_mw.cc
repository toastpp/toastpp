#include "stoastlib.h"
#include "solver_mw.h"
#include "pscaler.h"
#include "fwdsolver.h"
#include "of.h"
#include "solverlm_mw.h"
#include "solverlm2_mw.h"
#include "solverpcg_mw.h"
#include "solverlbfgs_mw.h"
//#include "solverlin.h"

using namespace std;

// ==========================================================================
// Solver virtual base class

Solver_MW *Solver_MW::Create (SOLVER solver)
{
    switch (solver) {
    case SOLVER_PCG:   return new SolverPCG_MW; break;
    case SOLVER_LM:    return new SolverLM_MW; break;
    case SOLVER_LBFGS: return new SolverLBFGS_MW; break;
    case SOLVER_LM_UNDER: return new SolverLM2_MW; break;
    default:           return 0;
    }
}

Solver_MW *Solver_MW::Create (ParamParser *_pp)
{
    char cbuf[256];
    Solver_MW *s = 0;

    if (_pp->GetString ("SOLVER", cbuf)) {
	if (!strcasecmp (cbuf, "PCG"))
	    s = new SolverPCG_MW (_pp);
	else if (!strcasecmp (cbuf, "LM"))
	    s = new SolverLM_MW (_pp);
	else if (!strcasecmp (cbuf, "LBFGS"))
	    s = new SolverLBFGS_MW (_pp);
	else if (!strcasecmp (cbuf, "LM2"))
	    s = new SolverLM2_MW (_pp);
	//else if (!strcasecmp (cbuf, "LINEAR"))
	//    s = new SolverLIN (_pp);
    }
    while (!s) {
	int cmd;
	cout << "\nSelect main solver:\n";
	cout << "(1) PCG (preconditioned conjugate gradient)\n";
	cout << "(2) LM (Levenberg-Marquardt)\n";
	cout << "(3) LM2 (LM with implicit chromophore blocks)\n";
	//cout << "(4) BFGS (Broyden-Fletcher-Goldfarb-Shanno)\n";
	cout << "(5) LBFGS (limited memory version of BFGS)\n";
	//cout << "(6) LINEAR (difference reconstructions only)\n";
	cout << "[1|2|3|4|5|6] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: s = new SolverPCG_MW (_pp); break;
	case 2: s = new SolverLM_MW (_pp); break;
	case 3: s = new SolverLM2_MW (_pp); break;
	case 5: s = new SolverLBFGS_MW (_pp); break;
	//case 6: s = new SolverLIN (_pp); break;
	}
    }
    s->ReadParams (*_pp);
    s->WriteParams (*_pp);
    return s;
}
