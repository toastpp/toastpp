// -*-C++-*-
// ==========================================================================
// Interface to SuperLU solver (singlecomplex version)
// ==========================================================================

#ifndef __FWDSOLVER_CSLU_H
#define __FWDSOLVER_CSLU_H

#include "fwdsolver.h"

class CSuperLU_engine;

class STOASTLIB CSuperLU {
public:
    CSuperLU ();
    ~CSuperLU ();
    void Reset (const SCCompRowMatrix *F);
    void CalcField (const SCVector &qvec, SCVector &phi,
	IterativeSolverResult *res) const;

private:
    const SCCompRowMatrix *A;
    CSuperLU_engine *engine;
};

#endif // !__FWDSOLVER_CSLU_H
