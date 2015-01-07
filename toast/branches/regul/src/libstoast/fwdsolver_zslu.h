// -*-C++-*-
// ==========================================================================
// Interface to SuperLU solver (doublecomplex version)
// ==========================================================================

#ifndef __FWDSOLVER_ZSLU_H
#define __FWDSOLVER_ZSLU_H

#include "fwdsolver.h"

class ZSuperLU_engine;

class STOASTLIB ZSuperLU {
public:
    ZSuperLU ();
    ~ZSuperLU ();
    void Reset (const CCompRowMatrix *F);
    void CalcField (const CVector &qvec, CVector &phi,
	IterativeSolverResult *res) const;
    void CalcFields (const CCompRowMatrix &qvec, CVector *phi,
	IterativeSolverResult *res) const;

private:
    const CCompRowMatrix *A;
    ZSuperLU_engine *engine;
};

#endif // !__FWDSOLVER_ZSLU_H
