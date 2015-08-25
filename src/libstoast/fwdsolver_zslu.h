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
    ZSuperLU (int n=1);
    ~ZSuperLU ();
    void Reset (const CCompRowMatrix *F);
    void CalcField (const CVector &qvec, CVector &phi,
	IterativeSolverResult *res, int en=0) const;
    void CalcFields (const CCompRowMatrix &qvec, CVector *phi,
	IterativeSolverResult *res, int en=0) const;

private:
    const CCompRowMatrix *A;
    ZSuperLU_engine **engine;
	int nengine;
};

#endif // !__FWDSOLVER_ZSLU_H
