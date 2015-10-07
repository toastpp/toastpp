// -*-C++-*-
// ==========================================================================
// Interface to SuperLU solver (doublecomplex version)
// ==========================================================================

#ifndef __FWDSOLVER_ZSLU_H
#define __FWDSOLVER_ZSLU_H

#include "fwdsolver.h"
#include "slu_zdefs.h"
#include "supermatrix.h"

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
    void AllocMatrix (const CCompRowMatrix *F);
    void Deallocate();
    void Solve (SuperMatrix *sB, SuperMatrix *sX) const;
    const CCompRowMatrix *A;
    int *perm_c;
    int *perm_r;
    int *etree;
    double *R, *C;
    bool allocated;
    mutable superlu_options_t options;
    mutable SuperMatrix sA, sL, sU;
};

#endif // !__FWDSOLVER_ZSLU_H
