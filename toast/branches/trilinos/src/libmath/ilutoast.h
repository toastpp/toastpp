#ifndef __ILUTOAST_H
#define __ILUTOAST_H

#ifdef HAVE_ILU

#include <ilupack.h>
#include "mathlib.h"

#define ILUPERM_NULL   0
#define ILUPERM_MMD    1
#define ILUPERM_RCM    2
#define ILUPERM_ND     3
#define ILUPERM_INDSET 4
#define ILUPERM_AMF    5
#define ILUPERM_AMD    6
#define ILUPERM_PQ     7

void ILUSetPermtype (int ptype);

//void ILUSolveDGNL (const RCompRowMatrix &A, const RVector &b, RVector &x);
int ILUSolveZGNL (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol=1e-10, double droptol=1e-3, int maxit=500);
int ILUSolveZSYM (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol=1e-10, double droptol=1e-3, int maxit=500);
void CreateZmat (const CCompRowMatrix &A, Zmat *mat);
void CreateDmat (const RCompRowMatrix &A, Dmat *mat);
void CreateCmat (const SCCompRowMatrix &A, Cmat *mat);
void CreateSmat (const FCompRowMatrix &A, Smat *mat);

#endif // HAVE_ILU

#endif // !__ILUTOAST_H
