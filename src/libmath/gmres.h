#ifndef __GMRES_H
#define __GMRES_H

#include "mathlib.h"

// =====================================================================
//    GMRES(restart) algorithm according to
//    
//    R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
//    V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Solvers
//    SIAM, Philadelphia, 1994
// =====================================================================

template<class MT>
int gmres (int restart, const TMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, TPreconditioner<MT> *precon, double &elim, int maxit, 
    void (*clbk)(void*) = 0);
// Solve Ax = b with generalised minimal residual method

template<class MT>
int gmres (int restart, TVector<MT> (*Av_clbk)(const TVector<MT> &v,
    void *context), const TVector<MT> &b, TVector<MT> &x,
    TPreconditioner<MT> *precon, double &elim, int maxit, void *context);
// This "matrix-less" version of gmres can be used whenever matrix A is
// not available in explicit form. The user provides a callback function
// (Av_clbk) which is called whenever the product of the matrix with a
// vector v is required.

#endif // !__GMRES_H
