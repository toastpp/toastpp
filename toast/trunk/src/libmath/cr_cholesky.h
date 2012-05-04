// ============================================================================
// TOAST v.15                                         (c) Martin Schweiger 1999
// cr_cholesky (module libmath)
//
// Cholesky factorisation routines for SCR (sparse compressed row) matrices
// ============================================================================

#ifndef __CR_CHOLESKY_H
#define __CR_CHOLESKY_H

MATHLIB int symbolic_cholesky_factor (int dim, idxtype *rowptr, idxtype *colidx,
    idxtype *&frowptr, idxtype *&fcolidx);
// For a symmetric positive definite matrix of dimension dim x dim, with
// sparse structure defined by rowptr and colidx (SCR format), generate the
// sparse structure of the lower triangle of its Cholesky factor (excluding
// diagonal) and return it in frowptr and fcolidx.
// frowptr and fcolidx should not be allocated before the call, and should
// be freed after use.

#endif // !__CR_CHOLESKY_H
