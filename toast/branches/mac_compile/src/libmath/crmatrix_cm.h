// -*-C++-*-
// ==========================================================================
// Module mathlib
// File crmatrix_cm.h
// Declaration of class SCCompRowMatrixMixed ("complex mixed")
//
// This is derived from TCompRowMatrix<scomplex> and adds some methods
// for double precision complex arguments.
//
// Inheritance:
// ------------
// TGenericSparseMatrix<scomplex>
//      ----> TCompRowMatrix<scomplex>
//            ----> SCCompRowMatrixMixed
// ==========================================================================

#ifndef __CRMATRIX_CM
#define __CRMATRIX_CM

#include "crmatrix.h"

// ==========================================================================
// Nonmember declarations

class SCCompRowMatrixMixed;
class SCPreconditionerMixed;

int GMRES (const SCCompRowMatrixMixed &A, const CVector &b, CVector &x,
    double &tol, const SCPreconditionerMixed *precon = 0, int restart = 10,
    int maxit = 0, void (*clbk)(void*) = 0);
// Specialised GMRES solver using a single-precision complex compressed row
// matrix, and double precision complex vectors.

// ==========================================================================
// class SCCompRowMatrixMixed

class SCCompRowMatrixMixed: public TCompRowMatrix<std::complex<float> > {
public:
    // constructors, destructor
    SCCompRowMatrixMixed ();
    SCCompRowMatrixMixed (int rows, int cols);
    SCCompRowMatrixMixed (int rows, int cols,
        const idxtype *_rowptr, const idxtype *_colidx,
        const std::complex<float> *data = 0);
    SCCompRowMatrixMixed (const SCCompRowMatrixMixed &m);
    ~SCCompRowMatrixMixed ();

    int SparseRow (int r, idxtype *ci, std::complex<double> *rv) const;

    // matrix x vector methods using double precision vectors
    void Ax (const CVector &x, CVector &b) const;
    void Ax (const CVector &x, CVector &b, int r1, int r2) const;
    void ATx (const CVector &x, CVector &b) const;

    // mixed precision operators
    inline CVector operator* (const CVector &x) const
    { CVector b; Ax (x, b); return b; }
    // matrix x vector multiplication
};

#endif // !__CRMATRIX_CM
