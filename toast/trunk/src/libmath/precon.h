// -*-C++-*-
// ==========================================================================
// Module mathlib
// File precon.h
// Declaration of tmplate class TPreconditioner
// Preconditioner for iterative matrix solution methods
// ==========================================================================

#ifndef __PRECON_H
#define __PRECON_H

#include "vector.h"
#include "matrix.h"
#include "dgmatrix.h"
#include "crmatrix.h"

typedef enum {
    PRECON_NULL,
    PRECON_DIAG,
    PRECON_ICH,
    PRECON_DILU,
    PRECON_CG_MULTIGRID
} PreconType;

// ==========================================================================
// class TPreconditioner

template<class MT> class TPreconditioner {
public:
    TPreconditioner() {}
    virtual ~TPreconditioner() {}

    virtual PreconType Type() = 0;

    virtual void Reset (const TMatrix<MT> *A) = 0;
    // Reset preconditioner for matrix A
    
    virtual void Apply (const TVector<MT> &r, TVector<MT> &s) const = 0;
    // Apply preconditioner to r and return the result in s
    // e.g. s = M^-1 r for a preconditioner matrix M

    static TPreconditioner *Create (PreconType type);
    // create a new preconditioner of type 'type' and return pointer to it
};

// ==========================================================================
// class TPrecon_Null
// Dummy preconditioner (M = I)

template<class MT> class TPrecon_Null: public TPreconditioner<MT> {
public: 
    TPrecon_Null() {}
    PreconType Type() { return PRECON_NULL; }
    void Reset (const TMatrix<MT>*) {}
    void Apply (const TVector<MT> &r, TVector<MT> &s) const { s = r; }
};

// ==========================================================================
// class TPrecon_Diag
// diagonal preconditioner (M = diag(A))

template<class MT> class TPrecon_Diag: public TPreconditioner<MT> {
public:
    TPrecon_Diag() {}
    PreconType Type() { return PRECON_DIAG; }
    void Reset (const TMatrix<MT> *A);
    void ResetFromDiagonal (const TVector<MT> &diag);
    void Apply (const TVector<MT> &r, TVector<MT> &s) const;

private:
    TVector<MT> idiag;
};

// ==========================================================================
// class TPrecon_IC
// Incomplete Cholesky preconditioner

template<class MT> class TPrecon_IC: public TPreconditioner<MT> {
public:
    TPrecon_IC() {}
    PreconType Type() { return PRECON_ICH; }
    void Reset (const TMatrix<MT> *A);
    void Apply (const TVector<MT> &r, TVector<MT> &s) const;

private:
    TCompRowMatrix<MT> L;
    TVector<MT>d;
};

// ==========================================================================
// class TPrecon_DILU
// D-ILU: Simple incomplete LU preconditioner where only diagonal elements
// are modified

template<class MT> class TPrecon_DILU: public TPreconditioner<MT> {
public:
    TPrecon_DILU() {}
    PreconType Type() { return PRECON_DILU; }
    void Reset (const TMatrix<MT> *);
    void Apply (const TVector<MT> &r, TVector<MT> &s) const;

private:
    int dim;                   // problem dimension
    TVector<MT>ipivot;         // inverse pivots
    const TMatrix<MT> *A;  // pointer to matrix
};

// ==========================================================================
// class TPrecon_CG_Multigrid
// Uses multigrid CG solution as preconditioner
// Suitable for SPD systems

template<class MT> class TPrecon_CG_Multigrid: public TPreconditioner<MT> {
public:
    TPrecon_CG_Multigrid (const TCompRowMatrix<MT> *AA,
        TCompRowMatrix<MT> **PP, TCompRowMatrix<MT> **RR,
        TPreconditioner<MT> **pprecon, int nlvl, int nngs = 2, int nnmg = 2) : 
      A(AA), P(PP), R(RR), precon(pprecon), maxlevel(nlvl), ngs(nngs),
      nmg(nnmg) {}

    PreconType Type() { return PRECON_CG_MULTIGRID; }

    void Reset (const TMatrix<MT> *);
    // _A points to an *array* of nlvl system matrices, from finest to
    // coarsest grid

    void Apply (const TVector<MT> &r, TVector<MT> &s) const;

private:
    const TCompRowMatrix<MT> *A; // pointer to system matrix array
    TCompRowMatrix<MT> **P;      // interpolator matrix array
    TCompRowMatrix<MT> **R;      // restrictor matrix array
    TPreconditioner<MT> **precon; // coarse-level CG preconditioner
    int maxlevel, ngs, nmg;
    double MGM (const TVector<MT> &r, TVector<MT> &x, int level) const;
    void Smoother (const TCompRowMatrix<MT> &A, const TVector<MT> &r,
        TVector<MT> &x, int itermax) const;
};

// ==========================================================================
// template typedefs

typedef TPreconditioner<double>         RPreconditioner;
typedef TPreconditioner<toast::complex> CPreconditioner;
typedef TPreconditioner<scomplex>       SCPreconditioner;

typedef TPrecon_Null<double>            RPrecon_Null;
typedef TPrecon_Null<toast::complex>    CPrecon_Null;
typedef TPrecon_Null<scomplex>          SCPrecon_Null;

typedef TPrecon_Diag<double>            RPrecon_Diag;
typedef TPrecon_Diag<toast::complex>    CPrecon_Diag;
typedef TPrecon_Diag<scomplex>          SCPrecon_Diag;

typedef TPrecon_IC<double>              RPrecon_IC;
typedef TPrecon_IC<toast::complex>      CPrecon_IC;
typedef TPrecon_IC<scomplex>            SCPrecon_IC;

typedef TPrecon_DILU<double>            RPrecon_DILU;
typedef TPrecon_DILU<toast::complex>    CPrecon_DILU;
typedef TPrecon_DILU<scomplex>          SCPrecon_DILU;

typedef TPrecon_CG_Multigrid<double>    RPrecon_CG_Multigrid;
typedef TPrecon_CG_Multigrid<toast::complex> CPrecon_CG_Multigrid;
typedef TPrecon_CG_Multigrid<scomplex>  SCPrecon_CG_Multigrid;


// ==========================================================================
// Specialised preconditioner: single-precision complex compressed-row
// matrix (requires separate implementation to work with double-precision
// vectors)

class SCCompRowMatrixMixed;

class SCPreconditionerMixed {
public:
    SCPreconditionerMixed() {}
    virtual ~SCPreconditionerMixed() {}

    virtual PreconType Type() = 0;

    virtual void Reset (const SCCompRowMatrixMixed *A) = 0;
    // Reset preconditioner for matrix A
    
    virtual void Apply (const CVector &r, CVector &s) const = 0;
    // Apply preconditioner to r and return the result in s
    // e.g. s = M^-1 r for a preconditioner matrix M

    static SCPreconditionerMixed *Create (PreconType type);
    // create a new preconditioner of type 'type' and return pointer to it
};

class SCPreconMixed_Null: public SCPreconditionerMixed {
public:
    SCPreconMixed_Null() {}
    PreconType Type() { return PRECON_NULL; }
    void Reset (const SCCompRowMatrixMixed*) {}
    void Apply (const CVector &r, CVector &s) const { s = r; }
};

class SCPreconMixed_Diag: public SCPreconditionerMixed {
public:
    SCPreconMixed_Diag() {}
    PreconType Type() { return PRECON_DIAG; }
    void Reset (const SCCompRowMatrixMixed *A);
    void Apply (const CVector &r, CVector &s) const;

private:
    CVector idiag;
};

class SCPreconMixed_DILU: public SCPreconditionerMixed {
public:
    SCPreconMixed_DILU() {}
    PreconType Type() { return PRECON_DILU; }
    void Reset (const SCCompRowMatrixMixed *);
    void Apply (const CVector &r, CVector &s) const;

private:
    int dim;                       // problem dimension
    CVector ipivot;                // inverse pivots
    const SCCompRowMatrixMixed *A; // pointer to matrix
};


#endif // !__PRECON_H
