// -*-C++-*-
// ==========================================================================
// Module mathlib
// File spmatrix.h
// Declaration of template class TSparseMatrix ('template sparse matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RSparseMatrix = TSparseMatrix<double>	('real')
//	FSparseMatrix = TSparseMatrix<float>	('float')
//	CSparseMatrix = TSparseMatrix<complex>	('complex')
//	ISparseMatrix = TSparseMatrix<int>	('integer')
//	SparseMatrix  = TSparseMatrix<double>	for backward compatibility
//
// Inheritance:
// ------------
// TRootMatrix ----> TSparseMatrix
// ==========================================================================

#ifndef __SPMATRIX_H
#define __SPMATRIX_H

// ==========================================================================
// class TSparseMatrix

template<class MT> class TSparseMatrix: public TRootMatrix<MT> {
public:
    TSparseMatrix ();
    TSparseMatrix (int r, int c);
    TSparseMatrix (const TSparseMatrix<MT> &mat);
    virtual ~TSparseMatrix ();
    operator TMatrix<MT> ();
    // cast to full matrix
  
    MatrixStorage StorageType() const { return MATRIX_COMPROW; }

    void Copy (const TSparseMatrix<MT> &mat);
    void New (int r, int c);

    void Initialise (int *pindex);
    // allocates pindex[i] entries for each row vector i

    void Initialise (int pindex);
    // allocates pindex entries for each row vector

    void Initialise (int *rowptr, int *colidx);
    // initialise matrix from the index lists of a compressed-row
    // type matrix. This allocates storage AND initialises the index
    // list.

    virtual MT Get (int r, int c) const;
    virtual void Put (int r, int c, MT val);
    virtual void Add (int r, int c, MT val);
    TVector<MT> Row (int r) const;
    TVector<MT> Col (int c) const;
    void Unlink ();
    void Shrink ();
    // removes all unused entries (i.e. those with index pointer -1) from
    // the matrix

    void GetFillIn (int *fillin);
    // Returns the physical size of the matrix, by filling array 'fillin' with
    // the pdim value of each row.
    // 'fillin' must be allocated to the correct size before the call
    // Note: this corresponds not necessarily to the number of nonzeros,
    // since allocated elements may have value zero

    void Clear ();
    // sets all entries to zero but without deallocating space

    TSparseMatrix<MT> operator= (const TSparseMatrix<MT> &mat);

    TSparseMatrix<MT> operator* (const MT &mt) const;
    // scalar post multiplication

    TVector<MT> operator* (const TVector<MT> &x) const;
    // matrix * vector operation

    virtual void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    // alternative function to operator*; returns the result in parameter b
    // and avoids local vector creation and copy

    TVector<MT> ATx (const TVector<MT> &x) const;
    // efficient coding of transpose(A) * x

    void ATx (const TVector<MT> &x, TVector<MT> &b) const;
    // alternative function to ATx(Vector&); returns the result in parameter
    // b and avoids local vector creation and copy

    TVector<MT> diag () const;  // returns matrix diagonal as vector

    void SparseOutput (ostream &os) const;

    // friends
    friend TSparseMatrix<MT> transpose (const TSparseMatrix<MT> &A);
	// transpose

    friend void CHdecompFill (const TSparseMatrix<MT> &A,
	TSparseMatrix<int> &L);
    // Calculates the fillin-structure of matrix A after Cholesky decomposition
    // and returns it in lower triangular matrix L. (nonzero entries in the
    // decomposed matrix are indicated by entries of 1 in L)
    // Note that diagonal elements are omitted

    friend bool CHdecomp (const TSparseMatrix<MT> &A, TSparseMatrix<MT> &L,
	TSparseMatrix<MT> &LT, TVector<MT> &d, bool reallocate,
	bool recover);
    // Returns Cholesky decomposition of sparse matrix A.
    // On exit, vector d contains the diagonal of the decomposition, matrix
    // L contains the lower triangle (without the diagonal) of the
    // decomposition, and LT is the transpose of L
    // If 'reallocate' is false then L and LT are not reset, i.e. are assumed
    // to have entries for exactly the required elements allocated on entry.

    friend bool IncompleteCHdecomp (const TSparseMatrix<MT> &A,
        TSparseMatrix<MT> &L, TSparseMatrix<MT> &LT, TVector<MT> &d,
	bool reallocate, bool recover);
    // Returns incomplete Cholesky decomposition of sparse matrix A.
    // On exit, vector d contains the diagonal of the decomposition, matrix
    // L contains the lower triangle (without the diagonal) of the
    // decomposition, and LT is the transpose of L. L and LT retain the
    // fill-in pattern of A. A is assumed symmetric, but must contain the full
    // matrix (not triangle only)
    // If 'reallocate' is false then L and LT are not reset, i.e. are assumed
    // to have entries for exactly the required elements allocated on entry.

    friend void LineCHdecomp (const TSparseMatrix<MT> &A, TMatrix<MT> &T,
	int p);
	// line-wise Cholesky decomposition of matrix A. A must be banded with
	// bandwidth b. T must be a bxb matrix whose contents must be retained
	// between consecutive calls to LineCHdecomp. For the first call, T
	// must be blank. p is the line number to be decomposed, and must be
	// incremented by 1 for each consecutive call, starting from 0.
	// On return, the first line of T contains the decomposition of line p.

    friend TVector<MT> CHsubst (const TSparseMatrix<MT> &L,
	const TSparseMatrix<MT> &LT, const TVector<MT> &d,
	const TVector<MT> &b);

    friend void CHsubst (const TSparseMatrix<MT> &L,
	const TSparseMatrix<MT> &LT, const TVector<MT> &d,
	const TVector<MT> &b, TVector<MT> &x);
	// alternative syntax; returns solution in x, avoiding local creation
	// and copy of vector

    friend TVector<MT> CGsolve (const TSparseMatrix<MT> &A,
	const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
	int *niter = 0);
	// Non-preconditioned conjugate gradient solver

    friend TVector<MT> PCGsolve (const TSparseMatrix<MT> &A,
	const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
	int *niter = 0);
	// Preconditioned conjugate gradient solution of Ax=b. A must be
	// symmetric and positive definite. xinit may be set to an initial
	// guess, otherwise x will be initialised to 0. If niter != 0 it
	// returns the number of CG iterations. This version is internally
	// preconditioned with the diagonal of A.

    friend TVector<MT> PCGsolve (const TSparseMatrix<MT> &A,
	TVector<MT> P, const TVector<MT> &b, TVector<MT> *xinit = 0,
	double err_limit = 1e-18, int *niter = 0);
	// as above, but expects the diagonal preconditioner passed as argument
	// P.

    friend TVector<MT> XCGsolve (const TSparseMatrix<MT> &A,
	const TVector<MT> &P, const TVector<MT> &b, TVector<MT> *xinit,
	double err_limit, int *niter);
	// as above, but expects the diagonal preconditioner passed as argument
	// P.

    friend TVector<MT> PCGsolve (const TSparseMatrix<MT> &A,
	const TSparseMatrix<MT> &Pr, const TSparseMatrix<MT> &PrT,
	const TVector<MT> &PrD, const TVector<MT> &b, TVector<MT> *xinit,
	double err_limit, int *niter);
	// as above, but uses a non-diagonal symmetric matrix as preconditioner:
	// on input, Pr and PrT are the lower triangle and transpose of the
	// preconditioner without the diagonal, and PrD is the diagonal

    friend TVector<MT> BiPCGsolve (const TSparseMatrix<MT> &A,
	const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
	int *niter = 0);

    friend void ComplexBiCGsolve (const TSparseMatrix<MT> &Are,
	const TSparseMatrix<MT> &Aim, const TVector<MT> &bre,
	const TVector<MT> &bim, TVector<MT> &xre, TVector<MT> &xim,
	double err_limit, int *niter);

#ifdef MATH_DEBUG
    TSparseVector<MT> &operator[] (int i) const;	// out-of-line
#else
    TSparseVector<MT> &operator[] (int i) const { return data[i]; }
#endif

protected:
    virtual void Allocate (int r, int c);
    TSparseVector<MT> *data;
};

// ==========================================================================
// inline functions

#ifndef MATH_DEBUG
template<class MT>
inline void TSparseMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    for (int r = 0; r < rows; r++) b[r] = data[r] & x;
}
#endif // !MATH_DEBUG

// ==========================================================================
// friend prototypes

#ifdef NEED_FRIEND_PT

template<class MT>
TSparseMatrix<MT> transpose (const TSparseMatrix<MT> &A);

template<class MT>
void CHdecompFill (const TSparseMatrix<MT> &A, TSparseMatrix<int> &L);

template<class MT>
bool CHdecomp (const TSparseMatrix<MT> &A, TSparseMatrix<MT> &L,
    TSparseMatrix<MT> &LT, TVector<MT> &d, bool reallocate, bool recover);

template<class MT>
bool IncompleteCHdecomp (const TSparseMatrix<MT> &A, TSparseMatrix<MT> &L,
    TSparseMatrix<MT> &LT, TVector<MT> &d, bool reallocate, bool recover);

template<class MT>
void LineCHdecomp (const TSparseMatrix<MT> &A, TMatrix<MT> &T, int p);

template<class MT>
TVector<MT> CHsubst (const TSparseMatrix<MT> &L, const TSparseMatrix<MT> &LT,
    const TVector<MT> &d, const TVector<MT> &b);

template<class MT>
void CHsubst (const TSparseMatrix<MT> &L, const TSparseMatrix<MT> &LT,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x);

template<class MT>
TVector<MT> CGsolve (const TSparseMatrix<MT> &A,
    const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
    int *niter = 0);

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A, TVector<MT> P,
    const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
    int *niter = 0);

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A,
    const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
    int *niter = 0);

template<class MT>
TVector<MT> XCGsolve (const TSparseMatrix<MT> &A,
    const TVector<MT> &P, const TVector<MT> &b, TVector<MT> *xinit = 0,
    double err_limit = 1e-18, int *niter = 0);

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A, const TSparseMatrix<MT> &Pr,
    const TSparseMatrix<MT> &PrT, const TVector<MT> &PrD, const TVector<MT> &b,
    TVector<MT> *xinit, double err_limit, int *niter);

template<class MT>
TVector<MT> BiPCGsolve (const TSparseMatrix<MT> &A,
    const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
    int *niter = 0);

#endif

// ==========================================================================
// typedefs for specific instances of `TMatrix'

typedef TSparseMatrix<double>	RSparseMatrix;	// 'real'
typedef TSparseMatrix<float>	FSparseMatrix;	// 'float'
typedef TSparseMatrix<complex>	CSparseMatrix;	// 'complex'
typedef TSparseMatrix<int>	ISparseMatrix;	// 'integer'
typedef TSparseMatrix<double>	SparseMatrix;	// for backward compatibility


#endif // !__SPMATRIX_H
