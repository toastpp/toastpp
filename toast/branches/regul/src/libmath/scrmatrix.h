// -*-C++-*-
// ==========================================================================
// Module mathlib
// File scrmatrix.h
// Declaration of template class TSymCompRowMatrix ('template symmetric
// compressed-row matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RSymCompRowMatrix = TSymCompRowMatrix<double>	  ('real')
//	FSymCompRowMatrix = TSymCompRowMatrix<float>	  ('float')
//	CSymCompRowMatrix = TSymCompRowMatrix<complex>    ('complex')
//	ISymCompRowMatrix = TSymCompRowMatrix<int>	  ('integer')
//
// Inheritance:
// ------------
// TGenericSparseMatrix ----> TSymCompRowMatrix
// ==========================================================================

#ifndef __SCRMATRIX_H
#define __SCRMATRIX_H

#include "gsmatrix.h"

// ==========================================================================
// Nonmember declarations

template<class MT>class TSymCompRowMatrix;

template<class MT>
bool CholeskyFactorize (const TSymCompRowMatrix<MT> &A, TCompRowMatrix<MT> &L,
    TVector<MT> &d, bool recover = false);

template<class MT>
std::istream &operator>> (std::istream &is, TSymCompRowMatrix<MT> &m);

template<class MT>
std::ostream &operator<< (std::ostream &os, const TSymCompRowMatrix<MT> &m);

// ==========================================================================
// class TSymCompRowMatrix

template<class MT> class TSymCompRowMatrix: public TGenericSparseMatrix<MT> {

public:


    void New (int nrows, int ncols);
   
    TSymCompRowMatrix ();
    TSymCompRowMatrix (int rows, int cols);

    TSymCompRowMatrix (int rows, int cols,
        const int *rptr, const int *cidx, const MT *data=0);
    // creates matrix of logical dimension rows x cols. fill-in structure
    // is given by _rowptr and _colidx (but entries in the upper triangle are
    // ignored. 

    //~TSymCompRowMatrix ();

    MatrixStorage StorageType() const { return MATRIX_SYMCOMPROW; }

    void Initialise (const int *_rowptr, const int *_colidx,
        const MT *data = 0);
    // Redefine sparse structure and assign data
    // If data are not provided, all entries are set to zero

    void AllowColIndexing (bool yes);
    // enable/disable column access index lists

    void SetColAccess (bool yes) const;
    // create/delete index arrays for column access

    MT &operator() (int r, int c);
    MT Get (int r, int c) const;

    void Add (int r, int c, const MT &val)
    { if (r >= c) (*this)(r,c) += val; }

    TVector<MT> Row (int r) const;
    // Returns row r as a vector

    TVector<MT> Col (int c) const
    { ERROR_UNDEF; return TVector<MT>(); }
    // Returns column c as a vector

    int SparseRow (int r, idxtype *ci, MT *rv) const;

    void RowScale (const TVector<MT> &scale)
    { ERROR_UNDEF; }
    // scales the rows with 'scale'

    void ColScale (const TVector<MT> &scale)
    { ERROR_UNDEF; }
    // scales the columns with 'scale'

    MT GetNext (int &r, int &c) const;

    void Ax (const TVector<MT> &x, TVector<MT> &b) const;

    void Ax (const TVector<MT> &x, TVector<MT> &b, int r1, int r2) const;

    void ATx (const TVector<MT> &x, TVector<MT> &b) const;

    int Shrink ();
    // Remove all entries with zero value. Return value is the number of
    // eliminated entries.

    void SymbolicCholeskyFactorize (idxtype *&frowptr, idxtype *&fcolidx) const;
    // Calculate the sparse fill-in structure of the lower triangle of
    // the Cholesky factorisation of the matrix (excluding diagonal)
    // and return in index lists `frowptr' and 'fcolidx'

    friend bool CholeskyFactorize<> (const TSymCompRowMatrix<MT> &A,
        TCompRowMatrix<MT> &L, TVector<MT> &d, bool recover);
    // Perform Cholesky factorisation of A and return the result in lower
    // triangle L and diagonal d (L does not contain diagonal elements)
    // L must have been initialised with the index lists returned from
    // A.CalculateCholeskyFillin()

    idxtype *rowptr;
    idxtype *colidx;

    mutable idxtype *colptr, *rowidx, *vofs;
    mutable bool col_access;
    // index arrays and flag for column access

    friend std::istream &operator>> <> (std::istream &is,
        TSymCompRowMatrix<MT> &m);
    friend std::ostream &operator<< <> (std::ostream &os,
        const TSymCompRowMatrix<MT> &m);
    // IO in generic format

private:
    int Get_index (int r, int c) const;
    // returns offset into data array of element at row r and column c
    // returns -1 if entry does not exist
    // This always returns -1 if r < c

    mutable int iterator_pos;

    bool allow_col_indexing;
};

// ==========================================================================
// friend prototypes

#ifdef NEED_FRIEND_PT

template<class MT>
bool CholeskyFactorize (const TSymCompRowMatrix<MT> &A, TCompRowMatrix<MT> &L,
    TVector<MT> &d, bool recover);

#endif // NEED_FRIEND_PT

// ==========================================================================
// typedefs for specific instances of `TCompRowMatrix'

typedef TSymCompRowMatrix<double>	RSymCompRowMatrix;	// 'real'
typedef TSymCompRowMatrix<float>	FSymCompRowMatrix;	// 'float'
typedef TSymCompRowMatrix<std::complex<double> > CSymCompRowMatrix;	// 'complex'
typedef TSymCompRowMatrix<std::complex<float> >	SCSymCompRowMatrix;	// 'complex'
typedef TSymCompRowMatrix<int>	        ISymCompRowMatrix;	// 'integer'

#endif // !__SCRMATRIX_H
