// -*-C++-*-
// ==========================================================================
// Module mathlib 				   Martin Schweiger - 23.5.96
// File symatrix.h
// Declaration of template class TSymMatrix ('template symmetric matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RSymMatrix = TSymMatrix<double>		('real')
//	FSymMatrix = TSymMatrix<float>		('float')
//	CSymMatrix = TSymMatrix<complex>	('complex')
//	ISymMatrix = TSymMatrix<int>		('integer')
//	SymMatrix  = TSymMatrix<double>		for backward compatibility
//
// Inheritance:
// ------------
// TMatrix ----> TSymMatrix
//
// Note:
// -----
// TSymMatrix is derived from TMatrix rather than TSquareMatrix because
// the data structure is different, so inheriting functionality from the
// TMatrix branch would be a problem.
// TSymMatrix only stores the lower left triangle of the matrix.
// ==========================================================================

#ifndef __SYMATRIX_H
#define __SYMATRIX_H

#include "matrix.h"

// ==========================================================================
// Nonmember declarations

template<class MT> class TSymMatrix;

template<class MT>
bool CHdecomp (TSymMatrix<MT> &a, bool recover = false);

template<class MT>
TVector<MT> CHsubst (const TSymMatrix<MT> &a, const TVector<MT> &b);

// ==========================================================================
// class TSymMatrix

template<class MT> class TSymMatrix: public TMatrix<MT> {
public:
    TSymMatrix (): TMatrix<MT> ()
    { nz = 0; }
    // Create a 0 x 0 matrix

    TSymMatrix (int r, int c): TMatrix<MT> (r, c) {
        dASSERT(r==c, "Symmetric matrix must be square");
	Alloc(r); Zero();
    }
    // Create r x c matrix filled with 0's
    // r == c is required

    TSymMatrix (int n): TMatrix<MT> (n, n)
    { Alloc(n); Zero(); }
    // Create n x n matrix filled with 0's

    TSymMatrix (const TSymMatrix<MT> &m): TMatrix<MT> (m)
    { Alloc (m.rows); memcpy (val, m.val, nz*sizeof(MT)); }
    // Create a copy of m

    TSymMatrix (int n, const char *valstr);
    // Create n x n symmetric matrix and initialise from string which contains
    // values for lower triangle in row order

    ~TSymMatrix () { Unlink (); }
    // Destructor

    MatrixStorage StorageType() const { return MATRIX_SYM; }
    // Matrix element storage method

    inline void New (int r, int c) {
        dASSERT(r==c, "Symmetric matrix must be square");
	New_dirty(r); Zero();
    }
    // Resize matrix and zero all elements
    // obsolete, use Zero(r,c) instead

    inline void Zero () { memset (val, 0, nz*sizeof(MT)); }
    // Zero all elements

    inline void Zero (int n) { New_dirty(n); Zero(); }
    // Resize to n x n square matrix and zero all elements

    inline void Zero (int r, int c) {
        dASSERT(r==c, "Symmetric matrix must be square");
	New_dirty(r); Zero();
    }
    // Resize to r x c matrix and zero all elements. r==c is required

    void Identity ();
    // Set diagonal elements to 1, all others to 0
    // only valid for template types which can cast (MT)1

    inline void Identity (int n) { New_dirty(n); Identity(); }
    // Resize to n x n square matrix and set to identity

    inline void Identity (int r, int c) {
        dASSERT(r==c, "Symmetric matrix must be square");
	New_dirty(r); Identity();
    }
    // Resize to r x c matrix and set to identity. r==c is required

    inline MT Get (int r, int c) const {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols,
		 "Index out of range");
	return val[Idx(r,c)];
    }
    // Retrieve value of an element

    inline const MT operator() (int r, int c) const {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols,
		 "Index out of range");
	return val[Idx(r,c)];
    }
    // Retrieve value of an element (alternative syntax to 'Get')

    inline MT &operator() (int r, int c) {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols,
		 "Index out of range");
	return val[Idx(r,c)];
    }
    // Retrieve reference to element

    inline TSymMatrix<MT> &operator*= (MT f)
    { for (int i = 0; i < nz; i++) val[i] *= f; return *this; }

    inline TSymMatrix<MT> &operator/= (MT f)
    { for (int i = 0; i < nz; i++) val[i] /= f; return *this; }

    TVector<MT> Row (int r) const;
    // Retrieve a row

    inline TVector<MT> Col (int c) const { return Row(c); }
    // Retrieve a column. Identical to Row(c) since symmetric

    int SparseRow (int r, idxtype *ci, MT *rv) const;
    // See TMatrix

    void ColScale (const TVector<MT> &scale);
    // scales the columns with 'scale'

    void RowScale (const TVector<MT> &scale);
    // scales the rows with 'scale'

    void AddDiag (const TVector<MT> &d);
    // adds vector d to diagonal of *this

    void AddDiag (const MT &d);
    // adds 'd' to all diagonal elements

    void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    // Matrix x Vector multiplication: b = Ax

    inline void ATx (const TVector<MT> &x, TVector<MT> &b) const
    { return Ax(x,b); }
    // b = transpose(A) * x = Ax (since symmetric)

    TSymMatrix<MT> &operator= (const TSymMatrix<MT> &m);
    // *this = m

    TSymMatrix<MT> &operator= (const MT &mt);
    // *this[i] = mt for all i

    TSymMatrix<MT> operator+ (const TSymMatrix<MT> &m) const;
    // *this + m

    TSymMatrix<MT> operator- (const TSymMatrix<MT> &m) const;
    // *this - m

    TSymMatrix<MT> &operator+= (const TSymMatrix<MT> &m) {
        dASSERT(this->rows == m.rows, "Matrices have different size");
	for (int i = 0; i < nz; i++) val[i] += m.val[i];
	return *this;
    }
    // *this = *this + m

    TSymMatrix<MT> &operator-= (const TSymMatrix<MT> &m) {
        dASSERT(this->rows == m.rows, "Matrices have different size");
	for (int i = 0; i < nz; i++) val[i] -= m.val[i];
	return *this;
    }
    // *this = *this - m

    TSymMatrix<MT> operator* (const MT &mt) const;
    // *this * mt

    TVector<MT> operator* (const TVector<MT> &x) const;
    // matrix * vector operation

    // **** friends ****

    friend bool CHdecomp<> (TSymMatrix<MT> &a, bool recover);
    // replaces `a' with its Choleski decomposition

    friend TVector<MT> CHsubst<> (const TSymMatrix<MT> &a,
        const TVector<MT> &b);

    MT *data_buffer() { return val; }
    const MT *data_buffer() const { return val; }

protected:
    MT *val; // data array (row vectors for lower triangle)
    int nz;  // length of data vector: nz = n*(n+1)/2

    inline int Idx (int r, int c) const
    { return (r >= c ? c + ((r*(r+1)) >> 1) : r + ((c*(c+1)) >> 1)); }
    // data array index for element (r,c)

private:
    inline void Alloc (int n) { if (n) val = new MT[nz = (n*(n+1))/2]; }
    // allocate data array

    inline void Unlink () { if (nz) delete[]val, nz = 0; }
    // deallocate data array

    void New_dirty (int n);
    // This version of new leaves the elements undefined and is slightly
    // faster than 'New'. It should only be used where initialisation
    // follows immediately
};


// ==========================================================================
// typedefs for specific instances of `TSymMatrix'

typedef TSymMatrix<double>   RSymMatrix;	// 'real'
typedef TSymMatrix<float>    FSymMatrix;	// 'float'
typedef TSymMatrix<std::complex<double> >  CSymMatrix;	// 'complex'
typedef TSymMatrix<std::complex<float> > SCSymMatrix;	// 'single complex'
typedef TSymMatrix<int>	     ISymMatrix;	// 'integer'

// ==========================================================================
// extern declarations of TSymMatrix (only required for VS)

#ifndef __SYMATRIX_CC
//extern template class MATHLIB TSymMatrix<double>;
//extern template class MATHLIB TSymMatrix<float>;
//extern template class MATHLIB TSymMatrix<toast::complex>;
//extern template class MATHLIB TSymMatrix<scomplex>;
//extern template class MATHLIB TSymMatrix<int>;
#endif // !__SYMATRIX_CC

#endif // !__SYMATRIX_H
