// -*-C++-*-
// ==========================================================================
// Module mathlib 				   Martin Schweiger - 30.5.96
// File bsmatrix.h
// Declaration of template class TBandSymMatrix ('template banded symmetric
// matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RBandSymMatrix = TBandSymMatrix<double>		('real')
//	FBandSymMatrix = TBandSymMatrix<float>		('float')
//	CBandSymMatrix = TBandSymMatrix<complex>	('complex')
//	IBandSymMatrix = TBandSymMatrix<int>		('integer')
//	BandSymMatrix  = TBandSymMatrix<double>		backward compatibility
//
// Inheritance:
// ------------
// TRootMatrix ----> TBandSymMatrix
//
// Storage convention:
// -------------------
// full matrix:	a b c 0 0 0	TBandSymMatrix:	0 0 a
//		b a b c 0 0			0 b a
//		c b a b c 0			c b a
//		0 c b a b c			c b a
//		0 0 c b a b			c b a
//		0 0 0 c b a			c b a
// ==========================================================================

#ifndef __BSMATRIX_H
#define __BSMATRIX_H

// ==========================================================================
// class TBandSymMatrix

template<class MT> class TBandSymMatrix: public TRootMatrix<MT> {
public:
    TBandSymMatrix ();
    TBandSymMatrix (int rc, int hb);
    TBandSymMatrix (const TBandSymMatrix<MT> &A);
    virtual ~TBandSymMatrix ();

    TVector<MT> &operator[] (int i) const;
    // line-vector reference

    MT Get (int r, int c) const;
    TVector<MT> Row (int r) const;
    TVector<MT> Col (int c) const;
    void New (int rc, int hb);
    void Copy (const TBandSymMatrix &mat);

    void Unlink ();
    // unlinks the matrix from its data block

    bool PIndex (int i, int j, int &r, int &c);
    // converts logical matrix indices i,j to physical indices r,c
    // returns TRUE if physical mapping is found, FALSE otherwise

    void LIndex (int r, int c, int &i, int &j);
    // converts physical matrix position r,c to logical position i,j

    TVector<MT> operator* (const TVector<MT> &x) const;
    // matrix * vector operation

    virtual void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    // alternative function to operator*; returns the result in parameter b
    // and avoids local vector creation and copy

    MT LineVecMul (int line, const TVector<MT> &x) const;
    // inner product of line `line' of the BandSymMatrix and vector `x'

    // **** friends ****

    friend TSquareMatrix<MT> ToSquare (TBandSymMatrix<MT> &A);
    // converts the BandSymMatrix into a SquareMatrix

    friend TSymMatrix<MT> ToSym (TBandSymMatrix<MT> &A);
    // converts the BandSymMatrix into a SymMatrix

    friend bool CHdecomp (TBandSymMatrix<MT> &a, bool recover = FALSE);
    // Replaces `a' with its Cholesky decomposition.
    // If recover=TRUE then the function will return even if the matrix is
    // not positive definite, otherwise the error handler is invoked.
    // If the function returns, the return value indicates success (i.e.
    // positive-definiteness). A function which calls CHdecomp with
    // recover=TRUE is responsible for handling failures.

    friend TVector<MT> CHsubst (const TBandSymMatrix<MT> &a,
	const TVector<MT> &b);
    // Cholesky backward and forward substitution

    friend TVector<MT> PCHsubst (const TBandSymMatrix<MT> &a,
	const TVector<MT> &b, const TVector<MT> &p);
    // preconditioned Cholesky substitution

protected:
    int hband;
    void Allocate (int rc, int hband);
    TVector<MT> *data;	// data are stored in row vectors over the band
};


// ==========================================================================
// inline functions

template<class MT>
inline TBandSymMatrix<MT>::TBandSymMatrix ()
: TRootMatrix<MT> ()
{
    data = 0;
    hband = 0;
}

template<class MT>
inline TBandSymMatrix<MT>::TBandSymMatrix (int rc, int hb)
: TRootMatrix<MT> (rc, hband)
{
    data = 0;
    Allocate (rc, hb);
}

template<class MT>
inline TBandSymMatrix<MT>::TBandSymMatrix (const TBandSymMatrix<MT> &A)
{
    data = 0;
    Allocate (A.rows, A.hband);
    Copy (A);
}

template<class MT>
inline TBandSymMatrix<MT>::~TBandSymMatrix ()
{
    Unlink ();
}

// Unlink current data block and allocate a new one
template<class MT>
inline void TBandSymMatrix<MT>::New (int rc, int hb)
{
    Unlink ();
    Allocate (rc, hb);
}

#ifndef MATH_DEBUG
template<class MT>
inline TVector<MT>& TBandSymMatrix<MT>::operator[] (int i) const
{
    return data[i];
}
#endif // !MATH_DEBUG

#ifndef MATH_DEBUG
template<class MT>
inline void TBandSymMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    int r, c, k;
    MT br;

    for (r = 0; r < rows; r++) {
	TVector<MT> &lr = data[r];
	for (c = hband-1, k = r, br = 0; c >= 0 && k >= 0; c--, k--)
	    br += lr[c] * x[k];
	for (c = hband-2, k = r+1; c >= 0 && k < rows; c--, k++)
	    br += data[k][c] * x[k];
	b[r] = br;
    }
}
#endif // !MATH_DEBUG

// ==========================================================================
// friend prototypes

#ifdef NEED_FRIEND_PT

template<class MT>
TSquareMatrix<MT> ToSquare (TBandSymMatrix<MT> &A);

template<class MT>
TSymMatrix<MT> ToSym (TBandSymMatrix<MT> &A);

template<class MT>
bool CHdecomp (TBandSymMatrix<MT> &A, bool recover = FALSE);

template<class MT>
TVector<MT> CHsubst (const TBandSymMatrix<MT> &A, const TVector<MT> &b);

template<class MT>
TVector<MT> PCHsubst (const TBandSymMatrix<MT> &A, const TVector<MT> &b,
    const TVector<MT> &p);

#endif

// ==========================================================================
// typedefs for specific instances of `TBandSymMatrix'

typedef TBandSymMatrix<double>  RBandSymMatrix;	// 'real'
typedef TBandSymMatrix<float>   FBandSymMatrix;	// 'float'
typedef TBandSymMatrix<complex> CBandSymMatrix;	// 'complex'
typedef TBandSymMatrix<int>	IBandSymMatrix;	// 'integer'
typedef TBandSymMatrix<double>  BandSymMatrix;	// for backward compatibility

#endif // !__BSMATRIX_H
