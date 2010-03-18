// -*-C++-*-
// ============================================================================
// General dense matrix class
// ============================================================================

#ifndef __DNSMATRIX_H
#define __DNSMATRIX_H

#include "toastdef.h"
#include <string.h>
#include "vector.h"
#include "matrix.h"

// ==========================================================================
// Nonmember declarations

template<class MT> class TSymMatrix;
template<class MT> class TDenseMatrix;

template<class MT>
MATHLIB TDenseMatrix<MT> transpose (const TDenseMatrix<MT> &A);

template<class MT>
MATHLIB TDenseMatrix<MT> cath (const TDenseMatrix<MT> &A,
    const TDenseMatrix<MT> &B);

template<class MT>
MATHLIB TDenseMatrix<MT> catv (const TDenseMatrix<MT> &A,
    const TDenseMatrix<MT> &B);

template<class MT>
MATHLIB TSymMatrix<MT> ATA (const TDenseMatrix<MT> &A);

template<class MT>
MATHLIB TSymMatrix<MT> AAT (const TDenseMatrix<MT> &A);

template<class MT>
MATHLIB TDenseMatrix<MT> kron (const TDenseMatrix<MT> &A, const TDenseMatrix<MT> &B);

template<class MT>
MATHLIB MT det (const TDenseMatrix<MT> &A, TDenseMatrix<MT> *Ai = 0);

template<class MT>
MATHLIB TDenseMatrix<MT> inverse (const TDenseMatrix<MT> &A);

template<class MT>
MATHLIB int QRFactorize (TDenseMatrix<MT> &A, TVector<MT> &c, TVector<MT> &d);

template<class MT>
MATHLIB void RSolve (const TDenseMatrix<MT> &A, const TVector<MT> &d,
    TVector<MT> &b);

template<class MT>
MATHLIB void QRSolve (const TDenseMatrix<MT> &A, const TVector<MT> &c,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x);

template<class MT>
MATHLIB void LUFactorize (TDenseMatrix<MT> &a, IVector &indx, double &d);

template<class MT>
MATHLIB void LUSolve (const TDenseMatrix<MT> &a, const IVector &indx, TVector<MT> &b);

template<class MT>
MATHLIB std::istream &operator>> (std::istream &is, TDenseMatrix<MT> &m);

template<class MT>
MATHLIB std::ostream &operator<< (std::ostream &os, const TDenseMatrix<MT> &m);

// ==========================================================================
// class TDenseMatrix
// ==========================================================================
/**
 * \brief Dense matrix class.
 *
 * A fully populated matrix type that stores its elements in row order.
 * Supports all operations inherited from the matrix base class interface,
 * such as row and column norms, row and column scaling, matrix x vector
 * operations, and various direct and iterative linear solvers.
 */

template<class MT> class TDenseMatrix: public TMatrix<MT> {

public:
    TDenseMatrix (): TMatrix<MT> ()
    { rc = 0; }
    // Create 0 x 0 matrix

    TDenseMatrix (int r, int c): TMatrix<MT> (r, c)
    { Alloc (r,c); Zero(); }
    // Create r x c matrix filled with 0's

    TDenseMatrix (int n): TMatrix<MT> (n, n)
    { Alloc (n,n); Zero(); }
    // Create a square n x n matrix filled with 0's

    TDenseMatrix (const TDenseMatrix<MT> &m): TMatrix<MT> (m)
    { Alloc (m.rows, m.cols); memcpy (val, m.val, rc*sizeof(MT)); }
    // Create copy of matrix m

    /**
     * \brief Construct a matrix by copying a sub-block of another matrix.
     * \param m source matrix
     * \param i0 low row index of sub-block
     * \param j0 low column index of sub-block
     * \param i1 high row index+1 of sub-block
     * \param j1 high column index+1 of sub-block
     * \note This constructor creates a dense matrix of dimensions
     *   i1-i0 x j1-j0 by copying a block of matrix m.
     * \note 0 <= i0 <= i1 <= m.nRows() and 0 <= j0 <= j1 <= m.nCols() is
     *   required.
     * \note This is not a reference constructor: the new matrix does not
     *   share its data block with the source matrix.
     */
    TDenseMatrix (const TDenseMatrix<MT> &m, int i0, int j0, int i1, int j1);

    TDenseMatrix (int r, int c, MT *values): TMatrix<MT> (r, c)
    { Alloc (r,c); memcpy (val, values, rc*sizeof(MT)); }
    // Create r x c matrix and initialise from 'values' array, which
    // contains row vectors and must be at least of dimension r*c

    TDenseMatrix (int r, int c, const char *valstr);
    // Create r x c matrix and initialise from a string which contains
    // at least r*c MT values in row order

    TDenseMatrix (const TSymMatrix<MT> &A);
    // Create a dense matrix from a symmetric matrix

    ~TDenseMatrix () { Unlink(); }
    // Destructor

    inline MatrixStorage StorageType () const { return MATRIX_DENSE; }
    // Matrix element storage method

    inline void New (int r, int c) { New_dirty(r,c); Zero(); }
    // Resize matrix and zero all elements
    // obsolete, use Zero(r,c) instead

    inline void Zero () { memset (val, 0, rc*sizeof(MT)); }
    // Zero all elements

    inline void Zero (int n) { New_dirty(n,n); Zero(); }
    // Resize to n x n square matrix and zero all elements

    inline void Zero (int r, int c) { New_dirty(r,c); Zero(); }
    // Resize to r x c matrix and zero all elements

    void Identity ();
    // Set diagonal elements to 1, all others to 0
    // only valid for template types which can cast (MT)1

    inline void Identity (int n) { New_dirty(n,n); Identity(); }
    // Resize to n x n square matrix and set to identity

    inline void Identity (int r, int c) { New_dirty(r,c); Identity(); }
    // Resize to r x c matrix and set to identity

    inline MT Get (int r, int c) const {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols, Index out of range);
	return val[r*this->cols + c];
    }
    // Retrieve value of an element

    inline void Set (int r, int c, MT v) {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols, Index out of range);
	val[r*this->cols + c] = v;
    }
    // Write value into element

    inline MT operator() (int r, int c) const { 
      return Get(r,c); }
    // Retrieve value of an element (alternative syntax to 'Get')

    inline MT &operator() (int r, int c) {
        dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols, Index out of range);
	return val[r*this->cols + c];
    }
    // Retrieve reference to element

    inline MT *valptr() { return val; }
    // Return pointer to data array

    inline TVector<MT> Row (int r) const {
        dASSERT (r >= 0 && r < this->rows, Index out of range);
	return TVector<MT>(this->cols, val+r*this->cols);
    }
    // Retrieve a row

    TVector<MT> Col (int c) const;
    // Retrieve a column

    int SparseRow (int r, int *ci, MT *rv) const;

    void SetRow (int r, const TVector<MT> &rval);
    // replace row 'r' of *this with 'rval'

    void ColScale (const TVector<MT> &scale);
    // scales the columns with 'scale'

    void RowScale (const TVector<MT> &scale);
    // scales the rows with 'scale'

    /**
     * \brief Returns a vector of row sums.
     * \return Row sum vector.
     * \note The size of the returned vector equals the number of rows of
     *   the matrix. Each entry i contains the sum of elements of row i.
     * \sa RowSumSq, ColSum, ColSumSq
     */
    TVector<MT> RowSum() const;

    /**
     * \brief Returns a vector of sum of squares for each row.
     * \return Vector of sum of squares.
     * \note The size of the returned vector equals the number of rows of
     *   the matrix. Each entry i contains the sum of the squares of elements
     *   of row i.
     * \sa RowSum, ColSum, ColSumSq
     */
    TVector<MT> RowSumSq() const;

    /**
     * \brief Returns a vector of column sums.
     * \return Column sum vector.
     * \note The size of the returned vector equals the number of columns
     *   of the matrix. Each entry j contains the sum of elements of column j.
     * \sa RowSum, RowSumSq, ColSumSq
     */
    TVector<MT> ColSum() const;

    /**
     * \brief Returns a vector of sum of squares for each column.
     * \return Vector of sum of squares.
     * \note The size of the returned vector equals the number of columns
     *   of the matrix. Each entry j contains the sum of the squares of
     *   elements of column j.
     * \sa ColSum, RowSum, RowSumSq
     */
    TVector<MT> ColSumSq() const;

    void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    // Matrix x Vector multiplication: b = Ax

    void ATx (const TVector<MT> &x, TVector<MT> &b) const;
    // b = transpose(A) * x

    void AB (const TDenseMatrix<MT> &A, const TDenseMatrix<MT> &B);
    // *this = AB
      
    MT colmult (int c1, int c2);
    // inner product of columns c1 and c2

    TDenseMatrix<MT> &operator= (const TDenseMatrix<MT> &m);
    // *this = m

    TDenseMatrix<MT> &operator= (MT v);
    // *this = v

    TDenseMatrix<MT> operator+ (const TDenseMatrix<MT> &m) const;
    // *this + m

    TDenseMatrix<MT> operator- (const TDenseMatrix<MT> &m) const;
    // *this - m

    TDenseMatrix<MT> operator* (const TDenseMatrix<MT> &m) const
    { TDenseMatrix<MT> tmp; tmp.AB(*this,m); return tmp; }
    // *this * m

    inline TVector<MT> operator* (const TVector<MT> &x) const
    { TVector<MT> b(this->rows); Ax (x, b); return b; }
    // *this * x

    inline TDenseMatrix<MT> operator* (const MT &mt) const {
        TDenseMatrix<MT> tmp(*this);
	for (int i = 0; i < rc; i++) tmp.val[i] *= mt;
	return tmp;
    }
    // scalar matrix multiplication

    inline TDenseMatrix<MT> operator/ (const MT &mt) const {
        TDenseMatrix<MT> tmp(*this);
	for (int i = 0; i < rc; i++) tmp.val[i] /= mt;
	return tmp;
    }
    // scalar matrix division

    inline TDenseMatrix<MT> &operator*= (const MT &mt) {
        for (int i = 0; i < rc; i++) val[i] *= mt;
	return *this;
    }
    // scalar matrix multiplication on *this

    inline TDenseMatrix<MT> &operator/= (const MT &mt) {
        return *this *= (MT)1 / mt;
    }
    // scalar matrix division on *this

    friend MATHLIB TDenseMatrix<MT> transpose<> (const TDenseMatrix<MT> &A);
    // returns transpose(A)

    /**
     * \brief Concatenate two matrices horizontally
     * \param A first matrix argument
     * \param B second matrix argument
     * \return Result of concatenation [A B]
     * \note A and B must have the same number of rows.
     * \note This method has the functionality of the MATLAB construct
     *   C = (A,B)
     */
    friend MATHLIB TDenseMatrix<MT> cath<> (const TDenseMatrix<MT> &A,
        const TDenseMatrix<MT> &B);

    /**
     * \brief Concatenate two matrices vertically
     * \param A first matrix argument
     * \param B second matrix argument
     * \return Result of concatenation
     *    \f[ C = \left[ \begin{array}{c} A\\B \end{array} \right] \f]
     * \note A and B must have the same number of columns
     * \note This method has the functionality of the MATLAB construct
     *   C = (A;B)
     */
    friend MATHLIB TDenseMatrix<MT> catv<> (const TDenseMatrix<MT> &A,
        const TDenseMatrix<MT> &B);

    friend MATHLIB TSymMatrix<MT> ATA<> (const TDenseMatrix<MT> &A);
    // returns transpose(A) * A

    friend MATHLIB TSymMatrix<MT> AAT<> (const TDenseMatrix<MT> &A);
    // returns A * transpose(A)

    friend MATHLIB TDenseMatrix<MT> kron<> (const TDenseMatrix<MT> &A,
        const TDenseMatrix<MT> &B);
    // Kronecker matrix product

    friend MATHLIB MT det<> (const TDenseMatrix<MT> &A, TDenseMatrix<MT> *Ai);
    // Returns the determinant of A. If Ai != 0 and order <= 3x3 then it
    // is set to the inverse of A; this is more efficient than a separate
    // call to 'inverse'

    friend MATHLIB TDenseMatrix<MT> inverse<> (const TDenseMatrix<MT> &A);
    // uses LU decomposition for matrices > 3x3

    friend MATHLIB int QRFactorize<> (TDenseMatrix<MT> &A, TVector<MT> &c,
        TVector<MT> &d);
    // QR decomposition. Return value != 0 indicates singularity

    friend MATHLIB void RSolve<> (const TDenseMatrix<MT> &A, const TVector<MT> &d,
	TVector<MT> &b);
    // Solves set of linear equations Rx=b where R is upper triangular matrix
    // stored in A and d. On return b is overwritten by the result x.

    friend MATHLIB void QRSolve<> (const TDenseMatrix<MT> &A, const TVector<MT> &c,
        const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x);
    // Solves set of linear equations Ax=b, where A, c and d are the result of
    // a preceeding call to QRFactorize

    friend MATHLIB void LUFactorize<MT> (TDenseMatrix<MT> &a, IVector &indx,
        double &d);
    // Replaces a with its LU factorisation. On exit, indx contains the
    // permutation vector from partial pivoting, and d = +/-1 depending on
    // even/odd number of interchanges

    friend MATHLIB void LUSolve<MT> (const TDenseMatrix<MT> &a,
        const IVector &indx, TVector<MT> &b);
  
    MT *data_buffer() { return val; }
    const MT *data_buffer() const { return val; }

    friend MATHLIB std::istream &operator>> <> (std::istream &is, TDenseMatrix<MT> &m);
    friend MATHLIB std::ostream &operator<< <> (std::ostream &os,
        const TDenseMatrix<MT> &m);
    // IO in generic format

protected:
    MT *val; // data array (row vectors)
    int rc;  // rows*cols = number of entries in data array

private:
    inline void Alloc (int r, int c) { if (rc = r*c) val = new MT[rc]; }
    // allocate data array of length r*c

    inline void Unlink () { if (rc) delete []val, rc = 0; }
    // deallocate data array

    void New_dirty (int r, int c);
    // This version of New leaves the elements undefined and is slightly
    // faster than 'New'. It should only be used where initialisation
    // follows immediately
};

// ==========================================================================
// typedefs for specific instances of `TDenseMatrix'

typedef TDenseMatrix<double>   RDenseMatrix;	// 'double real'
typedef TDenseMatrix<float>    FDenseMatrix;    // 'single real'
typedef TDenseMatrix<toast::complex>  CDenseMatrix;	// 'complex'
typedef TDenseMatrix<scomplex> SCDenseMatrix;   // 'single complex'
typedef TDenseMatrix<int>      IDenseMatrix;    // 'integer'

#ifndef MATHLIB_IMPLEMENTATION
extern template class MATHLIB TDenseMatrix<double>;
extern template class MATHLIB TDenseMatrix<float>;
extern template class MATHLIB TDenseMatrix<toast::complex>;
extern template class MATHLIB TDenseMatrix<scomplex>;
extern template class MATHLIB TDenseMatrix<int>;
#endif // MATHLIB_IMPLEMENTATION

#endif // !__DNSMATRIX_H