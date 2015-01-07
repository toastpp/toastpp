// -*-C++-*-
// ==========================================================================
// Module mathlib
// File dgmatrix.h
// Declaration of template class TDiagMatrix ('template diagonal matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RDiagMatrix = TDiagMatrix<double>	('real')
//	FDiagMatrix = TDiagMatrix<float>	('float')
//	CDiagMatrix = TDiagMatrix<complex>	('complex')
//	IDiagMatrix = TDiagMatrix<int>		('integer')
//	DiagMatrix  = TDiagMatrix<double>	for backward compatibility
//
// Inheritance:
// ------------
// TGenericSparseMatrix ----> TDiagMatrix
// ==========================================================================

#ifndef __DGMATRIX_H
#define __DGMATRIX_H

#include "gsmatrix.h"
#include "dnsmatrix.h"

// ==========================================================================
// class TDiagMatrix
// ==========================================================================
/**
 * \brief Diagonal matrix class
 *
 * A sparse matrix class which has nonzeros only on the diagonal.
 * Also supports non-square matrices.
 * Diagonal elements are always allocated, even if they are zero.
 * The following template types have been instantiated:
 * <table rows=2>
 * <tr><td>TDiagMatrix<double></td><td>RDiagMatrix</td></tr>
 * <tr><td>TDiagMatrix<float></td><td>FDiagMatrix</td></tr>
 * <tr><td>TDiagMatrix<toast::complex></td><td>CDiagMatrix</td></tr>
 * <tr><td>TDiagMatrix<int></td><td>IDiagMatrix</td></tr>
 * </table>
 */
template<class MT> class TDiagMatrix: public TGenericSparseMatrix<MT> {
public:
    /**
     * \brief Creates a diagonal matrix of dimension 0 x 0.
     */
    TDiagMatrix ();

    /**
     * \brief Creates a diagonal matrix of dimension r x c.
     * \param r number of rows
     * \param c number of columns
     * \param v diagonal value (default: zero)
     */
    TDiagMatrix (int r, int c, const MT v = (MT)0);

    /**
     * \brief Creates a diagonal matrix as a copy of another.
     * \param mat source matrix
     */
    TDiagMatrix (const TDiagMatrix<MT> &mat);

    /**
     * \brief Matrix destructor.
     */
    virtual ~TDiagMatrix ();

    /**
     * \brief Returns the matrix storage type.
     * \return MATRIX_DIAG
     */
    MatrixStorage StorageType() const { return MATRIX_DIAG; }

    /**
     * \brief Resizes and resets the matrix.
     * \param nrows new number of rows
     * \param ncols new number of columns
     * \note All diagonal elements are reset to zero.
     */
    void New (int nrows, int ncols);

    /**
     * \brief Retrieves a matrix element.
     * \param r matrix row (0 <= r < nRows())
     * \param c matrix column (0 <= c < nCols())
     * \return Matrix element (*this)<sub>r,c</sub>
     * \note For any off-diagonal elements, this returns zero.
     */
    virtual MT Get (int r, int c) const
	{ return (r == c) ? this->val[r] : (MT)0; }

    /**
     * \brief Returns a vector containing a copy of row 'r'.
     * \param r row index (0 <= r < nRows())
     * \return Vector containing row r (dimension nCols())
     * \note At most, element r of the returned vector (if it exists) is
     *    non-zero.
     */
    TVector<MT> Row (int r) const;

    /**
     * \brief Returns a vector containing a copy of column 'c'.
     * \param c column index (0 <= c < nCols())
     * \return Vector containing column c (dimension nRows())
     * \note At most, element c of the returned vector (if it exists) is
     *   non-zero.
     */
    TVector<MT> Col (int c) const;

    /**
     * \brief Returns a row of the matrix in sparse format.
     * \param r row index (0 <= r < nRows())
     * \param colidx column index array
     * \param val element value array
     * \return Number of allocated matrix entries
     * \note The index and value arrays must be allocated by the caller and be
     *   of size 1 or larger.
     * \note This function will return at most one element (the diagonal
     *   element of row r, if it exists)
     * \note colidx and val can therefore be pointers instead of arrays.
     */
    int SparseRow (int r, idxtype *colidx, MT *val) const;

    /**
     * \brief Scales the matrix columns.
     * \param scale scaling vector (dimension nCols())
     * \note Multiplies each column c of the matrix with scale[c]
     */
    void ColScale (const TVector<MT> &scale);

    /**
     * \brief Scales the matrix rows.
     * \param scale scaling vector (dimensioin nRows())
     * \note Multiplies each row r of the matrix with scale[r]
     */
    void RowScale (const TVector<MT> &scale);

    /**
     * \brief Matrix assignment
     * \param mat right-hand assignment argument
     * \return The new matrix (a copy of mat)
     * \note Sets *this = mat; resizes the matrix if required.
     */
    TDiagMatrix &operator= (const TDiagMatrix<MT> &mat);

    /**
     * \brief Constant diagonal assignment
     * \param v new diagonal value
     * \return The new matrix
     * \note Sets all diagonal elements of *this to v. Matrix dimension does
     *   not change.
     */
    TDiagMatrix &operator= (const MT &v);

    /**
     * \brief Matrix copy operation.
     * \param mat source matrix
     * \note Makes *this a copy of mat. Resizes the matrix if required.
     */
    void Copy (const TDiagMatrix<MT> &mat);

    /**
     * \brief Cast to full matrix.
     * \note Creates a full matrix from the diagonal matrix with zero
     *   off-diagonal elements.
     */
    operator TDenseMatrix<MT> ();

    /**
     * \brief Element access.
     * \param r row index (0 <= r < nRows())
     * \param c column index (0 <= c < nCols())
     * \return matrix element (r,c)
     * \note For all off-diagonal elements (r != c) this returns a reference
     *   to a static zero value.
     */
    MT &operator() (int r, int c);

    /**
     * \brief Matrix addition.
     * \param mat matrix operand
     * \return *this + mat
     * \note *this is unchanged. *this and mat must have the same dimensions.
     */
    TDiagMatrix operator+ (const TDiagMatrix &mat) const;

    /**
     * \brief Matrix subtraction.
     * \param mat matrix operand
     * \return *this - mat
     * \note *this is unchanged. *this and mat must have the same dimensions.
     */
    TDiagMatrix operator- (const TDiagMatrix &mat) const;

    /**
     * \brief Checks allocation of a matrix element.
     * \param r row index (0 <= r < nRows())
     * \param c column index (0 <= c < nCols())
     * \returns \e true for diagonal elements (r==c), \e false otherwise.
     */
    bool Exists (int r, int c) const;

    /**
     * \brief Returns data array index for an element.
     * \param r row index (0 <= r < nRows())
     * \param c column index (0 <= c < nCols())
     * \return index of element (r,c) in data array, or -1 if not allocated.
     */
    int Get_index (int r, int c) const;

    /**
     * \brief Element iterator.
     * \param r row index
     * \param c column index
     * \return matrix element
     * \note This function returns the nonzero elements of the matrix in
     *   sequence, together with their row and column index.
     * \note In the first call, set r to a value < 0 to retrieve the first
     *   nonzero element. On exit, r and c are set to the element's row and
     *   column indices.
     * \note In subsequent calls, set r and c to the results of the previous
     *   call to get the next element.
     * \note A return value r < 0 indicates the end of the nonzero elements.
     * \note Returns the diagonal elements in order, starting from (0,0).
     */
    MT GetNext (int &r, int &c) const;

    /**
     * \brief Matrix-vector multiplication
     * \param x vector operand (dimension nCols())
     * \param b result of the multiplication.
     * \note Returns *this * x in b.
     * \note b is resized if required.
     */
    void Ax (const TVector<MT> &x, TVector<MT> &b) const;

    /**
     * \brief Partial matrix-vector multiplication
     * \param x vector operand (dimension nCols())
     * \param b result of the multiplication
     * \param r1 first row to be processed
     * \param r2 last row+1 to be processed
     * \note For all rows r with r1 <= r < r2 of the matrix, this sets b[r]
     *   to the result of the inner product *this[r] * x.
     */
    void Ax (const TVector<MT> &x, TVector<MT> &b, int r1, int r2) const;

    /**
     * \brief Transpose matrix-vector multiplication
     * \param x vector operand (dimension nRows())
     * \param b result of the multiplication
     * \note Returns (*this)<sup>T</sup> * x in b.
     * \note b is resized if required.
     */
    void ATx (const TVector<MT> &x, TVector<MT> &b) const;
};

// ==========================================================================
// typedefs for specific instances of `TDiagMatrix'
// ==========================================================================
typedef TDiagMatrix<double>	RDiagMatrix;	// 'real'
typedef TDiagMatrix<float>	FDiagMatrix;	// 'float'
typedef TDiagMatrix<std::complex<double> >	CDiagMatrix;	// 'complex'
typedef TDiagMatrix<int>	IDiagMatrix;	// 'integer'
typedef TDiagMatrix<double>	DiagMatrix;	// for backward compatibility
// ==========================================================================

#endif // !__DGMATRIX_H
