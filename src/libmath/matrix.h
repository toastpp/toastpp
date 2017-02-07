// -*-C++-*-
// ==========================================================================
// Module mathlib
// File matrix.h
// Declaration of template class TMatrix
// This is the root for all matrix classes and defines generic matrix methods
//
// Inheritance:
// ------------
// TMatrix ----> ...
// ==========================================================================

#ifndef __MATRIX_H
#define __MATRIX_H

#include "vector.h"

/// \defgroup matrix_storage matrix storage type
//@{
enum MatrixStorage {
    MATRIX_DENSE,      ///< dense matrix storage
    MATRIX_SYM,        ///< store lower triangle + diagonal of symmetric matrix
    MATRIX_DIAG,       ///< store diagonal matrix as vector
    MATRIX_COORD,      ///< coordinate storage (sparse)
    MATRIX_COMPROW,    ///< compressed row storage (sparse)
    MATRIX_COMPCOL,    ///< compressed column storage (sparse)
    MATRIX_SYMCOMPROW  ///< symmetric compressed row storage (sparse)
};
//@}

template<class MT> class TSymMatrix;
template<class MT> class TPreconditioner;

struct IterativeSolverResult {
    int it_count;
    double rel_err;
};

// ==========================================================================
// Nonmember declarations
// ==========================================================================

template<class MT> class TMatrix;

template<class MT>
TVector<MT> Ax (const TMatrix<MT> &A, const TVector<MT> &x)
{ TVector<MT> b; A.Ax(x,b); return b; }
// Returns the product A x

template<class MT>
TVector<MT> ATx (const TMatrix<MT> &A, const TVector<MT> &x)
{ TVector<MT> b; A.ATx(x,b); return b; }
// Returns the product A^T x

template<class MT>
TSymMatrix<MT> ATA (const TMatrix<MT> &A);
// Returns symmetric matrix A^T A (dimension A.Cols() x A.Cols())

template<class MT>
TSymMatrix<MT> AAT (const TMatrix<MT> &A);
// Returns symmetric matrix A A^T (dimension A.Rows() x A.Rows())

template<class MT>
TVector<MT> ATA_diag (const TMatrix<MT> &A);
// Returns the diagonal of A^T A as a vector

template<class MT>
int PCG (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0);
// PCG (preconditioned conjugate gradient) linear solver
// Solves Ax = b for given tolerance and preconditioner

template<class MT>
void PCG (const TMatrix<MT> &A, const TVector<MT> *b, TVector<MT> *x,
   int nrhs, double tol, int maxit = 0, TPreconditioner<MT> *precon = 0,
   IterativeSolverResult *res = 0);
// PCG solver for multiple right-hand sides

template<class MT>
int PCG (TVector<MT> (*Mv_clbk)(const TVector<MT> &v,
    void *context), void *context, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon = 0, 
    int maxit = 0);

template<class MT>
int BiPCG (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0);
// BiPCG (preconditioned bi-conjugate gradient) linear solver
// Solves Ax = b for given tolerance and preconditioner

template<class MT>
int BiCGSTAB (const TMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
void BiCGSTAB (const TMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x, int nrhs, double tol, int maxit = 0,
    TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0);
// BiCGSTAB solver for multiple right-hand sides

template<class MT>
int BiCGSTAB (TVector<MT> (*Mv_clbk)(const TVector<MT> &v,
    void *context), void *context, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon = 0, 
    int maxit = 0);

template<class MT>
int BiCGSTAB (void (*MVM)(TVector<MT> &),  const TVector<MT> &b, 
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon  = 0, 
    int maxit = 0);

template<class MT>
int GMRES (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon = 0, int restart = 10,
    int maxit = 0, void (*clbk)(void*) = 0);
// GMRES (generalised minimum residual) linear solver
// Solves Ax = b for given tolerance and preconditioner

/**
 * \brief GMRES iterative solver for linear system Ax=b: operator version.
 * \param Mv_clbk Pointer to callback function that implements operator Av.
 *   for any given vector v. The callback function must have the following format:
 *   TVector<MT> Av (const TVector<MT>&v, void *context)
 *   where context is a pointer to user data passed as a parameter to GMRES.
 * \param context pointer to user-defined data which will be passed on to the
 *   operator callback function.
 * \param b right-hand side
 * \param x solution vector
 * \param tol tolerance limit (relative residual target)
 * \param precon pointer to preconditioner instance (or NULL for none)
 * \param restart restart interval
 * \param iter pointer to variable receiving actual number of iterations
 *   (or NULL if not required)
 * \param res pointer to variable receiving final residual (or NULL if not
 *   required)
 * \return Currently always 0.
 */
template<class MT>
int GMRES (TVector<MT> (*Mv_clbk)(const TVector<MT> &v, void *context),
    void *context, const TVector<MT> &b, TVector<MT> &x, double tol,
    TPreconditioner<MT> *precon = 0, int restart = 10, int maxit = 0,
    int *iter = 0, double *res = 0);

template<class MT>
std::ostream &operator<< (std::ostream &os, const TMatrix<MT> &mat);

// ==========================================================================
// class TMatrix
// ==========================================================================
/**
 * \brief Virtual base class for all matrix types (dense and sparse)
 *
 * The following template types have been instantiated:
 * - TMatrix<double> (RMatrix)
 * - TMatrix<float> (FMatrix)
 * - TMatrix<complex> (CMatrix)
 * - TMatrix<scomplex> (SCMatrix)
 * - TMatrix<int> (IMatrix>
 */
template<class MT> class TMatrix {
public:
    enum RC { ROW, COL };  // row/column identifier

    /**
     * \brief Create a matrix of size 0 x 0.
     */
    TMatrix ()
    { rows = cols = 0; }

    /**
     * \brief Create a matrix of logical size nrows x ncols
     * \param nrows number of rows
     * \param ncols number of columns
     * \remarks This constructor does not allocate memory for data storage.
     *   Storage details are left to derived classes.
     */
    TMatrix (int nrows, int ncols)
    { rows = nrows, cols = ncols; }

    /**
     * \brief Create a matrix as a copy of another matrix
     * \param m matrix to be copied
     */
    TMatrix (const TMatrix<MT> &m);

    /**
     * \brief Destroy the matrix
     * \remarks Derived classes deallocate their data blocks here.
     */
    virtual ~TMatrix() {}

    /**
     * \brief Return a matrix dimension
     * \param rc ROW or COL
     * \return size of the specified dimension
     * \sa nRows(), nCols()
     */
    int Dim (RC rc) const { return (rc==ROW ? rows : cols); }

    /**
     * \brief Return number of rows of the matrix
     * \return number of rows
     * \sa nCols(), Dim()
     */
    int nRows () const { return rows; }

    /**
     * \brief Return number of columns of the matrix
     * \return number of columns
     * \sa nRows(), Dim()
     */
    int nCols () const { return cols; }

    /**
     * \brief Return sparse storage flag
     * \return true for sparse, false for dense matrix classes
     */
    bool isSparse() const { return StorageType() != MATRIX_DENSE; }

    /// \brief Return dense storage flag
    /// \return true for dense, false for sparse matrix classes
    bool isFull() const { return StorageType() == MATRIX_DENSE; }

    /**
     * \brief Matrix storage class
     * \return storage class identifier
     */
    virtual MatrixStorage StorageType() const = 0;

    /**
     * \brief Resize and reset the matrix.
     *
     * Resets the logical size of the matrix to nrows x ncols. Derived classes
     * reallocate the data block and any index arrays, and reset the matrix
     * elements to zero.
     * \param nrows new number of matrix rows
     * \param ncols new number of matrix columns
     */
    virtual void New (int nrows, int ncols);

    /**
     * \brief Retrieve a matrix element
     * \param r matrix row (0 <= r < nRows())
     * \param c matrix column (0 <= c < nCols())
     * \return matrix element (*this)<sub>r,c</sub>
     * \note This is a read operation and returns the element value. For
     *   writing operations, use Put() or operator().
     * \sa operator()
     */
    virtual MT Get (int r, int c) const = 0;

    /**
     * \brief Matrix element access (read only)
     * \param r row index (0 <= r < nRows())
     * \param c column index (0 <= c < nCols())
     * \return matrix element (r,c)
     * \sa Get
     */
    inline MT operator () (int r, int c) const { return Get (r,c); }

    /**
     * \brief Returns a vector containing a copy of row `r'
     * \param r row index (>= 0)
     * \return vector containing row r
     * \note Sparse matrix types expand to a dense row, with missing entries
     *   filled with zeros, so this can be used as a "scatter" operation.
     * \sa SparseRow, SetRow
     */
    virtual TVector<MT> Row (int r) const = 0;

    /**
     * \brief Substitute a row of the matrix.
     *
     * Replaces row r by 'row'.
     * \param r row index (>= 0)
     * \param row vector containing new row values
     * \note For sparse matrix types, this function compresses the row vector
     *    by stripping zero entries, so this can be used as a "gather"
     *    operation.
     * \sa Row, SparseRow
     */
    virtual void SetRow (int r, const TVector<MT> &row)
    { ERROR_UNDEF; }

    /**
     * \brief Returns a row of the matrix in sparse format.
     *
     * Returns a list of column indices and values for all allocated
     * entries of row r. This is only really useful for sparse matrix types.
     * For dense matrices this simply returns the complete row, with indices
     * 0, 1, 2, ... n-1
     * \param r row index (>= 0)
     * \param colidx pointer to array of column indices
     * \param val pointer to array of element values
     * \return Actual number of allocated matrix entries in the row.
     * \note The arrays must be allocated by the caller and be of sufficient
     *   size.
     * \sa Row, SetRow
     */
    virtual int SparseRow (int r, idxtype *colidx, MT *val) const = 0;

    /**
     * \brief Returns a vector containing a copy of column 'c'
     * \param c column index (>= 0)
     * \return vector containing column c
     * \note Sparse matrix types expand to a dense column, with missing entries
     *   filled with zeros, so this can be used as a "scatter" operation.
     * \sa Row
     */
    virtual TVector<MT> Col (int c) const = 0;

    /**
     * \brief Returns the matrix diagonal as a vector
     * \return vector containing diagonal elements
     * \note The length of the returned vector is min(rows,cols)
     * \note The base class implementation uses the Get() method to retrieve
     *   the matrix elements. Derived matrix classes should overwrite this
     *   method to improve performance.
     */
    virtual TVector<MT> Diag () const;

    /**
     * \brief Returns vector of column norms
     * \return column L2-norm vector of dimension 'cols'
     * \note Returns vector \e c with
     *   \f[
     *   c_j = \sqrt{\sum_i M_{ij}^2}
     *   \f]
     * \note Uses the Col() method to retrieve column vectors. Derived classes
     *   which do not have an efficient Col() implementation should overload
     *   this method.
     */
    virtual TVector<MT> ColNorm () const;

    virtual void ColScale (const TVector<MT> &scale) = 0;
    // scales the columns with 'scale'

    virtual void RowScale (const TVector<MT> &scale) = 0;
    // scales the rows with 'scale'

    virtual void Unlink () = 0;
    // removes the matrix' link to its data block and deletes the data
    // block, if necessary

    virtual void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    // returns the result of Ax in b; this is an alternative to
    // operator* which avoids local vector creation and copy

    inline TVector<MT> operator* (const TVector<MT> &x) const
    { TVector<MT> b; Ax (x, b); return b; }
    // matrix x vector multiplication

    virtual void ATx (const TVector<MT> &x, TVector<MT> &b) const;
    // returns the result of A^T x in b

    inline TVector<MT> ATx (const TVector<MT> &x) const
    { TVector<MT> b; ATx (x, b); return b; }

    virtual void Transpone ()
    { ERROR_UNDEF; }
    // Replace matrix with its transpose

    virtual MT RowMult (int r, MT *x) const
    { ERROR_UNDEF; return 0; }
    // inner product of row r with full vector given by x
    // where dimension(x) >= ncols

    /**
     * \brief Write matrix to ASCII stream.
     *
     * Outputs the matrix in a generic format, with elements separated
     * by whitespace, and rows separated by line feeds.
     * \param os output stream
     */
    void Export (std::ostream &os) const;

    void Print (std::ostream &os = std::cout, int n = 80) const;
    // Pretty-print matrix in dense format to os, using a maximum of
    // n characters per line

    void PrintNzeroGraph (char *fname);
    // Generates a PPM bitmap file 'fname' with the image of the nonzero
    // elements of the matrix.

    // **** friends ****

    /// Return transp(*this) * *this as a symmetric matrix.
    friend TSymMatrix<MT> ATA<> (const TMatrix<MT> &A);
    friend TSymMatrix<MT> AAT<> (const TMatrix<MT> &A);
    friend TVector<MT>ATA_diag<> (const TMatrix<MT> &A);

    friend int PCG<> (const TMatrix<MT> &A, const TVector<MT> &b,
        TVector<MT> &x, double &tol, TPreconditioner<MT> *precon,
        int maxit);

    friend void PCG<> (const TMatrix<MT> &A, const TVector<MT> *b,
        TVector<MT> *x, int nrhs, double tol, int maxit,
        TPreconditioner<MT> *precon, IterativeSolverResult *res);

    friend int BiCGSTAB<> (const TMatrix<MT> &A,
        const TVector<MT> &b, TVector<MT> &x, double &tol,
        TPreconditioner<MT> *precon, int maxit);
    // biconjugate gradient stabilised method.

    friend void BiCGSTAB<> (const TMatrix<MT> &A,
        const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit,
        TPreconditioner<MT> *precon, IterativeSolverResult *res);

    // Explicit linear solvers

    virtual int pcg (const TVector<MT> &b, TVector<MT> &x,
        double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0)
        const;

    virtual void pcg (const TVector<MT> *b, TVector<MT> *x,
        int nrhs, double tol, int maxit = 0,
        TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0)
        const;

    virtual int bicgstab (const TVector<MT> &b, TVector<MT> &x,
        double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0)
        const;

    virtual void bicgstab (const TVector<MT> *b, TVector<MT> *x,
        int nrhs, double tol, int maxit = 0,
        TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0)
        const;

    friend std::ostream &operator<< <> (std::ostream &os,
        const TMatrix<MT> &mat);
    // stream output

protected:
    int rows, cols;	// number of rows and columns of the matrix
};

// ==========================================================================
// typedefs for specific instances of `TMatrix'

typedef TMatrix<double>   RMatrix;	// 'double real'
typedef TMatrix<float>    FMatrix;      // 'single real'
typedef TMatrix<std::complex<double> >  CMatrix;	// 'complex'
typedef TMatrix<std::complex<float> > SCMatrix;     // 'single complex'
typedef TMatrix<int>      IMatrix;      // 'integer'


// ==========================================================================
// ==========================================================================
// Member definitions

template<class MT>
TMatrix<MT>::TMatrix (const TMatrix<MT> &m)
{
    rows = m.rows;
    cols = m.cols;
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::New (int nrows, int ncols)
{
    rows = nrows;
    cols = ncols;
}

// --------------------------------------------------------------------------

template<class MT>
TVector<MT> TMatrix<MT>::Diag () const
{
    int i, n = (rows < cols ? rows : cols);
    TVector<MT> diag(n);
    for (i = 0; i < n; i++)
        diag[i] = Get(i,i);
    return diag;
}

// --------------------------------------------------------------------------

template<class MT>
TVector<MT> TMatrix<MT>::ColNorm () const
{
    int j;
    TVector<MT> cn(cols);
    for (j = 0; j < cols; j++) {
	cn[j] = (MT)l2norm (Col(j));
    }
    return cn;
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT (cols == x.Dim(), "Argument 1: vector has wrong size");
    if (b.Dim() != rows) b.New(rows);
    for (int i = 0; i < rows; i++) b[i] = Row(i) & x;
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT (rows == x.Dim(), "Argument 1: vector has wrong size");
    if (b.Dim() != cols) b.New(cols);
    for (int i = 0; i < cols; i++) b[i] = Col(i) & x;
}

// --------------------------------------------------------------------------

template<class MT>
TSymMatrix<MT> ATA (const TMatrix<MT> &A)
{
    int i, j, nc = A.nCols();
    TSymMatrix<MT> ata(nc);
    for (i = 0; i < nc; i++) {
        TVector<MT> col = A.Col(i);
	for (j = 0; j <= i; j++)
	    ata(i,j) = col & A.Col(j);
    }
    return ata;
}

// --------------------------------------------------------------------------

template<>
inline int TMatrix<float>::pcg (const FVector &b, FVector &x,
    double &tol, TPreconditioner<float> *precon, int maxit) const
{
    return PCG (*this, b, x, tol, precon, maxit);
}

template<>
inline int TMatrix<double>::pcg (const RVector &b, RVector &x,
    double &tol, TPreconditioner<double> *precon, int maxit) const
{
    return PCG (*this, b, x, tol, precon, maxit);
}

template<class MT>
int TMatrix<MT>::pcg (const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit) const
{
    ERROR_UNDEF;
    return 0;
}

// --------------------------------------------------------------------------

template<>
inline void TMatrix<float>::pcg (const FVector *b, FVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<float> *precon,
    IterativeSolverResult *res) const
{
    PCG (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<>
inline void TMatrix<double>::pcg (const RVector *b, RVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<double> *precon,
    IterativeSolverResult *res) const
{
    PCG (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<class MT>
void TMatrix<MT>::pcg (const TVector<MT> *b, TVector<MT> *x, int nrhs,
    double tol, int maxit, TPreconditioner<MT> *precon,
    IterativeSolverResult *res) const
{
    ERROR_UNDEF;
}

// --------------------------------------------------------------------------

template<>
inline int TMatrix<float>::bicgstab (const FVector &b, FVector &x,
    double &tol, TPreconditioner<float> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

template<>
inline int TMatrix<double>::bicgstab (const RVector &b, RVector &x,
    double &tol, TPreconditioner<double> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

template<class MT>
int TMatrix<MT>::bicgstab (const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit) const
{
    ERROR_UNDEF;
    return 0;
}

// --------------------------------------------------------------------------

template<>
inline void TMatrix<float>::bicgstab (const FVector *b, FVector *x, int nrhs,
    double tol, int maxit,TPreconditioner<float> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<>
inline void TMatrix<double>::bicgstab (const RVector *b, RVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<double> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<class MT>
void TMatrix<MT>::bicgstab (const TVector<MT> *b, TVector<MT> *x, int nrhs,
    double tol, int maxit, TPreconditioner<MT> *precon,
    IterativeSolverResult *res) const
{
    ERROR_UNDEF;
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::Export (std::ostream &os) const
{
    int r, c;
    for (r = 0; r < rows; r++) {
	for (c = 0; c < cols; c++) {
	    os << Get(r,c);
	    os << (c == cols-1 ? '\n' : ' ');
	}
    }
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::Print (std::ostream &os, int n) const
{
    int r, c, nc, minc, maxc, maxlen = 0;

    maxlen = 11;
    nc = n / (maxlen+1);
    if (cols <= nc) {
        for (r = 0; r < rows; r++) {
	    for (c = 0; c < cols; c++) {
	        os.width(maxlen);
		os.precision(5);
	        os << Get(r,c) << ' ';
	    }
	    os << std::endl;
	}
    } else {
        minc = 0;
	while (minc < cols) {
	    maxc = minc+nc;
	    if (maxc > cols) maxc = cols;
	    os << "*** Columns " << (minc+1) << " through " << maxc
	       << " ***\n";
	    for (r = 0; r < rows; r++) {
	        for (c = minc; c < maxc; c++) {
		    os.width(maxlen);
		    os.precision(5);
		    os << Get(r,c) << ' ';
		}
		os << std::endl;
	    }
	    minc = maxc;
	}
    }	
}

// --------------------------------------------------------------------------

template<class MT>
void TMatrix<MT>::PrintNzeroGraph (char *fname)
{
    int i, r, c;
    std::ofstream ofs (fname);
    ofs << "P1" << std::endl;
    ofs << "# CREATOR: TOAST TMatrix::PrintNzeroGraph" << std::endl;
    ofs << cols << ' ' << rows << std::endl;
    for (i = r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
	    ofs << (Get(r,c) != (MT)0 ? '1' : '0');
	    ofs << (++i % 35 ? ' ' : '\n');
	}
    }
}

// --------------------------------------------------------------------------

#endif // !__MATRIX_H
