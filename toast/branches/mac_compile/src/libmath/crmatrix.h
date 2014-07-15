// -*-C++-*-
// ==========================================================================
// Module mathlib
// File crmatrix.h
// Declaration of template class TCompRowMatrix ('template compressed-row
//  matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RCompRowMatrix = TCompRowMatrix<double>	  ('real')
//	FCompRowMatrix = TCompRowMatrix<float>	  ('float')
//	CCompRowMatrix = TCompRowMatrix<complex>  ('complex')
//	ICompRowMatrix = TCompRowMatrix<int>	  ('integer')
//
// Inheritance:
// ------------
// TGenericSparseMatrix ----> TCompRowMatrix
//
// Notes:
// ------
// * TCompRowMatrix does not support dynamic growth, i.e. a write operation
//   must not address a non-existing entry. If you need dynamic growth,
//   use a TCoordMatrix and convert it to a TCompRowMatrix once it is fully
//   assembled.
// ==========================================================================

#ifndef __CRMATRIX_H
#define __CRMATRIX_H

#include "gsmatrix.h"

// ==========================================================================
// Nonmember declarations

template<class MT>class TCompRowMatrix;
template<class MT>class TCoordMatrix;
template<class MT>class TPrecon_IC;
template<class MT>class TDenseMatrix;
template<class MT>class TDiagMatrix;

void BlockExpand (int *rowptr, int *colidx, int n,
		  int *&browptr, int *&bcolidx, int &bn,
		  int blockn, int blockm);
// expand the nonsparsity structure (rowptr,colidx) of an (n x m) matrix
// into a (n blockn x m blockm) matrix structure (browptr, bcolidx) by
// expanding every entry in the original matrix into an (blockn x blockm)
// block

template<class MT>
TCompRowMatrix<MT>transp (const TCompRowMatrix<MT> &m);

template<class MT>
double l2norm (const TCompRowMatrix<MT> &A);

template<class MT>
TCompRowMatrix<MT> kron (const TCompRowMatrix<MT> &A,
    const TCompRowMatrix<MT> &B);

template<class MT>
TCompRowMatrix<MT> cath (const TCompRowMatrix<MT> &A,
    const TCompRowMatrix<MT> &B);

template<class MT>
TCompRowMatrix<MT> catv (const TCompRowMatrix<MT> &A,
    const TCompRowMatrix<MT> &B);

template<class MT>
bool CholeskyFactorize (const TCompRowMatrix<MT> &A, TCompRowMatrix<MT> &L,
    TVector<MT> &d, bool recover = false);

template<class MT>
bool IncompleteCholeskyFactorize (const TCompRowMatrix<MT> &A,
    TCompRowMatrix<MT> &L, TVector<MT> &d, bool recover = false);
// Perform an incomplete Cholesky factorisation without creating additional
// fillin. Lower triangle (without diagonal) stored in L, diagonal stored in
// d. If recover==true, then non-SPD matrices will not cause a fatal error.

template<class MT>
void LUFactorize (TCompRowMatrix<MT> &A, bool LUrealloc = true);

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TVector<MT> &b, TVector<MT> &x);

template<class MT>
TVector<MT> CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TVector<MT> &b);

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TDenseMatrix<MT> &BT, TDenseMatrix<MT> &XT, int n = 0);

template<class MT>
void LUSolve (const TCompRowMatrix<MT> &LU, const TVector<MT> &b,
    TVector<MT> &x);

template<class MT>
void LU (TCompRowMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x);

template<class MT>
int ILUSolve (TCompRowMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double tol = 1e-10, double droptol = 1e-3, int maxit = 500);
// solve general system Ax=b with ILU solver

template<class MT>
int ILUSymSolve (TCompRowMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double tol = 1e-10, double droptol = 1e-3, int maxit = 500);
// solve symmetric system Ax=b with ILU solver

#ifdef USE_CUDA_FLOAT
template<class MT>
int PCG (const TCompRowMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
void PCG (const TCompRowMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x, int nrhs, double tol, int maxit = 0,
    TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0);

template<class MT>
int BiCGSTAB (const TCompRowMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
void BiCGSTAB (const TCompRowMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x, int nrhs, double tol, int maxit = 0,
    TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0);
#endif // USE_CUDA_FLOAT


#ifdef ML_INTERFACE

template<class MT>
int ML_matvec (ML_Operator *Amat, int in_length, double p[],
    int out_length, double ap[]);

template<class MT>
int ML_getrow (ML_Operator *Amat, int N_requested_rows,
    int requested_rows[], int allocated_space, int columns[],
    double values[], int row_lenghts[]);

#endif // ML_INTERFACE

template<class MT>
std::istream &operator>> (std::istream &is, TCompRowMatrix<MT> &m);

template<class MT>
std::ostream &operator<< (std::ostream &os,
    const TCompRowMatrix<MT> &m);

// ==========================================================================
// class TCompRowMatrix
// ==========================================================================
/**
 * \brief Compressed-row sparse matrix class
 *
 * This matrix class represents its data by a data array 'val' of length nz, a
 * column index array 'colidx' of length nz, and a row pointer array 'rowptr'
 * of length nr+1, where nz is the number of allocated (nonzero) elements, and
 * nr is the row dimension of the matrix.
 * \note Column indices are zero-based, and are stored sequentially (in
 *   ascending order) for each row
 * \note Row pointers are zero based, and point to the entry in the data
 *   array containing the first element of the row.
 * \note rowptr[0] is always 0; rowptr[nr+1] is always equal to nz.
 * \note The number of entries in row i (i >= 0) is rowptr[i+1]-rowptr[i]
 */
template<class MT> class TCompRowMatrix: public TGenericSparseMatrix<MT> {

    friend class TPrecon_IC<MT>;
    friend class CForwardSolver;
    friend class CForwardSolverPP;
    friend class SCCompRowMatrixMixed;

public:

    /**
     * \brief Creates a compressed-row matrix of dimension 0 x 0.
     */
    TCompRowMatrix ();

    /**
     * \brief Creates a compressed-row matrix of dimension rows x cols.
     * \param rows row dimension (>= 0)
     * \param cols column dimension (>= 0)
     * \note No space is allocated for matrix element storage. Use Initialise()
     *   to allocate data space.
     */
    TCompRowMatrix (int rows, int cols);

    /**
     * \brief Creates a compressed-row matrix of dimension rows x cols, and
     *   sets up data storage with specified fill structure.
     * \param rows row dimension (>= 0)
     * \param cols column dimension (>= 0)
     * \param _rowptr row pointer array of length rows+1
     * \param _colidx column index array of length nzero
     * \param data data array of length nzero
     * \note nzero is the number of allocated entries, and is equal to
     *   _rowptr[rows].
     * \note If no data array is provided, all allocated entries are
     *   initialised to zero.
     */
    TCompRowMatrix (int rows, int cols,
        const idxtype *_rowptr, const idxtype *_colidx,
	const MT *data = 0);

    TCompRowMatrix (int rows, int cols,
	idxtype *_rowptr, idxtype *_colidx,
	MT *data, CopyMode cmode);

    /**
     * \brief Creates a compressed-row matrix as a copy of matrix m
     * \param m matrix to create a copy from.
     * \note If m has allocated its column-access index arrays, then they are
     *   also allocated for *this.
     */
    TCompRowMatrix (const TCompRowMatrix<MT> &m);

    TCompRowMatrix (const TCoordMatrix<MT> &m);
    TCompRowMatrix (const TDenseMatrix<MT> &m);
    ~TCompRowMatrix ();

    MatrixStorage StorageType() const { return MATRIX_COMPROW; }

    TCompRowMatrix<MT> Submatrix (int r1, int r2, int c1, int c2) const;
    // return a submatrix of dimension r2-r1 x c2-c1 consisting of the block
    // (r1 ... r2-1) x (c1 ... c2-1)

    TCompRowMatrix<MT> Subrows (int r1, int r2) const
    { return Submatrix (r1, r2, 0, this->cols); }
    // return a submatrix of dimension r2-r1 x cols consisting of the block
    // (r1 ... r2-1) x (0 ... cols-1)

    TCompRowMatrix<MT> Subcols (int c1, int c2) const
    { return Submatrix (0, this->rows, c1, c2); }
    // return a submatrix of dimension rows x c2-c1 consisting of the block
    // (0 ... rows-1) x (c1 ... c2-1)

    friend TCompRowMatrix<std::complex<double> > cplx
    (const TCompRowMatrix<double> &A);
    // type conversion

    TCoordMatrix<MT> MakeCoordMatrix () const;
    // Convert *this into a coordinate storage matrix

    void New (int nrows, int ncols);
    // reset logical dimension to nrows x ncols
    // deallocate data and index lists, i.e. reset entries to zero
    // Disables column access (call SetColumnAccess to enable)

    void Identity (int n);
    // reset logical dimension to n x n, and set to identity matrix

    void DiagV (const TVector<MT> &x);
    // reset to diagonal matrix

    void Unlink ();

    void Initialise (const idxtype *_rowptr, const idxtype *_colidx,
        const MT *data = 0);
    // Redefine sparse structure and assign data
    // If data are not provided, all entries are set to zero

    int GetSparseStructure (const idxtype **_rowptr, const idxtype **_colidx) const
    { *_rowptr = rowptr, *_colidx = colidx;
      return TGenericSparseMatrix<MT>::nval; }

    TCompRowMatrix &operator= (const TCompRowMatrix<MT> &m);
    // assignment

    TCompRowMatrix &operator= (const TCoordMatrix<MT> &m);
    // assignment from coordinate storage matrix

    TCompRowMatrix<MT> operator- () const;
    // unary minus

    TCompRowMatrix<MT> operator* (MT f) const;
    // *this * f

    TCompRowMatrix<MT> operator* (const TDiagMatrix<MT> &D) const;
  //   TCompRowMatrix<MT> operator* ( RDiagMatrix &D) const;
    // this * diagonal D

    inline TVector<MT> operator* (const TVector<MT> &x) const
    { TVector<MT> b(TMatrix<MT>::rows); Ax(x,b); return b; }
    // Multiplies *this with x and returns the result

    inline TCompRowMatrix<MT> operator* (const TCompRowMatrix<MT> &m) const
    { TCompRowMatrix<MT> res; AB(m, res); return res; }
    // Matrix multiplication: returns *this * m

    TCompRowMatrix<MT> operator+ (const TCompRowMatrix<MT> &m) const;
    // *this + m

    TCompRowMatrix<MT> operator- (const TCompRowMatrix<MT> &m) const
    { return operator+(-m); }
    // *this - m

    TCompRowMatrix<MT> &operator+= (const TCompRowMatrix<MT> &m);

    TCompRowMatrix<MT> &operator-= (const TCompRowMatrix<MT> &m)
    { return operator+=(-m); }

    /**
     * \brief Concatenate two matrices horizontally
     * \param A first matrix argument
     * \param B second matrix argument
     * \return Result of concatenation [A B]
     * \note A and B must have the same number of rows.
     * \note This method has the functionality of the MATLAB construct
     *   C = (A,B)
     */
    friend TCompRowMatrix<MT> cath<> (const TCompRowMatrix<MT> &A,
	const TCompRowMatrix<MT> &B);

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
    friend TCompRowMatrix<MT> catv<> (const TCompRowMatrix<MT> &A,
	const TCompRowMatrix<MT> &B);

    TCompRowMatrix &merge (const TCompRowMatrix<MT> &m);
    // Merges (adds) 'm' into 'this'. m can have different dimensions and
    // different fill structure

    MT &operator() (int r, int c);
    MT Get (int r, int c) const;
    // access to element in row r, col c.
    // First version for write access, second for read access

    bool Exists (int r, int c) const;
    // true if an entry for the element is allocated

    TVector<MT> Row (int r) const;
    // Returns row r as a vector

    TVector<MT> Col (int c) const;
    // Returns column c as a vector

    int SparseRow (int r, idxtype *ci, MT *rv) const;
    // See TMatrix

    void SetRow (int r, const TVector<MT> &row);
    // Replace row r with 'row', removing zero elements

    void SetRows (int r0, const TCompRowMatrix<MT> &rows);
    // starting with row r0, replace matrix rows with 'rows'
 
    void RemoveRow(int c);
    //Remove row c;

    void ColScale (const TVector<MT> &scale);
    // scales the columns with 'scale'

    void RowScale (const TVector<MT> &scale);
    // scales the rows with 'scale'

    MT GetNext (int &r, int &c) const;

    void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    void Ax (const TVector<MT> &x, TVector<MT> &b, int r1, int r2) const;
    void ATx (const TVector<MT> &x, TVector<MT> &b) const;

    void Ax_cplx (const TVector<std::complex<double> > &x,
		  TVector<std::complex<double> > &b)
	const;
    // special case: multiply with complex vector (only instantiated for
    // real matrix)

    void ATx_cplx (const TVector<std::complex<double> > &x,
		   TVector<std::complex<double> > &b)const;
    // special case: multiply transpose with complex vector (only instantiated
    // for real matrix)

    void Ax (const TVector<MT> *x, TVector<MT> *b, int nrhs) const;

    void AB (const TCompRowMatrix<MT> &B, TCompRowMatrix<MT> &C) const;
    // sparse matrix x matrix multiplication: C = *this * B

    void Transpone ();
    // Replace matrix with its transpose

    friend TCompRowMatrix<MT> transp<> (const TCompRowMatrix<MT> &m);
    // returns transpose of m

    friend double l2norm<> (const TCompRowMatrix<MT> &A);

    friend TCompRowMatrix<MT> kron<> (const TCompRowMatrix<MT> &A,
        const TCompRowMatrix<MT> &B);
    // Kronecker matrix product

    MT RowMult (int r, MT *x) const;
    // inner product of row r with full vector given by x
    // where dimension(x) >= ncols

    void Sort() const;
    // sort row entries for ascending column index

    int Shrink ();
    // Remove all entries with zero value. Return value is the number of
    // eliminated entries.

    void SetColAccess (bool yes = true) const;
    // Initialise or remove index lists required for column access

    void SetDiagAccess (bool yes = true) const;
    // Initialise or remove index list required for diagonal access

    double LargestInRow (int r, int i = 0) const;
    // Returns magnitude of the largest entry in row r, starting from the
    // row's i-th entry

    double LargestInCol (int c, int i = 0) const;
    // Returns magnitude of the largest entry in column c, starting from the
    // column's i-th entry

    // Explicit linear solvers

    int pcg (const TVector<MT> &b, TVector<MT> &x,
        double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0)
        const;

    void pcg (const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol,
        int maxit = 0, TPreconditioner<MT> *precon = 0,
        IterativeSolverResult *res = 0) const;

    int bicgstab (const TVector<MT> &b, TVector<MT> &x,
        double &tol, TPreconditioner<MT> *precon = 0, int maxit = 0)
        const;
  
    void bicgstab (const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol,
        int maxit = 0, TPreconditioner<MT> *precon = 0,
        IterativeSolverResult *res = 0) const;

    void SymbolicCholeskyFactorize (idxtype *&frowptr, idxtype *&fcolidx)
        const;
    // Calculate the sparse fill-in structure of the lower triangle of
    // the Cholesky factorisation of the matrix (excluding diagonal)
    // and return in index lists `frowptr' and 'fcolidx'

    void CalculateIncompleteCholeskyFillin (idxtype *&frowptr,
        idxtype *&fcolidx) const;
    // Calculate sparse structure for the lower triangle of the the
    // incomplete Cholesky factorisation of the matrix. This is simply the
    // structure of the lower triangle of the original matrix, excluding the
    // diagonal. No additional fill-in is considered.

    friend bool CholeskyFactorize<> (const TCompRowMatrix<MT> &A,
        TCompRowMatrix<MT> &L, TVector<MT> &d, bool recover);
    // Perform Cholesky factorisation of A and return the result in lower
    // triangle L and diagonal d (L does not contain diagonal elements)
    // L must have been initialised with the index lists returned from
    // A.CalculateCholeskyFillin()

    friend bool IncompleteCholeskyFactorize<> (const TCompRowMatrix<MT> &A,
        TCompRowMatrix<MT> &L, TVector<MT> &d, bool recover);
    // As above, but does not create additional fill-in. The structure of L
    // must have been set according to A.CalculateIncompleteCholeskyFillin()

    friend void LUFactorize<> (TCompRowMatrix<MT> &A, bool LUrealloc);
    // Perform LU factorisation of A and return result in LU
    // Fill-in is calculated on the fly if data length of LU is 0
    // Otherwise LU is expected to have correct fillin structure (e.g. from
    // previous call to LUFactorize)

    friend void CholeskySolve<> (const TCompRowMatrix<MT> &L,
        const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x);
    // Returns solution of LL^T x = b in x
    // d contains the diagonal elements of L, as returned by CholeskyFactorize

    friend TVector<MT> CholeskySolve<> (const TCompRowMatrix<MT> &L,
        const TVector<MT> &d, const TVector<MT> &b);
    // As above, but returns solution directly as function result (involves
    // local vector construction)

    friend void CholeskySolve<> (const TCompRowMatrix<MT> &L,
        const TVector<MT> &d, const TDenseMatrix<MT> &BT, TDenseMatrix<MT> &XT,
	int n);
    // Computes the solution of LL^T XT^T = BT^T for up to n columns of BT^T
    // If n==0 then all columns of BT^T are used

    friend void LUSolve<> (const TCompRowMatrix<MT> &LU, const TVector<MT> &b,
        TVector<MT> &x);
    // Returns solution of LU x = b

    friend int ILUSolve<> (TCompRowMatrix<MT> &A, const TVector<MT> &b,
        TVector<MT> &x, double tol, double droptol, int maxit);
    // solve general system Ax=b with ILU solver

    friend int ILUSymSolve<> (TCompRowMatrix<MT> &A, const TVector<MT> &b,
        TVector<MT> &x, double tol, double droptol, int maxit);
    // solve symmetric system Ax=b with ILU solver

#ifdef USE_CUDA_FLOAT
    friend int PCG<> (const TCompRowMatrix<MT> &A,
        const TVector<MT> &b, TVector<MT> &x, double &tol,
        TPreconditioner<MT> *precon, int maxit);

    friend void PCG<> (const TCompRowMatrix<MT> &A,
        const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit,
        TPreconditioner<MT> *precon, IterativeSolverResult *res);

    friend int BiCGSTAB<> (const TCompRowMatrix<MT> &A,
        const TVector<MT> &b, TVector<MT> &x, double &tol,
        TPreconditioner<MT> *precon, int maxit);
    // biconjugate gradient stabilised method.

    friend void BiCGSTAB<> (const TCompRowMatrix<MT> &A,
        const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit,
        TPreconditioner<MT> *precon, IterativeSolverResult *res);

#endif

    void ExportHB (std::ostream &os);
    // output matrix in Harwell-Boeing format

    /**
     * \brief Write sparse matrix to ASCII output stream.
     *
     * This writes the nonzero elements of the matrix in a row-column-value
     * format
     * \code
     * <r_i> <c_i> <val_i>
     * \endcode
     * containing integer row and column index, and one (double) or two
     * (complex) floating point values, all white-space separated.
     * \note The row and column indices are 1-based.
     * \note This format can be read into MATLAB by using the 'load' command
     * and converted into a MATLAB sparse matrix by using the 'spconvert'
     * command.
     */
    void ExportRCV (std::ostream &os);

    void SplitExport (const char *rootname);
    // write row index, column index and data in separate files

    friend std::istream &operator>> <> (std::istream &is,
        TCompRowMatrix<MT> &m);
    friend std::ostream &operator<< <> (std::ostream &os,
        const TCompRowMatrix<MT> &m);
    // IO in generic format

protected:
    void ReplaceRow (int row, int nz, int *rcolidx, MT *rval = 0);
    // Replace 'row' with the nonzero structure defined by 'rcolindex' (nz
    // entries). Assign values if provided.
    // This method is expensive if nz is different from the current entry
    // count of 'row' (requires memory re-allocation)

    void ReplaceRow (int row, const TVector<MT>& vec);
    // Replace 'row' with the nonzero entries of the full row stored in
    // vector 'data'

    MT row_mult (int r1, int r2, int from, int to) const;
    // Returns dot product of rows r1 and r2, between columns `from' and 'to'
    // REQUIRES ROWS TO BE SORTED WITH ASCENDING COLUMN INDEX

    MT sparsevec_mult (int *idx1, MT *val1, int n1,
		       int *idx2, MT *val2, int n2) const;
    // Returns dot product of two sparse vectors given by index lists idx and
    // data array val (n is length of sparse array)
    // REQUIRES ROWS TO BE SORTED WITH ASCENDING COLUMN INDEX

public:
    idxtype *rowptr; 
    idxtype *colidx;

private:
    int Get_index (int r, int c) const;
    // returns offset into data array of element at row r and column c
    // returns -1 if entry does not exist

    MT Get_sorted (int r, int c) const;
    void Put_sorted (int r, int c, MT v);
    int Get_index_sorted (int r, int c) const;
    // get or set entries assuming the matrix is sorted

    void CholeskySubst (const TVector<MT> &d, const MT *b, MT *x) const;
    // Low-level driver routine for Cholesky forward and backward substitution

    mutable int iterator_pos;

    mutable bool sorted; // true if row entries are sorted

    mutable bool diag_access;
    // true if diagonal access has been initialised
    mutable int *diagptr;
    // val[diagptr[i]] retrieves the diagonal element of row i
    // if diagptr[i] < 0 then i has no diagonal entry

    mutable bool col_access; // true if column access has been initialised
    // the following lists are only valid if col_access == true
    mutable idxtype *colptr;
    // offset to the first entry for column i in rowidx
    mutable idxtype *rowidx;
    // rowidx[colptr[i]+j] contains the row number of the j-th entry of col i
    mutable idxtype *vofs;
    // vofs[colptr[i]+j] contains the offset into data array val of the j-th
    // entry of column i

#ifdef ML_INTERFACE
friend int ML_matvec<> (ML_Operator *Amat, int in_length, double p[],
    int out_length, double ap[]);
friend int ML_getrow<> (ML_Operator *Amat, int N_requested_rows,
    int requested_rows[], int allocated_space, int columns[],
    double values[], int row_lenghts[]);
#endif // ML_INTERFACE

};

// ==========================================================================
// typedefs for specific instances of `TCompRowMatrix'

typedef TCompRowMatrix<double>	 RCompRowMatrix;	// 'real'
typedef TCompRowMatrix<float>	 FCompRowMatrix;	// 'float'
typedef TCompRowMatrix<std::complex<double> > CCompRowMatrix;	// 'complex'
typedef TCompRowMatrix<std::complex<float> > SCCompRowMatrix;	// 'single complex'
typedef TCompRowMatrix<int>	 ICompRowMatrix;	// 'integer'

// ==========================================================================
// extern declarations of TCompRowMatrix (only required for VS)

#ifdef UNDEF
//#ifndef __CRMATRIX_CC
extern template class MATHLIB TCompRowMatrix<double>;
extern template class MATHLIB TCompRowMatrix<float>;
extern template class MATHLIB TCompRowMatrix<std::complex<double> >;
extern template class MATHLIB TCompRowMatrix<std::complex<float> >;
extern template class MATHLIB TCompRowMatrix<int>;
#endif // !__CRMATRIX_CC

#endif // !__CRMATRIX_H
