// -*-C++-*-
// ==========================================================================
// Module mathlib
// File crmatrix_mpi.h
// Declaration of template class TCompRowMatrixMPI
// Distributed (MPI) version of TCompRowMatrix (template sparse matrix in
// compressed row format)
// Each MPI process operates on a block of matrix rows (r0 <= r < r1)
// ==========================================================================

#ifndef __CRMATRIXMPI_H
#define __CRMATRIXMPI_H

#include "mathlib.h"

// ==========================================================================
// class TCompRowMatrixMPI
// ==========================================================================
/**
 * \brief Distributed compressed-row sparse matrix class
 *
 * A sparse matrix type that distributes itself across the available
 * processes. The class interface is designed to be compatible with a subset
 * of the \ref TCompRowMatrix interface, hiding the parallelism within the
 * class implementation.
 *
 * Upon data allocation, an instance of %TCompRowMatrixMPI assigns a block of
 * rows r0 <= r < r1 to each process, such that the number of nonzero elements
 * is approximately evenly distributed. Each block is allocated as a compressed
 * row sub- matrix representing part of the complete matrix. Matrix operations
 * such as Ax or ATx are performed in parallel.
 */
template<class MT> class TCompRowMatrixMPI: public TGenericSparseMatrix<MT> {

public:
    /**
     * \brief Creates a matrix of dimension 0 x 0.
     */
    TCompRowMatrixMPI ();

    /**
     * \brief Creates a matrix of dimension rows x cols.
     * \param rows row dimension (>= 0)
     * \param cols column dimension (>= 0)
     * \note No space is allocated for matrix element storage. Use
     *   \ref Initialise to allocate data space.
     */
    TCompRowMatrixMPI (int rows, int cols);

    TCompRowMatrixMPI (int rows, int cols,
	const int *_rowptr, const int *_colidx,
	int proc_nrows, const int *proc_rows);

    /**
     * \brief Creates a matrix of dimension rows x cols, and sets up data
     *   storage with the specified fill structure.
     * \param rows row dimension (>= 0)
     * \param cols column dimension (>= 0)
     * \param _rowptr row pointer array of length rows+1
     * \param _colidx column index array of length nzero
     * \param _data optional data array of length nzero
     * \note nzero is the number of allocated entries, and is equal to
     *   _rowptr[rows].
     * \note If no data array is provided, all allocated entries are
     *   initialised to zero.
     */
    TCompRowMatrixMPI (int rows, int cols,
        const int *_rowptr, const int *_colidx,
        const MT *_data = 0);

    /**
     * \brief Matrix destructor.
     */
    ~TCompRowMatrixMPI ();

    /**
     * \brief Storage class identifier
     */
    MatrixStorage StorageType() const
    { return MATRIX_COMPROW; }

    /**
     * \brief Re-allocate fill structure and assign values.
     * \param _rowptr row pointer array of length rows+1
     * \param _colidx column index array of length nzero
     * \param _data optional data array of length nzero
     * \note nzero is the number of allocated entries, and is equal to
     *   _rowptr[rows].
     * \note If no data array is provided, all allocated entries are
     *   initialised to zero.
     */
    void Initialise (const int *_rowptr, const int *_colidx,
        const MT *_data = 0);

    /**
     * \brief Returns the row range of the current process.
     * \param [out] _r0 low row index (>= 0)
     * \param [out] _r1 high row index+1
     * \note This method returns the range of rows the current process is
     *   operating on.
     */
    void MPIRange (int *_r0, int *_r1) const
    { *_r0 = r0; *_r1 = r1; }

    /**
     * \brief Zero all elements, but keep fill structure
     * \note This method resets all allocated elements to zero, but does
     *   not deallocate the data and index arrays.
     */
    void Zero ();

    /**
     * \brief Retrieve a matrix element
     * \param r matrix row (0 <= r < nRows())
     * \param c matrix column (0 <= c < nCols())
     * \return matrix element (*this)<sub>r,c</sub>
     * \note This is a read operation and returns the element value. For
     *   writing operations, use Put().
     */
    MT Get (int r, int c) const;

    /**
     * \brief Returns a vector containing a copy of row `r'
     * \param r row index (>= 0)
     * \return vector containing row r
     * \note The matrix row is expanded into a full vector, replacing
     *   un-allocated elements with zeros.
     */
    TVector<MT> Row (int r) const;

    /**
     * \brief Returns a row of the matrix in sparse format.
     *
     * Returns a list of column indices and values for all allocated
     * entries of row r.
     * \param[in] r row index (>= 0)
     * \param[out] colidx pointer to array of column indices
     * \param[out] val pointer to array of element values
     * \return Actual number of allocated matrix entries in the row.
     * \note The arrays must be allocated by the caller and be of sufficient
     *   size.
     * \sa Row, SetRow
     */
    int SparseRow (int r, int *colidx, MT *val) const;

    /**
     * \brief Returns a vector containing a copy of column 'c'
     * \param c column index (>= 0)
     * \return vector containing column c
     * \note Sparse matrix types expand to a dense column, with missing entries
     *   filled with zeros, so this can be used as a "scatter" operation.
     * \sa Row
     */
    TVector<MT> Col (int c) const;

    virtual void RowScale (const TVector<MT> &scale);
    // scales the rows with 'scale'

    virtual void ColScale (const TVector<MT> &scale);
    // scales the columns with 'scale'

    virtual void Unlink ();
    // removes the matrix' link to its data block and deletes the data
    // block, if necessary

    bool Exists (int r, int c) const;

    MT &operator() (int r, int c);

    int Get_index (int r, int c) const;

    MT GetNext (int &r, int &c) const;

    /**
     * \brief Matrix-vector product.
     * \param[in] x vector argument (length ncols)
     * \param[out] b result vector (length nrows)
     * \note Computes Ax = b
     */
    void Ax (const TVector<MT> &x, TVector<MT> &b) const;

    void Ax (const TVector<MT> &x, TVector<MT> &b, int i1, int i2) const;

    void ATx (const TVector<MT> &x, TVector<MT> &b) const
    { xERROR("Not implemented"); }

    // ======================================================================
    // Process-specific functions
    // These methods are invoked only for the current process

    /**
     * \brief Add a value to a matrix element (process-specific)
     * \param r row index (>= 0)
     * \param c column index (>= 0)
     * \param val added value
     * \note This method is invoked only for the current process. Therefore,
     *   r must be in the row range of the process (r0 <= r < r1)
     */
    void Add_proc (int r, int c, MT val);
  
  //private:
public:
    /**
     * \brief MPI data initialisation
     * \note Sets the sze and rnk parameters.
     */
    void MPIinit();

    int rnk;    ///< MPI process id (0 <= rnk < \ref sze)
    int sze;    ///< Number of MPI processes (>= 1)

    /**
     * \brief Array of row pointers
     * \note The rowptr array contains the indices of the first nonzero
     *   entry of each row in the data array. Note that for the MPI version,
     *   the indices refer to the process data sub-block, not to the global
     *   data array.
     * \note The dimension of rowptr is m+1, where m = r1-r0.
     */
    int *rowptr;

    /**
     * \brief Array of column indices
     * \note The colidx array contains the column indices (0-based) of each
     *   allocated data entry for the process sub-block.
     * \note The dimension of colidx is \ref nval.
     */
    int *colidx;

    int r0;     ///< Low row index for this process - OBSOLETE
    int r1;     ///< High row index + 1 for this process - OBSOLETE
    int my_nr;  ///< Number of rows managed by this process
    int *my_r;  ///< List of rows owned by this process

    int *mpi_r0; ///< array of row offsets for all processes
    int *mpi_nr; ///< array of row numbers for all processes

    MPI_Datatype mpitp; ///< MPI data type corresponding to template type
};


// ==========================================================================
// typedefs for specific instances of `TCompRowMatrixMPI'

typedef TCompRowMatrixMPI<double>   RCompRowMatrixMPI;
typedef TCompRowMatrixMPI<float>    FCompRowMatrixMPI;
typedef TCompRowMatrixMPI<toast::complex> CCompRowMatrixMPI;
typedef TCompRowMatrixMPI<scomplex> SCCompRowMatrixMPI;
typedef TCompRowMatrixMPI<int>      ICompRowMatrixMPI;

// ==========================================================================
// extern declarations of TCompRowMatrix (only required for VS)

#ifndef MATHLIB_IMPLEMENTATION
extern template class MATHLIB TCompRowMatrixMPI<double>;
extern template class MATHLIB TCompRowMatrixMPI<float>;
extern template class MATHLIB TCompRowMatrixMPI<toast::complex>;
extern template class MATHLIB TCompRowMatrixMPI<scomplex>;
extern template class MATHLIB TCompRowMatrixMPI<int>;
#endif // MATHLIB_IMPLEMENTATION

#endif // !__CRMATRIXMPI_H
