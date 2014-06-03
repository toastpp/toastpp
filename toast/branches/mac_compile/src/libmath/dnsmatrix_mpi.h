// -*-C++-*-
// ==========================================================================
// General dense matrix class
// Distributed (MPI) version
// Each MPI process operates on a block of matrix rows (r0 <= r < r1)
// ==========================================================================

#ifndef __DNSMATRIXMPI_H
#define __DNSMATRIXMPI_H

#include "mathlib.h"

// ==========================================================================
// Nonmember declarations
// ==========================================================================

// ==========================================================================
// class TDenseMatrixMPI
// ==========================================================================
/**
 * \brief Distributed dense matrix class.
 *
 * A dense matrix type that distributes itself across the available processes,
 * by assigning blocks of rows of equal size to each process.
 *
 * Each process only allocates the memory required to store the associated
 * sub-block. Methods such as matrix x vector (Ax) and transpose(matrix) x
 * vector (ATx) are automaticall performed in parallel, without the need
 * for any special code by the application developer.
 */

template<class MT> class TDenseMatrixMPI: public TMatrix<MT> {

public:
    /**
     * \brief Defines a matrix of size 0 x 0.
     */
    TDenseMatrixMPI ();

    /**
     * \brief Defines a matrix of size r x c, and allocates space for a
     *   sub-set of rows.
     * \param r number of rows
     * \param c number of columns
     * \note This constructor automatically subdivides the matrix according
     *   to the process rank and MPI size.
     */
    TDenseMatrixMPI (int r, int c);

    /**
     * \Destructor.
     */
    ~TDenseMatrixMPI ();

    inline MatrixStorage StorageType () const { return MATRIX_DENSE; }
    // Matrix element storage method

    /**
     * \brief Resize matrix to new dimensions.
     * \param r new row dimension
     * \param c new column dimension
     * \note All matrix entries are reset to zero by this method.
     * \note This is a blocking function: it must be called for all 
     *   processes.
     */
    void New (int r, int c);

    /**
     * \brief Returns the row range of the current process.
     * \param [out] _r0 low row index (>= 0)
     * \param [out] _r1 high row index+1
     * \note This method returns the range of rows the current process is
     *   operating on.
     */
    void MPIRange (int *_r0, int *_r1) const { *_r0 = r0; *_r1 = r1; }

    MT Get (int r, int c) const;

    MT &operator() (int r, int c);

    /**
     * \brief Returns a row of the matrix.
     * \param r row index (>= 0)
     * \return Row vector. The size of the returned vector equals the number
     *   of columns of the matrix.
     * \note This is a blocking method. The process responsible for row r
     *   broadcasts its value to all other processes. Each process then
     *   returns the same row vector.
     */
    TVector<MT> Row (int r) const;

    TVector<MT> Col (int c) const;

    int SparseRow (int r, int *ci, MT *rv) const;

    /**
     * \brief Replace entries in row r with the values provided in rval.
     * \param r row number (>= 0)
     * \param rval replacement row
     * \note The length of rval must be equal to the number of columns in the
     *   matrix.
     * \note This method can be called for all processes, but it affects only
     *   the process responsible for row r.
     */
    void SetRow (int r, const TVector<MT> &rval);

    void ColScale (const TVector<MT> &scale);

    void RowScale (const TVector<MT> &scale);

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

    /**
     * \brief Matrix x vector product.
     * \param [in] x Right-hand operand.
     * \param [out] b Result vector
     * \note The size of x must equal the number of columns of the matrix.
     * \note The size of b must equal the number of rows of the matrix.
     * \note This is a synchronising method: The result b is collected from
     *   all contributions, and re-distributed to all processes.
     * \sa ATx
     */
    void Ax (const TVector<MT> &x, TVector<MT> &b) const;

    /**
     * \brief Matrix-transpose x vector product.
     * \param [in] x Right-hand operand.
     * \param [out] b Result vector
     * \note The size of x must equal the number of rows of the matrix
     * \note The size of b must equal the number of columns of the matrix
     * \note This is a synchronising method: The result b is collected from
     *   all contributions, and re-distributed to all processes.
     * \sa Ax
     */
    void ATx (const TVector<MT> &x, TVector<MT> &b) const;

    /**
     * \brief Returns MPI type compatible with template type
     */
    MPI_Datatype MPIType () const;

private:
    void Unlink ();
    void Alloc (int r, int c);

    MT *val;    // data array (only contains relevant sub-block)
    int r0, r1; // row index range
    int nr;     // number of rows for process
    int nr_max; // max number of rows for any process
    int rnk, sze; // MPI rank of process, MPI size
};


// ==========================================================================
// typedefs for specific instances of `TDenseMatrixMPI'

typedef TDenseMatrixMPI<double>   RDenseMatrixMPI;	// 'double real'
typedef TDenseMatrixMPI<float>    FDenseMatrixMPI;    // 'single real'
typedef TDenseMatrixMPI<toast::complex>  CDenseMatrixMPI;	// 'complex'
typedef TDenseMatrixMPI<scomplex> SCDenseMatrixMPI;   // 'single complex'
typedef TDenseMatrixMPI<int>      IDenseMatrixMPI;    // 'integer'

#ifndef MATHLIB_IMPLEMENTATION
extern template class MATHLIB TDenseMatrixMPI<double>;
extern template class MATHLIB TDenseMatrixMPI<float>;
extern template class MATHLIB TDenseMatrixMPI<toast::complex>;
extern template class MATHLIB TDenseMatrixMPI<scomplex>;
extern template class MATHLIB TDenseMatrixMPI<int>;
#endif // MATHLIB_IMPLEMENTATION

#endif // !__DNSMATRIXMPI_H
