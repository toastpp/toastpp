// -*-C++-*-

#ifndef __CRMATRIXTRILINOS_H
#define __CRMATRIXTRILINOS_H

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
//#include <zoltan.h>
#include <iostream>
#include "mathlib.h"

class zoltanMesh;

template<class MT> class TCompRowMatrixTrilinos
: public TGenericSparseMatrix<MT>
{
public:
    typedef Kokkos::DefaultNode::DefaultNodeType node_type;
    typedef Tpetra::CrsMatrix<MT, int, int, node_type> matrix_type;
    typedef Tpetra::Map<int, int, node_type> map_type;
  
    /**
     * \brief Default constructor. Creates a matrix of dimension 0x0 without
     *   an associated map.
     */
    TCompRowMatrixTrilinos ();

    /**
     * \brief Construct a matrix based on a row distribution map and a global
     *   matrix structure.
     * \param rows global number of rows
     * \param cols global number of columns
     * \param _rowptr global row pointer array of matrix elements
     * \param _colidx global column index array of matrix elements
     * \param _map pointer to Tpetra::Map containing the row distribution.
     */
    TCompRowMatrixTrilinos (int rows, int cols, int *_rowptr, int *_colidx,
        Teuchos::RCP<const map_type>_map);

    virtual ~TCompRowMatrixTrilinos ();

    /**
     * \brief Storage class identifier
     */
    MatrixStorage StorageType() const
    { return MATRIX_COMPROW; }

    void SetPartition (Teuchos::RCP<const map_type> _map);

    /**
     * \brief Retrieve a matrix element
     * \param r matrix row (0 <= r < nRows())
     * \param c matrix column (0 <= c < nCols())
     * \return matrix element (*this)<sub>r,c</sub>
     * \note This is a read operation and returns the element value. For
     *   writing operations, use Put().
     */
    MT Get (int r, int c) const
    {
        xERROR("Not implemented");
        return (MT)0; 
    }

    /**
     * \brief Returns a vector containing a copy of row `r'
     * \param r row index (>= 0)
     * \return vector containing row r
     * \note The matrix row is expanded into a full vector, replacing
     *   un-allocated elements with zeros.
     */
    TVector<MT> Row (int r) const
    {
        xERROR("Not implemented");
	return TVector<MT>();
    }

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
    int SparseRow (int r, int *colidx, MT *val) const
    {
        xERROR("Not implemented");
        return 0;
    }

    /**
     * \brief Returns a vector containing a copy of column 'c'
     * \param c column index (>= 0)
     * \return vector containing column c
     * \note Sparse matrix types expand to a dense column, with missing entries
     *   filled with zeros, so this can be used as a "scatter" operation.
     * \sa Row
     */
    TVector<MT> Col (int c) const
    {
        xERROR("Not implemented");
	return TVector<MT>();
    }

    virtual void RowScale (const TVector<MT> &scale)
    {
        xERROR("Not implemented");
    }
    // scales the rows with 'scale'

    virtual void ColScale (const TVector<MT> &scale)
    {
        xERROR("Not implemented");
    }
    // scales the columns with 'scale'

    /**
     * \brief Matrix-vector product.
     * \param[in] x vector argument (length ncols)
     * \param[out] b result vector (length nrows)
     * \note Computes Ax = b
     */
    void Ax (const TVector<MT> &x, TVector<MT> &b) const
    { xERROR("Not implemented"); }

    void Ax (const TVector<MT> &x, TVector<MT> &b, int i1, int i2) const
    { xERROR("Not implemented"); }

    void ATx (const TVector<MT> &x, TVector<MT> &b) const
    { xERROR("Not implemented"); }

     bool Exists (int r, int c) const
    {
        xERROR("Not implemented");
        return false;
    }

    MT &operator() (int r, int c)
    {
        xERROR("Not implemented");
	static MT ret = (MT)0;
	return ret;
    }

    int Get_index (int r, int c) const
    {
        xERROR("Not implemented");
	return 0;
    }

    MT GetNext (int &r, int &c) const
    {
        xERROR("Not implemented");
	return (MT)0;
    }

     /**
     * \brief Return the Teuchos interface of the matrix.
     */
    Teuchos::RCP<matrix_type> GetMatrix ()
    { return teuchosMatrix; }

    void Assemble (zoltanMesh *mesh, void *data,
                   void(*clbkAssemble)(TCompRowMatrix<MT>&,zoltanMesh*,void*));

private:
    idxtype *rowptr;
    idxtype *colidx;
    int numMyElements;
    int numGlobElements;
    size_t numNzero;
    Teuchos::ArrayRCP<size_t>NumNz;
    Teuchos::RCP<const map_type> map;
    Teuchos::RCP<matrix_type> teuchosMatrix;
    Teuchos::ArrayView<const int> myGlobalElements;
    TCompRowMatrix<MT> toastMatrix;
};

#endif // !__TEST_TRILINOS_ZOLTAN_CLASS_H
