// ==========================================================================
// Module mathlib
// File crmatrix_cusp.cc
// Definition of template class cuspTCompRowMatrix ('CUSP template
// compressed-row matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "crmatrix_cusp.h"

// ==========================================================================

template<class MT>
cuspTCompRowMatrix<MT>::cuspTCompRowMatrix (): TCompRowMatrix<MT>()
{
    cusp::csr_matrix<int,MT,cusp::host_memory> A(0,0,0);
    dA = A;
}

// ==========================================================================

template<class MT>
cuspTCompRowMatrix<MT>::cuspTCompRowMatrix (int rows, int cols,
    const idxtype *_rowptr, const idxtype *_colidx, const MT *data)
  : TCompRowMatrix<MT>(rows, cols, _rowptr, _colidx, data)
{
    // Map matrix structure into CUSP format
    MapHostToDev();
}

// ==========================================================================

template<class MT>
void cuspTCompRowMatrix<MT>::Set (const TCompRowMatrix<MT> &M)
{
    TCompRowMatrix<MT>::operator= (M);
    MapHostToDev();
}

// ==========================================================================

template<class MT>
TCompRowMatrix<MT> &cuspTCompRowMatrix<MT>::Get () const
{
    MapDevToHost();
    return *(TCompRowMatrix<MT>*)this;
}

// ==========================================================================

template<class MT>
void cuspTCompRowMatrix<MT>::MapHostToDev ()
{
    int i, nz = this->rowptr[this->rows];
    cusp::csr_matrix<int,MT,cusp::host_memory> A(this->rows,this->cols,nz);
    for (i = 0; i <= this->rows; i++)
	A.row_offsets[i] = this->rowptr[i];
    for (i = 0; i < nz; i++)
	A.column_indices[i] = this->colidx[i];
    for (i = 0; i < nz; i++)
        A.values[i] = this->val[i];
    dA = A;
}

// ==========================================================================

template<class MT>
void cuspTCompRowMatrix<MT>::MapDevToHost () const
{
    int i, nz;
    cusp::csr_matrix<int,MT,cusp::host_memory> A = dA;
    for (i = 0; i <= this->rows; i++)
        this->rowptr[i] = A.row_offsets[i];
    nz = this->rowptr[this->rows];
    for (i = 0; i < nz; i++) {
        this->colidx[i] = A.column_indices[i];
	this->val[i] = A.values[i];
    }
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB cuspTCompRowMatrix<float>;

#endif // NEED_EXPLICIT_INSTANTIATION
