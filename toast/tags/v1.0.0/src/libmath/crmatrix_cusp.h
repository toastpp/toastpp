// -*-C++-*-
// ==========================================================================
// Module mathlib
// File crmatrix_cusp.h
// Declaration of template class cuspTCompRowMatrix ('CUSP template
// compressed-row matrix')
//
// This matrix format defines a compressed row matrix living on the GPU
// device memory and formatted for the CUSP GPU sparse solver library
//
// Currently only single-precision real and complex types are supported
//
// The following typedefs for specific template types have been defined
// for convenience:
//	cuspFCompRowMatrix = cuspTCompRowMatrix<float>	  ('float')
//	cuspSCCompRowMatrix = cuspTCompRowMatrix<scomplex>  ('single complex')
//
// Inheritance:
// ------------
// TCompRowMatrix ----> cuspTCompRowMatrix
// ==========================================================================

#ifndef __CRMATRIX_CUSP_H
#define __CRMATRIX_CUSP_H

#include "crmatrix.h"
#include <thrust/extrema.h>
#include <cusp/csr_matrix.h>

// ==========================================================================
// class cuspTCompRowMatrix
// ==========================================================================

template<class MT> class cuspTCompRowMatrix: public TCompRowMatrix<MT> {
public:
    cuspTCompRowMatrix ();

    cuspTCompRowMatrix (int rows, int cols,
        const idxtype *_rowptr, const idxtype *_colidx, const MT *data = 0);

    // Set from CPU matrix
    void Set (const TCompRowMatrix<MT> &M);

    // Copy matrix from device memory to CPU matrix
    TCompRowMatrix<MT> &Get () const;

protected:
    void MapHostToDev ();
    void MapDevToHost () const;
    cusp::csr_matrix<int,MT,cusp::device_memory> dA;
};

#endif // !__CRMATRIX_CUSP_H
