// -*-C++-*-
// ==========================================================================
// Module mathlib
// File vector_cusp.h
// Declaration of template class cuspTVector ('CUSP template vector')
// ==========================================================================

#ifndef __VECTOR_CUSP_H
#define __VECTOR_CUSP_H

#include "vector.h"
#include <thrust/extrema.h>
#include <cusp/csr_matrix.h>

// ==========================================================================
// Nonmember declarations

template<class VT> class cuspTVector;

// ==========================================================================
// class cuspTVector

template<class VT> class cuspTVector: public TVector<VT> {
public:
    cuspTVector ();

    cuspTVector (int dim, VT *values);

    // Copy from vector in host memory
    void Set (const TVector<VT> &v);

    // Return vector in host memory format
    TVector<VT> &Get ();

    cuspTVector<VT> &operator= (const cuspTVector<VT> &v);

protected:
    void MapHostToDev ();
    void MapDevToHost ();
    cusp::array1d<VT,cusp::device_memory> dv;
};

#endif // !__VECTOR_CUSP_H
