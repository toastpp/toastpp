#include "mathlib.h"
#include "vector_cusp.h"

// ==========================================================================

template<class VT>
cuspTVector<VT>::cuspTVector (): TVector<VT>()
{
    cusp::array1d<VT,cusp::host_memory>v(0);
    dv = v;
}

// ==========================================================================

template<class VT>
cuspTVector<VT>::cuspTVector (int dim, VT *values)
  : TVector<VT>(dim)
{
    // Map vector into CUSP format
    int i;
    cusp::array1d<VT,cusp::host_memory>v(dim);
    for (i = 0; i < dim; i++)
        v[i] = values[i];
    
    // Copy to device memory
    dv = v;
}

// ==========================================================================

template<class VT>
void cuspTVector<VT>::Set (const TVector<VT> &v)
{
    TVector<VT>::operator= (v);
    MapHostToDev();
}

// ==========================================================================

template<class VT>
TVector<VT> &cuspTVector<VT>::Get ()
{
    MapDevToHost();
    return *(TVector<VT>*)this;
}

// ==========================================================================

template<class VT>
cuspTVector<VT> &cuspTVector<VT>::operator= (const cuspTVector<VT> &v)
{
    dv = v.dv;
    return *this;
}

// ==========================================================================

template<class VT>
void cuspTVector<VT>::MapHostToDev ()
{
    int i;
    cusp::array1d<VT,cusp::host_memory> v(this->size);
    for (i = 0; i < this->size; i++)
        v[i] = this->data[i];
    dv = v;
}

// ==========================================================================

template<class VT>
void cuspTVector<VT>::MapDevToHost ()
{
    int i;
    cusp::array1d<VT,cusp::host_memory> v = dv;
    if (v.size() != this->size) New(v.size());
    for (i = 0; i < this->size; i++)
        this->data[i] = v[i];    
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB cuspTVector<float>;

#endif // NEED_EXPLICIT_INSTANTIATION
