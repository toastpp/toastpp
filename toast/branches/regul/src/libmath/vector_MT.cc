// ==========================================================================
// Module mathlib
// File vector_MT.cc
// Multi-threaded versions of TVector<> methods
// Uses pthreads library
// ==========================================================================

#include "mathlib.h"

#if THREAD_LEVEL==1

#define MATHLIB_IMPLEMENTATION

using namespace toast;

// ==========================================================================
// l2normsq

template<class VT>
void l2normsq_engine (void *context, int i0, int i1)
{
    typedef struct {
	const TVector<VT> *v;
	double sum;
    } THDATA;
    THDATA *thdata = (THDATA*)context;

    double term, sum = 0.0;
    const TVector<VT> &v = *thdata->v;
    for (int i = i0; i < i1; i++) {
	term = norm(v[i]);
	sum += term*term;
    }

    g_tpool->LockUserMutex();
    thdata->sum += sum;
    g_tpool->UnlockUserMutex();
}

#ifdef NEED_EXPLICIT_INSTANTIATION
template void l2normsq_engine<double> (void *context, int i0, int i1);
template void l2normsq_engine<float> (void *context, int i0, int i1);
template void l2normsq_engine<complex> (void *context, int i0, int i1);
template void l2normsq_engine<scomplex> (void *context, int i0, int i1);
template void l2normsq_engine<int> (void *context, int i0, int i1);
#endif



template<class VT>
MATHLIB double l2normsq (const TVector<VT> &v)
{
    dASSERT(g_tpool, ThreadPool not initialised);

    int grain = (v.Dim()+NUMTHREAD-1)/NUMTHREAD;
    static struct {
	const TVector<VT> *v;
	double sum;
    } thdata;

    int i;

    thdata.v = &v;
    thdata.sum = 0.0;

    g_tpool->ProcessSequence (l2normsq_engine<VT>, &thdata, 0, v.Dim(), grain);

    return thdata.sum;
}

#ifdef NEED_EXPLICIT_INSTANTIATION

template MATHLIB double   l2normsq (const RVector &vec);
template MATHLIB double   l2normsq (const FVector &vec);
template MATHLIB double   l2normsq (const CVector &vec);
template MATHLIB double   l2normsq (const SCVector &vec);
template MATHLIB double   l2normsq (const IVector &vec);

#endif // NEED_EXPLICIT_INSTANTIATION

#endif // THREAD_LEVEL==1
