#ifndef __TOASTMPI_H
#define __TOASTMPI_H

#include <mpi.h>
#include "toastdef.h"
#include "complex.h"

#ifdef TOAST_MPI
template<class T> class TMPI {
public:
    TMPI() {}
    static MPI_Datatype MPIType();
    static int Rank();
    static int Size();
};
#endif

#ifndef MATHLIB_IMPLEMENTATION
extern template class MATHLIB TMPI<double>;
extern template class MATHLIB TMPI<float>;
extern template class MATHLIB TMPI<toast::complex>;
extern template class MATHLIB TMPI<scomplex>;
extern template class MATHLIB TMPI<int>;
#endif // MATHLIB_IMPLEMENTATION

#endif // !__TOASTMPI_H
