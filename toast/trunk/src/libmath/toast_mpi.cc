#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "toast_mpi.h"

template <class T>
int TMPI<T>::Rank()
{
    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
    return rnk;
}

template <class T>
int TMPI<T>::Size()
{
    int sze;
    MPI_Comm_size(MPI_COMM_WORLD, &sze);
    return sze;
}

template <class T>
MPI_Datatype TMPI<T>::MPIType()
{
    return MPI_INT;
}

template<>
MPI_Datatype TMPI<double>::MPIType()
{
    return MPI_DOUBLE;
}

template<>
MPI_Datatype TMPI<toast::complex>::MPIType()
{
    return MPI_DOUBLE_COMPLEX;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TMPI<double>;
template class MATHLIB TMPI<float>;
template class MATHLIB TMPI<toast::complex>;
template class MATHLIB TMPI<scomplex>;
template class MATHLIB TMPI<int>;

#endif // NEED_EXPLICIT_INSTANTIATION
