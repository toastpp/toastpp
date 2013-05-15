#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "toast_mpi.h"

int TMPI::rank = -1;
int TMPI::size = -1;

int TMPI::Rank()
{
    if (rank < 0) 
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int TMPI::Size()
{
    if (size < 0)
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

template <class T>
MPI_Datatype TMPI::MPIType ()
{
    return MPI_INT;
}

template<>
MPI_Datatype TMPI::MPIType<double> ()
{
    return MPI_DOUBLE;
}

template<>
MPI_Datatype TMPI::MPIType<float> ()
{
    return MPI_FLOAT;
}

template<>
MPI_Datatype TMPI::MPIType<toast::complex> ()
{
    return MPI_DOUBLE_COMPLEX;
}

template<>
MPI_Datatype TMPI::MPIType<scomplex> ()
{
    return MPI_COMPLEX;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template MATHLIB MPI_Datatype TMPI::MPIType<int>();

#endif // NEED_EXPLICIT_INSTANTIATION
