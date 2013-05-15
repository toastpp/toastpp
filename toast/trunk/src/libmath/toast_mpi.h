#ifndef __TOASTMPI_H
#define __TOASTMPI_H

#include <mpi.h>
#include "toastdef.h"
#include "complex.h"

#ifdef TOAST_MPI
class TMPI {
public:
    TMPI() {}
    static int Rank();
    static int Size();
    template<class T> static MPI_Datatype MPIType ();

private:
    static int rank;
    static int size;
};
#endif

#endif // !__TOASTMPI_H
