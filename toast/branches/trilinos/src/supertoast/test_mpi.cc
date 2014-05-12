#include "mathlib.h"
#include "dnsmatrix_mpi.h"
#include "timing.h"
#include <iostream> 
#include <iomanip>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    const int n = 1024;

    int mpi_rank, mpi_size;
    int mpi_status = MPI_Init(&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);

    RDenseMatrixMPI Mpar(n,n);

    int i, j;
    RVector tmp(n);
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++)
	    tmp[j] = rand()-RAND_MAX*0.5;
	Mpar.SetRow(i, tmp);
    }

    RVector x(n);
    for (i = 0; i < n; i++)
	x[i] = rand()-RAND_MAX*0.5;

    RVector bpar(n);
    tic();
    for (i = 0; i < 100; i++)
	Mpar.ATx(x,bpar);
    double t_par = toc();

    if (!mpi_rank) {
	double nm = l2norm(bpar);
	cerr << "Norm of result: " << nm << endl;
	cerr << "Time:           " << setprecision(8) << t_par << endl;
    }

    MPI_Finalize();
    return 0;
}
