// -*-C++-*-

#include <stdio.h>
#include <iostream>
#include "spmv/spmv_device.h"
#include "spmv/sparse_conversions.h"

// Sparse matrix format. Options: CSR, ELL, HYB
#define MATFORMAT HYB

const int max_ell_nbr = 50; // entries-per-row limit for ELL format
const int hyb_ell_nbr = 30; // ell-entries-per-row limit for HYB format

template<class T>
void Ax_spmv(int nr, int nc, T *Av, int *Ar, int *Ac,
	     const T *x, T *b)
{
    csr_matrix<int,T> csr;
    csr.num_rows = nr;
    csr.num_cols = nc;
    csr.num_nonzeros = Ar[nr];
    csr.Ap = Ar;
    csr.Aj = Ac;
    csr.Ax = Av;

    T *x_loc = copy_array (x, nc, HOST_MEMORY, DEVICE_MEMORY);
    T *b_loc = copy_array (b, nr, HOST_MEMORY, DEVICE_MEMORY);

#if MATFORMAT == CSR

    // map from CSR to ELL format
    ell_matrix<int,T> ell = csr_to_ell(csr,max_ell_nbr);
    if (!ell.num_nonzeros) {
	std::cerr << "ELL format: max entries per row exceeded" << std::endl;
	return;
    }
    ell_matrix<int,T> ell_loc = copy_matrix_to_device (ell);
    spmv_ell_device<int,T>(ell_loc, x_loc, b_loc);
    delete_device_matrix(ell_loc);
    delete_host_matrix(ell);

#elif MATFORMAT == HYB

    // map from CSR to HYB format
    hyb_matrix<int,T> hyb = csr_to_hyb(csr,hyb_ell_nbr);
    hyb_matrix<int,T> hyb_loc = copy_matrix_to_device (hyb);
    spmv_hyb_device<int,T>(hyb_loc, x_loc, b_loc);
    delete_device_matrix(hyb_loc);
    delete_host_matrix(hyb);

#else

    csr_matrix<int,T> csr_loc = copy_matrix_to_device (csr);
    spmv_csr_scalar_device<int,T>(csr_loc, x_loc, b_loc);
    delete_device_matrix(csr_loc);

#endif

    memcpy_to_host (b, b_loc, nr);

    delete_array(x_loc, DEVICE_MEMORY);
    delete_array(b_loc, DEVICE_MEMORY);
}

template<class T>
void Ax_spmv(int nr, int nc, T *Av, int *Ar, int *Ac,
	     const T **x, T **b, int nrhs)
{
    csr_matrix<int,T> csr;
    csr.num_rows = nr;
    csr.num_cols = nc;
    csr.num_nonzeros = Ar[nr];
    csr.Ap = Ar;
    csr.Aj = Ac;
    csr.Ax = Av;

    T *x_loc = copy_array (x[0], nc, HOST_MEMORY, DEVICE_MEMORY);
    T *b_loc = copy_array (b[0], nr, HOST_MEMORY, DEVICE_MEMORY);

#if MATFORMAT == ELL

    // map from CSR to ELL format
    ell_matrix<int,T> ell = csr_to_ell(csr,max_ell_nbr);
    if (!ell.num_nonzeros) {
	std::cerr << "ELL format: max entries per row exceeded" << std::endl;
	return;
    }
    ell_matrix<int,T> ell_loc = copy_matrix_to_device (ell);

#elif MATFORMAT == HYB

    // map from CSR to HYB format
    hyb_matrix<int,T> hyb = csr_to_hyb(csr,hyb_ell_nbr);
    hyb_matrix<int,T> hyb_loc = copy_matrix_to_device (hyb);

#else

    csr_matrix<int,T> csr_loc = copy_matrix_to_device (csr);
 
#endif

    for (int i = 0; i < nrhs; i++) {
#if MATFORMAT == ELL
	spmv_ell_device<int,T>(ell_loc, x_loc, b_loc);
#elif MATFORMAT == HYB
	spmv_hyb_device<int,T>(hyb_loc, x_loc, b_loc);
#else
	spmv_csr_scalar_device<int,T>(csr_loc, x_loc, b_loc);
#endif
	memcpy_to_host (b[i], b_loc, nr);

	if (i < nrhs-1) {
	    memcpy_to_device (x_loc, x[i+1], nc);
	    memcpy_to_device (b_loc, b[i+1], nr);
	}
    }

#if MATFORMAT == ELL
    delete_device_matrix(ell_loc);
#elif MATFORMAT == HYB
    delete_device_matrix(hyb_loc);
#else
    delete_device_matrix(csr_loc);
#endif
    delete_array(x_loc, DEVICE_MEMORY);
    delete_array(b_loc, DEVICE_MEMORY);
}

// ===========================================================

//#ifdef NEED_EXPLICIT_INSTANTIATION

template void Ax_spmv (int nr, int nc, double *Av, int *Ar, int *Ac,
    const double *x, double *b);
template void Ax_spmv (int nr, int nc, float *Av, int *Ar, int *Ac,
    const float *x, float *b);

template void Ax_spmv (int nr, int nc, double *Av, int *Ar, int *Ac,
    const double **x, double **b, int nrhs);
template void Ax_spmv (int nr, int nc, float *Av, int *Ar, int *Ac,
    const float **x, float **b, int nrhs);

//#endif // NEED_EXPLICIT_INSTANTIATION
