#ifndef __TOASTCUDA_H
#define __TOASTCUDA_H

//void cuda_vecadd (float *a, float *b, float *c, int n);

bool cuda_SetDevice (int device);
void cuda_EchoDeviceProperties ();
void cuda_Init (int device);
 
struct SolverResult {
    int it_count;
    double rel_error;
};

typedef enum {
    CUSP_PRECON_IDENTITY,
    CUSP_PRECON_DIAGONAL,
    CUSP_PRECON_AINV,
    CUSP_PRECON_SMOOTHED_AGGREGATION
} CuspPreconType;

void cuda_SetCuspPrecon (CuspPreconType precontp);

// ===========================================================================
// single precision

void cuda_Ax (const float *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const float *x, float *b);

void cuda_Ax_cplx (const scomplex *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const scomplex *x, scomplex *b);

template<class T>
void cuda_CG (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T *b, T *x,
    T tol, int maxit, SolverResult *res = 0);

template<class T>
void cuda_CG (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b, T **x, int nrhs,
    T tol, int maxit, SolverResult *res = 0);

template<class T>
void cuda_BiCGSTAB (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T *b, T *x,
    T tol, int maxit, SolverResult *res = 0);

template<class T>
void cuda_BiCGSTAB (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b, T **x, int nrhs,
    T tol, int maxit, SolverResult *res = 0);

template<class T,class TR>
void cuda_BiCGSTAB_cplx (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T *b, T *x,
    TR tol, int maxit, SolverResult *res = 0);

template<class T,class TR>
void cuda_BiCGSTAB_cplx (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b, T **x,
    int nrhs, TR tol, int maxit, SolverResult *res = 0);

template<class T>
void Tstep_loop (int n, int nq, int nm,
    const T *K0_val, const int *K0_rowptr, const int *K0_colidx,
    const T *K1_val, const int *K1_rowptr, const int *K1_colidx,
    T **qvec_val, T **mvec_val, T *proj,
    T tol, int maxit, int nstep);

#endif // !__TOASTCUDA_H
