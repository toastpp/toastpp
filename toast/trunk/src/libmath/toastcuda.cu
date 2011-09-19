// -*-C++-*-

#include "toastdef.h"
#define MATHLIB DLLEXPORT
#include "complex.h"
#include "scomplex.h"
#include "toastcuda.h"

#include <cuda_runtime.h>
#include <iostream>
#include <cusp/csr_matrix.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/smoothed_aggregation.h>

#define MATHLIB DLLEXPORT

using namespace std;

static CuspPreconType g_precontp = CUSP_PRECON_IDENTITY;

// ===========================================================

void cuda_EchoDeviceProperties ()
{
    struct cudaDeviceProp prop;
    int dev, devno;
    cudaGetDevice (&dev);
    cudaGetDeviceCount (&devno);
    cout << "Device no: " << dev << " of " << devno << endl;
    cudaGetDeviceProperties (&prop, 0);
    cout << "Device: " << prop.name << endl;
    cout << "Global mem: " << prop.totalGlobalMem << endl;
    cout << "Shared mem per block: " << prop.sharedMemPerBlock << endl;
    cout << "Regs per block: " << prop.regsPerBlock << endl;
    cout << "Warp size: " << prop.warpSize << endl;
    cout << "Max mem pitch: " << prop.memPitch << endl;
    cout << "Max threads per block: " << prop.maxThreadsPerBlock << endl;
    cout << "Max threads dim: " << prop.maxThreadsDim[0] << ' '
	 << prop.maxThreadsDim[1] << ' ' << prop.maxThreadsDim[2] << endl;
    cout << "Max grid size: " << prop.maxGridSize[0] << ' '
	 << prop.maxGridSize[1] << ' ' << prop.maxGridSize[2] << endl;
    cout << "Clock rate: " << prop.clockRate << "kHz" << endl;
    cout << "Total const mem: " << prop.totalConstMem << endl;
    cout << "Compute capability: " << prop.major << "." << prop.minor << endl;
    cout << "Texture alignment: " << prop.textureAlignment << endl;
    cout << "Device overlap: " << prop.deviceOverlap << endl;
    cout << "Multiprocessors: " << prop.multiProcessorCount << endl;
    cout << "Kernel exec timeout enabled: " << prop.kernelExecTimeoutEnabled
	 << endl;
    cout << "Integrated device: " << prop.integrated << endl;
    cout << "Can map host memory: " << prop.canMapHostMemory << endl;
    cout << "Compute mode: " << prop.computeMode << endl;
}

// ===========================================================

bool cuda_SetDevice (int device)
{
    return (cudaSetDevice (device) == cudaSuccess);
}

// ===========================================================

void cuda_Init (int device)
{
    cuda_SetDevice (device);
    cuda_EchoDeviceProperties ();
}

// ===========================================================

void cuda_SetCuspPrecon (CuspPreconType precontp)
{
    g_precontp = precontp;
}

// ===========================================================

template<typename ValueType, typename MemorySpace, typename IndexType=int>
class dynamic_preconditioner
: public cusp::linear_operator<ValueType,MemorySpace,IndexType>
{
public:
    template<typename Matrix> dynamic_preconditioner (Matrix &matrix,
        CuspPreconType precontp) : precontp(precontp)
    {
        switch (precontp) {
	case CUSP_PRECON_IDENTITY:
	    precon_ident = new cusp::identity_operator<ValueType,MemorySpace>
	    (matrix.num_rows,matrix.num_cols);
	    break;
	case CUSP_PRECON_DIAGONAL:
	    precon_diag = new cusp::precond::diagonal<ValueType,MemorySpace>
	      (matrix);
	    break;
	case CUSP_PRECON_AINV:
	    precon_ainv = new cusp::precond::bridson_ainv<ValueType,MemorySpace>
	      (matrix, 0.01, 1000, false,1);
	    break;
	case CUSP_PRECON_SMOOTHED_AGGREGATION:
	    precon_smagg = new cusp::precond::smoothed_aggregation<IndexType,
	      ValueType,MemorySpace> (matrix);
	    break;
	}
    }

    ~dynamic_preconditioner ()
    {
        switch (precontp) {
	case CUSP_PRECON_IDENTITY: delete precon_ident; break;
	case CUSP_PRECON_DIAGONAL: delete precon_diag;  break;
	case CUSP_PRECON_AINV:     delete precon_ainv;  break;
	case CUSP_PRECON_SMOOTHED_AGGREGATION:
	    delete precon_smagg; break;
	}
    }

    template<typename VectorType1, typename VectorType2>
    void operator()(const VectorType1 &x, VectorType2 &y) const
    {
        switch (precontp) {
	case CUSP_PRECON_IDENTITY:
	    (*precon_ident)(x,y); break;
	case CUSP_PRECON_DIAGONAL:
	    (*precon_diag)(x,y); break;
	case CUSP_PRECON_AINV:
	    (*precon_ainv)(x,y); break;
	case CUSP_PRECON_SMOOTHED_AGGREGATION:
	    (*precon_smagg)(x,y); break;
	}
    }

private:
    CuspPreconType precontp;
    cusp::identity_operator<ValueType,MemorySpace> *precon_ident;
    cusp::precond::diagonal<ValueType,MemorySpace> *precon_diag;
    cusp::precond::bridson_ainv<ValueType,MemorySpace> *precon_ainv;
    cusp::precond::smoothed_aggregation<IndexType,ValueType,MemorySpace>
      *precon_smagg;
};

// ===========================================================

void cuda_Ax (const float *A_val, const int *A_rowptr, const int *A_colidx,
    int m, int n, const float *x_val, float *b_val)
{
    // Copy toast to cusp matrix format
    int i, nz = A_rowptr[m];
    cusp::csr_matrix<int,float,cusp::host_memory> A(m,n,nz);
    cusp::array1d<float,cusp::host_memory>x(n);

    for (i = 0; i <= m; i++)
	A.row_offsets[i] = A_rowptr[i];
    for (i = 0; i < nz; i++) {
	A.column_indices[i] = A_colidx[i];
	A.values[i] = A_val[i];
    }
    for (i = 0; i < n; i++)
	x[i] = x_val[i];

    // Copy to device memory
    cusp::csr_matrix<int,float,cusp::device_memory>dA = A;
    cusp::array1d<float,cusp::device_memory> dx = x;

    cusp::array1d<float,cusp::device_memory> db(m,0);
    cusp::multiply(dA,dx,db);

    // Copy result back to host memory
    cusp::array1d<float,cusp::host_memory>b = db;

    // Map back to toast format
    for (i = 0; i < m; i++)
	b_val[i] = b[i];
}

// ===========================================================

void cuda_Ax_cplx (const scomplex *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const scomplex *x_val, scomplex *b_val)
{
    // solve this by expanding the complex problem into a real one:
    // | A_re  -A_im | | x_re |     | b_re |
    // | A_im   A_re | | x_im |  =  | b_im |  

    // Copy toast to cusp matrix format
    int i, j, k, p0, p1, np, nz = A_rowptr[m];
    int nz_real = nz*4;
    int m_real = m*2;
    int n_real = n*2;

    cusp::csr_matrix<int,float,cusp::host_memory> A(m_real,n_real,nz_real);
    cusp::array1d<float,cusp::host_memory>x(n_real);

    for (i = 0; i <= m; i++) {
        A.row_offsets[i] = A_rowptr[i]*2;
        A.row_offsets[i+m] = nz*2 + A_rowptr[i]*2;
    }
    for (i = 0; i < m; i++) {
        k = A.row_offsets[i];
	p0 = A_rowptr[i];
	p1 = A_rowptr[i+1];
	np = p1-p0;
        for (j = 0; j < np; j++) {
	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].re;
	    A.values[k+j+np] = -A_val[p0+j].im;
	}
	k = A.row_offsets[i+m];
	for (j = 0; j < np; j++) {
  	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].im;
	    A.values[k+j+np] = A_val[p0+j].re;
	}
    }
    for (i = 0; i < n; i++) {
        x[i] = x_val[i].re;
	x[i+n] = x_val[i].im;
    }

    // Copy to device memory
    cusp::csr_matrix<int,float,cusp::device_memory>dA = A;
    cusp::array1d<float,cusp::device_memory> dx = x;

    cusp::array1d<float,cusp::device_memory> db(m_real,0);
    cusp::multiply(dA,dx,db);

    // Copy result back to host memory
    cusp::array1d<float,cusp::host_memory>b = db;

    // Map back to toast format
    for (i = 0; i < m; i++) {
	b_val[i].re = b[i];
	b_val[i].im = b[i+m];
    }
}

// ===========================================================

template<class T>
void cuda_CG (const T *A_val, const int *A_rowptr,
	      const int *A_colidx, int m, int n,
	      const T *b_val, T *x_val,
	      T tol, int maxit, SolverResult *res)
{
    cerr << "In CG_mono (" << (sizeof(T)==4 ? "single":"double")
	 << ")" << endl;

    if (!maxit) maxit = m+1;

    // Copy toast to cusp matrix format
    int i, nz = A_rowptr[m];
    cusp::csr_matrix<int,T,cusp::host_memory> A(m,n,nz);
    cusp::array1d<T,cusp::host_memory>b(m);

    for (i = 0; i <= m; i++)
	A.row_offsets[i] = A_rowptr[i];
    for (i = 0; i < nz; i++) {
	A.column_indices[i] = A_colidx[i];
	A.values[i] = A_val[i];
    }
    for (i = 0; i < m; i++)
	b[i] = b_val[i];

    // Copy to device memory
    cusp::csr_matrix<int,T,cusp::device_memory>dA = A;

    cusp::array1d<T,cusp::device_memory> dx(n,0);
    cusp::array1d<T,cusp::device_memory> db = b;

    cusp::default_monitor<T> monitor(db,maxit,tol);

    // set preconditioner
    dynamic_preconditioner<T,cusp::device_memory> *M =
        new dynamic_preconditioner<T, cusp::device_memory>(dA, g_precontp);

    cusp::krylov::cg(dA,dx,db,monitor,*M);

    delete M;

    // Copy result back to host memory
    cusp::array1d<T,cusp::host_memory>x = dx;

    // Map back to toast format
    for (i = 0; i < n; i++)
	x_val[i] = x[i];

    if (res) {
        if (monitor.converged()) {
	    res->it_count = monitor.iteration_count();
	} else {
	    res->it_count = maxit;
	}
	res->rel_error = monitor.residual_norm() * tol/monitor.tolerance();
    }
}

// ===========================================================

template<class T>
void cuda_CG (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b_val, T **x_val,
    int nrhs, T tol, int maxit, SolverResult *res)
{
    cerr << "In CG_multi (" << (sizeof(T)==4 ? "single":"double")
	 << ")" << endl;

    if (!maxit) maxit = m+1;
    if (res) res->it_count = 0;

    // Copy toast to cusp matrix format
    int i, iq, nz = A_rowptr[m];
    cusp::csr_matrix<int,T,cusp::host_memory> A(m,n,nz);
    cusp::array1d<T,cusp::host_memory>b(m);

    for (i = 0; i <= m; i++)
	A.row_offsets[i] = A_rowptr[i];
    for (i = 0; i < nz; i++) {
	A.column_indices[i] = A_colidx[i];
	A.values[i] = A_val[i];
    }
    // Copy to device memory
    cusp::csr_matrix<int,T,cusp::device_memory>dA = A;

    // set preconditioner
    dynamic_preconditioner<T,cusp::device_memory> *M =
        new dynamic_preconditioner<T, cusp::device_memory>(dA, g_precontp);

    for (iq = 0; iq < nrhs; iq++) {
        const T *bq_val = b_val[iq];
	T *xq_val = x_val[iq];

	for (i = 0; i < m; i++)
	    b[i] = bq_val[i];

	cusp::array1d<T,cusp::device_memory> dx(n,0);
	cusp::array1d<T,cusp::device_memory> db = b;

	cusp::default_monitor<T> monitor(db,maxit,tol);

	cusp::krylov::cg(dA,dx,db,monitor,*M);

	// Copy result back to host memory
	cusp::array1d<T,cusp::host_memory>x = dx;

	// Map back to toast format
	for (i = 0; i < n; i++)
	    xq_val[i] = x[i];

	if (res) {
  	    int it_count;
	    double rel_error;
	    if (monitor.converged()) {
	        it_count = monitor.iteration_count();
	    } else {
	        it_count = maxit;
	    }
	    rel_error = monitor.residual_norm() * tol/monitor.tolerance();
	    res->it_count += it_count;
	    if (!iq || rel_error > res->rel_error) res->rel_error = rel_error;
	}
    }

    delete M;
}

// ===========================================================

template<class T>
void cuda_BiCGSTAB (const T *A_val, const int *A_rowptr,
		    const int *A_colidx, int m, int n,
                    const T *b_val, T *x_val,
		    T tol, int maxit, SolverResult *res)
{
    cerr << "In BiCGSTAB_mono (" << (sizeof(T)==4 ? "single":"double")
	 << ")" << endl;

    if (!maxit) maxit = m+1;

    // Copy toast to cusp matrix format
    int i, nz = A_rowptr[m];
    cusp::csr_matrix<int,T,cusp::host_memory> A(m,n,nz);
    cusp::array1d<T,cusp::host_memory>b(m);

    for (i = 0; i <= m; i++)
	A.row_offsets[i] = A_rowptr[i];
    for (i = 0; i < nz; i++) {
	A.column_indices[i] = A_colidx[i];
	A.values[i] = A_val[i];
    }
    for (i = 0; i < m; i++)
	b[i] = b_val[i];

    // Copy to device memory
    cusp::csr_matrix<int,T,cusp::device_memory>dA = A;

    cusp::array1d<T,cusp::device_memory> dx(n,0);
    cusp::array1d<T,cusp::device_memory> db = b;

    cusp::default_monitor<T> monitor(db,maxit,tol);

    // set preconditioner
    dynamic_preconditioner<T,cusp::device_memory> *M =
        new dynamic_preconditioner<T, cusp::device_memory>(dA, g_precontp);

    cusp::krylov::bicgstab(dA,dx,db,monitor,*M);

    delete M;

    // Copy result back to host memory
    cusp::array1d<T,cusp::host_memory>x = dx;

    // Map back to toast format
    for (i = 0; i < n; i++)
	x_val[i] = x[i];

    if (res) {
        if (monitor.converged()) {
	    res->it_count = monitor.iteration_count();
	} else {
	    res->it_count = maxit;
	}
	res->rel_error = monitor.residual_norm() * tol/monitor.tolerance();
    }
}

// ===========================================================

template<class T>
void cuda_BiCGSTAB (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b_val, T **x_val,
    int nrhs, T tol, int maxit, SolverResult *res)
{
    cerr << "In BiCGSTAB_multi (" << (sizeof(T)==4 ? "single":"double")
	 << ")" << endl;

    if (!maxit) maxit = m+1;
    if (res) res->it_count = 0;

    // Copy toast to cusp matrix format
    int i, iq, nz = A_rowptr[m];
    cusp::csr_matrix<int,T,cusp::host_memory> A(m,n,nz);
    cusp::array1d<T,cusp::host_memory>b(m);

    for (i = 0; i <= m; i++)
	A.row_offsets[i] = A_rowptr[i];
    for (i = 0; i < nz; i++) {
	A.column_indices[i] = A_colidx[i];
	A.values[i] = A_val[i];
    }
    // Copy to device memory
    cusp::csr_matrix<int,T,cusp::device_memory>dA = A;

    // set preconditioner
    dynamic_preconditioner<T,cusp::device_memory> *M =
        new dynamic_preconditioner<T, cusp::device_memory>(dA, g_precontp);

    for (iq = 0; iq < nrhs; iq++) {
        const T *bq_val = b_val[iq];
	T *xq_val = x_val[iq];

	for (i = 0; i < m; i++)
	    b[i] = bq_val[i];

	cusp::array1d<T,cusp::device_memory> dx(n,0);
	cusp::array1d<T,cusp::device_memory> db = b;

	cusp::default_monitor<T> monitor(db,maxit,tol);

	cusp::krylov::bicgstab(dA,dx,db,monitor,*M);

	// Copy result back to host memory
	cusp::array1d<T,cusp::host_memory>x = dx;

	// Map back to toast format
	for (i = 0; i < n; i++)
	    xq_val[i] = x[i];

	if (res) {
  	    int it_count;
	    double rel_error;
	    if (monitor.converged()) {
	        it_count = monitor.iteration_count();
	    } else {
	        it_count = maxit;
	    }
	    rel_error = monitor.residual_norm() * tol/monitor.tolerance();
	    res->it_count += it_count;
	    if (!iq || rel_error > res->rel_error) res->rel_error = rel_error;
	}
    }

    delete M;
}

// ===========================================================

template<class T,class TR>
void cuda_BiCGSTAB_cplx (const T *A_val, const int *A_rowptr,
		    const int *A_colidx, int m, int n,
                    const T *b_val, T *x_val,
		    TR tol, int maxit, SolverResult *res)
{
    cerr << "In BiCGSTAB_cplx_mono (" << (sizeof(TR)==4 ? "single":"double")
	 << ")" << endl;

    // solve this by expanding the complex problem into a real one:
    // | A_re  -A_im | | x_re |     | b_re |
    // | A_im   A_re | | x_im |  =  | b_im |  

    if (!maxit) maxit = m+1;

    int i, j, k, p0, p1, np;
    int nz = A_rowptr[m];
    int nz_real = nz*4;
    int m_real  = m*2;
    int n_real  = n*2;

    cusp::csr_matrix<int,TR,cusp::host_memory> A(m_real,n_real,nz_real);
    cusp::array1d<TR,cusp::host_memory>b(m_real);

    for (i = 0; i <= m; i++) {
        A.row_offsets[i] = A_rowptr[i]*2;
        A.row_offsets[i+m] = nz*2 + A_rowptr[i]*2;
    }
    for (i = 0; i < m; i++) {
        k = A.row_offsets[i];
	p0 = A_rowptr[i];
	p1 = A_rowptr[i+1];
	np = p1-p0;
        for (j = 0; j < np; j++) {
	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].re;
	    A.values[k+j+np] = -A_val[p0+j].im;
	}
	k = A.row_offsets[i+m];
	for (j = 0; j < np; j++) {
  	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].im;
	    A.values[k+j+np] = A_val[p0+j].re;
	}
    }
    for (i = 0; i < m; i++) {
        b[i] = b_val[i].re;
	b[i+m] = b_val[i].im;
    }

    // Copy to device memory
    cusp::csr_matrix<int,TR,cusp::device_memory>dA = A;

    cusp::array1d<TR,cusp::device_memory> dx(m_real,0);
    cusp::array1d<TR,cusp::device_memory> db = b;

    cusp::default_monitor<TR> monitor(db,maxit,tol);

    // set preconditioner
    dynamic_preconditioner<TR,cusp::device_memory> *M =
        new dynamic_preconditioner<TR, cusp::device_memory>(dA, g_precontp);

    cusp::krylov::bicgstab(dA,dx,db,monitor,*M);

    delete M;

    cusp::array1d<TR,cusp::host_memory>x = dx;
    for (i = 0; i < n; i++) {
	x_val[i].re = x[i];
	x_val[i].im = x[i+n];
    }

    if (res) {
        if (monitor.converged()) {
	    res->it_count = monitor.iteration_count();
	} else {
	    res->it_count = maxit;
	}
	res->rel_error = monitor.residual_norm() * tol/monitor.tolerance();
    }
}

// ===========================================================

template<class T,class TR>
void cuda_BiCGSTAB_cplx (const T *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const T **b_val,
    T **x_val, int nrhs, TR tol, int maxit, SolverResult *res)
{
    cerr << "In BiCGSTAB_cplx_multi (" << (sizeof(TR)==4 ? "single":"double")
	 << ")" << endl;

    // solve this by expanding the complex problem into a real one:
    // | A_re  -A_im | | x_re |     | b_re |
    // | A_im   A_re | | x_im |  =  | b_im |  

    if (!maxit) maxit = m+1;
    if (res) res->it_count = 0;

    int i, j, iq, k, p0, p1, np;
    int nz = A_rowptr[m];
    int nz_real = nz*4;
    int m_real  = m*2;
    int n_real  = n*2;

    cusp::csr_matrix<int,TR,cusp::host_memory> A(m_real,n_real,nz_real);
    cusp::array1d<TR,cusp::host_memory>b(m_real);

    for (i = 0; i <= m; i++) {
        A.row_offsets[i] = A_rowptr[i]*2;
        A.row_offsets[i+m] = nz*2 + A_rowptr[i]*2;
    }
    for (i = 0; i < m; i++) {
        k = A.row_offsets[i];
	p0 = A_rowptr[i];
	p1 = A_rowptr[i+1];
	np = p1-p0;
        for (j = 0; j < np; j++) {
	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].re;
	    A.values[k+j+np] = -A_val[p0+j].im;
	}
	k = A.row_offsets[i+m];
	for (j = 0; j < np; j++) {
  	    A.column_indices[k+j] = A_colidx[p0+j];
	    A.column_indices[k+j+np] = A_colidx[p0+j]+n;
	    A.values[k+j] = A_val[p0+j].im;
	    A.values[k+j+np] = A_val[p0+j].re;
	}
    }

    // Copy to device memory
    cusp::csr_matrix<int,TR,cusp::device_memory>dA = A;

    // set preconditioner
    dynamic_preconditioner<TR,cusp::device_memory> *M =
        new dynamic_preconditioner<TR, cusp::device_memory>(dA, g_precontp);

    for (iq = 0; iq < nrhs; iq++) {
        const T *bq_val = b_val[iq];
	T *xq_val = x_val[iq];

	for (i = 0; i < m; i++) {
	    b[i] = bq_val[i].re;
	    b[i+m] = bq_val[i].im;
	}

	cusp::array1d<TR,cusp::device_memory> dx(m_real,0);
	cusp::array1d<TR,cusp::device_memory> db = b;

	cusp::default_monitor<TR> monitor(db,maxit,tol);

	cusp::krylov::bicgstab(dA,dx,db,monitor,*M);

	// Copy result back to host memory
	cusp::array1d<TR,cusp::host_memory>x = dx;

	// Map back to toast format
	for (i = 0; i < n; i++) {
	    xq_val[i].re = x[i];
	    xq_val[i].im = x[i+n];
	}

	if (res) {
  	    int it_count;
	    double rel_error;
	    if (monitor.converged()) {
	        it_count = monitor.iteration_count();
	    } else {
	        it_count = maxit;
	    }
	    rel_error = monitor.residual_norm() * tol/monitor.tolerance();
	    res->it_count += it_count;
	    if (!iq || rel_error > res->rel_error) res->rel_error = rel_error;
	}
    }

    delete M;
}

// ===========================================================

template<class T>
void Tstep_loop (int n, int nq, int nm,
    const T *K0_val, const int *K0_rowptr, const int *K0_colidx,
    const T *K1_val, const int *K1_rowptr, const int *K1_colidx,
    T **qvec_val, T **mvec_val, T *proj,
    T tol, int maxit, int nstep)
{
    if (!maxit) maxit = n+1;

    // copy toast to cusp matrix format
    int i, j, step, idx, nz = K0_rowptr[n];
    cusp::csr_matrix<int,T,cusp::host_memory> K0(n,n,nz);
    for (i = 0; i <= n; i++)
	K0.row_offsets[i] = K0_rowptr[i];
    for (i = 0; i < nz; i++) {
	K0.column_indices[i] = K0_colidx[i];
	K0.values[i] = K0_val[i];
    }
    
    nz = K1_rowptr[n];
    cusp::csr_matrix<int,T,cusp::host_memory> K1(n,n,nz);
    for (i = 0; i <= n; i++)
	K1.row_offsets[i] = K1_rowptr[i];
    for (i = 0; i < nz; i++) {
	K1.column_indices[i] = K1_colidx[i];
	K1.values[i] = K1_val[i];
    }

    // copy to device memory
    cusp::csr_matrix<int,T,cusp::device_memory>dK0 = K0;
    cusp::csr_matrix<int,T,cusp::device_memory>dK1 = K1;

    // preconditioner
    //cusp::precond::diagonal<T, cusp::device_memory> M(dK1);
    cusp::precond::bridson_ainv<T, cusp::device_memory> M(dK1);

    // copy field vectors
    cusp::array1d<T,cusp::device_memory> ddphi[nq];
    for (i = 0; i < nq; i++) {
	cusp::array1d<T,cusp::host_memory> dphi_tmp(n,0);
	ddphi[i] = dphi_tmp;  // is there a better way to set the size?
    }

    // copy measurement vectors
    cusp::array1d<T,cusp::device_memory> dmvec[nm];
    for (i = 0; i < nm; i++) {
	cusp::array1d<T,cusp::host_memory> mvec_tmp(n);
	for (j = 0; j < n; j++)
	    mvec_tmp[j] = mvec_val[i][j];
	dmvec[i] = mvec_tmp;
    }

    // intermediate source vector
    cusp::array1d<T,cusp::device_memory> qi(n,0);
    idx = 0;

    // initial fields
    for (i = 0; i < nq; i++) {
	cusp::array1d<T,cusp::host_memory>qvec_tmp(n);
	for (j = 0; j < n; j++)
	    qvec_tmp[j] = qvec_val[i][j];
	qi = qvec_tmp;
	cusp::default_monitor<T> monitor(qi,maxit,tol);
	cusp::krylov::bicgstab(dK1,ddphi[i],qi,monitor,M);
	// projection at first time step
	for (j = 0; j < nm; j++)
	    proj[idx++] = cusp::blas::dot (ddphi[i],dmvec[j]);
    }

    for (step = 1; step < nstep; step++) {
	std::cerr << "step " << step << std::endl;
	for (i = 0; i < nq; i++) {
	    cusp::multiply (dK0, ddphi[i], qi);
	    cusp::default_monitor<T> monitor(qi,maxit,tol);
	    cusp::krylov::bicgstab(dK1,ddphi[i],qi,monitor,M);
	    for (j = 0; j < nm; j++)
	    	proj[idx++] = cusp::blas::dot (ddphi[i],dmvec[j]);
	}
    }
}

// ===========================================================

//#ifdef NEED_EXPLICIT_INSTANTIATION

// Single-precision instantiations
template void cuda_CG (const float *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const float *b, float *x,
    float tol, int maxit, SolverResult *res);

template void cuda_CG (const float *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const float **b, float **x, int nrhs,
    float tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB (const float *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const float *b, float *x,
    float tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB (const float *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const float **b, float **x, int nrhs,
    float tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB_cplx (const scomplex *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const scomplex *b, scomplex *x,
    float tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB_cplx (const scomplex *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const scomplex **b, scomplex **x,
    int nrhs, float tol, int maxit, SolverResult *res);

// Double-precision instantiations
// NOTE: This works only with compute capability 1.3 or higher
// nvcc -arch=sm_13
// nvcc -arch=sm_20

template void cuda_CG (const double *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const double *b, double *x,
    double tol, int maxit, SolverResult *res);

template void cuda_CG (const double *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const double **b, double **x, int nrhs,
    double tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB (const double *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const double *b, double *x,
    double tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB (const double *A_val, const int *A_rowptr,
    const int *A_colidx, int m, int n, const double **b, double **x, int nrhs,
    double tol, int maxit, SolverResult *res);

template void cuda_BiCGSTAB_cplx (const toast::complex *A_val,
    const int *A_rowptr, const int *A_colidx, int m, int n,
    const toast::complex *b, toast::complex *x, double tol, int maxit,
    SolverResult *res);

template void cuda_BiCGSTAB_cplx (const toast::complex *A_val,
    const int *A_rowptr, const int *A_colidx, int m, int n,
    const toast::complex **b, toast::complex **x, int nrhs, double tol,
    int maxit, SolverResult *res);

template void Tstep_loop (int n, int nq, int nm,
    const double *K0_val, const int *K0_rowptr, const int *K0_colidx,
    const double *K1_val, const int *K1_rowptr, const int *K1_colidx,
    double **dphi_val, double **mvec, double *proj,
    double tol, int maxit, int nstep);

//#endif // NEED_EXPLICIT_INSTANTIATION
