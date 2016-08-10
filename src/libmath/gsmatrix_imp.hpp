// ==========================================================================
// Module mathlib
// File gsmatrix.cc
// Definition of template class TGenericSparseMatrix
// ==========================================================================

#define __GSMATRIX_CC
#define MATHLIB_IMPLEMENTATION
#define DBG_TIMING

#include "mathlib.h"
#include "timing.h"

#ifdef TOAST_PARALLEL
#include "task.h"
#endif

using namespace std;

// ==========================================================================
// QR factorisation and solver
// ==========================================================================

template<class MT>
int QRFactorize (TGenericSparseMatrix<MT> &A, TVector<MT> &x,
    TVector<MT> &d)
{
    int n = A.nRows();
    int m = A.nCols();
    int i, j, k;
    double sum;
    TVector<MT> Ak(n), Aj(n);
    TCompRowMatrix<MT> *A_cr;
    MT b, c;
    if (A.StorageType() == MATRIX_COMPROW)
      A_cr = (TCompRowMatrix<MT>*)&A;
    else A_cr = 0;

    A.Transpone(); // can be skipped if A symmetric

    for (k = 0; k < m; k++) {
        Ak = A.Row(k); // scatter
        for (sum = 0, i = k; i < n; i++)
	    sum += Ak[i]*Ak[i];
	d[k] = (Ak[k] < 0 ? -sqrt(sum) : sqrt(sum));
	b = sqrt (2.0*d[k] * (Ak[k] + d[k]));

	Ak[k] = (Ak[k]+d[k])/b;
	for (i = k+1; i < n; i++) Ak[i] /= b;
	for (j = k+1; j < m; j++) {
	    Aj = A.Row(j); // scatter
	    if (A_cr)
	        sum = SparseDotp (Ak, A_cr->colidx+A_cr->rowptr[k],
				  A_cr->rowptr[k+1]-A_cr->rowptr[k],
				  Aj, A_cr->colidx+A_cr->rowptr[j],
				  A_cr->rowptr[j+1]-A_cr->rowptr[j],
				  k, n-1);
	    else
	        for (sum = 0, i = k; i < n; i++)
		    sum += Ak[i]*Aj[i];
	    c = 2.0*sum;
	    for (i = k; i < n; i++)
	        Aj[i] -= c*Ak[i];
	    A.SetRow (j, Aj); // gather
	}
	A.SetRow (k, Ak); // gather
    }
    return 0;
}

template<class MT>
void RSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &d,
    TVector<MT> &b)
{
    int i, j;
    int m = A.nRows(); // note that we work with the transpose of A_factor
    double sum;
    b[m-1] /= -d[m-1];
    for (i = m-2; i >= 0; i--) {
        for (sum = 0.0, j = i+1; j < m; j++)
	    sum += A.Get(j,i) * b[j];
	b[i] = (b[i]-sum) / -d[i];
    }
}

template<class MT>
void QRSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &c,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x)
{
    int n = A.nCols(); // note that we work with the transpose of A_factor
    int m = A.nRows();
    int i, j;
    double sum;
    MT f;

    x = b;
    for (j = 0; j < m; j++) {
        for (sum = 0.0, i = j; i < n; i++)
	    sum += A.Get(j,i) * x[i];
	f = 2.0*sum;
	for (i = j; i < n; i++)
	    x[i] -= f * A.Get(j,i);
    }
    RSolve (A, d, x);
}

// ==========================================================================
// Conjugate gradient (CG) solver
// ==========================================================================

#ifndef CG_PARALLEL // **** serial implementation ****

template<class MT>
int CG (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    TIC; // timing

    dASSERT(A.nRows() == A.nCols(), "Matrix not square");
    dASSERT(b.Dim() == A.nRows(), "Dimension mismatch");
    dASSERT(x.Dim() == A.nRows(), "Dimension mismatch");

    double dnew, dold, alpha, beta;
    double err, bnorm;
    int niter, dim = x.Dim();
    if (!maxit) maxit = dim+1;
    TVector<MT> r(dim), d(dim), q(dim);
    r = b - (A*x);
    if (precon) precon->Apply (r, d);
    else d = r;
    dnew = r & d;
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        A.Ax (d, q);
	alpha = dnew / (d & q);
	x += d * alpha;
	r -= q * alpha;
	dold = dnew;
	if (precon) {
	    precon->Apply (r, q);
	    dnew = r & q;
	    beta = dnew/dold;
	    d *= beta;
	    d += q;
	} else {
	    dnew = r & r;
	    beta = dnew/dold;
	    d *= beta;
	    d += r;
	}
	// check convergence
	err = l2norm (r);
	//cerr << niter << " " << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;

    TOCADD(cgtime);

    IterCount += niter;
    return niter;
}

#else // **** threaded implementation **** : not currently preconditioned

template<class MT>
struct cgdata1 {
    const TGenericSparseMatrix<MT> *A;
    TVector<MT> *d;
    TVector<MT> *q;
    double *dq;
    pthread_mutex_t *mutexp;
};

template<class MT>
struct cgdata2 {
    TVector<MT> *x;
    TVector<MT> *r;
    TVector<MT> *d;
    TVector<MT> *q;
    double *alpha;
    double *dnew;
    pthread_mutex_t *mutexp;
};

template<class MT>
struct cgdata3 {
    TVector<MT> *d;
    TVector<MT> *r;
    double *beta;
};

template<class MT>
void TGenericSparseMatrix<MT>::cg_loop1 (void *arg, int i1, int i2)
{
    register int i;
    struct cgdata1<MT> *data = (struct cgdata1<MT>*)arg;
    TVector<MT> &d = *data->d;
    TVector<MT> &q = *data->q;
    data->A->Ax(d, q, i1, i2);
    double local_dq = 0.0;
    for (i = i1; i < i2; i++)
        local_dq += d[i] * q[i];
    pthread_mutex_lock (data->mutexp);
    *data->dq += local_dq;
    pthread_mutex_unlock (data->mutexp);
}

template<class MT>
void TGenericSparseMatrix<MT>::cg_loop2 (void *arg, int i1, int i2)
{
    struct cgdata2<MT> *data = (struct cgdata2<MT>*)arg;
    register int i;
    double local_dnew = 0.0;
    TVector<MT> &x = *data->x;
    TVector<MT> &r = *data->r;
    TVector<MT> &d = *data->d;
    TVector<MT> &q = *data->q;
    double alpha = *data->alpha;
    for (i = i1; i < i2; i++) {
        x[i] += d[i] * alpha;
        r[i] -= q[i] * alpha;
	local_dnew += r[i] * r[i];
    }
    pthread_mutex_lock (data->mutexp);
    *data->dnew += local_dnew;
    pthread_mutex_unlock (data->mutexp);
}

template<class MT>
void TGenericSparseMatrix<MT>::cg_loop3 (void *arg, int i1, int i2)
{
    struct cgdata3<MT> *data = (struct cgdata3<MT>*)arg;
    register int i;
    TVector<MT> &d = *data->d;
    TVector<MT> &r = *data->r;
    double beta = *data->beta;
    for (i = i1; i < i2; i++) {
        d[i] *= beta;
	d[i] += r[i];
    }
}

template<class MT>
int CG (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    dASSERT(A.nRows() == A.nCols(), "Matrix not square");
    dASSERT(b.Dim() == A.nRows(), "Dimension mismatch");
    dASSERT(x.Dim() == A.nRows(), "Dimension mismatch");

    double dnew, dold, d0, alpha, beta = 0.0, dq;
    int niter, dim = x.Dim();
    if (!maxit) maxit = dim+1;
    TVector<MT> r(dim), d(dim), q(dim);
    r = b - (A*x);
    d = r;
    dnew = r & d;
    d0 = tol*tol * dnew;

    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    int grain = (dim+Task::GetThreadCount()-1)/Task::GetThreadCount();
    struct cgdata1<MT> data1 = {&A, &d, &q, &dq, &mutex};
    struct cgdata2<MT> data2 = {&x, &r, &d, &q, &alpha, &dnew, &mutex};
    struct cgdata3<MT> data3 = {&d, &r, &beta};

    for (niter = 0; niter < maxit && dnew > d0; niter++) {
        dq = 0.0;
	g_tpool->ProcessSequence (TGenericSparseMatrix<MT>::cg_loop1,
				  &data1, 0, dim, grain);
	alpha = dnew/dq;
        dold  = dnew;
	dnew  = 0.0;
	cerr << "Starting loop2" << endl;
	g_tpool->ProcessSequence (TGenericSparseMatrix<MT>::cg_loop2,
				  &data2, 0, dim, grain);
	cerr << "Finished loop2" << endl;
	beta = dnew/dold;
	g_tpool->ProcessSequence (TGenericSparseMatrix<MT>::cg_loop3,
				  &data3, 0, dim, grain);
    }
    tol *= sqrt (dnew/d0);
    return niter;
}

#endif // CG_PARALLEL

template<> // specialisation: complex
inline int CG<std::complex<double> > (const CGenericSparseMatrix &A,
    const CVector &b, CVector &x, double &tol, CPreconditioner *cprecon,
    int maxit)
{
    // NOT IMPLEMENTED YET
    return 0;
}

template<> // specialisation: complex
inline int CG<std::complex<float> > (const SCGenericSparseMatrix &A,
    const SCVector &b, SCVector &x, double &tol, SCPreconditioner *cprecon,
    int maxit)
{
    // NOT IMPLEMENTED YET
    return 0;
}

// ==========================================================================
// Preconditioned bi-conjugate gradient (BiCG) solver
// from: netlib (http://www.netlib.org/linalg/html_templates/node32.html)
// ==========================================================================

template<class MT>
int BiCG (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    double rho=0, alpha, aden, beta, bden, err, bnorm;
    TVector<MT> r(dim), rd(dim), z(dim), zd(dim), p(dim), pd(dim), q(dim), qd(dim);
    if (!maxit) maxit = dim+1;
    r = rd = b - (A*x);
    bnorm = l2norm(b);

    for (niter = 0; niter < maxit; niter++) {
        if (precon) {
	    precon->Apply (r, z);
	    precon->Apply (rd, zd);
	    // note this assumes symmetric preconditioner!
	    // generally we need to solve M^T zd = rd, rather than M zd = rd
	} else {
	    z = r;
	    zd = rd;
	}

	bden = rho;
	rho = z & rd;
	xASSERT (rho != 0.0, "BiCG solver fails");
	if (!niter) {
	    p = z;
	    pd = zd;
	} else {
	    beta = rho/bden;
	    for (i = 0; i < dim; i++) {
	        p[i]  *= beta; p[i]  += z[i];
		pd[i] *= beta; pd[i] += zd[i];
	    }
	}
	A.Ax (p, q);
	A.ATx (pd, qd);
	aden = pd & q;
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    x[i]  += alpha * p[i];
	    r[i]  -= alpha * q[i];
	    rd[i] -= alpha * qd[i];
	}

	// check convergence
	err = l2norm(r);
	//cerr << "BiCG tol=" << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

template<> // specialisation: complex
inline int BiCG<std::complex<double> > (const CGenericSparseMatrix &A,
    const CVector &b, CVector &x, double &tol, CPreconditioner *cprecon,
    int maxit)
{
    int i, niter, dim = x.Dim();
    double rho=0, alpha, aden, beta, bden, err, bnorm;
    if (!maxit) maxit = dim+1;
    CVector r(dim), rd(dim), z(dim), zd(dim), p(dim), pd(dim), q(dim), qd(dim);
    r = rd = b - (A*x);
    bnorm = l2norm(b);

    for (niter = 0; niter < maxit; niter++) {
        if (cprecon) {
	    cprecon->Apply (r, z);
	    cprecon->Apply (rd, zd);
	    // note this assumes symmetric preconditioner!
	    // generally we need to solve M^T zd = rd, rather than M zd = rd
	} else {
	    z = r;
	    zd = rd;
	}
	bden = rho;
	for (rho = 0.0, i = 0; i < dim; i++)
  	    rho += z[i].real()*rd[i].real() + z[i].imag()*rd[i].imag();
	xASSERT(rho != 0.0, "BiCG solver fails");
	if (!niter) {
	    p = z;
	    pd = zd;
	} else {
	    beta = rho/bden;
	    for (i = 0; i < dim; i++) {
	        p[i]  *= beta;  p[i]  += z[i];
		pd[i] *= beta;  pd[i] += zd[i];
	    }
	}
	A.Ax (p, q);
	A.ATx (pd, qd);
	for (aden = 0.0, i = 0; i < dim; i++)
	    aden += pd[i].real()*q[i].real() + pd[i].imag()*q[i].imag();
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    x[i]  += p[i]*alpha;
	    r[i]  -= q[i]*alpha;
	    rd[i] -= qd[i]*alpha;
	}

	// check convergence
	for (err = 0.0, i = 0; i < dim; i++)
	    err += std::norm(r[i]);
	err = sqrt (err);
	//cerr << "BiCG tol=" << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

template<class MT>
int GaussSeidel (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, int maxit)
{
    int i, j, k, nz, niter;
    idxtype *colidx = new idxtype[A.nCols()];
    MT *val = new MT[A.nCols()];
    MT aii;
    double bnorm = l2norm(b);

    for (niter = 0; niter < maxit; niter++) {
        TVector<MT> xp(x);
        for (i = 0; i < A.nRows(); i++) {
	    x[i] = b[i];
	    aii = (MT)0;
	    nz = A.SparseRow (i, colidx, val);
	    for (k = 0; k < nz; k++) {
	        j = colidx[k];
		if (j == i) aii = val[k];
		else x[i] -= val[k]*x[j];
	    }
	    if (aii != (MT)0) x[i] /= aii;
	    else xERROR("Zero diagonal element in matrix");
	}
	// stopping criteria
	if (l2norm (x-xp) / l2norm(x) < tol) break; // no improvement
	if (l2norm (b-A*x) / bnorm < tol) break;   // convergence
    }
    return niter;
}

// ==========================================================================
// iterative linear sparse solver
// ==========================================================================

template<class MT>
int IterativeSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int niter;

    switch (itmethod_general) {
    case ITMETHOD_CG:
        niter = A.pcg (b, x, tol, precon, maxit);
        //niter = CG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICG:
        niter = BiCG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICGSTAB:
        niter = A.bicgstab (b, x, tol, precon, maxit);
        //niter = BiCGSTAB (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_GMRES:
        niter = GMRES (A, b, x, tol, precon, 20, maxit, 0);
	break;
    case ITMETHOD_GAUSSSEIDEL:
        niter = GaussSeidel (A, b, x, tol, maxit);
	break;
    default:
        ERROR_UNDEF;
    }
    return niter;
}

template<> // specialisation for single complex case
inline int IterativeSolve (const SCGenericSparseMatrix &A, const SCVector &b,
    SCVector &x, double &tol, SCPreconditioner *precon, int maxit)
{
    int niter = 0;

    switch (itmethod_complex) {
    case ITMETHOD_CG:
        niter = A.pcg (b, x, tol, precon, maxit);
        //niter = CG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICG:
        ERROR_UNDEF;
        //niter = BiCG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICGSTAB:
        niter = A.bicgstab (b, x, tol, precon, maxit);
	break;
    case ITMETHOD_GMRES:
        ERROR_UNDEF;
        //niter = GMRES (A, b, x, tol, precon, 10, 0);
	break;
    case ITMETHOD_GAUSSSEIDEL:
        ERROR_UNDEF;
        //niter = GaussSeidel (A, b, x, tol, maxit);
	break;
    default:
        ERROR_UNDEF;
    }
    return niter;
}

template<> // specialisation for complex case
inline int IterativeSolve<std::complex<double> > (const CGenericSparseMatrix &A,
    const CVector &b, CVector &x, double &tol, CPreconditioner *precon,
    int maxit)
{
    int niter;

    switch (itmethod_complex) {
    case ITMETHOD_CG:
        niter = CG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICG:
        niter = BiCG (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_BICGSTAB:
        niter = BiCGSTAB (A, b, x, tol, precon, maxit);
	break;
    case ITMETHOD_GMRES:
        niter = GMRES (A, b, x, tol, precon, 20, maxit, 0);
	break;
    default:
        ERROR_UNDEF;
    }
    return niter;
}

// ==========================================================================
// multiple right-hand sides

template<class MT>
void IterativeSolve (const TGenericSparseMatrix<MT> &A,
    const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit,
    TPreconditioner<MT> *precon, IterativeSolverResult *res)
{
    switch (itmethod_general) {
    case ITMETHOD_CG:
        A.pcg (b, x, nrhs, tol, maxit, precon, res);
	break;
    case ITMETHOD_BICGSTAB:
        A.bicgstab (b, x, nrhs, tol, maxit, precon, res);
	break;
    default:
        ERROR_UNDEF;
    }
}

template<> // specialisation: single complex
inline void IterativeSolve (const SCGenericSparseMatrix &A, const SCVector *b,
    SCVector *x, int nrhs, double tol, int maxit,
    SCPreconditioner *precon, IterativeSolverResult *res)
{
    switch (itmethod_complex) {
    case ITMETHOD_BICGSTAB:
        A.bicgstab (b, x, nrhs, tol, maxit, precon, res);
	break;
    default:
        ERROR_UNDEF;
    }
}

template<> // specialisation: complex
inline void IterativeSolve (const CGenericSparseMatrix &A, const CVector *b,
    CVector *x, int nrhs, double tol, int maxit, CPreconditioner *precon,
    IterativeSolverResult *res)
{
    switch (itmethod_complex) {
    case ITMETHOD_BICGSTAB:
        A.bicgstab (b, x, nrhs, tol, maxit, precon, res);
	break;
    case ITMETHOD_GMRES:
        for (int i = 0; i < nrhs; i++) {
  	    double err = tol;
	    GMRES (A, b[i], x[i], err, precon, 20, maxit, 0);
	}
	break;
    default:
        ERROR_UNDEF;
    }
}

// ==========================================================================
// OBSOLETE

template<class MT>
int ComplexBiCGSolve (const TGenericSparseMatrix<MT>& Are,
    const TGenericSparseMatrix<MT>& Aim, const TVector<MT>& bre,
    const TVector<MT>& bim, TVector<MT>& xre, TVector<MT>& xim,
    double &tol, int maxit)
{
    dASSERT(Are.nRows() == Aim.nRows() && Are.nCols() == Aim.nCols(),
	   "Matrix dimensions differ.");
    dASSERT(bre.Dim() == bim.Dim() && xre.Dim() == xim.Dim(),
	   "Vector dimensions differ.");
    dASSERT(Are.nRows() == bre.Dim(), "Matrix and vector incompatible.");

    const double EPS = 1e-10;
    int i, k = 0, dim = bre.Dim();
    double bknum, bkden, akden, alpha, err;
    double bnorm = sqrt ((bre & bre) + (bim & bim));
    TVector<MT> pre(dim), pim(dim);
    TVector<MT> pdre(dim), pdim(dim);
    if (!maxit) maxit = dim+1;
    xre = 0.0; xim = 0.0;
    TVector<MT> rre = bre - (Are*xre - Aim*xim);
    TVector<MT> rim = bim - (Aim*xre + Are*xim);
    TVector<MT> rdre = rre;
    TVector<MT> rdim = rim;
    TVector<MT> preconre = Are.Diag();
    TVector<MT> preconim = Are.Diag();
    for (i=0; i<dim; i++) {
	preconre[i] = (fabs (preconre[i]) < EPS ? 1.0 : 1.0 / preconre[i]);
	preconim[i] = (fabs (preconim[i]) < EPS ? 1.0 : 1.0 / preconim[i]);
    }
    TVector<MT> zre = rre * preconre;
    TVector<MT> zim = rim * preconim;
    TVector<MT> zdre = zre;
    TVector<MT> zdim = zim;
    do {
	bknum = (zre & rdre) + (zim & rdim);
	pre = (k ? pre * (bknum/bkden) + zre : zre);
	pim = (k ? pim * (bknum/bkden) + zim : zim);
	pdre = (k ? pdre * (bknum/bkden) + zdre : zdre);
	pdim = (k ? pdim * (bknum/bkden) + zdim : zdim);
	bkden = bknum;
	zre =  Are*pre - Aim*pim;	// z = A * p
	zim =  Aim*pre + Are*pim;
	zdre = Are*pdre + Aim*pdim;	// zd = A^T * p
	zdim = Are*pdim - Aim*pdre;
	akden = (zre & pdre) + (zim & pdim);
	alpha = bknum / akden;
	for (i=0; i<dim; i++) {
	    xre[i] += pre[i] * alpha;
	    xim[i] += pim[i] * alpha;
	    rre[i] -= zre[i] * alpha;
	    rim[i] -= zim[i] * alpha;
	    zre[i]  = rre[i] * preconre[i];
	    zim[i]  = rim[i] * preconim[i];
	    rdre[i] -= zdre[i] * alpha;
	    rdim[i] -= zdim[i] * alpha;
	    zdre[i]  = rdre[i] * preconre[i];
	    zdim[i]  = rdim[i] * preconim[i];
	}
	err = sqrt ((rre & rre) + (rim & rim)) / bnorm;
    } while (++k < maxit && err > tol);
    return k;
}

// ==========================================================================
// class and friend instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TGenericSparseMatrix<double>;
template class MATHLIB TGenericSparseMatrix<float>;
template class MATHLIB TGenericSparseMatrix<toast::complex>;
template class MATHLIB TGenericSparseMatrix<scomplex>;
template class MATHLIB TGenericSparseMatrix<int>;

template int QRFactorize (TGenericSparseMatrix<double> &A, TVector<double> &c,
    TVector<double> &d);
template void RSolve (const TGenericSparseMatrix<double> &A,
    const TVector<double> &d, TVector<double> &b);
template void QRSolve (const TGenericSparseMatrix<double> &A,
    const TVector<double> &c, const TVector<double> &d,
    const TVector<double> &b, TVector<double> &x);

template MATHLIB int IterativeSolve (const RGenericSparseMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit);
template MATHLIB int IterativeSolve (const FGenericSparseMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit);
template MATHLIB void IterativeSolve (const RGenericSparseMatrix &A,
    const RVector *b, RVector *x, int nrhs, double tol, int maxit,
    RPreconditioner *precon, IterativeSolverResult *res);
template MATHLIB void IterativeSolve (const FGenericSparseMatrix &A,
    const FVector *b, FVector *x, int nrhs, double tol, int maxit,
    FPreconditioner *precon, IterativeSolverResult *res);


template MATHLIB int CG (const RGenericSparseMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit);
template MATHLIB int CG (const FGenericSparseMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit);

template MATHLIB int BiCG (const FGenericSparseMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit);
template MATHLIB int BiCG (const RGenericSparseMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit);

template MATHLIB int GaussSeidel (const FGenericSparseMatrix &A,
    const FVector &b, FVector &x, double &tol, int maxit);
template MATHLIB int GaussSeidel (const RGenericSparseMatrix &A,
    const RVector &b, RVector &x, double &tol, int maxit);
template MATHLIB int GaussSeidel (const CGenericSparseMatrix &A,
    const CVector &b, CVector &x, double &tol, int maxit);
template MATHLIB int GaussSeidel (const SCGenericSparseMatrix &A,
    const SCVector &b, SCVector &x, double &tol, int maxit);

// Note that specialisations are not explicitly instantiated

template int ComplexBiCGSolve (const RGenericSparseMatrix &Are,
    const RGenericSparseMatrix &Aim, const RVector &bre,
    const RVector &bim, RVector &xre, RVector &xim, double &tol, int maxit);

#endif // NEED_EXPLICIT_INSTANTIATION
