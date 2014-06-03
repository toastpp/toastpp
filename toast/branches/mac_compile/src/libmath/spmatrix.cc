// ==========================================================================
// Module mathlib
// File spmatrix.cc
// Definition of template class TSparseMatrix ('template sparse matrix')
// ==========================================================================

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "mathdef.h"
#include "complex.h"
#include "vector.h"
#include "spvector.h"
#include "rtmatrix.h"
#include "matrix.h"
#include "sqmatrix.h"
#include "symatrix.h"
#include "bsmatrix.h"
#include "spmatrix.h"

// ==========================================================================
// member definitions

template<class MT>
TSparseMatrix<MT>::TSparseMatrix ()
: TRootMatrix<MT> ()
{
    data = 0;
}

template<class MT>
TSparseMatrix<MT>::TSparseMatrix (int r, int c)
: TRootMatrix<MT> ()
{
    data = 0;
    Allocate (r, c);
}

template<class MT>
TSparseMatrix<MT>::TSparseMatrix (const TSparseMatrix<MT> &mat)
: TRootMatrix<MT> ()
{
    data = 0;
    Allocate (mat.rows, mat.cols);
    Copy (mat);
}

template<class MT>
TSparseMatrix<MT>::~TSparseMatrix ()
{
    Unlink ();
}

template<class MT>
TSparseMatrix<MT>::operator TMatrix<MT> ()
{
    TMatrix<MT> tmp (rows, cols);
    for (int i = 0; i < rows; i++) tmp[i] = data[i];
    return tmp;
}

template<class MT>
void TSparseMatrix<MT>::Copy (const TSparseMatrix<MT> &mat)
{
    dASSERT(rows == mat.rows && cols == mat.cols,
	Matrices have different size.);
    for (int i = 0; i < rows; i++) data[i].Copy (mat.data[i]);
}

template<class MT>
void TSparseMatrix<MT>::New (int r, int c)
{
    Unlink ();
    Allocate (r, c);
}

template<class MT>
void TSparseMatrix<MT>::Unlink ()
{
    if (data) delete[]data;
    data = 0;
    rows = cols = 0;
}

template<class MT>
void TSparseMatrix<MT>::Shrink ()
{
    for (int r = 0; r < rows; r++) data[r].Shrink ();
}

template<class MT>
void TSparseMatrix<MT>::GetFillIn (int *fillin)
{
    for (int r = 0; r < rows; r++)
        fillin[r] = data[r].PDim();
}

template<class MT>
void TSparseMatrix<MT>::Clear ()
{
    for (int r = 0; r < rows; r++) data[r].Clear ();
}

template<class MT>
void TSparseMatrix<MT>::Allocate (int r, int c)
{
    dASSERT(!data, Data block present. Use Unlink first.);
    if (r) {
	data = new TSparseVector<MT>[r];
	dASSERT (data, Memory allocation failed.);
	for (int i = 0; i < r; i++) data[i].New (c);
	rows = r, cols = c;
    }
}

template<class MT>
void TSparseMatrix<MT>::Initialise (int *pindex)
{
    for (int r = 0; r < rows; r++)
	data[r].Allocate (pindex[r]);
}

template<class MT>
void TSparseMatrix<MT>::Initialise (int pindex)
{
    for (int r = 0; r < rows; r++)
	data[r].Allocate (pindex);
}

template<class MT>
void TSparseMatrix<MT>::Initialise (int *rowptr, int *colidx)
{
    for (int r = 0; r < rows; r++) {
        int nc = rowptr[r+1] - rowptr[r];
	data[r].Initialise (nc, colidx);
	colidx += nc;
    }	  
}

template<class MT>
MT TSparseMatrix<MT>::Get (int r, int c) const
{
    dASSERT(r >= 0 && r < rows, Index out of range);
    return data[r].Get (c);
}

template<class MT>
void TSparseMatrix<MT>::Put (int r, int c, MT val)
{
    dASSERT(r >= 0 && r < rows, Index out of range);
    data[r].Put (c, val);
}

template<class MT>
void TSparseMatrix<MT>::Add (int r, int c, MT val)
{
    dASSERT(r >= 0 && r < rows, Index out of range);
    data[r].Add (c, val);
}

template<class MT>
TVector<MT> TSparseMatrix<MT>::Row (int r) const
{
    return ToVector (data[r]);
}

template<class MT>
TVector<MT> TSparseMatrix<MT>::Col (int c) const
{
    int ind;
    TVector<MT> col(rows);
    for (int i = 0; i < rows; i++)
	if ((ind = data[i].PIndex(c)) >= 0) col[i] = data[i][ind];
    return col;
}

template<class MT>
TSparseMatrix<MT> TSparseMatrix<MT>::operator= (const TSparseMatrix<MT> &mat)
{
    Copy (mat);
    return *this;
}

template<class MT>
TSparseMatrix<MT> TSparseMatrix<MT>::operator* (const MT &mt) const
{
    TSparseMatrix<MT> tmp(*this);
    for (int i = 0; i < rows; i++) tmp[i] *= mt;
    return tmp;
}

template<class MT>
TVector<MT> TSparseMatrix<MT>::operator* (const TVector<MT> &x) const
{
    TVector<MT> b(rows);
    Ax (x, b);
    return b;
}

#ifdef MATH_DEBUG // otherwise inline
template<class MT>
void TSparseMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == cols, Invalid size - vector x);
    dASSERT(b.Dim() == rows, Invalid size - vector b);
    for (int r = 0; r < rows; r++) b[r] = data[r] & x;
}
#endif // MATH_DEBUG

template<class MT>
TVector<MT> TSparseMatrix<MT>::ATx (const TVector<MT> &x) const
{
    dASSERT(rows == x.Dim(), Matrix and vector not compatible.);
    TVector<MT> y(rows);
    for (int r = 0; r < rows; r++)
	for (int c = 0; c < data[r].PDim(); c++)
	    y[data[r].LIndex(c)] += data[r][c] * x[r];
    return y;
}

template<class MT>
void TSparseMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == rows, Vector x wrong size.);
    dASSERT(b.Dim() == cols, Vector b wrong size.);
    int r, c;
    for (c = 0; c < cols; c++) b[c] = (MT)0;
    for (r = 0; r < rows; r++)
	for (c = 0; c < data[r].PDim(); c++)
	    b[data[r].LIndex(c)] += data[r][c] * x[r];
}

template<class MT>
TVector<MT> TSparseMatrix<MT>::diag () const
{
    int n = min (rows, cols);
    TVector<MT> d(n);
    for (int i = 0; i < n; i++) d[i] = Get (i, i);
    return d;
}

template<class MT>
void TSparseMatrix<MT>::SparseOutput (ostream &os) const
{
    for (int i = 0; i < rows; i++) {
	os << "row " << i << ": " ;
	data[i].SparseOutput (os);
	os << endl;
    }
}

#ifdef MATH_DEBUG // otherwise inline
template<class MT>
TSparseVector<MT> &TSparseMatrix<MT>::operator[] (int i) const
{
    dASSERT(i >= 0 && i < rows, Index out of range);
    return data[i];
}
#endif // MATH_DEBUG


// ==========================================================================
// friend definitions

// transpose

template<class MT>
TSparseMatrix<MT> transpose (const TSparseMatrix<MT> &A)
{
    int row, col, i, j, *pindex;
    TSparseMatrix<MT> AT;

    // pass 1: create fill-in list for transpose matrix
    row = A.Dim(ROW), col = A.Dim(COL);
    pindex = new int[col];
    for (i = 0; i < col; i++) pindex[i] = 0;
    for (i = 0; i < row; i++)
	for (j = 0; j < A[i].PDim(); j++) pindex[A[i].LIndex(j)]++;
    AT.New (col, row);
    AT.Initialise (pindex);

    // pass 2: fill transpose
    for (i = 0; i < row; i++)
	for (j = 0; j < A[i].PDim(); j++)
	    AT.Put (A[i].LIndex(j), i, A[i][j]);

    return AT;
}

template<class MT>
bool CHdecomp (const TSparseMatrix<MT> &A, TSparseMatrix<MT> &L,
	       TSparseMatrix<MT> &LT, TVector<MT> &d, bool reallocate,
	       bool recover)
{
    LOGOUT2_ENTER;
    const double EPS = 1e-10;
    bool ok = TRUE;
    int i, j;
    int n = A.Dim(ROW);
    MT x, diag, val;

    if (reallocate) {
      L.New (n,n);
      LT.New (n,n);
      d.New(n);
    } else {
      L.Clear();
      LT.Clear();
    }

    LOGOUT2_INIT_PROGRESSBAR ("CHdecomp", 50, n);
    for (i = 0; i < n; i++) {
        diag = A.Get (i,i);
	TSparseVector<MT> &Li = L[i];
	for (x = (MT)0, j = 0; j < L[i].PDim(); j++)
	    x += Li[j] * Li[j];
	if (diag > x) {
	    d[i] = sqrt (diag - x);
	} else {
	    if (!recover) xERROR(Matrix not positive definite);
	    ok = FALSE;
	    d[i] = EPS;  // hack it positive
	}
	for (j = i+1; j < n; j++) {
	    if (i) {
	        x = iprod_sorted (Li, L[j], 0, i-1);
		// Li and Lj have been built with ascending index lists
		// so we can use the more efficient version of iprod here
		val = (A.Get(j,i)-x)/d[i];
	    } else {
	        val = A.Get(j,i)/d[i];
	    }
	    if (val) L.Put (j, i, val);
	}
	LOGOUT2_PROGRESS(i);
    }
    if (reallocate) L.Shrink();  // remove unused entries
    LT = transpose (L);
    LOGOUT2_EXIT;
    return ok;
}

template<class MT>
void CHdecompFill (const TSparseMatrix<MT> &A, TSparseMatrix<int> &L)
{
    bool nzero;
    int i, j;
    int n = A.Dim(ROW);

    L.New (n,n);

    for (i = 0; i < n; i++) {
	TSparseVector<int> &Li = L[i];
	for (j = i+1; j < n; j++) {
	    nzero = (A.Get(j,i) != (MT)0);
	    if (!nzero && i)
	        nzero = overlap_sorted (Li, L[j], 0, i-1);
	    if (nzero)
	        L.Put (j, i, 1);
	}

	if (!(i%1000)) cerr << i << endl;
    }
    L.Shrink();  // remove unused entries
}

template<class MT>
bool IncompleteCHdecomp (const TSparseMatrix<MT> &A, TSparseMatrix<MT> &L,
	       TSparseMatrix<MT> &LT, TVector<MT> &d, bool reallocate,
	       bool recover)
{
    const double EPS = 1e-10;
    bool ok = TRUE;
    int i, j, k;
    int n = A.Dim(ROW);
    MT x, diag, val, Aj;

    if (reallocate) {
      L.New (n,n);
      LT.New (n,n);
      d.New(n);
    } else {
      L.Clear();
      LT.Clear();
    }

    for (i = 0; i < n; i++) {
        diag = A.Get (i,i);
	TSparseVector<MT> &Li = L[i];
	for (x = (MT)0, j = 0; j < L[i].PDim(); j++)
	    x += Li[j] * Li[j];
	if (diag > x) {
	    d[i] = sqrt (diag - x);
	} else {
	    if (!recover) xERROR(Matrix not positive definite);
	    ok = FALSE;
	    d[i] = EPS;  // hack it positive
	}
	for (k = 0; k < A[i].PDim(); k++) {
	    j = A[i].LIndex(k);
	    if (j <= i) continue;
	    Aj = A[i][k];
	    if (i) {
	        x = iprod_sorted (Li, L[j], 0, i-1);
		// Li and Lj have been built with ascending index lists
		// so we can use the more efficient version of iprod here
		val = (Aj-x)/d[i];
	    } else {
	        val = Aj/d[i];
	    }
	    L.Put (j, i, val);
	}
    }
    if (reallocate) L.Shrink();  // remove unused entries
    LT = transpose (L);
    return ok;
}

// Cholesky decomposition of a single line of a sparse banded matrix

template<class MT>
void LineCHdecomp (const TSparseMatrix<MT> &A, TMatrix<MT> &T, int p)
{
    int i, j, hbnd = T.Dim(ROW), intdof = A.Dim(ROW);
    MT x;

    // shift frontal mask one step
//    for (i = 1; i < hbnd; i++)
//    	for (j = i; j < hbnd; j++) T[i-1][j-1] = T[i][j];
    for (i = 1; i < hbnd; i++) T[i-1].Copy (T[i], i-1, i, hbnd-i);

    for (i = p; i < p+hbnd; i++)
	if (i < intdof && (j = A[i].PIndex(p)) >= 0) T[i-p][hbnd-1] = A[i][j];
	else T[i-p][hbnd-1] = 0.0;

    for (x = T[0][hbnd-1], j = hbnd-2; j >= 0; j--) x -= T[0][j] * T[0][j];
    // need to check positive definiteness: x > 0
    T[0][hbnd-1] = sqrt (x);
    
    for (i = 1; i < hbnd; i++) {
	for (x = T[i][hbnd-1], j = hbnd-2; j >= i; j--) x -= T[0][j] * T[i][j];
	T[i][hbnd-1] = x / T[0][hbnd-1];
    }
}

// Cholesky substitution
// Vector d contains the diagonal elements, and Matrix L the lower triangle
// (without the diagonal) of the CH decomposition. LT is the transpose of L.

template<class MT>
TVector<MT> CHsubst (const TSparseMatrix<MT> &L, const TSparseMatrix<MT> &LT,
    const TVector<MT> &d, const TVector<MT> &b)
{
    TVector<MT> x(L.Dim(ROW));
    CHsubst (L, LT, d, b, x);
    return x;
}

template<class MT>
void CHsubst (const TSparseMatrix<MT> &L, const TSparseMatrix<MT> &LT,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x)
{
    int i, k, pd, n = L.Dim(ROW);
    MT sum;

    for (i = 0; i < n; i++) {
	TSparseVector<MT> &Li = L[i];
	pd = Li.PDim();
	//bd = Li.base->data;
	//bi = Li.base->index;
	for (sum = b[i], k = 0; k < pd; k++)
	    sum -= Li[k] * x[Li.LIndex(k)];
	x[i] = sum / d[i];
    }
    for (i = n-1; i >= 0; i--) {
	TSparseVector<MT> &LTi = LT[i];
	pd = LTi.PDim();
	//bd = LTi.base->data;
	//bi = LTi.base->index;
	for (sum = x[i], k = 0; k < pd; k++)
	    sum -= LTi[k] * x[LTi.LIndex(k)];
	x[i] = sum / d[i];
    }
}

// ***************************************************************************
// The following functions CGsolve and PCGsolve are variants of the Conjugate
// Gradient solver Ax=b (without/with preconditioning) for symmetric positive
// definite sparse matrices A.
// For all versions, xinit is the initial guess to the solution (should be
// set to NULL if not available), err_limit is the convergence criterion.
// The number of iterations performed is returned in niter.

// ***************************************************************************
// Non-preconditioned version : Ax=b

template<class MT>
TVector<MT> CGsolve (const TSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> *xinit, double err_limit, int *niter)
{
    double dnew, dold, d0, alpha, beta;
    int i, k, dim = b.Dim();
    TVector<MT> x(dim);
    TVector<MT> Ap(dim);
    if (xinit) x = *xinit;
    TVector<MT> r = b - (A * x);
    TVector<MT> p = r;
    dnew = r & p;
    d0 = dnew * err_limit*err_limit;
    for (k = 0; k < dim && dnew > d0; k++) {
	A.Ax (p, Ap);
	alpha = dnew / (p & Ap);
	for (i = 0; i < dim; i++) {
	    x[i] += p[i] * alpha;
	    r[i] -= Ap[i] * alpha;
	}
	dold = dnew;
	dnew = r & r;
	beta = dnew / dold;
	for (i = 0; i < dim; i++)
	    p[i] = p[i] * beta + r[i];
    }
    if (niter) *niter = k;
    return x;
}

// ***************************************************************************
// Preconditioned version: (P^-1 A)x = P^-1 b
// internally uses diagonal of A as preconditioner P

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> *xinit, double err_limit, int *niter)
{
    TVector<MT> PrD = A.diag();
    return PCGsolve (A, PrD, b, xinit, err_limit, niter);
}

// ***************************************************************************
// This version uses diagonal preconditioner, (P^-1 A) x = P^-1 b

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A, TVector<MT> P,
    const TVector<MT> &b, TVector<MT> *xinit, double err_limit,
    int *niter)
{
    const double EPS = 1e-10;
    double dnew, dold, d0, alpha, beta;
    int i, k, dim = b.Dim();
    dASSERT(dim == P.Dim(), Invalid size of preconditioner vector.);
    TVector<MT> x(dim);
    TVector<MT> z(dim);
    if (xinit) x = *xinit;
    TVector<MT> r = b - (A * x);
    for (i = 0; i < dim; i++)	// invert preconditioner, since z = P^-1 r
	P[i] = (norm (P[i]) < EPS ? (MT)1 : (MT)1/P[i]);
    TVector<MT> p = P * r;
    dnew = r & p;
    d0 = dnew * err_limit*err_limit;
    for (k = 0; k < dim && dnew > d0; k++) {
	A.Ax (p, z);
	alpha = dnew / (p & z);
	for (i = 0; i < dim; i++) {
	    x[i] += p[i] * alpha;
	    r[i] -= z[i] * alpha;
	    z[i]  = r[i] * P[i];
	}
	dold = dnew;
	dnew = r & z;
	beta = dnew / dold;
	for (i = 0; i < dim; i++)
	    p[i] = p[i] * beta + z[i];
    }
    if (niter) *niter = k;
    return x;
}

// ***************************************************************************
// This version uses diagonal preconditioner, (A P^-1) (P x) = b

template<class MT>
TVector<MT> XCGsolve (const TSparseMatrix<MT> &A, const TVector<MT> &P,
    const TVector<MT> &b, TVector<MT> *xinit = 0, double err_limit = 1e-18,
    int *niter = 0)
{
    const double EPS = 1e-10;
    double dnew, dold, d0, alpha, beta;
    int i, k, dim = b.Dim();
    TVector<MT> x(dim);
    TVector<MT> Ap(dim);
    TVector<MT> bnorm(dim);
    if (xinit) x = *xinit;
    TVector<MT> Pi(dim);
    for (i = 0; i < dim; i++) {
	bnorm[i] = (MT)1;
	Pi[i] = (norm(P[i]) < EPS ? (MT)1 : (MT)1/P[i]);
    }
    TVector<MT> r = bnorm - (A * x);
    TVector<MT> p = r;
    dnew = r & p;
    d0 = dnew * err_limit*err_limit;
    for (k = 0; k < dim && dnew > d0; k++) {
	A.Ax (p, Ap);
	alpha = dnew / (p & Ap);
	for (i = 0; i < dim; i++) {
	    x[i] += p[i] * alpha;
	    r[i] -= Ap[i] * alpha;
	}
	dold = dnew;
	dnew = r & r;
	beta = dnew / dold;
	for (i = 0; i < dim; i++)
	    p[i] = p[i] * beta + r[i];
    }
    if (niter) *niter = k;
    return Pi * x;
}

void ComplexBiCGsolve (const SparseMatrix& Are,
    const SparseMatrix& Aim, const RVector& bre, const RVector& bim,
    RVector& xre, RVector& xim, double err_limit, int *niter)
{
    dASSERT(Are.Dim(ROW) == Aim.Dim(ROW) && Are.Dim(COL) == Aim.Dim(COL),
	   Matrix dimensions differ.);
    dASSERT(bre.Dim() == bim.Dim() && xre.Dim() == xim.Dim(),
	   Vector dimensions differ.);
    dASSERT(Are.Dim(ROW) == bre.Dim(), Matrix and vector incompatible.);

    const double EPS = 1e-10;
    int i, k = 0, dim = bre.Dim();
    double bknum, bkden, akden, alpha, err;
    double bnorm = sqrt ((bre & bre) + (bim & bim));
    RVector pre(dim), pim(dim);
    RVector pdre(dim), pdim(dim);
    xre = 0.0; xim = 0.0;
    RVector rre = bre - (Are*xre - Aim*xim);
    RVector rim = bim - (Aim*xre + Are*xim);
    RVector rdre = rre;
    RVector rdim = rim;
    RVector preconre = Are.diag();
    RVector preconim = Are.diag();
    for (i=0; i<dim; i++) {
	preconre[i] = (fabs (preconre[i]) < EPS ? 1.0 : 1.0 / preconre[i]);
	preconim[i] = (fabs (preconim[i]) < EPS ? 1.0 : 1.0 / preconim[i]);
    }
    RVector zre = rre * preconre;
    RVector zim = rim * preconim;
    RVector zdre = zre;
    RVector zdim = zim;
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
    } while (++k < dim && err > err_limit);
    if (niter != NULL) *niter = k;
}

// ***************************************************************************
// this version uses an incomplete Cholesky decomposition as preconditioner:
// PrD is the diagonal, Pr and PrT are the lower triangle and the transpose
// of the decomposition (without the diagonal).

template<class MT>
TVector<MT> PCGsolve (const TSparseMatrix<MT> &A, const TSparseMatrix<MT> &Pr,
    const TSparseMatrix<MT> &PrT, const TVector<MT> &PrD,
    const TVector<MT> &b, TVector<MT> *xinit, double err_limit, int *niter)
{
    double dnew, dold, d0, alpha, beta;
    int i, k, dim = b.Dim();
    TVector<MT> x(dim), z(dim), p(dim);
    if (xinit) x = *xinit;
    TVector<MT> r = b - (A * x);
    CHsubst (Pr, PrT, PrD, r, p);
    dnew = r & p;
    d0 = dnew * err_limit*err_limit;
    for (k = 0; k < dim && dnew > d0; k++) {
	A.Ax (p, z);
	alpha = dnew / (p & z);
	for (i = 0; i < dim; i++) {
	    x[i] += p[i] * alpha;
	    r[i] -= z[i] * alpha;
	}
	CHsubst (Pr, PrT, PrD, r, z);
	dold = dnew;
	dnew = r & z;
	beta = dnew / dold;
	for (i = 0; i < dim; i++)
	    p[i] = p[i] * beta + z[i];
    }
    if (niter) *niter = k;
    return x;
}

// ***************************************************************************
// complex versions not implemented (yet), so the following are just dummies

TVector<complex> CGsolve (const CSparseMatrix &, const CVector &,
    const CVector &, CVector*, double, int*)
{
    return CVector();
}

TVector<complex> PCGsolve (const CSparseMatrix &, const CVector &,
    const CVector &, CVector*, double, int*)
{
    return CVector();
}

TVector<complex> PCGsolve (const CSparseMatrix &, const CVector &,
    CVector*, double, int*)
{
    return CVector();
}

TVector<complex> XCGsolve (const CSparseMatrix &, const CVector &,
    const CVector &, CVector*, double, int*)
{
    return CVector();
}

TVector<complex> PCGsolve (const CSparseMatrix &, const CSparseMatrix &,
    const CSparseMatrix &, const CVector &, const CVector &, CVector *,
    double, int *)
{
    return CVector();
}

// ***************************************************************************
// Biconjugate gradient solver; this version is valid for nonsymmetric,
// non-positive definite matrices

template<class MT>
TVector<MT> BiPCGsolve (const TSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> *xinit, double err_limit, int *niter)
{
    dASSERT(A.Dim(ROW) == b.Dim(), Matrix and vector are incompatible.);

    int i, k = 0, dim = b.Dim();
    const double EPS = 1e-10;
    MT bknum, bkden, bk, akden, alpha;
    double err, bnorm = l2norm (b);
    TVector<MT> p(dim);
    TVector<MT> pd(dim);
    TVector<MT> x(dim);
    if (xinit) x = *xinit;	// initialise solution
    TVector<MT> r = b - (A * x);
    TVector<MT> rd = r;
    TVector<MT> precon = A.diag ();
    for (i = 0; i < dim; i++)
	precon[i] = (norm (precon[i]) < EPS ? (MT)1 : (MT)1/precon[i]);
    TVector<MT> z = r * precon;
    TVector<MT> zd = z;
    do {
	bknum = z & rd;
	if (k) {
	    bk = bknum/bkden;
	    for (i = 0; i < dim; i++) {
		p[i]  = p[i] * bk + z[i];
		pd[i] = pd[i] * bk + zd[i];
	    }
	} else {
	    p  = z;
	    pd = zd;
	}
	bkden = bknum;
	A.Ax (p, z);		// z = A * p
	A.ATx (pd, zd);		// zd = transpose(A) * pd
	akden = z & pd;
	alpha = bknum / akden;
	for (i = 0; i < dim; i++) {
	    x[i]  += p[i] * alpha;
	    r[i]  -= z[i] * alpha;
	    z[i]   = r[i] * precon[i];
	    rd[i] -= zd[i] * alpha;
	    zd[i]  = rd[i] * precon[i];
	}
	err = l2norm (r) / bnorm;
    } while (++k < dim && err > err_limit);
    if (niter) *niter = k;
    return x;
}


// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TSparseMatrix<double>;
template class TSparseMatrix<float>;
template class TSparseMatrix<complex>;
template class TSparseMatrix<int>;

template RSparseMatrix transpose (const RSparseMatrix &A);

template void CHdecompFill (const RSparseMatrix &A, ISparseMatrix &L);

template bool CHdecomp (const RSparseMatrix &A, RSparseMatrix &L,
    RSparseMatrix &LT, RVector &d, bool reallocate, bool recover);

template bool IncompleteCHdecomp (const RSparseMatrix &A,
    RSparseMatrix &L, RSparseMatrix &LT, RVector &d,
    bool reallocate, bool recover);

template void LineCHdecomp (const RSparseMatrix &A, RMatrix &T, int p);

template RVector CHsubst (const RSparseMatrix &L,
    const RSparseMatrix &LT, const RVector &d, const RVector &b);

template void CHsubst (const RSparseMatrix &L, const RSparseMatrix &LT,
    const RVector &d, const RVector &b, RVector &x);

template RVector CGsolve (const RSparseMatrix &A, const RVector &b,
    RVector *xinit, double err_limit, int *niter);
template FVector CGsolve (const FSparseMatrix &A, const FVector &b,
    FVector *xinit, double err_limit, int *niter);

template RVector PCGsolve (const RSparseMatrix &A, const RVector &b,
    RVector *xinit, double err_limit, int *niter);
template FVector PCGsolve (const FSparseMatrix &A, const FVector &b,
    FVector *xinit, double err_limit, int *niter);

template RVector PCGsolve (const RSparseMatrix &A, RVector P,
    const RVector &b, RVector *xinit, double err_limit, int *niter);
template FVector PCGsolve (const FSparseMatrix &A, FVector P,
    const FVector &b, FVector *xinit, double err_limit, int *niter);

template RVector XCGsolve (const RSparseMatrix &A, const RVector &P,
    const RVector &b, RVector *xinit, double err_limit, int *niter);
template FVector XCGsolve (const FSparseMatrix &A, const FVector &P,
    const FVector &b, FVector *xinit, double err_limit, int *niter);

template RVector PCGsolve (const RSparseMatrix &A, const RSparseMatrix &Pr,
    const RSparseMatrix &PrT, const RVector &PrD,
    const RVector &b, RVector *xinit, double err_limit, int *niter);

template RVector BiPCGsolve (const RSparseMatrix &A, const RVector &b,
    RVector *xinit, double err_limit, int *niter);
template FVector BiPCGsolve (const FSparseMatrix &A, const FVector &b,
    FVector *xinit, double err_limit, int *niter);
template CVector BiPCGsolve (const CSparseMatrix &A, const CVector &b,
    CVector *xinit, double err_limit, int *niter);

#endif // NEED_EXPLICIT_INSTANTIATION

