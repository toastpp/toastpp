// ==========================================================================
// Module mathlib
// File bsmatrix.cc
// Definition of template class TBandSymMatrix ('template banded symmetric
// matrix')
// ==========================================================================

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "mathdef.h"
#include "complex.h"
#include "vector.h"
#include "rtmatrix.h"
#include "matrix.h"
#include "sqmatrix.h"
#include "symatrix.h"
#include "bsmatrix.h"

// ==========================================================================
// member definitions

template<class MT>
void TBandSymMatrix<MT>::Copy (const TBandSymMatrix<MT> &mat)
{
    dASSERT(rows == mat.rows && hband == mat.hband,
	Matrices have different size.);
    for (int i = 0; i < rows; i++) data[i].Copy (mat.data[i]);
}

template<class MT>
void TBandSymMatrix<MT>::Unlink ()
{
    if (data) delete []data;
    data = 0;
    rows = cols = hband = 0;
}

template<class MT>
void TBandSymMatrix<MT>::Allocate (int rc, int hb)
{
    dASSERT (!data, Data block present. Use Unlink first.);
    data = new TVector<MT>[rc];
    dASSERT(data, Memory allocation failed.);
    for (int i = 0; i < rc; i++) data[i].New (hb);
    rows = cols = rc, hband = hb;
}

template<class MT>
MT TBandSymMatrix<MT>::Get (int r, int c) const
{
    dASSERT(r >= 0 && r < rows && c >= 0 && c < rows, Index out of range.);
    if (r < c) { int tmp = r; r = c; c = tmp; }
    int j = hband-1 - (r-c);
    if (j < 0) return (MT)0;
    return data[r][j];
}

template<class MT>
TVector<MT> TBandSymMatrix<MT>::Row (int r) const
{
    int c;
    TVector<MT> row(cols);
    for (c = 0; c < cols; c++) row[c] = Get (r, c);
    return row;
}

template<class MT>
TVector<MT> TBandSymMatrix<MT>::Col (int c) const
{
    TVector<MT> col(rows);
    for (int r = 0; r < rows; r++) col[r] = Get (r, c);
    return col;
}

template<class MT>
bool TBandSymMatrix<MT>::PIndex (int i, int j, int &r, int &c)
{
    r = i;
    c = hband + j - i - 1;
    return (r >= 0 && r < rows && c >= 0 && c < hband);
}

template<class MT>
void TBandSymMatrix<MT>::LIndex (int r, int c, int &i, int &j)
{
    i = r;
    j = c - hband + i + 1;
}

template<class MT>
TVector<MT> TBandSymMatrix<MT>::operator* (const TVector<MT> &x) const
{
    TVector<MT> b(rows);
    Ax (x, b);
    return b;
}

#ifdef MATH_DEBUG // otherwise inline
template<class MT>
void TBandSymMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    int r, c, k;
    MT br;

    dASSERT(rows == x.Dim(), Vector x wrong size);
    dASSERT(rows == b.Dim(), Vector b wrong size);
    for (r = 0; r < rows; r++) {
	TVector<MT> &lr = data[r];
	for (c = hband-1, k = r, br = 0; c >= 0 && k >= 0; c--, k--)
	    br += lr[c] * x[k];
	for (c = hband-2, k = r+1; c >= 0 && k < rows; c--, k++)
	    br += data[k][c] * x[k];
	b[r] = br;
    }
}
#endif // MATH_DEBUG

template<class MT>
MT TBandSymMatrix<MT>::LineVecMul (int line, const TVector<MT> &x) const
{
    int i, j;
    MT sum = (MT)0;

    RANGE_CHECK(line >= 0 && line < rows);
    dASSERT(x.Dim() == rows, Vector wrong size);

    for (i = max (0, hband-line-1); i < hband; i++)
        sum += data[line][i] * x[i+line-hband+1];

    for (i = hband-2; i >= 0; i--) {
        if ((j = line-i+hband-1) >= rows) break;
	sum += data[j][i] * x[j];
    }
    return sum;
}

#ifdef MATH_DEBUG // otherwise inline
template<class MT>
TVector<MT> &TBandSymMatrix<MT>::operator[] (int i) const
{
    dASSERT(i >= 0 && i < rows, Index out of range);
    return data[i];
}
#endif // MATH_DEBUG


// ==========================================================================
// friend definitions

template<class MT>
TSquareMatrix<MT> ToSquare (TBandSymMatrix<MT> &A)
{
    int i, j, jj;
    TSquareMatrix<MT> tmp(A.Dim(ROW), A.Dim(COL));
    for (i = 0; i < A.Dim(ROW); i++)
	for (j = 0; j < A.hband; j++)
	    if ((jj = i+j-A.hband+1) >= 0)
		tmp[i][jj] = tmp[jj][i] = A[i][j];
    return tmp;
}

template<class MT>
TSymMatrix<MT> ToSym (TBandSymMatrix<MT> &A)
{
    int i, j, jj;
    TSymMatrix<MT> tmp(A.Dim(ROW));
    for (i = 0; i < A.Dim(ROW); i++)
	for (j = 0; j < A.hband; j++)
	    if ((jj = i+j-A.hband+1) >= 0)
		tmp[i][jj] = A[i][j];
    return tmp;
}
/*
template<class MT>
void CHdecomp (TBandSymMatrix<MT> &A)
{
    int n = A.Dim(ROW);
    int hband = A.hband;
    int w = hband - 1;
    int i, j, k, la, lb, l;
    MT x;

    for (i = 0; i < n; i++) {
	for (x = (MT)0, j = 0; j < w; j++)
	    x += A[i][j] * A[i][j];
	xASSERT(A[i][w] > x, "Matrix not positive definite.");
	A[i][w] = sqrt (A[i][w] - x);
	for (k = 1; k <= w && i+k < n; k++) {
	    x = (MT)0;
	    la = i+k, lb = w-k;
	    for (l = lb-1; l >= 0; l--) x += A[la][l] * A[i][l+k];
	    A[la][lb] = (A[la][lb] - x) / A[i][w];
	}
    }
}
*/
template<class MT>
bool CHdecomp (TBandSymMatrix<MT> &A, bool recover)
{
    const double EPS = 1e-10;
    bool ok = TRUE;
    int i, j, k, wj;
    int n = A.Dim(ROW);
    int w = A.hband - 1;
    MT x;

    for (i = 0; i < n; i++) {
	TVector<MT> &Ai = A[i];
	for (x = (MT)0, j = 0; j < w; j++) x += Ai[j] * Ai[j];
	if (Ai[w] > x) {
	    Ai[w] = sqrt (Ai[w] - x);
	} else {
	    if (!recover) xERROR(Matrix not positive definite);
	    ok = FALSE;
	    Ai[w] = EPS;  // hack it positive
	}
	for (j = 1; j <= w && i+j < n; j++) {
	    x = (MT)0;
	    wj = w-j;
	    TVector<MT> &Aij = A[i+j];
	    for (k = wj-1; k >= 0; k--) x += Aij[k] * Ai[k+j];
	    Aij[wj] = (Aij[wj] - x) / Ai[w];
	}
    }
    return ok;
}

bool CHdecomp (TBandSymMatrix<complex> &A, bool recover)
{
    int n = A.Dim(ROW);
    int hband = A.hband;
    int w = hband - 1;
    int i, j, k, la, lb, l;
    complex x;

    for (i = 0; i < n; i++) {
	for (x = complex(0,0), j = 0; j < w; j++)
	    x += A[i][j] * A[i][j];
	A[i][w] = sqrt (A[i][w] - x);
	for (k = 1; k <= w && i+k < n; k++) {
	    x = complex(0,0);
	    la = i+k, lb = w-k;
	    for (l = lb-1; l >= 0; l--) x += A[la][l] * A[i][l+k];
	    A[la][lb] = (A[la][lb] - x) / A[i][w];
	}
    }
    return TRUE;
}

template<class MT>
TVector<MT> CHsubst (const TBandSymMatrix<MT> &A, const TVector<MT> &b)
{
    int n = A.Dim(ROW);
    int hband = A.hband;
    int i, j, k;
    MT sum;
    TVector<MT> x(n);

    dASSERT(n == b.Dim(), Matrix and vector not compatible.);

    for (i = 0; i < n; i++) {
	TVector<MT> &ai = A[i];
	for (sum = b[i], k = i-1, j = hband-2; k >= 0 && j >= 0; j--, k--)
	    sum -= ai[j] * x[k];
	x[i] = sum / ai[hband-1];
    }
    for (i = n-1; i >= 0; i--) {
	for (sum = x[i], k = i+1, j = hband-2; k < n && j >= 0; j--, k++)
	    sum -= A[k][j] * x[k];
	x[i] = sum / A[i][hband-1];
    }
    return x;
}

template<class MT>
TVector<MT> PCHsubst (const TBandSymMatrix<MT> &A, const TVector<MT> &b,
    const TVector<MT> &p)
{
    int n = A.Dim(ROW);
    int hband = A.hband;
    int i, j, k;
    MT sum;
    TVector<MT> x(n);

    dASSERT(n == b.Dim(), Matrix and vector not compatible.);

    for (i = 0; i < n; i++) {
	for (sum = b[i], k = i-1, j = hband-2; k >= 0 && j >= 0; k--, j--)
	    sum -= A[i][j] * x[k];
	x[i] = sum / A[i][hband-1];
    }
    for (i = n-1; i >= 0; i--) {
	for (sum = x[i], k = i+1, j = hband-2; k < n && j >= 0; k++, j--) {
	    sum -= A[k][j] * p[k] * x[k];
	}
	x[i] = sum / (A[i][hband-1] * p[i]);
    }
    return x*p;

}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TBandSymMatrix<double>;
template class TBandSymMatrix<float>;
template class TBandSymMatrix<complex>;
template class TBandSymMatrix<int>;

template RSquareMatrix ToSquare (RBandSymMatrix &A);
template FSquareMatrix ToSquare (FBandSymMatrix &A);
template CSquareMatrix ToSquare (CBandSymMatrix &A);
template ISquareMatrix ToSquare (IBandSymMatrix &A);

template RSymMatrix ToSym (RBandSymMatrix &A);
template FSymMatrix ToSym (FBandSymMatrix &A);
template CSymMatrix ToSym (CBandSymMatrix &A);
template ISymMatrix ToSym (IBandSymMatrix &A);

template bool CHdecomp (RBandSymMatrix &A, bool recover);
template bool CHdecomp (FBandSymMatrix &A, bool recover);
// complex version is specially instantiated

template RVector CHsubst (const RBandSymMatrix &A, const RVector &b);
template FVector CHsubst (const FBandSymMatrix &A, const FVector &b);
template CVector CHsubst (const CBandSymMatrix &A, const CVector &b);
template IVector CHsubst (const IBandSymMatrix &A, const IVector &b);

template RVector PCHsubst (const RBandSymMatrix &A, const RVector &b,
    const RVector &p);
template FVector PCHsubst (const FBandSymMatrix &A, const FVector &b,
    const FVector &p);
template CVector PCHsubst (const CBandSymMatrix &A, const CVector &b,
    const CVector &p);
template IVector PCHsubst (const IBandSymMatrix &A, const IVector &b,
    const IVector &p);

#endif // NEED_EXPLICIT_INSTANTIATION

