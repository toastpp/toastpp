// ==========================================================================
// Module mathlib
// File symatrix.cc
// Definition of template class TSymMatrix ('template symmetric matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __SYMATRIX_CC

#include <iostream>
#include <sstream>
#include "mathlib.h"

// ==========================================================================
// member definitions

template<class MT>
TSymMatrix<MT>::TSymMatrix (int n, const char *valstr): TMatrix<MT> (n, n)
{
    Alloc (n);
    std::istringstream iss (valstr);
    for (int i = 0; i < nz; i++) iss >> val[i];
}

template<class MT>
void TSymMatrix<MT>::New_dirty (int n)
{
    if (n != this->rows) {            // realloc only if size differs
        TMatrix<MT>::New (n, n);      // set nrows and ncols
        Unlink ();                    // dealloc current data array
	Alloc (n);                    // and re-allocate with new size
    }
}

template<class MT>
void TSymMatrix<MT>::Identity ()
{
    Zero();
    int i, incr = 0, n = (this->rows < this->cols ? this->rows : this->cols);
    MT *v = val;
    for (i = 0; i < n; i++) {
        *(v+incr) = (MT)1;
        v += ++incr;
    }
}

template<class MT>
TVector<MT> TSymMatrix<MT>::Row (int r) const
{
    TVector<MT> row(this->cols);
    MT *v = val + Idx(r,0);
    int i;
    for (i = 0; i < r; i++) row[i] = *v++;
    for (; i < this->cols; i++) {
	row[i] = *v;
	v += ++r; // jump to next row
    }
    return row;
}

template<class MT>
int TSymMatrix<MT>::SparseRow (int r, idxtype *ci, MT *rv) const
{
    int i;
    MT *v = val + Idx(r,0);
    for (i = 0; i < r; i++) {
        ci[i] = i;
        rv[i] = *v++;
    }
    for (; i < this->cols; i++) {
      ci[i] = i;
        rv[i] = *v;
	v += ++r;
    }
    return this->cols;
}

template<class MT>
void TSymMatrix<MT>::ColScale (const TVector<MT> &scale)
{
    xERROR ("Operation not permitted for symmetric matrix");
}

template<class MT>
void TSymMatrix<MT>::RowScale (const TVector<MT> &scale)
{
    xERROR ("Operation not permitted for symmetric matrix");
}

template<class MT>
void TSymMatrix<MT>::AddDiag (const TVector<MT> &d)
{
    dASSERT(this->rows == d.Dim(),
        "Argument 1 wrong size (expected %d, actual %d)", this->rows, d.Dim());

    int i, incr = 0;

    MT *v = val;
    for (i = 0; i < this->rows; i++) {
        *(v+incr) += d[i];
	v += ++incr;
    }
}

template<class MT>
void TSymMatrix<MT>::AddDiag (const MT &d)
{
    int i, incr = 0;

    MT *v = val;
    for (i = 0; i < this->rows; i++) {
        *(v+incr) += d;
	v += ++incr;
    }
}

template<class MT>
void TSymMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == x.Dim(), "Argument 1: vector has wrong dimension");

    int r, c, incr;
    MT *row;
    if (b.Dim() != this->rows) b.New (this->rows); // resize
    for (r = 0; r < this->rows; r++) {
        MT &br = b[r];
	row = val + Idx(r,0);
	for (c = 0, br = (MT)0; c <= r; c++)
	    br += *row++ * x[c];
	for (--row, incr = r; c < this->cols; c++) {
	    row += ++incr;
	    br += *row * x[c];
	}
    }
}

template<class MT>
TSymMatrix<MT> &TSymMatrix<MT>::operator= (const TSymMatrix<MT> &m)
{
    if (this->rows != m.rows) New_dirty (m.rows);  // resize
    memcpy (val, m.val, nz*sizeof(MT));            // copy
    return *this;
}

template<class MT>
TSymMatrix<MT> &TSymMatrix<MT>::operator= (const MT &mt)
{
    for (int i = 0; i < nz; i++) val[i] = mt;
    return *this;
}

template<class MT>
TSymMatrix<MT> TSymMatrix<MT>::operator+ (const TSymMatrix<MT> &m) const
{
    dASSERT(this->rows == m.rows, "Matrices have different size.");
    TSymMatrix<MT> tmp(*this);
    for (int i = 0; i < nz; i++) tmp.val[i] += m.val[i];
    return tmp;
}

template<class MT>
TSymMatrix<MT> TSymMatrix<MT>::operator- (const TSymMatrix<MT> &m) const
{
    dASSERT(this->rows == m.rows, "Matrices have different size.");
    TSymMatrix<MT> tmp(*this);
    for (int i = 0; i < nz; i++) tmp.val[i] -= m.val[i];
    return tmp;
}

template<class MT>
TSymMatrix<MT> TSymMatrix<MT>::operator* (const MT &mt) const
{
    TSymMatrix<MT> tmp(*this);
    for (int i = 0; i < nz; i++) tmp.val[i] *= mt;
    return tmp;
}

template<class MT>
TVector<MT> TSymMatrix<MT>::operator* (const TVector<MT> &x) const
{
    TVector<MT> b(this->rows);
    Ax (x, b);
    return b;
}

// ==========================================================================
// friend definitions

template<>
inline bool CHdecomp<double> (TSymMatrix<double> &a, bool recover)
{
    // *** Cholesky decomposition ***
    // Specialisation for type double (checks for positive definiteness)

    const double EPS = 1e-10;
    double sum, *Ai, *Aj, iAi;
    bool ok = true;
    int i, j, k, n = a.nRows();

    Ai = a.val;          // beginning of first row
    for (i = 0; i < n; i++) {
        Aj = Ai;

	// diagonal element
	for (sum = Ai[i], k = 0; k < i; k++) sum -= *Ai++ * *Aj++;
	if (sum > 0.0) {
	    *Ai = sqrt (sum);
	} else {
	    if (!recover) xERROR("Matrix not positive definite");
	    ok = false;  // set 'bad' flag
	    *Ai = EPS;   // force diagonal element positive
	}
	iAi = 1.0 / *Ai; // store 1/diag for later

	// column below diagonal
	for (j = i+1; j < n; j++) {
	    Ai -= i;     // return to beginning of row i
	    Aj += j-i;   // forward to beginning of next row j
	    for (sum = Aj[i], k = 0; k < i; k++) sum -= *Ai++ * *Aj++;
	    *Aj = sum * iAi;
	}

	Ai++;            // forward to beginning of next row i
    }
    return ok;
}

template<class MT>
bool CHdecomp (TSymMatrix<MT> &a, bool recover)
{
    // *** Cholesky decomposition ***
  
    MT sum, *Ai, *Aj;
    int i, j, k, n = a.nRows();

    Ai = a.val;          // beginning of first row
    for (i = 0; i < n; i++) {
        Aj = Ai;

	// diagonal element
	for (sum = Ai[i], k = 0; k < i; k++) sum -= *Ai++ * *Aj++;
	*Ai = sqrt (sum);

	// column below diagonal
	for (j = i+1; j < n; j++) {
	    Ai -= i;     // return to beginning of row i
	    Aj += j-i;   // forward to beginning of next row j
	    for (sum = Aj[i], k = 0; k < i; k++) sum -= *Ai++ * *Aj++;
	    *Aj = sum / *Ai;
	}

	Ai++;            // forward to beginning of next row i
    }
    return true;
}

template<class MT>
TVector<MT> CHsubst (const TSymMatrix<MT> &a, const TVector<MT> &b)
{
    int i, k, n = a.Dim(TMatrix<MT>::ROW);
    MT sum, *Ai, *Ak;
    TVector<MT> x(n);

    for (i = 0; i < n; i++) {
        Ai = a.val + a.Idx(i,0);
	for (sum = b[i], k = i-1; k >= 0; k--)
	    sum -= *(Ai+k) * x[k];
	x[i] = sum / *(Ai+i);
    }
    for (i = n-1; i >= 0; i--) {
        Ai = Ak = a.val + a.Idx(i,i);
	for (sum = x[i], k = i+1; k < n; k++)
	    sum -= *(Ak += k) * x[k];
	x[i] = sum / *Ai;
    }
    return x;
}


// ==========================================================================
// class and friend instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TSymMatrix<double>;
template class MATHLIB TSymMatrix<float>;
template class MATHLIB TSymMatrix<toast::complex>;
template class MATHLIB TSymMatrix<scomplex>;
template class MATHLIB TSymMatrix<int>;

//template MATHLIB bool CHdecomp (RSymMatrix &a, bool recover); // instantiated by specialisation
template MATHLIB bool CHdecomp (FSymMatrix &a, bool recover);
template MATHLIB bool CHdecomp (CSymMatrix &a, bool recover);
//template MATHLIB bool CHdecomp (ISymMatrix &a, bool recover);

template MATHLIB RVector CHsubst (const RSymMatrix &a, const RVector &b);
template MATHLIB FVector CHsubst (const FSymMatrix &a, const FVector &b);
template MATHLIB CVector CHsubst (const CSymMatrix &a, const CVector &b);
//template MATHLIB IVector CHsubst (const ISymMatrix &a, const IVector &b);

#endif // NEED_EXPLICIT_INSTANTIATION

