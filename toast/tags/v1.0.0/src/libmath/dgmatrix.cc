// ==========================================================================
// Module mathlib
// File dgmatrix.cc
// Definition of template class TDiagMatrix ('template diagonal matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include "mathlib.h"

using namespace toast;

// ==========================================================================
// member definitions

template<class MT>
TDiagMatrix<MT>::TDiagMatrix (): TGenericSparseMatrix<MT> ()
{}

template<class MT>
TDiagMatrix<MT>::TDiagMatrix (int r, int c, const MT v):
    TGenericSparseMatrix<MT> (r, c)
{
    this->nbuf = this->nval = min(r,c);
    if (this->nbuf) {
	this->val = new MT[this->nval];
	for (int i = 0; i < this->nval; i++) this->val[i] = v;
    }
}

template<class MT>
TDiagMatrix<MT>::TDiagMatrix (const TDiagMatrix<MT> &mat):
    TGenericSparseMatrix<MT> (mat.nRows(), mat.nCols())
{
    this->nbuf = this->nval = min(this->rows,this->cols);
    if (this->nbuf) {
	this->val = new MT[this->nval];
	for (int i = 0; i < this->nval; i++) this->val[i] = mat.val[i];
    }
}

template<class MT>
TDiagMatrix<MT>::~TDiagMatrix ()
{}

template<class MT>
void TDiagMatrix<MT>::New (int nrows, int ncols)
{
    TGenericSparseMatrix<MT>::New (nrows, ncols);
    this->nbuf = this->nval = min (nrows, ncols);
    this->val = new MT[this->nval];
    for (int i = 0; i < this->nval; i++) this->val[i] = (MT)0;
}

template<class MT>
TVector<MT> TDiagMatrix<MT>::Row (int r) const
{
    TVector<MT> row(this->cols);
    row[r] = this->val[r];
    return row;
}

template<class MT>
TVector<MT> TDiagMatrix<MT>::Col (int c) const
{
    TVector<MT> col(this->rows);
    col[c] = this->val[c];
    return col;
}

template<class MT>
int TDiagMatrix<MT>::SparseRow (int r, idxtype *colidx, MT *v) const
{
    if (r < this->nval) {
	colidx[0] = r;
	v[0] = this->val[r];
	return 1;
    } else {
	return 0;
    }
}

template<class MT>
void TDiagMatrix<MT>::ColScale (const TVector<MT> &scale)
{
    dASSERT(scale.Dim() == this->cols, "Argument 1: wrong size");
    int nz = min (this->cols, this->nval);
    for (int i = 0; i < nz; i++)
	this->val[i] *= scale[i];
}

template<class MT>
void TDiagMatrix<MT>::RowScale (const TVector<MT> &scale)
{
    dASSERT(scale.Dim() == this->rows, "Argument 1: wrong size");
    int nz = min (this->rows, this->nval);
    for (int i = 0; i < nz; i++)
	this->val[i] *= scale[i];
}

template<class MT>
TDiagMatrix<MT> &TDiagMatrix<MT>::operator= (const TDiagMatrix<MT> &mat)
{
    Copy (mat);
    return *this;
}

template<class MT>
TDiagMatrix<MT> &TDiagMatrix<MT>::operator= (const MT &v)
{
    for (int i = 0; i < this->nval; i++)
	this->val[i] = v;
    return *this;
}

template<class MT>
void TDiagMatrix<MT>::Copy (const TDiagMatrix<MT> &mat)
{
    int nz = min (mat.nRows(), mat.nCols());
    this->rows = mat.nRows();
    this->cols = mat.nCols();
    if (nz != this->nval) {
	this->Unlink();
	this->nbuf = this->nval = nz;
	if (this->nval) this->val = new MT[this->nval];
    }
    for (int i = 0; i < this->nval; i++)
	this->val[i] = mat.val[i];
}

template<class MT>
TDiagMatrix<MT>::operator TDenseMatrix<MT> ()
{
    TDenseMatrix<MT> tmp (this->rows, this->cols);
    for (int i = 0; i < this->nval; i++)
	tmp(i,i) = this->val[i];
    return tmp;
}

template<class MT>
MT &TDiagMatrix<MT>::operator() (int r, int c)
{
    static MT zero = (MT)0;
    dASSERT(r >= 0 && r < this->rows, "Argument 1 out of range");
    dASSERT(c >= 0 && c < this->cols, "Argument 2 out of range");
    if (r == c) return this->val[r];
    else return zero;
}

template<class MT>
TDiagMatrix<MT> TDiagMatrix<MT>::operator+ (const TDiagMatrix<MT> &mat) const
{
    dASSERT(this->rows == mat.rows && this->cols == mat.cols,
	    "Matrices have different size.");
    TDiagMatrix<MT> tmp(this->rows, this->cols);
    for (int i = 0; i < this->nval; i++)
	tmp.val[i] = this->val[i] + mat.val[i];
    return tmp;
}

template<class MT>
TDiagMatrix<MT> TDiagMatrix<MT>::operator- (const TDiagMatrix<MT> &mat) const
{
    dASSERT(this->rows == mat.rows && this->cols == mat.cols,
	    "Matrices have different size.");
    TDiagMatrix<MT> tmp(this->rows, this->cols);
    for (int i = 0; i < this->nval; i++)
	tmp.val[i] = this->val[i] - mat.val[i];
    return tmp;
}

template<class MT>
bool TDiagMatrix<MT>::Exists (int r, int c) const
{
    if (r >= 0 && r < this->rows &&
        c >= 0 && c < this->cols &&
	r == c) return true;
    else return false;
}

template<class MT>
int TDiagMatrix<MT>::Get_index (int r, int c) const
{
    dASSERT(r >= 0 && r < this->rows, "Argument 1 out of range");
    dASSERT(c >= 0 && c < this->cols, "Argument 2 out of range");
    if (r == c) return r;
    else return -1;
}

template<class MT>
MT TDiagMatrix<MT>::GetNext (int &r, int &c) const
{
    if (r < 0) { // first nonzero
	if (this->nval) {
	    r = c = 0;
	    return this->val[0];
	} else return (MT)0;
    } else { // next nonzero
	dASSERT(r == c, "Inconsistent index combination");
	if (r < this->nval-1) {
	    r++, c++;
	    return this->val[r];
	} else {
	    r = -1;
	    return (MT)0;
	}
    }
}

template<class MT>
void TDiagMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->cols,
	    "Parameter 1 invalid size (expected %d, actual %d)",
	    this->cols, x.Dim());

    if (b.Dim() != this->rows) b.New (this->rows);

    int nz = min (this->nval, this->rows);
    for (int i = 0; i < nz; i++)
	b[i] = this->val[i] * x[i];
}

template<class MT>
void TDiagMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b, int r1, int r2)
    const
{
    dASSERT(x.Dim() == this->cols,
	    "Parameter 1 invalid size (expected %d, actual %d)",
	    this->cols, x.Dim());

    if (b.Dim() != this->rows) b.New (this->rows);
    
    r2 = min (r2, min (this->rows, this->nval));
    r1 = min (r1, r2);
    for (int r = r1; r < r2; r++)
	b[r] = this->val[r] * x[r];
}

template<class MT>
void TDiagMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->rows,
	    "Parameter 1 invalid size (expected %d, actual %d",
	    this->rows, x.Dim());

    if (b.Dim() != this->cols) b.New (this->cols);

    int nz = min (this->nval, this->cols);
    for (int i = 0; i < nz; i++)
	b[i] = this->val[i] * x[i];
}

// ==========================================================================
// friend definitions

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TDiagMatrix<double>;
template class TDiagMatrix<float>;
template class TDiagMatrix<toast::complex>;
template class TDiagMatrix<int>;

#endif // NEED_EXPLICIT_INSTANTIATION

