// ==========================================================================
// Module mathlib
// File cdmatrix.cc
// Definition of template class TCoordMatrix ('template coordinate storage
// matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __CDMATRIX_CC

#include "mathlib.h"
#include <iostream>
#include <string.h>

using namespace std;

// ==========================================================================
// member definitions

template<class MT>
TCoordMatrix<MT>::TCoordMatrix ()
  : TGenericSparseMatrix<MT> ()
{}

template<class MT>
TCoordMatrix<MT>::TCoordMatrix (int rows, int cols)
  : TGenericSparseMatrix<MT> (rows, cols)
{}

template<class MT>
TCoordMatrix<MT>::TCoordMatrix (int rows, int cols, int vals,
    idxtype *_rowidx, idxtype *_colidx, MT *_data)
  : TGenericSparseMatrix<MT> (rows, cols, vals, _data)
{
    int i;

    if (vals) {
        rowidx = new idxtype[vals];
	colidx = new idxtype[vals];
	for (i = 0; i < vals; i++) rowidx[i] = _rowidx[i];
	for (i = 0; i < vals; i++) colidx[i] = _colidx[i];
    }
}

template<class MT>
TCoordMatrix<MT>::TCoordMatrix (const TCoordMatrix<MT> &m)
  : TGenericSparseMatrix<MT> (m)
{
    int i;

    if (this->nbuf) {
        rowidx = new idxtype[this->nbuf];
	colidx = new idxtype[this->nbuf];
	for (i = 0; i < this->nval; i++) rowidx[i] = m.rowidx[i];
	for (i = 0; i < this->nval; i++) colidx[i] = m.colidx[i];
    }
}

template<class MT>
TCoordMatrix<MT>::~TCoordMatrix ()
{
    if (this->nbuf) {
        delete []rowidx;
	delete []colidx;
    }
}

template<class MT>
void TCoordMatrix<MT>::New (int rows, int cols)
{
    Unlink ();
    TGenericSparseMatrix<MT>::New (rows, cols);
}

template<class MT>
void TCoordMatrix<MT>::Unlink ()
{
    if (this->nbuf) {
        delete []rowidx;
	delete []colidx;
    }
    TGenericSparseMatrix<MT>::Unlink();
}

template<class MT>
void TCoordMatrix<MT>::Initialise (int nv, idxtype *_rowidx, idxtype *_colidx,
    MT *data)
{
    int i;

    if (nv != this->nval) {
        if (this->nbuf) {
	    delete []rowidx;
	    delete []colidx;
	}
	if (nv) {
	    rowidx = new idxtype[nv];
	    colidx = new idxtype[nv];
	}
    }
    for (i = 0; i < nv; i++) rowidx[i] = _rowidx[i];
    for (i = 0; i < nv; i++) colidx[i] = _colidx[i];
    TGenericSparseMatrix<MT>::Initialise (nv, data);
}

template<class MT>
TCoordMatrix<MT> &TCoordMatrix<MT>::operator= (const TCoordMatrix<MT> &m)
{
    this->rows = m.rows;
    this->cols = m.cols;
    Initialise (m.nval, m.rowidx, m.colidx, m.val);
    return *this;
}

template<class MT>
MT &TCoordMatrix<MT>::operator() (int r, int c)
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    // assumes unsorted lists
    for (int i = 0; i < this->nval; i++)
        if (rowidx[i] == r && colidx[i] == c) return this->val[i];
    int idx = Insert(r,c);
    return this->val[idx];
}

template<class MT>
MT TCoordMatrix<MT>::Get (int r, int c) const
{
    const static MT zero = (MT)0;

    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    for (int i = 0; i < this->nval; i++)
        if (rowidx[i] == r && colidx[i] == c) return this->val[i];
    return zero;
}

template<class MT>
bool TCoordMatrix<MT>::Exists (int r, int c) const
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    for (int i = 0; i < this->nval; i++)
        if (rowidx[i] == r && colidx[i] == c) return true;
    return false;
}

template<class MT>
TVector<MT> TCoordMatrix<MT>::Row (int r) const
{
    ERROR_UNDEF;
    return TVector<MT>(); // dummy
}

template<class MT>
TVector<MT> TCoordMatrix<MT>::Col (int c) const
{
    ERROR_UNDEF;
    return TVector<MT>(); // dummy
}

template<class MT>
int TCoordMatrix<MT>::SparseRow (int r, idxtype *ci, MT *rv) const
{
    int i, j = 0;
    for (i = 0; i < this->nval; i++) {
        if (rowidx[i] == r) {
	    ci[j] = colidx[i];
	    rv[j++] = this->val[i];
	}
    }
    return j;
}

template<class MT>
int TCoordMatrix<MT>::Get_index (int r, int c) const
{
    for (int i = 0; i < this->nval; i++)
        if (rowidx[i] == r && colidx[i] == c) return i;
    return -1;
}

template<class MT>
void TCoordMatrix<MT>::ColScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->cols, "Argument 1: wrong size");
    for (int i = 0; i < this->nval; i++)
        this->val[i] *= scale[colidx[i]];
}

template<class MT>
void TCoordMatrix<MT>::RowScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->cols, "Argument 1: wrong size");
    for (int i = 0; i < this->nval; i++)
	this->val[i] *= scale[rowidx[i]];
}

template<class MT>
MT TCoordMatrix<MT>::GetNext (int &r, int &c) const
{
    if (r >= 0) {
        if (++iterator_pos >= this->nval) { // end reached
	    r = -1;
	    return (MT)0; // dummy
	}
    } else {
        if (!this->nval) return (MT)0;     // empty matrix
	iterator_pos = 0;
    }
    r = rowidx[iterator_pos];
    c = colidx[iterator_pos];
    return this->val[iterator_pos];
}

template<class MT>
void TCoordMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    b.Clear();
    for (int i = 0; i < this->nval; i++)
        b[rowidx[i]] += this->val[i] * x[colidx[i]];
}

template<class MT>
void TCoordMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b,
    int r1, int r2) const
{
    // this implementation is not efficient since matrix is not assumed
    // to be sorted

    int i;
    for (i = r1; i < r2; i++) b[i] = 0;
    
    for (int i = 0; i < this->nval; i++) {
        int ri = rowidx[i];
	if (ri >= r1 && ri < r2)
	    b[ri] += this->val[i] * x[colidx[i]];
    }
}

template<class MT>
void TCoordMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->rows, "Invalid size - vector x");
    dASSERT(b.Dim() == this->cols, "Invalid size - vector b");

    b.Clear();
    for (int i = 0; i < this->nval; i++)
        b[colidx[i]] += this->val[i] * x[rowidx[i]];
}

template<class MT>
void TCoordMatrix<MT>::Transpone ()
{
    // all we have to do is swap the row and column index lists
    idxtype *tmp_idx = rowidx;
    rowidx = colidx;
    colidx = tmp_idx;
    int tmp = this->rows;
    this->rows = this->cols;
    this->cols = tmp;
}

template<class MT>
void TCoordMatrix<MT>::Sort (bool roworder) const
{
    int i, ir, j, l, rr1, rr2;
	idxtype *r1, *r2;
    int n = this->nval;
    MT v, *vl = this->val-1;
    if (roworder) r1 = rowidx-1, r2 = colidx-1;
    else          r1 = colidx-1, r2 = rowidx-1;
    
    if (n < 2) return;
    l = (n >> 1)+1;
    ir = n;
    for (;;) {
        if (l > 1) {
	    --l;
	    rr1 = r1[l], rr2 = r2[l], v = vl[l];
	} else {
	    rr1 = r1[ir], rr2 = r2[ir]; v = vl[ir];
	    r1[ir] = r1[1], r2[ir] = r2[1], vl[ir] = vl[1];
	    if (--ir == 1) {
	        r1[1] = rr1, r2[1] = rr2, vl[1] = v;
		break;
	    }
	}
	i = l;
	j = l+l;
	while (j <= ir) {
	    if (j < ir && (r1[j] < r1[j+1] ||
			   (r1[j]==r1[j+1] && r2[j] < r2[j+1])))
	        j++;
	    if (rr1 < r1[j] || (rr1==r1[j] && rr2 < r2[j])) {
	        r1[i] = r1[j], r2[i] = r2[j], vl[i] = vl[j];
		i = j;
		j <<= 1;
	    } else j = ir+1;
	}
	r1[i] = rr1, r2[i] = rr2, vl[i] = v;
    }
}

template<class MT>
int TCoordMatrix<MT>::Insert (int r, int c, MT v)
{
    int i;

    if (this->nval == this->nbuf) { // no space left - reallocation required
        idxtype *tmp_rowidx = new idxtype[this->nbuf+BUFFER_CHUNK_SIZE];
	idxtype *tmp_colidx = new idxtype[this->nbuf+BUFFER_CHUNK_SIZE];
	for (i = 0; i < this->nval; i++) tmp_rowidx[i] = rowidx[i];
	for (i = 0; i < this->nval; i++) tmp_colidx[i] = colidx[i];
	if (this->nbuf) {
	    delete []rowidx;
	    delete []colidx;
	}
	rowidx = tmp_rowidx;
	colidx = tmp_colidx;
    }

    // add new entry to the end of the list
    rowidx[this->nval] = r;
    colidx[this->nval] = c;
    TGenericSparseMatrix<MT>::Append (v);
    return this->nval-1;
}

template<class MT>
TCoordMatrix<MT> cath (const TCoordMatrix<MT> &A,
		       const TCoordMatrix<MT> &B)
{
    // Concatenates matrices A and B horizontally
    // Note: this creates an unsorted matrix. May be better to preserve
    // sorting status if A and B are already sorted

    int i, nr, nc, jofs, nzA, nzB, nz;
    nr = A.nRows();
    nc = A.nCols()+B.nCols();
    jofs = A.nCols();
    xASSERT (nr == B.nRows(), "Matrix row dimensions do not match");

    const MT *Aval = A.ValPtr();
    const MT *Bval = B.ValPtr();
    nzA = A.nVal();
    nzB = B.nVal();
    nz  = nzA+nzB;

    idxtype *ri = new idxtype[nz];
    idxtype *ci = new idxtype[nz];
    MT  *v  = new MT[nz];

    // copy A
    for (i = 0; i < nzA; i++) {
	ri[i] = A.rowidx[i];
	ci[i] = A.colidx[i];
	v[i]  = Aval[i];
    }
    // copy B with column offset
    for (i = 0; i < nzB; i++) {
	ri[i+nzA] = B.rowidx[i];
	ci[i+nzA] = B.colidx[i]+jofs;
	v[i+nzA]  = Bval[i];
    }
    TCoordMatrix<MT> C(nr, nc, nz, ri, ci, v);
    delete []ri;
    delete []ci;
    delete []v;
    return C;
}

template<class MT>
TCoordMatrix<MT> catv (const TCoordMatrix<MT> &A,
		       const TCoordMatrix<MT> &B)
{
    // Concatenates matrices A and B vertically
    // If A and B are already row-major-sorted, result should be sorted
    // as well

    int i, nrA, nrB, nr, nc, nzA, nzB, nz;
    nrA = A.nRows();
    nrB = B.nRows();
    nr  = nrA+nrB;
    nc  = A.nCols();
    xASSERT (nc == B.nCols(), "Matrix column dimensions do not match");

    const MT *Aval = A.ValPtr();
    const MT *Bval = B.ValPtr();
    nzA = A.nVal();
    nzB = B.nVal();
    nz  = nzA+nzB;

    idxtype *ri = new idxtype[nz];
    idxtype *ci = new idxtype[nz];
    MT  *v  = new MT[nz];

    // copy A
    for (i = 0; i < nzA; i++) {
	ri[i] = A.rowidx[i];
	ci[i] = A.colidx[i];
	v[i]  = Aval[i];
    }
    // copy B with row offset
    for (i = 0; i < nzB; i++) {
	ri[i+nzA] = B.rowidx[i]+nrA;
	ci[i+nzA] = B.colidx[i];
	v[i+nzA]  = Bval[i];
    }
    TCoordMatrix<MT> C(nr, nc, nz, ri, ci, v);
    delete []ri;
    delete []ci;
    delete []v;
    return C;
}

template<class MT>
istream &operator>> (istream &is, TCoordMatrix<MT> &m)
{
    char cbuf[256], idstr[100];
    int i, j, k, nr, nc, nz;
    MT v;

    if (!is.getline(cbuf, 256)) return is;
    if (sscanf (cbuf, "%s%d%d%d", idstr, &nr, &nc, &nz) != 4 ||
	strcmp (idstr, "CoordMat")) return is;
    idxtype *tmp_rowidx = new idxtype[nz];
    idxtype *tmp_colidx = new idxtype[nz];
    MT  *tmp_data   = new MT[nz];
    for (k = 0; k < nz; k++) {
        is >> i >> j >> v;
	is.getline (cbuf, 256); // skip to eol
	tmp_rowidx[k] = i;
	tmp_colidx[k] = j;
	tmp_data[k] = v;
    }
    m.New (nr, nc);
    m.Initialise (nz, tmp_rowidx, tmp_colidx, tmp_data);
    return is;
}

template<class MT>
ostream &operator<< (ostream &os, const TCoordMatrix<MT> &m)
{
    os << "CoordMat " << m.rows << ' ' << m.cols << ' ' << m.nval << endl;
    for (int k = 0; k < m.nval; k++)
        os << m.rowidx[k] << ' ' << m.colidx[k] << ' ' << m.val[k] << endl;
    return os;
}

// ==========================================================================
// class and friend instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TCoordMatrix<double>;
template class MATHLIB TCoordMatrix<float>;
template class MATHLIB TCoordMatrix<toast::complex>;
template class MATHLIB TCoordMatrix<scomplex>;
template class MATHLIB TCoordMatrix<int>;

template MATHLIB TCoordMatrix<double> cath (const TCoordMatrix<double> &A,
    const TCoordMatrix<double> &B);
template MATHLIB TCoordMatrix<toast::complex> cath (const TCoordMatrix<toast::complex> &A,
    const TCoordMatrix<toast::complex> &B);

template MATHLIB TCoordMatrix<double> catv (const TCoordMatrix<double> &A,
	const TCoordMatrix<double> &B);
template MATHLIB TCoordMatrix<toast::complex> catv (const TCoordMatrix<toast::complex> &A,
	const TCoordMatrix<toast::complex> &B);

template MATHLIB istream &operator>> (istream &is, RCoordMatrix &m);
template MATHLIB istream &operator>> (istream &is, FCoordMatrix &m);
template MATHLIB istream &operator>> (istream &is, CCoordMatrix &m);
template MATHLIB istream &operator>> (istream &is, ICoordMatrix &m);

template MATHLIB ostream &operator<< (ostream &os, const RCoordMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const FCoordMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const CCoordMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const ICoordMatrix &m);

#endif // NEED_EXPLICIT_INSTANTIATION
