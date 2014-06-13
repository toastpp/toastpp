// ==========================================================================
// Module mathlib
// File scrmatrix.cc
// Definition of template class TSymCompRowMatrix ('template symmetric
// compressed-row matrix')
//
// Notes:
// This class only stores the lower triangle including diagonal.
// It is efficient for matrix x vector operations but not very
// efficient to extract rows or columns.
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include <iomanip>
#include <string.h>
#include "mathlib.h"

using namespace std;

template<class MT>
void TSymCompRowMatrix<MT>::New (int nrows, int ncols)
{
    if (this->nval) delete []colidx;
    if (this->rows != nrows) {
        delete []rowptr;
	rowptr = new idxtype[nrows+1];
    }
    for (int i = 0; i <= nrows; i++) rowptr[i] = 0;
    TGenericSparseMatrix<MT>::New (nrows, ncols);
}

template<class MT>
TSymCompRowMatrix<MT>::TSymCompRowMatrix ()
    : TGenericSparseMatrix<MT> ()
{
    rowptr = new idxtype[1];
    rowptr[0] = 0;

    // column access is off on startup
    allow_col_indexing = false;
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;
}

template<class MT>
TSymCompRowMatrix<MT>::TSymCompRowMatrix (int rows, int cols)
    : TGenericSparseMatrix<MT> (rows, cols)
{
    xASSERT(rows == cols, "Invalid arguments");
    rowptr = new idxtype[rows+1];
    for (int i = 0; i <= rows; i++) rowptr[i] = 0;

    // column access is off on startup
    allow_col_indexing = false;
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;
}

template<class MT>
TSymCompRowMatrix<MT>::TSymCompRowMatrix (int rows, int cols,
    const int *rptr, const int *cidx, const MT *data)
    : TGenericSparseMatrix<MT> (rows, cols)
{
    xASSERT(rows == cols, "Invalid arguments");
    rowptr = new idxtype[rows+1];

    // column access is off on startup
    allow_col_indexing = false;
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    Initialise (rptr, cidx, data);
}

#ifdef UNDEF
template<class MT>
TSymCompRowMatrix<MT>::~TSymCompRowMatrix()
{
    SetColAccess (false);
    delete []rowptr;
    if (nval) delete []colidx;
}
#endif

template<class MT>
void TSymCompRowMatrix<MT>::Initialise (const int *_rowptr, const int *_colidx,
    const MT *data)
{
    SetColAccess (false);

    int r, c, i, j, nz;
    int *nv = new int[this->rows];
    
    for (i = 0; i < this->rows; i++) nv[i] = 0;

    // calculate nonzeros in lower triangle
    for (r = 0; r < this->rows; r++) {
	for (i = _rowptr[r]; i < _rowptr[r+1]; i++) {
	    c = _colidx[i];
	    if (c <= r) nv[r]++; // ignore upper triangle
	}
    }

    rowptr[0] = 0;
    for (i = 0; i < this->rows; i++) rowptr[i+1] = rowptr[i] + nv[i];
    nz = rowptr[this->rows];
    if (nz != this->nval) {
	if (this->nval) delete []colidx;
	if (nz) colidx = new idxtype[nz];
    }
    TGenericSparseMatrix<MT>::Initialise (nz, 0);

    for (r = j = 0; r < this->rows; r++) {
	for (i = _rowptr[r]; i < _rowptr[r+1]; i++) {
	    c = _colidx[i];
	    if (c <= r) {
		colidx[j] = c;
		if (data) this->val[j] = data[i];
		j++;
	    }
	}
    }
    delete []nv;
}

template<class MT>
void TSymCompRowMatrix<MT>::AllowColIndexing (bool yes)
{
    if (yes == allow_col_indexing) return; // nothing to do
    allow_col_indexing = yes;
    if (!allow_col_indexing) SetColAccess (false);
}


template<class MT>
MT &TSymCompRowMatrix<MT>::operator() (int r, int c)
{
    static MT dummy;
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    if (r < c) { int tmp = r; r = c; c = tmp; }
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return this->val[rp];
    xERROR("Attempt to access non-existing entry");
    return dummy;
}

template<class MT>
MT TSymCompRowMatrix<MT>::Get (int r, int c) const
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    if (r < c) { int tmp = r; r = c; c = tmp; }
    const static MT zero = (MT)0;
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return this->val[rp];
    return zero;
}

template<class MT>
int TSymCompRowMatrix<MT>::Get_index (int r, int c) const
{
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return rp;
    return -1;
}

template<class MT>
TVector<MT> TSymCompRowMatrix<MT>::Row (int r) const
{
    SetColAccess (true);
    // make sure column index arrays are valid

    int i;
    TVector<MT> tmp (this->cols);

    // lower triangle
    for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	tmp[colidx[i]] = this->val[i];
    }

    if (col_access) {

	// upper triangle
	for (i = colptr[r]; i < colptr[r+1]; i++) {
	    tmp[rowidx[i]] = this->val[vofs[i]];
	}

    } else { // need to do it the hard way

	// now we need to search the upper triangular segment of the row
	int r2;
	for (r2 = r+1; r2 < this->rows; r2++) {
	    for (i = rowptr[r2]; i < rowptr[r2+1]; i++) {
		if (colidx[i] == r) {
		    tmp[r2] = this->val[i];
		    break;
		}
	    }
	}
    }

    return tmp;
}

template<class MT>
int TSymCompRowMatrix<MT>::SparseRow (int r, idxtype *ci, MT *rv) const
{
    TVector<MT> row = Row(r);
    int i, nz;
    for (i = nz = 0; i < this->cols; i++)
	if (row[i] != (MT)0) {
	    ci[nz] = i;
	    rv[nz] = row[i];
	    nz++;
	}
    return nz;
}

// ==========================================================================

template<class MT>
MT TSymCompRowMatrix<MT>::GetNext (int &r, int &c) const
{
    if (r >= 0) {
        if (++iterator_pos >= this->nval) { // end reached
	    r = -1;
	    return (MT)0; // dummy
	}
    } else {
        if (!this->nval) return (MT)0;    // empty matrix
        iterator_pos = r = 0;
    }
    while (iterator_pos >= rowptr[r+1]) r++;
    c = colidx[iterator_pos];
    return this->val[iterator_pos];
}

template<class MT>
void TSymCompRowMatrix<MT>::SetColAccess (bool yes) const
{
    if (yes) {
	if (!allow_col_indexing) return; // column indexing is disabled
	if (col_access) return; // already initialised - nothing to do
	int i, c, r, ri;
	
	// column entry count
	int *nval_col = new int[this->cols];
	for (c = 0; c < this->cols; c++) nval_col[c] = 0;
	for (i = 0; i < this->nval; i++) nval_col[colidx[i]]++;

	// init column pointer array
	colptr = new idxtype[this->cols+1];
	for (c = 0, colptr[0] = 0; c < this->cols; c++)
	    colptr[c+1] = colptr[c] + nval_col[c];
	
	// init row index offset arrays
	rowidx = new idxtype[this->nval];
	vofs = new idxtype[this->nval];
	for (c = 0; c < this->cols; c++) nval_col[c] = 0;
	for (r = 0; r < this->rows; r++) {
	    for (ri = rowptr[r]; ri < rowptr[r+1]; ri++) {
		c = colidx[ri];
		rowidx[colptr[c]+nval_col[c]] = r;
		vofs[colptr[c]+nval_col[c]] = ri;
		nval_col[c]++;
	    }
	}
	delete []nval_col;
	col_access = true;
    } else {
	if (colptr) delete []colptr;
	colptr = 0;
	if (rowidx) delete []rowidx;
	rowidx = 0;
	if (vofs) delete []vofs;
	vofs = 0;
	col_access = false;
    }
}

// ==========================================================================

template<class MT>
void TSymCompRowMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->cols,
	"Parameter 1 invalid size %d (expected %d)", x.Dim(), this->cols);
    if (b.Dim() != this->rows) b.New (this->rows);
    else                 b.Clear();

    int r, c, i;

    for (r = 0; r < this->rows; r++) {
	for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	    c = colidx[i];
	    b[r] += this->val[i] * x[c];
	    if (r > c) b[c] += this->val[i] * x[r];
	}
    }
}

template<class MT>
void TSymCompRowMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b,
    int r1, int r2) const
{
    dASSERT(x.Dim() == this->cols,
	"Parameter 1 invalid size %d (expected %d)", x.Dim(), this->cols);
    if (b.Dim() != this->rows) b.New (this->rows);
    else b.Clear();

    int r, c, i;

    for (r = r1; r < r2; r++) {
	for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	    c = colidx[i];
	    b[r] += this->val[i] * x[c];
	    if (r > c) b[c] += this->val[i] * x[r];
	}
    }
}

template<class MT>
void TSymCompRowMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    Ax (x, b);
}

template<>
void TSymCompRowMatrix<std::complex<double> >::ATx (
    const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    dASSERT(x.Dim() == cols,
        "Parameter 1 invalid size %d (expected %d)", x.Dim(), cols);
    if (b.Dim() != rows) b.New (rows);
    else                 b.Clear();

    int r, c, i;

    for (r = 0; r < rows; r++) {
	for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	    c = colidx[i];
	    b[r] += conj(val[i]) * x[c];
	    if (r > c) b[c] += conj(val[i]) * x[r];
	}
    }
}

template<class MT>
int TSymCompRowMatrix<MT>::Shrink ()
{
    int i, r, idx, nz = 0;
	idxtype *rp, *ci;
    MT *v;
    for (i = 0; i < this->nval; i++)
        if (this->val[i] != (MT)0) nz++;
    rp = new idxtype[this->rows+1];
    ci = new idxtype[nz];
    v  = new MT[nz];
    rp[0] = 0;
    for (r = idx = 0; r < this->rows; r++) {
        for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	    if (this->val[i] != (MT)0) {
	        v[idx]  = this->val[i];
		ci[idx] = colidx[i];
		idx++;
	    }
	}
	rp[r+1] = idx;
    }
    delete []rowptr;
    delete []colidx;
    delete []this->val;
    rowptr = rp;
    colidx = ci;
    this->val = v;
    i      = this->nval - nz;
    this->nval = nz;
    return i;
}

template<class MT>
void TSymCompRowMatrix<MT>::SymbolicCholeskyFactorize (idxtype *&frowptr,
    idxtype *&fcolidx) const
{
    symbolic_cholesky_factor (this->rows, rowptr, colidx, frowptr, fcolidx);
    // implemented in cr_cholesky.cc
}

template<class MT>
bool CholeskyFactorize (const TSymCompRowMatrix<MT> &A, TCompRowMatrix<MT> &L,
    TVector<MT> &d, bool recover)
{
    // This assumes that entries of L are sorted COLUMN-wise, which should
    // be the case if L has been initialised via SymbolicCholeskyFactorize
    // Do NOT call L.Sort since this will sort row-wise!

    return false; // finish this

#ifdef UNDEF

    const double EPS = 1e-10;
    int i, k, c, r, cc, term, n = A.nRows();
    MT idiag, ajk, *fullcol = new MT[n];
    bool ok = true;
    //if (!A.col_access) A.SetColAccess();
    if (!L.col_access) L.SetColAccess();
    // some shortcuts to improve performance in the inner loops
    int    *Lrowidx = L.rowidx;
    int    *Lvofs   = L.vofs;
    MT     *Lval    = L.val;
    for (c = 0; c < n; c++) {
	// expand current column
        memset (fullcol, 0, n*sizeof(double));
	for (i = A.rowptr[c]; i < A.rowptr[c+1]; i++)
	    fullcol[r] = A.val[i];

	// loop over left columns
	for (k = L.rowptr[c]; k < L.rowptr[c+1]; k++) {
	    if ((cc = L.colidx[k]) >= c) continue;
#ifdef FEM_DEBUG
	    bool ajk_set = false;
#endif
	    for (i = L.colptr[cc], term = L.colptr[cc+1]; i < term; i++) {
		if ((r = Lrowidx[i]) < c) continue;
		if (r == c) {
		    ajk = Lval[Lvofs[i]];
#ifdef FEM_DEBUG
		    ajk_set = true;
#endif
		}
		dASSERT(ajk_set, "Cholesky factor not column-sorted");
		fullcol[r] -= ajk * Lval[Lvofs[i]];
	    }
	}
	// diagonal element
	if (fullcol[c] > MT(0)) {   /* problem here, if using complex */
	    d[c] = sqrt (fullcol[c]);
	} else {
	    if (!recover) xERROR("Matrix not positive definite");
	    ok = false;
	    d[c] = EPS;
	}
	idiag = MT(1)/d[c];
	// scale column with diagonal
	for (i = L.colptr[c], term = L.colptr[c+1]; i < term; i++) {
	    r = Lrowidx[i];
	    Lval[Lvofs[i]] = fullcol[r] * idiag;
	}
    }
    delete []fullcol;
    return ok;
#endif
}

template<class MT>
istream &operator>> (istream &is, TSymCompRowMatrix<MT> &m)
{
    char cbuf[256];
    int i, nr, nc, nz;

    is.getline (cbuf, 256);
    if (strncmp (cbuf, "TSymCompRow", 11)) return is; // should set error flag
    sscanf (cbuf+11, "%d%d%d", &nr, &nc, &nz);
    if (nr == m.rows && nc == m.cols && nz == m.nval) { // no need to realloc
        for (i = 0; i <= nr; i++)
	    is >> m.rowptr[i];
	for (i = 0; i < nz; i++)
	    is >> m.colidx[i];
	for (i = 0; i < nz; i++)
	    is >> m.val[i];
    } else {
        int *rowptr = new int[nr+1];
	int *colidx = new int[nz];
	MT  *val    = new MT[nz];
        for (i = 0; i <= nr; i++)
	    is >> rowptr[i];
	for (i = 0; i < nz; i++)
	    is >> colidx[i];
	for (i = 0; i < nz; i++)
	    is >> val[i];
	m.New (nr, nc);
	m.Initialise (rowptr, colidx, val);
	delete []rowptr;
	delete []colidx;
	delete []val;
    }
    return is;
}

template<class MT>
ostream &operator<< (ostream &os, const TSymCompRowMatrix<MT> &m)
{
    int i;

    // header
    os << "TSymCompRow " << m.rows << ' ' << m.cols << ' ' << m.nval << endl;
    // row pointer list
    for (i = 0; i <= m.rows; i++)
        os << m.rowptr[i] << ' ';
    os << endl;
    // col index list
    for (i = 0; i < m.nval; i++)
        os << m.colidx[i] << ' ';
    os << endl;
    // data list
    for (i = 0; i < m.nval; i++)
        os << m.val[i] << ' ';
    os << endl;
    return os;
}


// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TSymCompRowMatrix<double>;
template class TSymCompRowMatrix<float>;
template class TSymCompRowMatrix<std::complex<double> >;
template class TSymCompRowMatrix<std::complex<float> >;
template class TSymCompRowMatrix<int>;

template bool CholeskyFactorize (const RSymCompRowMatrix &A, RCompRowMatrix &L,
    RVector &d, bool recover);

template istream &operator>> (istream &is, RSymCompRowMatrix &m);
template istream &operator>> (istream &is, CSymCompRowMatrix &m);
template istream &operator>> (istream &is, SCSymCompRowMatrix &m);
template ostream &operator<< (ostream &os, const RSymCompRowMatrix &m);
template ostream &operator<< (ostream &os, const CSymCompRowMatrix &m);
template ostream &operator<< (ostream &os, const SCSymCompRowMatrix &m);

#endif // NEED_EXPLICIT_INSTANTIATION
