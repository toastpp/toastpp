// ==========================================================================
// Module mathlib
// File crmatrix.cc
// Definition of template class TCompRowMatrix ('template compressed-row
//  matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __CRMATRIX_CC

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mathlib.h"
#include "ilutoast.h"

#ifdef USE_CUDA_FLOAT
#include "toastcuda.h"
#include "toastspmv.h"
#endif

#ifdef ML_INTERFACE
#include "ml_defs.h"
#include "ml_operator.h"
#endif // ML_INTERFACE

using namespace std;

// ==========================================================================
// member definitions

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix ()
  : TGenericSparseMatrix<MT> ()
{
    rowptr = new idxtype[1];
    rowptr[0] = 0;

    // diagonal access is off on startup
    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    // we regard a 0x0 matrix as sorted by definition
    sorted = true;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (int rows, int cols)
  : TGenericSparseMatrix<MT> (rows, cols)
{
    rowptr = new idxtype[rows+1];
    for (int i = 0; i <= rows; i++)
        rowptr[i] = 0;

    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    // we regard a matrix without entries as sorted by definition
    sorted = true;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (int rows, int cols,
    const idxtype *_rowptr, const idxtype *_colidx, const MT *data)
  : TGenericSparseMatrix<MT> (rows, cols, _rowptr[rows], data)
{
    int i;

    rowptr = new idxtype[rows+1];
    for (i = 0; i <= rows; i++) rowptr[i] = _rowptr[i];
    if (this->nval) {
        colidx = new idxtype[this->nval];
	for (i = 0; i < this->nval; i++) colidx[i] = _colidx[i];
    }

    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    // we don't know the sort status, so assume unsorted
    sorted = false;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (int rows, int cols,
    idxtype *_rowptr, idxtype *_colidx, MT *data,
    CopyMode cmode)
    : TGenericSparseMatrix<MT> (rows, cols, _rowptr[rows], data, cmode)
{
    if (cmode == DEEP_COPY) {
	int i;

	rowptr = new idxtype[rows+1];
	for (i = 0; i <= rows; i++) rowptr[i] = _rowptr[i];
	if (this->nval) {
	    colidx = new idxtype[this->nval];
	    for (i = 0; i < this->nval; i++) colidx[i] = _colidx[i];
	}
    } else {
	rowptr = _rowptr;
	colidx = _colidx;
    }

    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    // we don't know the sort status, so assume unsorted
    sorted = false;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (const TCompRowMatrix<MT> &m)
  : TGenericSparseMatrix<MT> (m)
{
    int i;

    rowptr = new idxtype[this->rows+1];
    for (i = 0; i <= this->rows; i++) rowptr[i] = m.rowptr[i];
    if (this->nval) {
        colidx = new idxtype[this->nval];
	for (i = 0; i < this->nval; i++) colidx[i] = m.colidx[i];
    }

    if ((diag_access = m.diag_access)) {
        diagptr = new int[this->rows];
	for (i = 0; i < this->rows; i++) diagptr[i] = m.diagptr[i];
    } else {
        diagptr = 0;
    }

    if ((col_access = m.col_access)) {
        colptr = new idxtype[this->cols+1];
	for (i = 0; i <= this->cols; i++) colptr[i] = m.colptr[i];
	rowidx = new idxtype[this->nval];
	vofs   = new idxtype[this->nval];
	for (i = 0; i < this->nval; i++) rowidx[i] = m.rowidx[i];
	for (i = 0; i < this->nval; i++) vofs[i]   = m.vofs[i];
    } else {
        colptr = 0;
	rowidx = 0;
	vofs   = 0;
    }
    sorted = m.sorted;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (const TCoordMatrix<MT> &m)
  : TGenericSparseMatrix<MT> (m)
{
    int i, r, *row_nval = new int[this->rows];
    for (r = 0; r < this->rows; r++) row_nval[r] = 0;
    for (i = 0; i < this->nval; i++) row_nval[m.RowIndex(i)]++;

    rowptr = new idxtype[this->rows+1];
    colidx = new idxtype[this->nval];

    for (r = 0, rowptr[0] = 0; r < this->rows; r++)
        rowptr[r+1] = rowptr[r] + row_nval[r];

    for (r = 0; r < this->rows; r++) row_nval[r] = 0;
    for (i = 0; i < this->nval; i++) {
        r = m.RowIndex(i);
	colidx[rowptr[r]+row_nval[r]] = m.ColIndex(i);
	this->val[rowptr[r]+row_nval[r]] = m.Val(i);
	row_nval[r]++;
    }
    delete []row_nval;

    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    // this should really be sorted = m.sorted
    sorted = false;
}

template<class MT>
TCompRowMatrix<MT>::TCompRowMatrix (const TDenseMatrix<MT> &m)
  : TGenericSparseMatrix<MT> (m.nRows(), m.nCols())
{
    int i, j;
    this->nval = this->rows*this->cols;
    rowptr = new idxtype[this->rows+1];
    colidx = new idxtype[this->nval];
    this->val = new MT[this->nval];
    rowptr[0] = 0;
    for (i = 0; i < this->rows; i++) {
        rowptr[i+1] = rowptr[i]+this->cols;
	for (j = 0; j < this->cols; j++) {
	    colidx[i*this->cols+j] = j;
	    this->val[i*this->cols+j] = m(i,j);
	}
    }

    diag_access = false;
    diagptr = 0;

    // column access is off on startup
    col_access = false;
    colptr = 0;
    rowidx = 0;
    vofs   = 0;

    sorted = true;
}

template<class MT>
TCompRowMatrix<MT>::~TCompRowMatrix ()
{
    SetColAccess (false); // deallocate column access index lists
    if (this->nbuf) { // otherwise we are using external buffers
	delete []rowptr;
	if (this->nval) delete []colidx;
    }
}

template<class MT>
void TCompRowMatrix<MT>::New (int nrows, int ncols)
{
    SetColAccess (false); // deallocate column access index lists
    if (this->nbuf) {
	if (this->nval) delete []colidx;
	if (this->rows != nrows) {
	    delete []rowptr;
	    rowptr = new idxtype[nrows+1];
	}
    } else {
	rowptr = new idxtype[nrows+1];
    }
    for (int i = 0; i <= nrows; i++) rowptr[i] = 0;
    TGenericSparseMatrix<MT>::New (nrows, ncols);
    sorted = true;
}

template<class MT>
void TCompRowMatrix<MT>::Identity (int n)
{
    int i;
    idxtype *rp = new idxtype[n+1];
    idxtype *ci = new idxtype[n];
    MT  *v  = new MT[n];
    for (i = 0; i <= n; i++)
	rp[i] = i;
    for (i = 0; i < n; i++) {
	ci[i] = i;
	v[i] = (MT)1;
    }
    New (n, n);
    Initialise (rp, ci, v);
    delete []rp;
    delete []ci;
    delete []v;
}

template<class MT>
void TCompRowMatrix<MT>::DiagV (const TVector<MT> &x)
{
    int i,n = x.Dim();
    idxtype *rp = new idxtype[n+1];
    idxtype *ci = new idxtype[n];
    MT  *v  = new MT[n];
    for (i = 0; i <= n; i++)
	rp[i] = i;
    for (i = 0; i < n; i++) {
	ci[i] = i;
	v[i] = x[i];
    }
    New (n, n);
    Initialise (rp, ci, v);
    delete []rp;
    delete []ci;
    delete []v;
}

template<class MT>
TCoordMatrix<MT> TCompRowMatrix<MT>::MakeCoordMatrix () const
{
    int i, j;
    idxtype *cd_rowidx = new idxtype[this->nval];
    idxtype *cd_colidx = new idxtype[this->nval];
    MT  *cd_val    = new MT[this->nval];
    for (i = 0; i < this->nval; i++) {
        cd_colidx[i] = colidx[i];
	cd_val[i]    = this->val[i];
    }
    for (i = 0; i < this->rows; i++) {
        for (j = rowptr[i]; j < rowptr[i+1]; j++)
	    cd_rowidx[j] = i;
    }
    return TCoordMatrix<MT> (this->rows, this->cols, this->nval,
        cd_rowidx, cd_colidx, cd_val);
}

template<class MT>
void TCompRowMatrix<MT>::Unlink ()
{
    SetColAccess (false); // deallocate column access index lists
    delete []rowptr;
    if (this->nval) delete []colidx;
    TGenericSparseMatrix<MT>::Unlink();
}

inline TCompRowMatrix<std::complex<double> > cplx (
    const TCompRowMatrix<double> &A)
{
    TCompRowMatrix<std::complex<double> > C(A.nRows(), A.nCols(), A.rowptr,
					    A.colidx);
    for (int i = 0; i < A.nval; i++)
	C.val[i] = A.val[i];
    return C;
}

template<class MT>
void TCompRowMatrix<MT>::Initialise (const idxtype *_rowptr, const idxtype *_colidx,
    const MT *data)
{
    int i, nv = _rowptr[this->rows];

    SetColAccess (false);
    if (nv != this->nval) {
        if (this->nval) delete []colidx;
	if (nv) colidx = new idxtype[nv];
    }
    for (i = 0; i < nv; i++) colidx[i] = _colidx[i];
    for (i = 0; i <= this->rows; i++) rowptr[i] = _rowptr[i];
    TGenericSparseMatrix<MT>::Initialise (nv, data);
    sorted = false;
}

template<class MT>
TCompRowMatrix<MT> TCompRowMatrix<MT>::Submatrix (int r1, int r2,
    int c1, int c2) const
{
    // pass 1: scan non-zeros
    int i, r, nv = 0;
    for (i = rowptr[r1]; i < rowptr[r2]; i++)
	if (colidx[i] >= c1 && colidx[i] < c2) nv++;

    // pass 2: create index and value lists
    idxtype *rp = new idxtype[r2-r1+1];
    idxtype *ci = new idxtype[nv];
    MT *v = new MT[nv];
    rp[0] = 0;
    for (r = r1, nv = 0; r < r2; r++) {
	rp[r-r1+1] = rp[r-r1];
	for (i = rowptr[r]; i < rowptr[r+1]; i++) {
	    if (colidx[i] >= c1 && colidx[i] < c2) {
		ci[nv] = colidx[i]-c1;
		v[nv] = this->val[i];
		rp[r-r1+1]++;
		nv++;
	    }
	}
    }

    // create matrix
    TCompRowMatrix A(r2-r1, c2-c1, rp, ci, v);
    
    // clean-up
    delete []rp;
    delete []ci;
    delete []v;

    return A;
}

template<class MT>
TCompRowMatrix<MT> cath (const TCompRowMatrix<MT> &A,
			 const TCompRowMatrix<MT> &B)
{
    // Concatenates matrices A and B horizontally

    int r, i, j, nr, nc, jofs, ncrA, ncrB, ncr, nzA, nzB, nz;
    nr = A.nRows();
    nc = A.nCols()+B.nCols();
    jofs = A.nCols();
    xASSERT (nr == B.nRows(), "Matrix row dimensions do not match");

    nzA = A.nVal();
    nzB = B.nVal();
    nz  = nzA+nzB;
    const MT *Aval = A.ValPtr();
    const MT *Bval = B.ValPtr();

    idxtype *rp = new idxtype[nr+1];
    idxtype *ci = new idxtype[nz];
    MT  *v  = new MT[nz];

    rp[0] = 0;
    for (r = i = 0; r < nr; r++) {
	ncrA = A.rowptr[r+1]-A.rowptr[r];
	ncrB = B.rowptr[r+1]-B.rowptr[r];
	ncr = ncrA+ncrB;
	rp[r+1] = rp[r] + ncr;
	for (j = A.rowptr[r]; j < A.rowptr[r+1]; j++) {
	    ci[i] = A.colidx[j];
	    v[i]  = Aval[j];
	    i++;
	}
	for (j = B.rowptr[r]; j < B.rowptr[r+1]; j++) {
	    ci[i] = B.colidx[j]+jofs;
	    v[i]  = Bval[j];
	    i++;
	}
    }
    TCompRowMatrix<MT> C(nr, nc, rp, ci, v);
    delete []rp;
    delete []ci;
    delete []v;
    return C;
}

template<class MT>
TCompRowMatrix<MT> catv (const TCompRowMatrix<MT> &A,
			 const TCompRowMatrix<MT> &B)
{
    // Concatenates matrices A and B vertically

    int nrA, nrB, nr, nc, nzA, nzB, nz, r, i;
    nrA = A.nRows();
    nrB = B.nRows();
    nr = nrA+nrB;
    nc = A.nCols();
    xASSERT (nc == B.nCols(), "Matrix column dimensions do not match");
    
    nzA = A.nVal();
    nzB = B.nVal();
    nz  = nzA+nzB;
    const MT *Aval = A.ValPtr();
    const MT *Bval = B.ValPtr();

    idxtype *rp = new idxtype[nr+1];
    idxtype *ci = new idxtype[nz];
    MT  *v  = new MT[nz];

    for (r = 0; r <= nrA; r++)
	rp[r] = A.rowptr[r];
    for (r = 1; r <= nrB; r++)
	rp[r+nrA] = B.rowptr[r]+nzA;

    for (i = 0; i < nzA; i++) {
	ci[i] = A.colidx[i];
	v[i]  = Aval[i];
    }
    for (i = 0; i < nzB; i++) {
	ci[i+nzA] = B.colidx[i];
	v[i+nzA]  = Bval[i];
    }
    TCompRowMatrix<MT> C(nr, nc, rp, ci, v);
    delete []rp;
    delete []ci;
    delete []v;
    return C;
}

template<class MT>
TCompRowMatrix<MT> &TCompRowMatrix<MT>::merge (const TCompRowMatrix<MT> &m)
{
    // first sort index lists
    if (!sorted) Sort();
    if (!m.sorted) m.Sort();

    int i, j, r, c1a, c2a, c1b, c2b, cola, colb;

    int max_nval = this->nval + m.nval; // max length of merged data array
    int tmp_rows = this->rows;
    if (m.rows > tmp_rows) tmp_rows = m.rows;
    idxtype *tmp_rowptr = new idxtype[tmp_rows+1];
    idxtype *tmp_colidx = new idxtype[max_nval];
    int tmp_nval = 0;
    MT *tmp_val = new MT[max_nval];
    tmp_rowptr[0] = 0;

    for (r = 0; r < tmp_rows; r++) {
        c1a = (r < this->rows ? rowptr[r] : rowptr[this->rows]);
	c2a = (r < this->rows ? rowptr[r+1] : rowptr[this->rows]);
	c1b = (r < m.rows ? m.rowptr[r] : m.rowptr[m.rows]);
	c2b = (r < m.rows ? m.rowptr[r+1] : m.rowptr[m.rows]);

	for (i = c1a, j = c1b; i < c2a && j < c2b;) {
	    cola = colidx[i];
	    colb = m.colidx[j];
	    if (cola < colb) {
	        tmp_colidx[tmp_nval] = cola;
		tmp_val[tmp_nval++] = this->val[i++];
	    } else if (cola > colb) {
	        tmp_colidx[tmp_nval] = colb;
		tmp_val[tmp_nval++] = m.val[j++];
	    } else {
	        tmp_colidx[tmp_nval] = cola;
		tmp_val[tmp_nval++] = this->val[i++] + m.val[j++];
	    }
	}
	for (; i < c2a; i++) {
	    tmp_colidx[tmp_nval] = colidx[i];
	    tmp_val[tmp_nval++] = this->val[i];
	}
	for (; j < c2b; j++) {
	    tmp_colidx[tmp_nval] = m.colidx[j];
	    tmp_val[tmp_nval++] = m.val[j];
	}
	tmp_rowptr[r+1] = tmp_nval;
    }
    Initialise (tmp_rowptr, tmp_colidx, tmp_val);
    sorted = true; // due to the assembly mechanism
    delete []tmp_rowptr;
    delete []tmp_colidx;
    delete []tmp_val;
    return *this;
}

template<class MT>
TCompRowMatrix<MT> &TCompRowMatrix<MT>::operator= (const TCompRowMatrix<MT> &m)
{
    this->cols = m.cols;
    if (this->rows != m.rows) {
        delete []rowptr;
	rowptr = new idxtype[(this->rows=m.rows)+1];
    }
    Initialise (m.rowptr, m.colidx, m.val);
    sorted = m.sorted;
    return *this;
}

template<class MT>
TCompRowMatrix<MT> &TCompRowMatrix<MT>::operator= (const TCoordMatrix<MT> &m)
{
    m.Sort (true);
    int r, nc, i;

    this->cols = m.cols;
    if (this->rows != m.rows) {
        delete []rowptr;
	rowptr = new idxtype[(this->rows=m.rows)+1];
    }
    rowptr[0] = 0;
    for (r = i = 0; r < m.rows; r++) {
        nc = 0;
	while (i < m.nval && m.rowidx[i] == r) nc++, i++;
	rowptr[r+1] = rowptr[r]+nc;
    }
    Initialise (rowptr, m.ColIndex(), m.val);
    sorted = true;
    return *this;
}

template<class MT>
TCompRowMatrix<MT> TCompRowMatrix<MT>::operator- () const
{
    TCompRowMatrix<MT> res(*this);
    for (int i = 0; i < this->nval; i++) res.val[i] = -this->val[i];
    return res;
}

template<class MT>
TCompRowMatrix<MT> TCompRowMatrix<MT>::operator* (MT f) const
{
    TCompRowMatrix<MT> res(*this);
    for (int i = 0; i < this->nval; i++)
	res.val[i] *= f;
    return res;
}

template<class MT> TCompRowMatrix<MT>  TCompRowMatrix<MT>::operator*
(const TDiagMatrix<MT> &D) const
{
    // check that matrices are compatible
    dASSERT (D.nRows() == this->cols, "Incompatible operators");

    TCompRowMatrix<MT> res(*this);
    for (int i = 0; i < this->rows; i++) {
	for(int j = rowptr[i]; j< rowptr[i+1] ; j++) {
	    int indx = colidx[j];
	    res.val[indx] *= D.Get(indx,indx);
	}
    }
    return res;
}

template<class MT>
TCompRowMatrix<MT> TCompRowMatrix<MT>::operator+ (const TCompRowMatrix<MT> &m)
    const
{
    // check that matrices are compatible
    dASSERT (m.rows == this->rows && m.cols == this->cols,
	     "Incompatible operators");

    // create the index lists of the result
    idxtype *rrowptr = new idxtype[this->rows+1];
    idxtype **rrow = new idxtype*[this->rows];
    MT  **rval = new MT*[this->rows];
    idxtype *r = new idxtype[this->cols];
    MT *v = new MT[this->cols];
    int *nzi = new int[this->rows];
    int i, j, idx, idxm, nz, rnz, col, colm;
    rrowptr[0] = 0;
    rnz = 0;
    for (i = 0; i < this->rows; i++) {
	idx  = rowptr[i];
	idxm = m.rowptr[i];
	nz = 0;
	while (idx < rowptr[i+1] || idxm < m.rowptr[i+1]) {
	    col  = (idx  < rowptr[i+1]   ? colidx[idx]    : this->cols);
	    colm = (idxm < m.rowptr[i+1] ? m.colidx[idxm] : this->cols);
	    if (col < colm) {
		r[nz] = col;
		v[nz++] = this->val[idx++];
	    } else if (colm < col) {
		r[nz] = colm;
		v[nz++] = m.val[idxm++];
	    } else {
		r[nz] = col;
		v[nz++] = this->val[idx++] + m.val[idxm++];
	    }
	}
	rrowptr[i+1] = rrowptr[i]+nz;
	rrow[i] = new idxtype[nz];
	rval[i] = new MT[nz];
	for (j = 0; j < nz; j++) {
	    rrow[i][j] = r[j];
	    rval[i][j] = v[j];
	}
	rnz += (nzi[i] = nz);
    }

    // now copy column indices and values into linear arrays
    idxtype *rcolidx = new idxtype[rnz];
    MT *rdata = new MT[rnz];
    for (i = idx = 0; i < this->rows; i++) {
	for (j = 0; j < nzi[i]; j++) {
	    rcolidx[idx] = rrow[i][j];
	    rdata[idx] = rval[i][j];
	    idx++;
	}
	delete []rrow[i];
	delete []rval[i];
    }
    delete []rrow;
    delete []rval;
    delete []nzi;
    delete []r;
    delete []v;

    // create result matrix
    TCompRowMatrix<MT> res (this->rows, this->cols, rrowptr, rcolidx, rdata);

    delete []rrowptr;
    delete []rcolidx;
    delete []rdata;
    return res;
}

template<class MT>
TCompRowMatrix<MT> &TCompRowMatrix<MT>::operator+= (
    const TCompRowMatrix<MT> &m)
{
    // check that matrices are compatible
    dASSERT (m.rows == this->rows && m.cols == this->cols,
	     "Incompatible operators");

    // create the index lists of the result
    idxtype *rrowptr = new idxtype[this->rows+1];
    idxtype **rrow = new idxtype*[this->rows];
    MT  **rval = new MT*[this->rows];
    idxtype *r = new idxtype[this->cols];
    MT *v = new MT[this->cols];
    idxtype *nzi = new idxtype[this->rows];
    int i, j, idx, idxm, nz, rnz, col, colm;
    rrowptr[0] = 0;
    rnz = 0;
    for (i = 0; i < this->rows; i++) {
	idx  = rowptr[i];
	idxm = m.rowptr[i];
	nz = 0;
	while (idx < rowptr[i+1] || idxm < m.rowptr[i+1]) {
	    col  = (idx  < rowptr[i+1]   ? colidx[idx]    : this->cols);
	    colm = (idxm < m.rowptr[i+1] ? m.colidx[idxm] : this->cols);
	    if (col < colm) {
		r[nz] = col;
		v[nz++] = this->val[idx++];
	    } else if (colm < col) {
		r[nz] = colm;
		v[nz++] = m.val[idxm++];
	    } else {
		r[nz] = col;
		v[nz++] = this->val[idx++] + m.val[idxm++];
	    }
	}
	rrowptr[i+1] = rrowptr[i]+nz;
	rrow[i] = new idxtype[nz];
	rval[i] = new MT[nz];
	for (j = 0; j < nz; j++) {
	    rrow[i][j] = r[j];
	    rval[i][j] = v[j];
	}
	rnz += (nzi[i] = nz);
    }

    // now copy column indices and values into linear arrays
    idxtype *rcolidx = new idxtype[rnz];
    MT *rdata = new MT[rnz];
    for (i = idx = 0; i < this->rows; i++) {
	for (j = 0; j < nzi[i]; j++) {
	    rcolidx[idx] = rrow[i][j];
	    rdata[idx] = rval[i][j];
	    idx++;
	}
	delete []rrow[i];
	delete []rval[i];
    }
    delete []rrow;
    delete []rval;
    delete []nzi;
    delete []r;
    delete []v;

    // update matrix
    Initialise (rrowptr, rcolidx, rdata);

    delete []rrowptr;
    delete []rcolidx;
    delete []rdata;
    return *this;
}

template<class MT>
MT &TCompRowMatrix<MT>::operator() (int r, int c)
{
    static MT dummy;
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return this->val[rp];
    xERROR("Attempt to access non-existing entry");
    return dummy;
}

template<class MT>
MT TCompRowMatrix<MT>::Get (int r, int c) const
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    const static MT zero = (MT)0;
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return this->val[rp];
    return zero;
}

template<class MT>
bool TCompRowMatrix<MT>::Exists (int r, int c) const
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return true;
    return false;
}

template<class MT>
int TCompRowMatrix<MT>::Get_index (int r, int c) const
{
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return rp;
    return -1;
}

template<class MT>
MT TCompRowMatrix<MT>::Get_sorted (int r, int c) const
{
    const static MT zero = (MT)0;
    int i0 = rowptr[r];
    int i1 = rowptr[r+1];
    if (i0 == i1) return zero;
    int im  = (i0+i1)/2;
    int cim = colidx[im];
    if (cim == c) return this->val[im];
    while (i1-i0 > 1) {
	if (cim < c) i0 = im;
	else i1 = im;
	im = (i0+i1)/2;
	if ((cim = colidx[im]) == c) return this->val[im];
    } 
    return zero;
}

template<class MT>
int TCompRowMatrix<MT>::Get_index_sorted (int r, int c) const
{
    int i0 = rowptr[r];
    int i1 = rowptr[r+1];
    if (i0 == i1) return -1;
    int im  = (i0+i1)/2;
    int cim = colidx[im];
    if (cim == c) return im;
    while (i1-i0 > 1) {
	if (cim < c) i0 = im;
	else i1 = im;
	im = (i0+i1)/2;
	if ((cim = colidx[im]) == c) return im;
    } 
    return -1;
}

template<class MT>
void TCompRowMatrix<MT>::Put_sorted (int r, int c, MT v)
{
    int i0 = rowptr[r];
    int i1 = rowptr[r+1];
    dASSERT (i0 != i1, "Entry does not exist");
    int im  = (i0+i1)/2;
    int cim = colidx[im];
    if (cim == c) { this->val[im] = v; return; }
    while (i1-i0 > 1) {
	if (cim < c) i0 = im;
	else i1 = im;
	im = (i0+i1)/2;
	if ((cim = colidx[im]) == c) { this->val[im] = v; return; }
    } 
    xERROR ("Entry does not exist");
}

template<class MT>
TVector<MT> TCompRowMatrix<MT>::Row (int r) const
{
    int i, i1 = rowptr[r], i2 = rowptr[r+1];
    TVector<MT> tmp (this->cols);
    for (i = i1; i < i2; i++)
        tmp[colidx[i]] = this->val[i];
    return tmp;
}

template<class MT>
TVector<MT> TCompRowMatrix<MT>::Col (int c) const
{
    if (!col_access) SetColAccess();
    TVector<MT> tmp (this->rows);
    int i, i1 = colptr[c], i2 = colptr[c+1];
    for (i = i1; i < i2; i++)
        tmp[rowidx[i]] = this->val[vofs[i]];
    return tmp;
}

template<class MT>
int TCompRowMatrix<MT>::SparseRow (int r, idxtype *ci, MT *rv) const
{
    int i, r0 = rowptr[r], nz = rowptr[r+1]-r0;
    for (i = 0; i < nz; i++) {
        ci[i] = colidx[r0+i];
	rv[i] = this->val[r0+i];
    }
    return nz;
}

template<class MT>
void TCompRowMatrix<MT>::SetRow (int r, const TVector<MT> &row)
{
    dASSERT(r >= 0 && r < this->rows, "Argument 1 out of range");
    dASSERT(row.Dim() == this->cols, "Argument 2 invalid dimension");

    // deflate row
    MT *rval = new MT[this->cols];
    idxtype *rcolidx = new idxtype[this->cols];
    int i, j, rnz = 0;

    for (i = 0; i < this->cols; i++) {
	if (!iszero(row[i])) {
	    rval[rnz] = row[i];
	    rcolidx[rnz] = i;
	    rnz++;
	}
    }

    // reallocate arrays if required
    int prnz = rowptr[r+1] - rowptr[r];
    int nzero = this->nval-prnz+rnz;
    if (nzero != this->nval) {
        idxtype *tmp_colidx = new idxtype[nzero];
	MT  *tmp_val    = new MT[nzero];
	if (this->nval) {
	    memcpy (tmp_colidx, colidx, rowptr[r]*sizeof(idxtype));
	    memcpy (tmp_colidx+rowptr[r]+rnz, colidx+rowptr[r+1],
		    (this->nval-rowptr[r+1])*sizeof(idxtype));
	    delete []colidx;
	    memcpy (tmp_val, this->val, rowptr[r]*sizeof(MT));
	    memcpy (tmp_val+rowptr[r]+rnz, this->val+rowptr[r+1],
		    (this->nval-rowptr[r+1])*sizeof(MT));
	    delete []this->val;
	}
	colidx = tmp_colidx;
	this->val = tmp_val;
    }

    // insert row
    for (i = 0, j = rowptr[r]; i < rnz; i++) {
        colidx[j+i] = rcolidx[i];
	this->val[j+i] = rval[i];
    }

    // re-align row pointers
    if (rnz != prnz) {
        int d = nzero-this->nval;
	for (i = r+1; i <= this->rows; i++) rowptr[i] += d;
	this->nval = nzero;
    }

    delete []rval;
    delete []rcolidx;
}

template<class MT>
void TCompRowMatrix<MT>::SetRows (int r0, const TCompRowMatrix<MT> &rws)
{
    int rr = min (rws.rows, this->rows-r0);   // number of rows to replace
    if (rr <= 0) return;

    int i, ncpy, ofs;
    int nnz_insert = rws.rowptr[rr];         // number of nonzeros to insert
    int nnz_remove = rowptr[r0+rr]-rowptr[r0]; // number of nonzeros to remove
    int nnz_new = this->nval+nnz_insert-nnz_remove;// nonzeros in updated matrix

    int *rp = new int[this->rows+1];
    int *ci = new int[nnz_new];
    MT *v   = new MT[nnz_new];

    for (i = 0; i <= r0; i++)
	rp[i] = rowptr[i];
    for (; i <= r0+rr; i++)
	rp[i] = rp[i-1] + rws.rowptr[i-r0]-rws.rowptr[i-r0-1];
    for (; i <= this->rows; i++)
	rp[i] = rowptr[i] + nnz_insert-nnz_remove;
    
    ofs = 0;
    ncpy = rowptr[r0];
    if (ncpy) {
	memcpy(ci+ofs, colidx, ncpy*sizeof(int));
	memcpy(v+ofs, this->val, ncpy*sizeof(MT));
	ofs += ncpy;
    }
    ncpy = rws.rowptr[rr];
    if (ncpy) {
	memcpy(ci+ofs, rws.colidx, ncpy*sizeof(int));
	memcpy(v+ofs, rws.val, ncpy*sizeof(MT));
	ofs += ncpy;
    }
    ncpy = rowptr[this->rows] - rowptr[r0+rr];
    if (ncpy) {
	memcpy(ci+ofs, colidx+rowptr[r0+rr], ncpy*sizeof(int));
	memcpy(v+ofs, this->val+rowptr[r0+rr], ncpy*sizeof(MT));
	ofs += ncpy;
    }
    Initialise(rp, ci, v);
    delete []rp;
    delete []ci;
    delete []v;

    xASSERT(ofs == this->nval, "SetRows: inconsistency");
}

template<class MT>
void TCompRowMatrix<MT>::RemoveRow (int nR)
{
 SetColAccess (false);

 int lm, j, k, g, s, f, r;

   lm = rowptr[nR+1]-rowptr[nR];
   //cerr<<"lm "<<lm<<endl;
   
   idxtype *lrowptr = new idxtype[this->rows];
   idxtype *lcolidx = new idxtype[this->nval-lm];
   MT  *lval    = new MT[this->nval-lm];

   
   for(j = 0; j < nR; j++)lrowptr[j] = rowptr[j];
   for(g = 0; g < rowptr[nR]; g++) lcolidx[g] = colidx[g];
   for(f = 0; f <  rowptr[nR]; f++) lval[f] = this->val[f];
   for(s = rowptr[nR]; s < this->nval-lm; s++) lcolidx[s] = colidx[s+lm];
   for(k = nR; k < this->rows; k++)  lrowptr[k] = rowptr[k+1] - lm;
   for(r = rowptr[nR]; r < this->nval-lm; r++) lval[r] = this->val[r+lm];

   delete []rowptr;
   delete []colidx;
   delete []this->val;
   rowptr = lrowptr;
   colidx = lcolidx;
   this->val    = lval;
   this->rows=this->rows-1;
   this->nval=this->nval-lm;

   //cerr<<"-row "<<nR<<" Removed"<<endl;
}

template<class MT>
void TCompRowMatrix<MT>::ColScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->cols, "Argument 1: wrong size");
    for (int i = 0; i < this->nval; i++)
        this->val[i] *= scale[colidx[i]];
}

template<class MT>
void TCompRowMatrix<MT>::RowScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->rows, "Argument 1: wrong size");
    for (int r = 0; r < this->rows; r++)
	for (int i = rowptr[r]; i < rowptr[r+1]; i++)
	    this->val[i] *= scale[r];
}

template<class MT>
MT TCompRowMatrix<MT>::GetNext (int &r, int &c) const
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

// ==========================================================================

#if THREAD_LEVEL==1

#ifdef OLD_AX_ENGINE
template<class MT>
void *Ax_engine (void *context)
{
    typedef struct {
        pthread_t ht;
	const TCompRowMatrix<MT> *A;
	const TVector<MT> *x;
	TVector<MT> *b;
        int ith;
    } THDATA;
    THDATA *thdata = (THDATA*)context;
    const TCompRowMatrix<MT> *A = thdata->A;
    const TVector<MT> *x = thdata->x;
    TVector<MT> *b = thdata->b;
    int r0 = (thdata->ith*b->Dim())/NUMTHREAD;
    int r1 = ((thdata->ith+1)*b->Dim())/NUMTHREAD;
    MT br;
    const MT *A_val = A->ValPtr();
    const MT *x_val = x->data_buffer();
    MT *b_val = b->data_buffer();
    int r, i = A->rowptr[r0], i2;
    int *colidx = A->colidx;

    for (r = r0; r < r1;) {
	i2 = A->rowptr[r+1];
	for (br = (MT)0; i < i2; i++)
	    br += A_val[i] * x_val[colidx[i]];
	b_val[r++] = br;
    }
    return NULL;
}
#else
template<class MT>
void Ax_engine (int ith, void *context)
{
    typedef struct {
	const TCompRowMatrix<MT> *A;
	const TVector<MT> *x;
	TVector<MT> *b;
    } THDATA;
    THDATA *thdata = (THDATA*)context;
    const TCompRowMatrix<MT> *A = thdata->A;
    const TVector<MT> *x = thdata->x;
    TVector<MT> *b = thdata->b;
    TVector<MT> b_local(b->Dim());

    MT br;
    const MT *A_val = A->ValPtr();
    const MT *x_val = x->data_buffer();
    MT *b_val = b->data_buffer();
    MT *b_local_val = b_local.data_buffer();
    int r, i = ith, i2;
    int stride = ThreadPool2::Pool()->NumThread();
    int *colidx = A->colidx;

    for (r = 0; r < A->nRows();) {
	i2 = A->rowptr[r+1];
	for (br = (MT)0; i < i2; i += stride)
	    br += A_val[i] * x_val[colidx[i]];
	b_local_val[r++] = br;
    }

    // add contribution to global b
    ThreadPool2::Pool()->MutexLock();
    for (r = 0; r < A->nRows(); r++)
        b_val[r] += b_local_val[r];
    ThreadPool2::Pool()->MutexUnlock();

    //ThreadPool2::Pool()->MutexLock();
    //cerr << "Ax done, thread " << ith << endl;
    //ThreadPool2::Pool()->MutexUnlock();
}
#endif

#ifdef NEED_EXPLICIT_INSTANTIATION
template void Ax_engine<double> (int ith, void *context);
template void Ax_engine<float> (int ith, void *context);
template void Ax_engine<complex> (int ith, void *context);
template void Ax_engine<scomplex> (int ith, void *context);
template void Ax_engine<int> (int ith, void *context);
#endif

template<class MT>
void TCompRowMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->cols,
	"Parameter 1 invalid size (expected %d, actual %d)",
        this->cols, x.Dim());
    if (b.Dim() != this->rows) b.New(this->rows);

    static struct {
	const TCompRowMatrix<MT> *A;
	const TVector<MT> *x;
	TVector<MT> *b;
    } thdata;

    thdata.A = this;
    thdata.x = &x;
    thdata.b = &b;
    
    b.Clear();
    ThreadPool2::Pool()->Invoke (Ax_engine<MT>, &thdata);
}

#else

template<class MT>
void TCompRowMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->cols,
	"Parameter 1 invalid size (expected %d, actual %d)",
        this->cols, x.Dim());
    if (b.Dim() != this->rows) b.New(this->rows);

    int r, i, i2;
    MT br;

    for (r = i = 0; r < this->rows;) {
	i2 = rowptr[r+1];
	for (br = (MT)0; i < i2; i++)
	    br += this->val[i] * x[colidx[i]];
	b[r++] = br;
    }
}
#endif // THREAD_LEVEL==1

#ifdef UNDEF// USE_SPBLAS
// warning: SPBLAS appears to have a size limit to the matrix and fails with
// an arithmetic exception if this is exceeded.
template<>
void TCompRowMatrix<double>::Ax (const TVector<double> &x, TVector<double> &b)
    const
{
    // NOT THREADSAFE: DO NOT USE THIS WITH THREADED CODE
    static double *dwork;
    static int ndwork = 0;

    dASSERT(x.Dim() == cols,
	"Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());
    if (b.Dim() != rows) b.New(rows);

    static int transa = 0;
    static int ccols = 1;
    static double alpha = 1.0, beta = 0.0;
    static int descra[9] = {0,0,0,0,1};
    if (ndwork < rows) {
        if (ndwork) delete []dwork;
	dwork = new double[ndwork=rows];
    }

    dcsrmmz_(transa, (int&)rows, ccols, (int&)cols, alpha, descra, val, colidx,
	     rowptr, rowptr+1, x.data_buffer(), (int&)cols, beta,
	     b.data_buffer(), (int&)rows, dwork, ndwork);
}
#endif // USE_SPBLAS

#ifdef USE_CUDA_FLOAT
// Specialisation: double precision
template<>
void TCompRowMatrix<double>::Ax (const TVector<double> &x, TVector<double> &b)
    const
{
    dASSERT(x.Dim() == cols,
	"Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());
    if (b.Dim() != rows) b.New(rows);

    Ax_spmv (rows, cols, val, rowptr, colidx, x.data_buffer(),
	     b.data_buffer());
}

// Specialisation: single precision
template<>
void TCompRowMatrix<float>::Ax (const TVector<float> &x, TVector<float> &b)
    const
{
    dASSERT(x.Dim() == cols,
	"Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());
    if (b.Dim() != rows) b.New(rows);

    //cuda_Ax (val, rowptr, colidx, rows, cols, x.data_buffer(),
    //    b.data_buffer());
    Ax_spmv (rows, cols, val, rowptr, colidx, x.data_buffer(),
	     b.data_buffer());
}

// Specialisation: single precision complex
template<>
void TCompRowMatrix<scomplex>::Ax (const TVector<scomplex> &x,
    TVector<scomplex> &b) const
{
    dASSERT(x.Dim() == cols,
	"Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());
    if (b.Dim() != rows) b.New(rows);

    cuda_Ax_cplx (val, rowptr, colidx, rows, cols, x.data_buffer(),
        b.data_buffer());
}

#endif // USE_CUDA_FLOAT

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b,
    int r1, int r2) const
{
    dASSERT(x.Dim() == this->cols,
        "Parameter 1 invalid size (expected %d, actual %d)",
	this->cols, x.Dim());

    int r, i2;
    register int i = rowptr[r1];
    register idxtype *pcolidx = colidx+i;
    register MT br, *pval = this->val+i;

    if (b.Dim() != this->rows) b.New (this->rows);

    for (r = r1; r < r2;) {
	i2 = rowptr[r+1];
	for (br = (MT)0; i < i2; i++)
	    br += *pval++ * x[*pcolidx++];
	b[r++] = br;
    }
}

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::Ax_cplx (const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    // Specialisation MT=double only (see below)
    xERROR("Method Ax_cplx not defined for return type");
}

template<>
inline void TCompRowMatrix<double>::Ax_cplx (
    const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    dASSERT(x.Dim() == cols, "Invalid size - vector x");

    if (b.Dim() != rows) b.New (rows);

    int r, i2;
    register int i;
	register idxtype *pcolidx = colidx;
    register double *pval = val;
    register double b_re, b_im;

    for (r = i = 0; r < rows;) {
	i2 = rowptr[r+1];
	b_re = b_im = 0.0;
	for (; i < i2; i++) {
	    b_re += x[*pcolidx].real() * *pval;
	    b_im += x[*pcolidx++].imag() * *pval++;
	}
	b[r++] = std::complex<double>(b_re, b_im);
    }
    INC_FLOPS_ADD (rowptr[rows]*2);
    INC_FLOPS_MUL (rowptr[rows]*2);
}

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(x.Dim() == this->rows,
        "Parameter 1 invalid size (expected %d, actual %d)",
		 this->rows, x.Dim());

    if (b.Dim() != this->cols) b.New (this->cols);
    else b.Clear();

    if (!col_access) SetColAccess();

    int i, c;
    for (i = c = 0; i < this->nval; i++) {
        while (colptr[c+1] <= i) c++;
	b[c] += this->val[vofs[i]] * x[rowidx[i]];
    }
}

template<> // specialisation:: complex
inline void TCompRowMatrix<std::complex<double> >::ATx (
    const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    dASSERT(x.Dim() == rows, "Invalid size - vector x");
    dASSERT(b.Dim() == cols, "Invalid size - vector b");

    if (!col_access) SetColAccess(); // should not be necessary!
    int i, c;
    for (c = 0; c < cols; c++) b[c] = 0;
    for (i = c = 0; i < nval; i++) {
        while (colptr[c+1] <= i) c++;
	b[c] += conj(val[vofs[i]]) * x[rowidx[i]];
    }
}

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::ATx_cplx (const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    // Specialisation MT=double only (see below)
    xERROR("Method ATx_cplx not defined for return type");
}

template<>
inline void TCompRowMatrix<double>::ATx_cplx (
    const TVector<std::complex<double> > &x,
    TVector<std::complex<double> > &b) const
{
    dASSERT(x.Dim() == rows, "Invalid size - vector x");
    dASSERT(b.Dim() == cols, "Invalid size - vector b");

    if (!col_access) SetColAccess(); // should not be necessary!
    int i, c;

    for (c = 0; c < cols; c++)
        b[c] = std::complex<double>(0,0);
    for (i = c = 0; i < nval; i++) {
	while (colptr[c+1] <= i) c++;
	b[c] += val[vofs[i]] * x[rowidx[i]];
    }
}

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::Ax (const TVector<MT> *x, TVector<MT> *b, int nrhs)
    const
{
    for (int i = 0; i < nrhs; i++)
	Ax (x[i], b[i]);
}

#ifdef USE_CUDA_FLOAT

template<>
void TCompRowMatrix<float>::Ax (const TVector<float> *x, TVector<float> *b,
    int nrhs) const
{
    static int bufsize = 16;
    static const float **xbuf = new const float*[bufsize];
    static float **bbuf = new float*[bufsize];
    if (nrhs > bufsize) {
	bufsize = nrhs;
	delete []xbuf;    xbuf = new const float*[bufsize];
	delete []bbuf;    bbuf = new float*[bufsize];
    }
    for (int i = 0; i < nrhs; i++) {
	xbuf[i] = x[i].data_buffer();
	bbuf[i] = b[i].data_buffer();
    }

    Ax_spmv (rows, cols, val, rowptr, colidx, xbuf, bbuf, nrhs);
}

template<>
void TCompRowMatrix<double>::Ax (const TVector<double> *x, TVector<double> *b,
    int nrhs) const
{
    static int bufsize = 16;
    static const double **xbuf = new const double*[bufsize];
    static double **bbuf = new double*[bufsize];
    if (nrhs > bufsize) {
	bufsize = nrhs;
	delete []xbuf;    xbuf = new const double*[bufsize];
	delete []bbuf;    bbuf = new double*[bufsize];
    }
    for (int i = 0; i < nrhs; i++) {
	xbuf[i] = x[i].data_buffer();
	bbuf[i] = b[i].data_buffer();
    }

    Ax_spmv (rows, cols, val, rowptr, colidx, xbuf, bbuf, nrhs);
}

#endif // USE_CUDA_FLOAT

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::AB (const TCompRowMatrix<MT> &B,
    TCompRowMatrix<MT> &C) const
{
    // Sparse matrix product C = AB

    int i, j, k, m, ra, ra1, ra2, rb, rb1, rb2, nzero;
    int nr = this->rows;
    int nc = B.cols;

    bool *fillrow = new bool[nc];
    idxtype *Crowptr  = new idxtype[nr+1];
    idxtype *Crowptr1 = Crowptr+1; // for faster access
    for (i = 0; i <= nr; i++) Crowptr[i] = 0;

    // pass 1: determine sparsity pattern of product matrix C
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) fillrow[j] = false;
        ra1 = rowptr[i];
	ra2 = rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = colidx[ra];
	    rb1 = B.rowptr[k];
	    rb2 = B.rowptr[k+1];
	    for (rb = rb1; rb < rb2; rb++) {
	        j = B.colidx[rb];
		if (!fillrow[j]) {
		    fillrow[j] = true;
		    Crowptr1[i]++;
		}
	    }
	}
    }
    for (i = 1; i <= nr; i++) Crowptr[i] += Crowptr[i-1];

    nzero = Crowptr[nr];
    idxtype *Ccolidx = new idxtype[nzero];
    MT *Cval     = new MT[nzero];
    MT *rowval   = new MT[nc];

    // row buffer
    int *vecidx  = new int[nc];
    int *vcolidx = new int[nc];
    MT *vval     = new MT[nc];
    int nnvec;
    for (j = 0; j < nc; j++) vecidx[j] = -1;

    // pass 2: create product matrix C
    for (i = 0; i < nr; i++) {
	ra1 = rowptr[i];
	ra2 = rowptr[i+1];
	nnvec = 0;
	for (ra = ra1; ra < ra2; ra++) {
	    k = colidx[ra];
	    rb1 = B.rowptr[k];
	    rb2 = B.rowptr[k+1];
	    for (rb = rb1; rb < rb2; rb++) {
	        j = B.colidx[rb];
		if (vecidx[j] >= 0) {
		    m = vecidx[j];
		    vval[m] += this->val[ra] * B.val[rb];
		} else {
		    vval[nnvec] = this->val[ra] * B.val[rb];
		    vcolidx[nnvec] = j;
		    vecidx[j] = nnvec++;
		}
	    }
	}
	for (j = 0, ra = Crowptr[i]; j < nnvec; j++) {
	    Cval[ra] = vval[j];
	    Ccolidx[ra] = vcolidx[j];
	    vecidx[vcolidx[j]] = -1;
	    ra++;
	}
    }

    C.New (nr, nc);
    C.Initialise (Crowptr, Ccolidx, Cval);

    // cleanup
    delete []vecidx;
    delete []vcolidx;
    delete []vval;

    delete []Crowptr;
    delete []Ccolidx;
    delete []Cval;
    delete []rowval;
    delete []fillrow;
}

#ifdef UNDEF // old version - slow but tested!
template<class MT>
void TCompRowMatrix<MT>::AB (const TCompRowMatrix<MT> &B,
    TCompRowMatrix<MT> &C) const
{
    // Sparse matrix product C = AB

    int i, j, k, ra, ra1, ra2, rb, rb1, rb2, nzero;
    int nr = this->rows;
    int nc = B.cols;

    bool *fillrow = new bool[nc];
    int *Crowptr  = new int[nr+1];
    int *Crowptr1 = Crowptr+1; // for faster access
    for (i = 0; i <= nr; i++) Crowptr[i] = 0;

    // pass 1: determine sparsity pattern of product matrix C
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) fillrow[j] = false;
        ra1 = rowptr[i];
	ra2 = rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = colidx[ra];
	    rb1 = B.rowptr[k];
	    rb2 = B.rowptr[k+1];
	    for (rb = rb1; rb < rb2; rb++) {
	        j = B.colidx[rb];
		if (!fillrow[j]) {
		    fillrow[j] = true;
		    Crowptr1[i]++;
		}
	    }
	}
    }
    for (i = 1; i <= nr; i++) Crowptr[i] += Crowptr[i-1];

    nzero = Crowptr[nr];
    int *Ccolidx = new int[nzero];
    MT *Cval     = new MT[nzero];
    MT *rowval   = new MT[nc];

    // pass 2: create product matrix C
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
	    fillrow[j] = false;
	    rowval[j] = (MT)0;
	}
	ra1 = rowptr[i];
	ra2 = rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = colidx[ra];
	    rb1 = B.rowptr[k];
	    rb2 = B.rowptr[k+1];
	    for (rb = rb1; rb < rb2; rb++) {
	        j = B.colidx[rb];
		fillrow[j] = true;
		rowval[j] += this->val[ra] * B.val[rb];
	    }
	}
	for (j = k = 0; j < nc; j++) {
	    if (fillrow[j]) {
	        ra = Crowptr[i]+k;
	        Cval[ra] = rowval[j];
		Ccolidx[ra] = j;
		k++;
	    }
	}
    }

    C.New (nr, nc);
    C.Initialise (Crowptr, Ccolidx, Cval);

    // cleanup
    delete []Crowptr;
    delete []Ccolidx;
    delete []Cval;
    delete []rowval;
    delete []fillrow;
}
#endif

template<class MT>
void TCompRowMatrix<MT>::Transpone ()
{
    SetColAccess (false);

    int i, j, k, idx;
    int *rcount = new int[this->cols];
    for (i = 0; i < this->cols; i++) rcount[i] = 0;
    for (i = 0; i < this->nval; i++) rcount[colidx[i]]++;
    idxtype *trowptr = new idxtype[this->cols+1];
    idxtype *tcolidx = new idxtype[this->nval];
    MT  *tval    = new MT[this->nval];
    trowptr[0] = 0;
    for (i = 0; i < this->cols; i++)
        trowptr[i+1] = trowptr[i]+rcount[i];
    for (i = 0; i < this->cols; i++) rcount[i] = 0;
    for (i = 0; i < this->rows; i++)
        for (k = rowptr[i]; k < rowptr[i+1]; k++) {
	    j = colidx[k];
	    idx = trowptr[j]+rcount[j];
	    tcolidx[idx] = i;
	    tval[idx] = this->val[k];
	    rcount[j]++;
	}
    delete []rcount;
    delete []rowptr;
    delete []colidx;
    delete []this->val;
    rowptr = trowptr;
    colidx = tcolidx;
    this->val = tval;
    i = this->rows, this->rows = this->cols, this->cols = i;
}

template<class MT>
TCompRowMatrix<MT> transp (const TCompRowMatrix<MT> &m)
{
    TCompRowMatrix<MT> mt(m);
    mt.Transpone();
    return mt;
}

template<class MT>
double l2norm (const TCompRowMatrix<MT> &A)
{
    double term, sum = 0.0;
    for (int i = 0; i < A.nval; i++) {
	term = norm(A.val[i]);
	sum += term*term;
    }
    return sqrt (sum);
}

template<class MT>
TCompRowMatrix<MT> kron (const TCompRowMatrix<MT> &A,
    const TCompRowMatrix<MT> &B)
{
    int ia, ib, ka, kb, ja, jb, i, idx;
    int va_i, vb_i, v_i;
    int na = A.nRows(), ma = A.nCols(), va = A.nVal();
    int nb = B.nRows(), mb = B.nCols(), vb = B.nVal();
    int n = na*nb, m = ma*mb, v = va*vb;
    MT a_ij;
    MT b_ij;
    const MT *Aval = A.ValPtr();
    const MT *Bval = B.ValPtr();

    idxtype *rowptr = new idxtype[n+1];
    idxtype *colidx = new idxtype[v];
    MT *val = new MT[v];

    rowptr[0] = 0;
    i = idx = 0;
    for (ia = 0; ia < na; ia++) {
	va_i = A.rowptr[ia+1] - A.rowptr[ia]; // nonzeros in row ia of A
	for (ib = 0; ib < nb; ib++) {
	    vb_i = B.rowptr[ib+1] - B.rowptr[ib]; // nonzeros in row ib of B
	    v_i = va_i * vb_i;
	    rowptr[i+1] = rowptr[i] + v_i;
	    for (ka = A.rowptr[ia]; ka < A.rowptr[ia+1]; ka++) {
		ja = A.colidx[ka];
		a_ij = Aval[ka];
		for (kb = B.rowptr[ib]; kb < B.rowptr[ib+1]; kb++) {
		    jb = B.colidx[kb];
		    b_ij = Bval[kb];
		    colidx[idx] = ja*mb + jb;
		    val[idx] = a_ij * b_ij;
		    idx++;
		}
	    }
	    i++;
	}
    }
    TCompRowMatrix<MT> C (n, m, rowptr, colidx, val);
    delete []rowptr;
    delete []colidx;
    delete []val;
    return C;
}

template<class MT>
MT TCompRowMatrix<MT>::RowMult (int r, MT *x) const
{
    dASSERT(r >= 0 && r < this->rows, "Argument 1 out of range");

    MT sum = (MT)0;
    int imax = rowptr[r+1];
    for (int i = rowptr[r]; i < imax; i++)
        sum += this->val[i] * x[colidx[i]];
    return sum;
}

template<class MT>
void TCompRowMatrix<MT>::Sort () const
{
    // This is just a quick and dirty insertion sort using Shell's method
    // (NR p331). Should be replaced if entries/row > 50

    if (sorted) return;
    int r, c, i1, i2, n, i, j, inc;
    MT v;

    for (r = 0; r < this->rows; r++) {
        i1  = rowptr[r];
	i2  = rowptr[r+1];
	n   = i2-i1;
	inc = 1;
	do { inc *= 3; inc++; } while (inc <= n);
	do {
	    inc /= 3;
	    for (i = inc+i1; i < i2; i++) {
		c = colidx[i];
	        v = this->val[i];
		j = i;
		while (colidx[j-inc] > c) {
		    colidx[j] = colidx[j-inc];
		    this->val[j] = this->val[j-inc];
		    j -= inc;
		    if (j-i1 < inc) break;
		}
		colidx[j] = c;
		this->val[j] = v;
	    }
	} while (inc > 1);
    }
    sorted = true;
    col_access = false;
}

template<class MT>
int TCompRowMatrix<MT>::Shrink ()
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
    i = this->nval - nz;
    this->nval = nz;
    return i;
}

template<class MT>
void TCompRowMatrix<MT>::SetColAccess (bool yes) const
{
    if (yes) {
        if (col_access) return; // already initialised - nothing to do
	SetColAccess (false);   // make sure all lists are de-allocated
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

template<class MT>
void TCompRowMatrix<MT>::SetDiagAccess (bool yes) const
{
    if (yes) {
        if (diag_access) return;
	if (!sorted) Sort();
	SetDiagAccess (false);
	diagptr = new int[this->rows];
	for (int i = 0; i < this->rows; i++)
	    diagptr[i] = Get_index_sorted (i,i);
	diag_access = true;
    } else {
        if (diagptr) delete []diagptr;
	diagptr = 0;
	diag_access = false;
    }
}

template<class MT>
double TCompRowMatrix<MT>::LargestInRow (int r, int i) const
{
    double magnitude, largest = 0.0;
    for (int j = rowptr[r]+i; j < rowptr[r+1]; j++)
        if ((magnitude = norm(this->val[j])) > largest)
	    largest = magnitude;
    return largest;
}

template<class MT>
double TCompRowMatrix<MT>::LargestInCol (int c, int i) const
{
    double magnitude, largest = 0.0;
    if (!col_access) SetColAccess();
    for (int j = colptr[c]+i; j < colptr[c+1]; j++)
        if ((magnitude = norm(this->val[vofs[j]])) > largest)
	    largest = magnitude;
    return largest;
}

template<class MT>
void TCompRowMatrix<MT>::ReplaceRow (int row, int nz, int *rcolidx, MT *rval)
{
    dASSERT(row >= 0 && row < this->rows, "Row index out of range");
    dASSERT(nz >= 0 && nz <= this->cols, "Nonzero count out of range");

    int i, rp0, rp1, nz_old, dnz;
    rp0 = rowptr[row];
    rp1 = rowptr[row+1];
    nz_old = rp1 - rp0;
    dnz = nz - nz_old;
    if (dnz) { // reallocation required
        int nval_new = this->nval + dnz;
	// new column index list
	idxtype *colidx_new = new idxtype[nval_new];
	memcpy (colidx_new, colidx, rp0 * sizeof(idxtype));
	memcpy (colidx_new+rp0+nz, colidx+rp1, (this->nval-rp1) * sizeof(idxtype));
	if (this->nval) delete []colidx;
	colidx = colidx_new;
	// new data block
	MT *val_new = new MT[nval_new];
	memcpy (val_new, this->val, rp0 * sizeof(MT));
	memcpy (val_new+rp0+nz, this->val+rp1, (this->nval-rp1) * sizeof(MT));
	if (this->nval) delete []this->val;
	this->val = val_new;
	// adjust row pointers
	for (i = row+1; i <= this->rows; i++) rowptr[i] += dnz;
	this->nval += dnz;
    }
    // insert new row
    memcpy (colidx+rp0, rcolidx, nz * sizeof(int));
    if (rval) memcpy (this->val+rp0, rval, nz * sizeof(MT));
    else      for (i = 0; i < nz; i++) this->val[rp0+i] = 0;
    col_access = false; // column index lists no longer valid
}

template<class MT>
void TCompRowMatrix<MT>::ReplaceRow (int row, const TVector<MT>& vec)
{
    dASSERT(vec.Dim() == this->rows, "Invalid row vector size");

    int i, nz;
    int *rcolidx = new int[this->rows];
    MT  *rval    = new MT[this->rows];
    for (i = nz = 0; i < this->rows; i++) {
        if (vec[i] != (MT)0) {
	    rcolidx[nz] = i;
	    rval[nz++]  = vec[i];
	}
    }
    ReplaceRow (row, nz, rcolidx, rval);
    delete []rcolidx;
    delete []rval;
}

inline void BlockExpand (int *rowptr, int *colidx, int n,
		  int *&browptr, int *&bcolidx, int &bn,
		  int blockn, int blockm)
{
    // expand an (n x m) matrix into an (n blockn x m blockm) matrix by
    // expanding every entry into a (blockn x blockm) block
    int nzero = rowptr[n];
    int bnzero = nzero * blockn * blockm;
    int i, ii, j, jj, idx, ncol, bncol;

    bn = n*blockn;
    browptr = new int[bn+1];
    bcolidx = new int[bnzero];
    browptr[0] = 0;
    for (i = idx = 0; i < n; i++) {
        ncol  = rowptr[i+1] - rowptr[i];
	bncol = ncol*blockm;
	for (ii = 0; ii < blockn; ii++) {
	    browptr[i*blockn+ii+1] = browptr[i*blockn+ii] + bncol;
	    for (j = 0; j < ncol; j++) {
	        for (jj = 0; jj < blockm; jj++)
		    bcolidx[idx++] = colidx[rowptr[i]+j]*blockm+jj;
	    }
	}
    }
}

template<class MT>
MT TCompRowMatrix<MT>::row_mult (int r1, int r2, int from, int to) const
{
    dASSERT(sorted, "Requires sorted CompRowMatrix");

    MT sum = (MT)0;
    int i1, i2, c1, c2;
    int end1 = rowptr[r1+1], end2 = rowptr[r2+1];

    if ((i1 = rowptr[r1]) == end1) return sum; // no entries in r1
    if ((i2 = rowptr[r2]) == end2) return sum; // no entries in r2
    if ((c1 = colidx[i1]) > to) return sum; // no entries within range in r1
    if ((c2 = colidx[i2]) > to) return sum; // no entries within range in r2

    for (;;) {
        if (c1 < c2) {
	    do {
	        if (++i1 == end1) return sum;
	    } while ((c1 = colidx[i1]) < c2);
	    if (c1 > to) return sum;
	} else if (c1 > c2) {
	    do {
	        if (++i2 == end2) return sum;
	    } while ((c2 = colidx[i2]) < c1);
	    if (c2 > to) return sum;
	}
	if (c1 == c2 && c1 >= from) {
	    sum += this->val[i1++] * this->val[i2++];
	    if (i1 == end1 || i2 == end2) return sum;
	    c1 = colidx[i1];
	    c2 = colidx[i2];
	}
    }
}

template<class MT>
MT TCompRowMatrix<MT>::sparsevec_mult (int *idx1, MT *val1, int n1,
				       int *idx2, MT *val2, int n2) const
{
    MT sum = (MT)0;
    if (!n1 || !n2) return sum;

    int i1, i2, c1, c2;

    i1 = 0; c1 = idx1[i1];
    i2 = 0; c2 = idx2[i2];
    
    for (;;) {
        if (c1 < c2) {
	    do {
	        if (++i1 == n1) return sum;
	    } while ((c1 = idx1[i1]) < c2);
	} else if (c1 > c2) {
	    do {
	        if (++i2 == n2) return sum;
	    } while ((c2 = idx2[i2]) < c1);
	}
	if (c1 == c2) {
	    sum += val1[i1++] * val2[i2++];
	    if (i1 == n1 || i2 == n2) return sum;
	    c1 = idx1[i1];
	    c2 = idx2[i2];
	}
    }
    return sum;
}

#ifdef UNDEF
bool overlap (int i, int j, int *col, int *next)
{
    // requires col to be sorted

    int lastcol = i-1;
    int c1, c2;

    if ((c1 = col[i]) < 0) return false; // no entries in row i
    if ((c2 = col[j]) < 0) return false; // no entries in row j

    if (c1 > lastcol) return false;
    if (c2 > lastcol) return false;

    for (;;) {
        if (c1 < c2) {
	    do {
	        if ((i = next[i]) < 0) return false; // end of row i
	    } while ((c1 = col[i]) < c2);
	    if (c1 > lastcol) return false;
	} else if (c2 < c1) {
	    do {
	        if ((j = next[j]) < 0) return false; // end of row j
	    } while ((c2 = col[j]) < c1);
	    if (c2 > lastcol) return false;
	}
	if (c1 == c2) return true;
    }
}
#endif

template<class MT>
void TCompRowMatrix<MT>::SymbolicCholeskyFactorize (idxtype *&frowptr,
    idxtype *&fcolidx) const
{
    dASSERT(this->rows == this->cols, "Requires square matrix");
    symbolic_cholesky_factor (this->rows, rowptr, colidx, frowptr, fcolidx);
    // implemented in cr_cholesky.cc
}

template<class MT>
void TCompRowMatrix<MT>::CalculateIncompleteCholeskyFillin (idxtype *&frowptr,
    idxtype *&fcolidx) const
{
    int i, j, nc;
    if (!sorted) Sort();

    // pass 1: construct row pointer list
    frowptr = new idxtype[this->rows+1];
    frowptr[0] = 0;
    for (i = 0; i < this->rows; i++) {
        for (nc = 0, j = rowptr[i]; j < rowptr[i+1]; j++) {
	    if (colidx[j] < i) nc++;
	    else break;
	}
	frowptr[i+1] = frowptr[i]+nc;
    }
    // pass 2: construct column index list
    fcolidx = new idxtype[frowptr[this->rows]];
    for (i = nc = 0; i < this->rows; i++) {
        for (j = rowptr[i]; j < rowptr[i+1]; j++) {
	    if (colidx[j] < i) fcolidx[nc++] = colidx[j];
	    else break;
	}
    }
}

template<class MT>
bool CholeskyFactorize (const TCompRowMatrix<MT> &A, TCompRowMatrix<MT> &L,
    TVector<MT> &d, bool recover)
{
    // This assumes that entries of L are sorted COLUMN-wise, which should
    // be the case if L has been initialised via SymbolicCholeskyFactorize
    // Do NOT call L.Sort since this will sort row-wise!

    const double EPS = 1e-10;
    int i, k, c, r, cc, term, n = A.nRows();
    MT idiag, ajk, *fullcol = new MT[n];
    bool ok = true;
    if (!A.col_access) A.SetColAccess();
    if (!L.col_access) L.SetColAccess();
    // some shortcuts to improve performance in the inner loops
    idxtype *Lrowidx = L.rowidx;
    idxtype *Lvofs   = L.vofs;
    MT     *Lval    = L.val;
    for (c = 0; c < n; c++) {
	// expand current column
        memset (fullcol, 0, n*sizeof(double));
	for (i = A.colptr[c]; i < A.colptr[c+1]; i++)
	    if ((r = A.rowidx[i]) >= c)
	        fullcol[r] = A.val[A.vofs[i]];

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
}

template<class MT>
bool IncompleteCholeskyFactorize (const TCompRowMatrix<MT> &A,
    TCompRowMatrix<MT> &L, TVector<MT> &d, bool recover)
{
    const double EPS = 1e-10;
    int i, j, k, n = A.nRows();
    bool ok = true;
    MT x, v, diag, Aj;

    L.Zero();
    if (!L.sorted) L.Sort();
    for (i = 0; i < n; i++) {
        diag = A.Get(i,i);
	x = (MT)0;
	for (j = L.rowptr[i]; j < L.rowptr[i+1]; j++)
	    x += L.val[j] * L.val[j];
	if (diag > x) {       /* problem here, if using complex */
	    d[i] = sqrt (diag-x);
	} else {
	    if (!recover) xERROR("Matrix not positive definite");
	    ok = false;
	    d[i] = EPS; // force positive
	}
	for (k = A.rowptr[i]; k < A.rowptr[i+1]; k++) {
	    j = A.colidx[k];
	    if (j <= i) continue;
	    Aj = A.val[k];
	    if (i) {
	        x = L.row_mult (i, j, 0, i-1);
		v = (Aj-x)/d[i];
	    } else {
	        v = Aj/d[i];
	    }
	    L(j,i) = v;
	}
    }
    return ok;
}

template<class MT>
void LUFactorize (TCompRowMatrix<MT> &A, bool LUrealloc)
{
    // implements Doolittle LU factorisation A=LU
    // with L unit lower triangular and U upper triangular
    // l_ij = (a_ij - Sum_{p=1}^{j-1} l_ip u_pj)/u_jj    (for i > j)
    // u_ij = (a_ij - Sum_{p=1}^{i-1} l_ip u_pj)         (for i <= j)
    
    int n = A.rows;
    int i, j, jj, cj, k, nn, rpk0, rpk1;
    TVector<MT> fullrow(n);
    int *fullrow_colidx = new int[n];
    MT  *fullrow_val    = new MT[n];
    
    for (k = 1; k < n; k++) {
        rpk0 = A.rowptr[k];
	rpk1 = A.rowptr[k+1];
	// scatter
	for (i = rpk0; i < rpk1; i++)
	    fullrow[A.colidx[i]] = A.val[i];

	for (i = rpk0; i < rpk1; i++) {
	    j = A.colidx[i];
	    if (j >= k) continue;
	    MT alpha = fullrow[j] / A.Get(j,j);
	    //A.val[i] = alpha;
	    fullrow[j] = alpha;
	    for (jj = A.rowptr[j]; jj < A.rowptr[j+1]; jj++) {
	        cj = A.colidx[jj];
	        if (cj <= j) continue;
		fullrow[cj] -= alpha * A.val[jj];
	    }
	}
	// gather
	if (LUrealloc) { // calculate fillin
	    for (nn = 0, i = rpk0; i < rpk1; i++) {
	        j = A.colidx[i];
		if (j <= k) {
		    fullrow_colidx[nn] = j;
		    fullrow_val[nn++] = fullrow[j];
		}
	    }
	    for (i = k+1; i < n; i++) {
	        if (fullrow[i] != (MT)0) {
		    fullrow_colidx[nn] = i;
		    fullrow_val[nn++] = fullrow[i];
		}
	    }
	    A.ReplaceRow (k, nn, fullrow_colidx, fullrow_val);
	} else { // assume LU fillin is already allocated
	    for (i = rpk0; i < rpk1; i++) {
	        j = A.colidx[i];
		/*if (j >= k)*/ A.val[i] = fullrow[j];
	    }
	}
    }

    // cleanup
    delete []fullrow_colidx;
    delete []fullrow_val;
}

template<class MT>
void TCompRowMatrix<MT>::CholeskySubst (const TVector<MT> &d,
    const MT *b, MT *x) const
{
    if (!col_access) SetColAccess();
    MT *Lval, sum;
    int i, k, k0, k1, n = this->rows;

    k = rowptr[0];
    Lval = this->val;
    idxtype *Lcolidx = colidx;
    for (i = 0; i < n; i++) {
	k1 = rowptr[i+1];
        for (sum = b[i]; k < k1; k++) {
	    sum -= *Lval++ * x[*Lcolidx++];
	}
	x[i] = sum / d[i];
    }
    k = colptr[n]-1;
    Lval = this->val;
    idxtype *Lrowidx = rowidx+k;
    idxtype *Lvofs   = vofs+k;
    for (i = n-1; i >= 0; i--) {
        k0 = colptr[i];
        for (sum = x[i]; k >= k0; k--)
	    sum -= Lval[*Lvofs--] * x[*Lrowidx--];
	x[i] = sum / d[i];
    }
}

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x)
{
    dASSERT(d.Dim() == L.nCols(), "Incompatible dimensions");
    dASSERT(d.Dim() == b.Dim(), "Incompatible dimensions");

    if (x.Dim() != d.Dim()) x.New (d.Dim());
    L.CholeskySubst (d, b.data_buffer(), x.data_buffer());
}

template<class MT>
TVector<MT> CholeskySolve (const TCompRowMatrix<MT> &L,
    const TVector<MT> &d, const TVector<MT> &b)
{
    TVector<MT> x(L.nRows()); 
    CholeskySolve (L, d, b, x);
    return x;
}

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TDenseMatrix<MT> &BT, TDenseMatrix<MT> &XT, int n)
{
    int i, m = L.nCols();
    dASSERT(m == BT.nCols(), "Incompatible dimensions");
    dASSERT(m == XT.nCols(), "Incompatible dimensions");

    if (n) {
        if (n > BT.nRows()) n = BT.nRows();
    } else {
        n = BT.nRows();
    }
    if (n > XT.nRows()) n = XT.nRows();
    const MT *pb = BT.data_buffer();
    MT *px = XT.data_buffer();
    for (i = 0; i < n; i++) {
        L.CholeskySubst (d, pb, px);
	pb += m;
	px += m;
    }
}

template<class MT>
void LUSolve (const TCompRowMatrix<MT> &LU, const TVector<MT> &b,
    TVector<MT> &x)
{
    int i, k, c, n = LU.nRows();
    MT diag, sum;

    for (i = 0; i < n; i++) {
        for (sum = b[i], k = LU.rowptr[i]; k < LU.rowptr[i+1]; k++)
	    if ((c = LU.colidx[k]) < i) sum -= LU.val[k] * x[c];
	x[i] = sum;
    }
    for (i = n-1; i >= 0; i--) {
        for (sum = x[i], k = LU.rowptr[i]; k < LU.rowptr[i+1]; k++) {
	    c = LU.colidx[k];
	    if (c > i) sum -= LU.val[k] * x[c];
	    else if (c == i) diag = LU.val[k];
	}
	x[i] = sum / diag;
    }
}

// ==========================================================================

template<class MT>
int ILUSolve (TCompRowMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double tol, double droptol, int maxit)
{
    xERROR ("ILUSolve: called for unsupported data type");
    // only implemented in data specialisations
    return 0;
}

template<>
inline int ILUSolve (TCompRowMatrix<std::complex<double> > &A,
    const TVector<std::complex<double> > &b,
    TVector<std::complex<double> > &x, double tol, double droptol, int maxit)
{
#ifdef HAVE_ILU
    return ILUSolveZGNL (A, b, x, tol, droptol, maxit);
#else
    xERROR("ILUPACK support not configured");
	return 0;
#endif
}

// ==========================================================================

template<class MT>
int ILUSymSolve (TCompRowMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double tol, double droptol, int maxit)
{
    xERROR ("ILUSymSolve: called for unsupported data type");
    // only implemented in data specialisations
    return 0;
}

template<>
inline int ILUSymSolve (TCompRowMatrix<std::complex<double> > &A,
    const TVector<std::complex<double> > &b,
    TVector<std::complex<double> > &x, double tol, double droptol, int maxit)
{
#ifdef HAVE_ILU
    return ILUSolveZSYM (A, b, x, tol, droptol, maxit);
#else
    xERROR("ILUPACK support not configured");
	return 0;
#endif
}

// ==========================================================================
// If we are using CUDA, the CompRowMatrix PCG and BiCGSTAB solvers are
// overloaded to make use of the CUSP implementations

#ifdef USE_CUDA_FLOAT


template<class MT>
int PCG (const TCompRowMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    return PCG((const TMatrix<MT>&)A, b, x, tol, precon, maxit);
}

template<> // specialisation: single precision
inline int PCG<float> (const FCompRowMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit)
{
    const float *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    float *x_val = x.data_buffer();
    const float *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    float ftol = (float)tol;
    cuda_CG (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, ftol, maxit,
		   &res);
    tol = (double)res.rel_error;
    return res.it_count;
}

template<> // specialisation: double precision
inline int PCG<double> (const RCompRowMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit)
{
    const double *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    double *x_val = x.data_buffer();
    const double *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    cuda_CG (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, tol, maxit,
		   &res);
    tol = res.rel_error;
    return res.it_count;
}

// ==========================================================================
// These versions use multiple right-hand sides

template<class MT>
void PCG (const TCompRowMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x, int nrhs, double tol, int maxit,
    TPreconditioner<MT> *precon, IterativeSolverResult *res)
{
    PCG ((const TMatrix<MT>&)A, b, x, nrhs, tol, maxit, precon, res);
}

template<> // specialisation: single precision
inline void PCG<float> (const FCompRowMatrix &A, const FVector *b, FVector *x,
    int nrhs, double tol, int maxit, FPreconditioner *precon,
    IterativeSolverResult *res)
{
    static int nbuf = 32;
    static float **x_val = new float*[nbuf];
    static const float **b_val = new const float*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new float*[nbuf];
	b_val = new const float*[nbuf];
    }
    
    const float *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    float ftol = (float)tol;
    cuda_CG (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs, ftol, maxit,
		   &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

template<> // specialisation: double precision
inline void PCG<double> (const RCompRowMatrix &A, const RVector *b,
    RVector *x, int nrhs, double tol, int maxit, RPreconditioner *precon,
    IterativeSolverResult *res)
{
    static int nbuf = 32;
    static double **x_val = new double*[nbuf];
    static const double **b_val = new const double*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new double*[nbuf];
	b_val = new const double*[nbuf];
    }
    
    const double *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    cuda_CG (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs, tol, maxit,
		   &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

// ==========================================================================

template<class MT>
int BiCGSTAB (const TCompRowMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    return BiCGSTAB((const TMatrix<MT>&)A, b, x, tol, precon, maxit);
    // by default, use the TMatrix implementation
}

template<> // specialisation: single precision
inline int BiCGSTAB<float> (const FCompRowMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit)
{
    const float *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    float *x_val = x.data_buffer();
    const float *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    float ftol = (float)tol;
    cuda_BiCGSTAB (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, ftol, maxit,
		   &res);
    tol = (double)res.rel_error;
    return res.it_count;
}

template<> // specialisation: single complex
inline int BiCGSTAB<scomplex> (const SCCompRowMatrix &A, const SCVector &b,
    SCVector &x, double &tol, SCPreconditioner *precon, int maxit)
{
    const scomplex *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    scomplex *x_val = x.data_buffer();
    const scomplex *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    float ftol = (float)tol;
    cuda_BiCGSTAB_cplx (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, ftol,
			maxit, &res);
    tol = (double)res.rel_error;
    return res.it_count;
}

template<> // specialisation: double precision
inline int BiCGSTAB<double> (const RCompRowMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit)
{
    const double *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    double *x_val = x.data_buffer();
    const double *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    cuda_BiCGSTAB (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, tol, maxit,
		   &res);
    tol = (double)res.rel_error;
    return res.it_count;
}

template<> // specialisation: double complex
inline int BiCGSTAB<complex> (const CCompRowMatrix &A, const CVector &b,
    CVector &x, double &tol, CPreconditioner *precon, int maxit)
{
    const complex *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    complex *x_val = x.data_buffer();
    const complex *b_val = b.data_buffer();
    int m = A.nRows();
    int n = A.nCols();
    SolverResult res;
    cuda_BiCGSTAB_cplx (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, tol,
			maxit, &res);
    tol = (double)res.rel_error;
    return res.it_count;
}

// ==========================================================================
// These versions use multiple right-hand sides

template<class MT>
void BiCGSTAB (const TCompRowMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x, int nrhs, double tol, int maxit,
    TPreconditioner<MT> *precon, IterativeSolverResult *res)
{
    BiCGSTAB ((const TMatrix<MT>&)A, b, x, nrhs, tol, maxit, precon, res);
    // by default, use the TMatrix implementation
}

template<> // specialisation: single precision
inline void BiCGSTAB<float> (const FCompRowMatrix &A, const FVector *b,
    FVector *x, int nrhs, double tol, int maxit, FPreconditioner *precon,
    IterativeSolverResult *res)
{
    static int nbuf = 32;
    static float **x_val = new float*[nbuf];
    static const float **b_val = new const float*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new float*[nbuf];
	b_val = new const float*[nbuf];
    }
    
    const float *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    float ftol = (float)tol;
    cuda_BiCGSTAB (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs, ftol,
        maxit, &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

template<> // specialisation: single complex
inline void BiCGSTAB<scomplex> (const SCCompRowMatrix &A, const SCVector *b,
    SCVector *x, int nrhs, double tol, int maxit,
    SCPreconditioner *precon, IterativeSolverResult *res)
{
    static int nbuf = 32;
    static scomplex **x_val = new scomplex*[nbuf];
    static const scomplex **b_val = new const scomplex*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new scomplex*[nbuf];
	b_val = new const scomplex*[nbuf];
    }
    
    const scomplex *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    float ftol = (float)tol;
    cuda_BiCGSTAB_cplx (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs,
        ftol, maxit, &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

template<> // specialisation: double precision
inline void BiCGSTAB<double> (const RCompRowMatrix &A, const RVector *b,
    RVector *x, int nrhs, double tol, int maxit, RPreconditioner *precon,
    IterativeSolverResult *res)
{
    static int nbuf = 32;
    static double **x_val = new double*[nbuf];
    static const double **b_val = new const double*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new double*[nbuf];
	b_val = new const double*[nbuf];
    }
    
    const double *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    cuda_BiCGSTAB (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs, tol,
        maxit, &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

template<> // specialisation: double complex
inline void BiCGSTAB<complex> (const CCompRowMatrix &A, const CVector *b,
    CVector *x, int nrhs, double tol, int maxit,
    CPreconditioner *precon, IterativeSolverResult *res)
{
    static int nbuf = 32;
    static complex **x_val = new complex*[nbuf];
    static const complex **b_val = new const complex*[nbuf];
    if (nrhs > nbuf) {
        nbuf = nrhs;
        delete []x_val;
	delete []b_val;
	x_val = new complex*[nbuf];
	b_val = new const complex*[nbuf];
    }
    
    const complex *A_val = A.ValPtr();
    const int *A_rowptr = A.rowptr;
    const int *A_colidx = A.colidx;
    for (int i = 0; i < nrhs; i++) {
        x_val[i] = x[i].data_buffer();
	b_val[i] = b[i].data_buffer();
    }
    int m = A.nRows();
    int n = A.nCols();
    SolverResult cuda_res;
    cuda_BiCGSTAB_cplx (A_val, A_rowptr, A_colidx, m, n, b_val, x_val, nrhs,
        tol, maxit, &cuda_res);
    if (res) {
        res->it_count = cuda_res.it_count;
	res->rel_err  = cuda_res.rel_error;
    }
}

#endif // CUDA_FLOAT

// ==========================================================================

template<class MT>
int TCompRowMatrix<MT>::pcg (const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit) const
{
    return TMatrix<MT>::pcg (b, x, tol, precon, maxit);
}

template<>
inline int TCompRowMatrix<float>::pcg (const FVector &b, FVector &x,
    double &tol, TPreconditioner<float> *precon, int maxit) const
{
    return PCG (*this, b, x, tol, precon, maxit);
}

template<>
inline int TCompRowMatrix<double>::pcg (const RVector &b, RVector &x,
    double &tol, TPreconditioner<double> *precon, int maxit) const
{
    return PCG (*this, b, x, tol, precon, maxit);
}

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::pcg (const TVector<MT> *b, TVector<MT> *x, int nrhs,
    double tol, int maxit, TPreconditioner<MT> *precon,
    IterativeSolverResult *res) const
{
    TMatrix<MT>::pcg (b, x, nrhs, tol, maxit, precon, res);
}

template<>
inline void TCompRowMatrix<float>::pcg (const FVector *b, FVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<float> *precon,
    IterativeSolverResult *res) const
{
    PCG (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<>
inline void TCompRowMatrix<double>::pcg (const RVector *b, RVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<double> *precon,
    IterativeSolverResult *res) const
{
    PCG (*this, b, x, nrhs, tol, maxit, precon, res);
}

// ==========================================================================

template<class MT>
int TCompRowMatrix<MT>::bicgstab (const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit) const
{
    //return TMatrix<MT>::bicgstab (b, x, tol, precon, maxit);
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

#ifdef UNDEF
template<>
int TCompRowMatrix<float>::bicgstab (const FVector &b, FVector &x,
    double &tol, TPreconditioner<float> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

template<>
int TCompRowMatrix<scomplex>::bicgstab (const SCVector &b, SCVector &x,
    double &tol, TPreconditioner<scomplex> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

template<>
int TCompRowMatrix<double>::bicgstab (const RVector &b, RVector &x,
    double &tol, TPreconditioner<double> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}

template<>
int TCompRowMatrix<complex>::bicgstab (const CVector &b, CVector &x,
    double &tol, TPreconditioner<complex> *precon, int maxit) const
{
    return BiCGSTAB (*this, b, x, tol, precon, maxit);
}
#endif

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::bicgstab (const TVector<MT> *b, TVector<MT> *x,
    int nrhs, double tol, int maxit, TPreconditioner<MT> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
    //TMatrix<MT>::bicgstab (b, x, nrhs, tol, maxit, precon, res);
}

#ifdef UNDEF
template<>
void TCompRowMatrix<float>::bicgstab (const FVector *b, FVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<float> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<>
void TCompRowMatrix<double>::bicgstab (const RVector *b, RVector *x, int nrhs,
    double tol, int maxit, TPreconditioner<double> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
}

template<>
void TCompRowMatrix<scomplex>::bicgstab (const SCVector *b, SCVector *x,
    int nrhs, double tol, int maxit, TPreconditioner<scomplex> *precon,
    IterativeSolverResult *res) const
{
    BiCGSTAB (*this, b, x, nrhs, tol, maxit, precon, res);
}
#endif

// ==========================================================================

template<class MT>
void TCompRowMatrix<MT>::ExportHB (ostream &os)
{
    int i, clines, ilines, vlines, tlines;

    if (!col_access) SetColAccess();  // storage in column-compressed format
    clines = (this->cols+10)/10;
    ilines = (this->nval+11)/12;
    vlines = (this->nval+4)/5;
    tlines = clines + ilines + vlines + 3;

    os.fill (' ');
    os.setf (ios::left, ios::adjustfield);
    os.width (72);
    os << "TOAST sparse matrix"; // title
    os.width (8);
    os << "";                    // key
    os << endl;
    os.setf (ios::right, ios::adjustfield);
    os.width (14);
    os << tlines;                // total number of lines ex header
    os.width (14);
    os << clines;                // number of lines for pointers
    os.width (14);
    os << ilines;                // number of lines for indices
    os.width (14);
    os << vlines;                // number of lines for values
    os.width (14);
    os << 0;                     // number of right hand sides
    os << endl;
    os.width (14);
    os.setf (ios::left, ios::adjustfield);
    os << "RUA";                 // matrix type (real symmetric)
    os.width (14);
    os.setf (ios::right, ios::adjustfield);
    os << this->rows;            // number of rows
    os.width (14);
    os << this->cols;            // number of columns
    os.width (14);
    os << this->nval;            // number of entries
    os.width (14);
    os << 0;                     // number of element matrix ents.
    os << endl;
    os.width (16);
    os.setf (ios::left, ios::adjustfield);
    os << "(10I8)";                          // format for pointers
    os.width (16);
    os << "(12I6)";                          // format for indices
    os.width (20);
    os << "(5E16.8)";                        // format for values
    os.width (20);
    os << "";                                // format for rhs
    os << endl;
    os.setf (ios::right, ios::adjustfield);

    for (i = 0; i <= this->cols; i++) {      // column pointers
        os.width (8);
        os << colptr[i]+1;
	if (!((i+1)%10) || i == this->cols)
	    os << endl;
    }
    for (i = 0; i < this->nval; i++) {       // row indices
        os.width (6);
        os << rowidx[i]+1;
	if (!((i+1)%12) || i == this->nval-1)
	    os << endl;
    }
    os.setf (ios::scientific, ios::floatfield);
    os.precision (8);
    for (i = 0; i < this->nval; i++) {
        os.width (16);
        os << this->val[i];
	if (!((i+1)%5) || i == this->nval-1)
	    os << endl;
    }
}

template<class MT>
void TCompRowMatrix<MT>::ExportRCV (ostream &os)
{
    int r, c, j;
    os << std::scientific << std::setprecision(10);
    for (r = 0; r < this->rows; r++) {
	for (j = rowptr[r]; j < rowptr[r+1]; j++) {
	    c = colidx[j];
	    os << r+1 << '\t' << c+1 << '\t' << this->val[j] << endl;
	}
    }
}

template<> // specialisation: complex
inline void TCompRowMatrix<std::complex<double> >::ExportRCV (ostream &os)
{
    int r, c, j;
    for (r = 0; r < this->rows; r++) {
	for (j = rowptr[r]; j < rowptr[r+1]; j++) {
	    c = colidx[j];
	    os << r+1 << '\t' << c+1 << '\t'
	       << this->val[j].real() << '\t' << this->val[j].imag() << endl;
	}
    }
}


template<class MT>
void TCompRowMatrix<MT>::SplitExport (const char *rootname)
{
    char cbuf[256];
    int i, j;

    // row index file
    strcpy (cbuf, rootname);
    strcat (cbuf, "_row.dat");
    ofstream ofs1 (cbuf);
    for (i = 0; i < this->rows; i++)
        for (j = rowptr[i]; j < rowptr[i+1]; j++)
	    ofs1 << i << endl;

    // column index file
    strcpy (cbuf, rootname);
    strcat (cbuf, "_col.dat");
    ofstream ofs2 (cbuf);
    for (i = 0; i < this->nval; i++)
        ofs2 << colidx[i] << endl;

    // data file
    strcpy (cbuf, rootname);
    strcat (cbuf, "_val.dat");
    ofstream ofs3 (cbuf);
    for (i = 0; i < this->nval; i++)
        ofs3 << this->val[i] << endl;
}

template<class MT>
istream &operator>> (istream &is, TCompRowMatrix<MT> &m)
{
    char cbuf[256];
    int i, nr, nc, nz;

    is.getline (cbuf, 256);
    if (strncmp (cbuf, "TCompRow", 8)) return is; // should set error flag
    sscanf (cbuf+8, "%d%d%d", &nr, &nc, &nz);
    if (nr == m.rows && nc == m.cols && nz == m.nval) { // no need to realloc
        for (i = 0; i <= nr; i++)
	    is >> m.rowptr[i];
	for (i = 0; i < nz; i++)
	    is >> m.colidx[i];
	for (i = 0; i < nz; i++)
	    is >> m.val[i];
    } else {
        idxtype *rowptr = new idxtype[nr+1];
	idxtype *colidx = new idxtype[nz];
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
ostream &operator<< (ostream &os, const TCompRowMatrix<MT> &m)
{
    int i;

    // header
    os << "TCompRow " << m.rows << ' ' << m.cols << ' ' << m.nval << endl;
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
// Interface between toast and the ML multigrid preconditioning package.

#ifdef ML_INTERFACE

template<class MT>
int ML_matvec (ML_Operator *Amat, int in_length, double p[],
    int out_length, double ap[])
{
    ERROR_UNDEF;
    return 0;
}

template<> // specialisation: double
int ML_matvec<double> (ML_Operator *Amat, int in_length, double p[],
    int out_length, double ap[])
{
    RCompRowMatrix *A = (RCompRowMatrix*)ML_Get_MyMatvecData (Amat);
    RVector b (in_length, p);
    RVector x (out_length);
    x = Ax (*A, b);
    memcpy (ap, x.data_buffer(), out_length*sizeof(double));
    // would be better if x could be constructed so as to directly use ap
    // as its data buffer
    return 1;
}

template<class MT>
int ML_getrow (ML_Operator *Amat, int N_requested_rows,
    int requested_rows[], int allocated_space, int columns[],
    double values[], int row_lenghts[])
{
    ERROR_UNDEF;
    return 0;
}

template<> // specialisation: double
int ML_getrow<double> (ML_Operator *Amat, int N_requested_rows,
    int requested_rows[], int allocated_space, int columns[],
    double values[], int row_lenghts[])
{
    int i, j, r, c;
    RCompRowMatrix *A = (RCompRowMatrix*)ML_Get_MyGetrowData (Amat);
    
    for (i = j = 0; i < N_requested_rows; i++) {
	r = requested_rows[i];
	for (c = A->rowptr[r]; c < A->rowptr[r+1]; c++) {
	    if (j == allocated_space) return 0; // out of space
	    columns[j] = A->colidx[c];
	    values[j] = A->val[c];
	    j++;
	}
	row_lengths[i] = A->rowptr[r+1] - A->rowptr[r];
    }
    return 0;
}

#endif // ML_INTERFACE

// ==========================================================================
// class and friend instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TCompRowMatrix<double>;
template class MATHLIB TCompRowMatrix<float>;
template class MATHLIB TCompRowMatrix<toast::complex>;
template class MATHLIB TCompRowMatrix<scomplex>;
template class MATHLIB TCompRowMatrix<int>;

template MATHLIB TCompRowMatrix<double> transp (const TCompRowMatrix<double> &m);
template MATHLIB TCompRowMatrix<toast::complex> transp (const TCompRowMatrix<toast::complex> &m);

template MATHLIB TCompRowMatrix<double> cath (const TCompRowMatrix<double> &A,
    const TCompRowMatrix<double> &B);
template MATHLIB TCompRowMatrix<toast::complex> cath (const TCompRowMatrix<toast::complex> &A,
    const TCompRowMatrix<toast::complex> &B);

template MATHLIB TCompRowMatrix<double> catv (const TCompRowMatrix<double> &A,
    const TCompRowMatrix<double> &B);
template MATHLIB TCompRowMatrix<toast::complex> catv (const TCompRowMatrix<toast::complex> &A,
    const TCompRowMatrix<toast::complex> &B);

template MATHLIB double l2norm (const TCompRowMatrix<double> &A);
template MATHLIB double l2norm (const TCompRowMatrix<toast::complex> &A);
template MATHLIB double l2norm (const TCompRowMatrix<scomplex> &A);

template MATHLIB TCompRowMatrix<double> kron (const TCompRowMatrix<double> &A,
    const TCompRowMatrix<double> &B);
template MATHLIB TCompRowMatrix<toast::complex> kron (const TCompRowMatrix<toast::complex> &A,
    const TCompRowMatrix<toast::complex> &B);
template MATHLIB TCompRowMatrix<scomplex> kron (const TCompRowMatrix<scomplex> &A,
    const TCompRowMatrix<scomplex> &B);
template MATHLIB TCompRowMatrix<int> kron (const TCompRowMatrix<int> &A,
    const TCompRowMatrix<int> &B);

template MATHLIB bool CholeskyFactorize (const FCompRowMatrix &A,
    FCompRowMatrix &L, FVector &d, bool recover);
template MATHLIB bool CholeskyFactorize (const RCompRowMatrix &A,
    RCompRowMatrix &L, RVector &d, bool recover);

// complex ones probably not allowed

template bool CholeskyFactorize (const CCompRowMatrix &A, CCompRowMatrix &L,
    CVector &d, bool recover);
template bool CholeskyFactorize (const SCCompRowMatrix &A, SCCompRowMatrix &L,
    SCVector &d, bool recover);

template MATHLIB bool IncompleteCholeskyFactorize (const FCompRowMatrix &A,
    FCompRowMatrix &L, FVector &d, bool recover);
template MATHLIB bool IncompleteCholeskyFactorize (const RCompRowMatrix &A,
    RCompRowMatrix &L, RVector &d, bool recover);
template MATHLIB bool IncompleteCholeskyFactorize (const CCompRowMatrix &A,
    CCompRowMatrix &L, CVector &d, bool recover);
template MATHLIB bool IncompleteCholeskyFactorize (const SCCompRowMatrix &A,
    SCCompRowMatrix &L, SCVector &d, bool recover);

template MATHLIB void LUFactorize (RCompRowMatrix &A, bool LUrealloc);
template MATHLIB void LUFactorize (CCompRowMatrix &A, bool LUrealloc);
template MATHLIB void LUFactorize (SCCompRowMatrix &A, bool LUrealloc);

template MATHLIB void CholeskySolve (const FCompRowMatrix &L,
    const FVector &d, const FVector &b, FVector &x);
template MATHLIB void CholeskySolve (const RCompRowMatrix &L,
    const RVector &d, const RVector &b, RVector &x);

template MATHLIB void CholeskySolve (const RCompRowMatrix &L,
    const RVector &d, const RDenseMatrix &BT, RDenseMatrix&XT, int n);
template MATHLIB void CholeskySolve (const CCompRowMatrix &L,
    const CVector &d, const CVector &b, CVector &x);
template MATHLIB void CholeskySolve (const SCCompRowMatrix &L, const SCVector &d,
    const SCVector &b, SCVector &x);

template MATHLIB RVector CholeskySolve (const RCompRowMatrix &L,
    const RVector &d, const RVector &b);
template MATHLIB CVector CholeskySolve (const CCompRowMatrix &L,
    const CVector &d, const CVector &b);
template MATHLIB SCVector CholeskySolve (const SCCompRowMatrix &L,
    const SCVector &d, const SCVector &b);

template MATHLIB void LUSolve (const RCompRowMatrix &LU, const RVector &b,
    RVector &x);
template MATHLIB void LUSolve (const CCompRowMatrix &LU, const CVector &b,
    CVector &x);
template MATHLIB void LUSolve (const SCCompRowMatrix &LU, const SCVector &b,
    SCVector &x);

#ifdef USE_CUDA_FLOAT
//template MATHLIB int PCG (const RCompRowMatrix &A, const RVector &b,
//    RVector &x, double &tol, const RPreconditioner *precon, int maxit);

//template MATHLIB void PCG (const RCompRowMatrix &A, const RVector *b,
//    RVector *x, int nrhs, double tol, int maxit, const RPreconditioner *precon,
//    IterativeSolverResult *res);

//template MATHLIB int BiCGSTAB (const RCompRowMatrix &A, const RVector &b,
//    RVector &x, double &tol, const RPreconditioner *precon, int maxit);
//template MATHLIB int BiCGSTAB (const CCompRowMatrix &A, const CVector &b,
//    CVector &x, double &tol, const CPreconditioner *precon, int maxit);

//template MATHLIB void BiCGSTAB (const RCompRowMatrix &A, const RVector *b,
//    RVector *x, int nrhs, double tol, int maxit, const RPreconditioner *precon,
//    IterativeSolverResult *res);
//template MATHLIB void BiCGSTAB (const CCompRowMatrix &A, const CVector *b,
//    CVector *x, int nrhs, double tol, int maxit, const CPreconditioner *precon,
//    IterativeSolverResult *res);
template MATHLIB int BiCGSTAB (const ICompRowMatrix &A, const IVector &b,
    IVector &x, double &tol, IPreconditioner *precon, int maxit);
template MATHLIB void BiCGSTAB (const ICompRowMatrix &A, const IVector *b,
    IVector *x, int nrhs, double tol, int maxit, IPreconditioner *precon,
    IterativeSolverResult *res);
#endif

template MATHLIB istream &operator>> (istream &is, RCompRowMatrix &m);
template MATHLIB istream &operator>> (istream &is, CCompRowMatrix &m);
template MATHLIB istream &operator>> (istream &is, SCCompRowMatrix &m);
template MATHLIB istream &operator>> (istream &is, ICompRowMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const RCompRowMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const CCompRowMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const SCCompRowMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const ICompRowMatrix &m);

#endif // NEED_EXPLICIT_INSTANTIATION
