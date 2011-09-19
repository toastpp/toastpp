// ==========================================================================
// Module mathlib
// File dnsmatrix.cc
// Definition of template class TDenseMatrix ('template dense matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <sstream>

#include "mathlib.h"
#include "fblas.h"

using namespace std;
using namespace toast;

template<class MT>
TDenseMatrix<MT>::TDenseMatrix (const TDenseMatrix<MT> &m,
    int i0, int j0, int i1, int j1): TMatrix<MT> ()
{
    if (i1 == END) i1 = m.nRows();
    if (j1 == END) j1 = m.nCols();

    dASSERT(i0 >= 0, "Argument 2 out of range");
    dASSERT(i1 >= i0, "Argument 4 >= argument 2 required");
    dASSERT(m.nRows() >= i1, "Argument 4 out of range");
    dASSERT(j0 >= 0, "Argument 3 out of range");
    dASSERT(j1 >= j0, "Argument 5 >= argument 3 required");
    dASSERT(m.nCols() >= j1, "Argument 5 out of range");

    int i, j;
    this->rows = i1-i0;
    this->cols = j1-j0;
    Alloc (this->rows, this->cols);

    MT *valt = val;
    MT *vals = m.val + (i0*m.nCols()+j0);
    for (i = 0; i < this->rows; i++) {
	memcpy (valt, vals, this->cols*sizeof(MT));
	valt += this->cols;
	vals += m.nCols();
    }
}

template<class MT>
TDenseMatrix<MT>::TDenseMatrix (int r, int c, const char *valstr)
: TMatrix<MT> (r, c)
{
    Alloc (r, c);
    std::istringstream iss (valstr);
    for (int i = 0; i < rc; i++) iss >> val[i];
}

template<class MT>
TDenseMatrix<MT>::TDenseMatrix (const TSymMatrix<MT> &A)
: TMatrix<MT> (A.nRows(), A.nCols())
{
    int i, j, n = A.nRows();
    Alloc (n, n);
    MT *v = val;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) *v++ = A(i,j);
}

template<class MT>
TDenseMatrix<MT>::TDenseMatrix (const TCompRowMatrix<MT> &A)
{
    int r, c, i;
    int m = A.nRows();
    int n = A.nCols();
    Alloc (m, n);
    const MT *Aval = A.ValPtr();

    for (r = 0; r < m; r++) {
	for (i = A.rowptr[r]; i < A.rowptr[r+1]; i++) {
	    c = A.colidx[i];
	    val[r*n+c] = Aval[i];
	}
    }
}

template<class MT>
void TDenseMatrix<MT>::New_dirty (int r, int c)
{
    TMatrix<MT>::New (r, c); // set nrows and ncols
    if (rc != r*c) {             // realloc only if size differs
        Unlink ();               // dealloc current data array
        Alloc (r, c);            // and re-allocate with new size
    }
}

template<class MT>
TVector<MT> TDenseMatrix<MT>::Col (int c) const
{
    TVector<MT> colvec(this->rows);
    for (int r = 0; r < this->rows; r++)
        colvec[r] = val[r*this->cols+c];
    return colvec;
}

template<class MT>
int TDenseMatrix<MT>::SparseRow (int r, idxtype *ci, MT *rv) const
{
    MT *rv0 = val+r*this->cols;
    for (int i = 0; i < this->cols; i++) {
        ci[i] = i;
	rv[i] = rv0[i];
    }
    return this->cols;
}

template<class MT>
void TDenseMatrix<MT>::SetRow (int r, const TVector<MT> &rval)
{
    RANGE_CHECK(r >= 0 && r < this->rows);
    dASSERT(rval.Dim() == this->cols, "Argument 2: wrong size");
    memcpy (val + (r*this->cols), rval.data_buffer(),
	    this->cols*sizeof(MT));
}

template<class MT>
void TDenseMatrix<MT>::ColScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->cols, "Argument 1: wrong size");
    for (int r = 0; r < this->rows; r++)
        for (int c = 0; c < this->cols; c++)
	    val[r*this->cols+c] *= scale[c];
}

template<class MT>
void TDenseMatrix<MT>::RowScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->rows, "Argument 1: wrong size");
    for (int r = 0; r < this->rows; r++)
	for (int c = 0; c < this->cols; c++)
	    val[r*this->cols+c] *= scale[r];
}

template<class MT>
TVector<MT> TDenseMatrix<MT>::RowSum () const
{
    int r, c, nr = this->rows, nc = this->cols;

    TVector<MT> res(nr);
    MT sum, *idx = val;
    for (r = 0; r < nr; r++) {
	sum = (MT)0;
	for (c = 0; c < nc; c++)
	    sum += *idx++;
	res[r] = sum;
    }
    return res;
}

template<class MT>
TVector<MT> TDenseMatrix<MT>::RowSumSq () const
{
    int r, c, nr = this->rows, nc = this->cols;

    TVector<MT> res(nr);
    MT sumsq, *idx = val;
    for (r = 0; r < nr; r++) {
	sumsq = (MT)0;
	for (c = 0; c < nc; c++) {
	    sumsq += (*idx)*(*idx);
	    idx++;
	}
	res[r] = sumsq;
    }
    return res;
}

template<class MT>
TVector<MT> TDenseMatrix<MT>::ColSum () const
{
    int r, c, nr = this->rows, nc = this->cols;

    TVector<MT> res(nc);
    MT *resptr = res.data_buffer();
    MT *idx = val;
    for (r = 0; r < nr; r++) {
	for (c = 0; c < nc; c++)
	    resptr[c] += *idx++;
    }
    return res;
}

template<class MT>
TVector<MT> TDenseMatrix<MT>::ColSumSq () const
{
    int r, c, nr = this->rows, nc = this->cols;

    TVector<MT> res(nc);
    MT *resptr = res.data_buffer();
    MT *idx = val;
    for (r = 0; r < nr; r++) {
	for (c = 0; c < nc; c++) {
	    resptr[c] += (*idx)*(*idx);
	    idx++;
	}
    }
    return res;
}

template<class MT>
void TDenseMatrix<MT>::Identity ()
{
    Zero();
    int i = (this->rows < this->cols ? this->rows : this->cols) - 1;
    for (; i >= 0; i--)
        val[i*this->cols + i] = (MT)1;
}

// ===========================================================================
// Ax(): Compute matrix x vector product

template<class MT>
void TDenseMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != this->rows) b.New (this->rows);  // resize

    int r, c;
    MT *row;
    for (r = 0; r < this->rows; r++) {
        MT &br = b[r];
	row = val + (r*this->cols);
	for (c = 0, br = (MT)0; c < this->cols; c++)
	    br += row[c] * x[c];
    }
}
#ifdef USE_BLAS_LEVEL2 // interface to BLAS level 2 xGEMV functions
template<>
MATHLIB void TDenseMatrix<double>::Ax (const TVector<double> &x, TVector<double> &b)
    const
{
    dASSERT(cols == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != rows) b.New (rows);  // resize

    // we need to transpose and flip the row and column dimensions
    // to convert from c++ row storage to fortran column storage
    static char trans = 'T';
    static double alpha = 1.0, beta = 0.0;
    static int incr = 1;

    dgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
template<>
void TDenseMatrix<float>::Ax (const TVector<float> &x, TVector<float> &b) const
{
    dASSERT(cols == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != rows) b.New (rows);  // resize;

    // we need to transpose and flip the row and column dimensions
    // to convert from c++ row storage to fortran column storage
    static char trans = 'T';
    static float alpha = 1.0f, beta = 0.0f;
    static int incr = 1;

    sgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
template<>
void TDenseMatrix<toast::complex>::Ax (const TVector<toast::complex> &x,
    TVector<toast::complex> &b) const
{
    dASSERT(cols == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != rows) b.New (rows);  // resize;

    // we need to transpose and flip the row and column dimensions
    // to convert from c++ row storage to fortran column storage
    static char trans = 'T';
    static toast::complex alpha(1,0), beta(0,0);
    static int incr = 1;

    zgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
#endif // USE_BLAS_LEVEL2

// ===========================================================================
// ATx(): Compute transpose(matrix) x vector product

template<class MT>
void TDenseMatrix<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->rows == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != this->cols) b.New (this->cols);  // resize
    
    int r, c;
    MT *col;
    for (c = 0; c < this->cols; c++) {
	col = val+c;
	MT &bc = b[c];
	for (r = 0, bc = (MT)0; r < this->rows; r++)
	    bc += col[r*this->cols] * x[r];
    }
}
#ifdef USE_BLAS_LEVEL2 // interface to BLAS level 2 xGEMV functions
template<>
MATHLIB void TDenseMatrix<double>::ATx (const TVector<double> &x, TVector<double> &b)
    const
{
    dASSERT(rows == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != cols) b.New (cols);  // resize
    
    static char trans = 'N';
    static double alpha = 1.0, beta = 0.0;
    static int incr = 1;

    dgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
template<>
void TDenseMatrix<float>::ATx (const TVector<float> &x, TVector<float> &b)
    const
{
    dASSERT(rows == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != cols) b.New (cols);  // resize
    
    static char trans = 'N';
    static float alpha = 1.0, beta = 0.0;
    static int incr = 1;

    sgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
template<>
void TDenseMatrix<toast::complex>::ATx (const TVector<toast::complex> &x,
    TVector<toast::complex> &b) const
{
    dASSERT(rows == x.Dim(), "Argument 1: vector has wrong dimension");
    if (b.Dim() != cols) b.New (cols);  // resize
    
    static char trans = 'N';
    static toast::complex alpha(1,0), beta(0,0);
    static int incr = 1;

    zgemv_(trans, (int&)cols, (int&)rows, alpha, val, (int&)cols,
           x.data_buffer(), incr, beta, b.data_buffer(), incr);
}
#endif // USE_BLAS_LEVEL2

// ===========================================================================
// AB(): Compute matrix x matrix product

template<class MT>
void TDenseMatrix<MT>::AB (const TDenseMatrix<MT> &A,
    const TDenseMatrix<MT> &B)
{
    dASSERT(A.cols == B.rows, "Matrix sizes do not match for this operation.");

    int r, c, i, rAcols, Acols = A.cols, Bcols = B.cols;
    if (this->rows != A.rows || this->cols != B.cols)
        New_dirty (A.rows, B.cols);      // resize
    for (r = 0; r < this->rows; r++) {
        rAcols = r*Acols;
        for (c = 0; c < this->cols; c++) {
	    MT &v = val[r*this->cols+c];
	    for (i = 0, v = 0; i < Acols; i++)
		v += A.val[rAcols+i] * B.val[i*Bcols+c];
	}
    }
}
#ifdef USE_BLAS_LEVEL3 // interface to BLAS level 3 xGEMM functions
template<>
MATHLIB void TDenseMatrix<double>::AB (const TDenseMatrix<double> &A,
    const TDenseMatrix<double> &B)
{
    dASSERT(A.cols == B.rows, "Matrix sizes do not match for this operation.");

    int Acols = A.cols, Bcols = B.cols;
    if (rows != A.rows || cols != B.cols)
        New_dirty (A.rows, B.cols);      // resize

    // Note: Column storage in Fortran means that matrix A and B are
    // interpreted as A = Af^T, B = Bf^T, and we need to calculate
    // Cf = (AB)^T = (Af^T Bf^T)^T = Bf Af
    // i.e. we need to flip the order of the operands

    static char transa = 'N', transb = 'N';
    static double alpha = 1.0, beta = 0.0;

    dgemm_(transa, transb, (int&)B.cols, (int&)A.rows, (int&)A.cols, alpha,
	   B.val, (int&)B.cols, A.val, (int&)A.cols, beta, val, (int&)B.cols);
}
template<>
MATHLIB void TDenseMatrix<float>::AB (const TDenseMatrix<float> &A,
    const TDenseMatrix<float> &B)
{
    dASSERT(A.cols == B.rows, "Matrix sizes do not match for this operation.");

    int Acols = A.cols, Bcols = B.cols;
    if (rows != A.rows || cols != B.cols)
        New_dirty (A.rows, B.cols);      // resize

    // Note: Column storage in Fortran means that matrix A and B are
    // interpreted as A = Af^T, B = Bf^T, and we need to calculate
    // Cf = (AB)^T = (Af^T Bf^T)^T = Bf Af
    // i.e. we need to flip the order of the operands

    static char transa = 'N', transb = 'N';
    static float alpha = 1.0f, beta = 0.0f;

    sgemm_(transa, transb, (int&)B.cols, (int&)A.rows, (int&)A.cols, alpha,
	   B.val, (int&)B.cols, A.val, (int&)A.cols, beta, val, (int&)B.cols);
}
template<>
MATHLIB void TDenseMatrix<toast::complex>::AB (const TDenseMatrix<toast::complex> &A,
    const TDenseMatrix<toast::complex> &B)
{
    dASSERT(A.cols == B.rows, "Matrix sizes do not match for this operation.");

    int Acols = A.cols, Bcols = B.cols;
    if (rows != A.rows || cols != B.cols)
        New_dirty (A.rows, B.cols);      // resize

    // Note: Column storage in Fortran means that matrix A and B are
    // interpreted as A = Af^T, B = Bf^T, and we need to calculate
    // Cf = (AB)^T = (Af^T Bf^T)^T = Bf Af
    // i.e. we need to flip the order of the operands

    static char transa = 'N', transb = 'N';
    static toast::complex alpha(1,0), beta = (0,0);

    zgemm_(transa, transb, (int&)B.cols, (int&)A.rows, (int&)A.cols, alpha,
	   B.val, (int&)B.cols, A.val, (int&)A.cols, beta, val, (int&)B.cols);
}
#endif // USE_BLAS_LEVEL3

// ===========================================================================
// ATA(): Compute A^T A matrix x matrix product

template<class MT>
MATHLIB TSymMatrix<MT> ATA (const TDenseMatrix<MT> &A)
{
    int r, c, k, nr = A.nRows(), nc = A.nCols();
    MT *Ar, *Ac, v;
    TSymMatrix<MT> ata(nc);
    for (r = 0; r < nc; r++) {
        Ar = A.val + r;
	for (c = 0; c <= r; c++) {
	    Ar = A.val + r;
	    Ac = A.val + c;
	    v = (MT)0;
	    //MT &ata_rc = ata(r,c);
	    for (k = 0; k < nr; k++) {
	        v += *Ar * *Ac;
		Ar += nc;
		Ac += nc;
	        //ata_rc += Ar[k*nc] * Ac[k*nc];
	    }
	    ata(r,c) = v;
	}
    }
    return ata;
}

#ifdef USE_BLAS_LEVEL3 // interface to BLAS level 3 xSYRK functions
template<>
MATHLIB TSymMatrix<double> ATA (const TDenseMatrix<double> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'N'; // 'N' actually indicates AAT, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request ATA
    static double alpha = 1.0;
    static double beta = 0.0;

    int n = A.nCols(), k = A.nRows();
    double *c = new double[n*n];

    dsyrk_(uplo, trans, n, k, alpha, A.val, n, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<double> ata(n);
    double *ata_buf = ata.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    ata_buf[idx++] = c[i*n+j];
    delete []c;
    return ata;
}
template<>
MATHLIB TSymMatrix<float> ATA (const TDenseMatrix<float> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'N'; // 'N' actually indicates AAT, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request ATA
    static float alpha = 1.0f;
    static float beta = 0.0f;

    int n = A.nCols(), k = A.nRows();
    float *c = new float[n*n];

    ssyrk_(uplo, trans, n, k, alpha, A.val, n, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<float> ata(n);
    float *ata_buf = ata.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    ata_buf[idx++] = c[i*n+j];
    delete []c;
    return ata;
}
template<>
MATHLIB TSymMatrix<toast::complex> ATA (const TDenseMatrix<toast::complex> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'N'; // 'N' actually indicates AAT, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request ATA
    static toast::complex alpha = toast::complex (1.0, 0.0);
    static toast::complex beta = toast::complex (0.0, 0.0);

    int n = A.nCols(), k = A.nRows();
    toast::complex *c = new toast::complex[n*n];

    zsyrk_(uplo, trans, n, k, alpha, A.val, n, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<toast::complex> ata(n);
    toast::complex *ata_buf = ata.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    ata_buf[idx++] = c[i*n+j];
    delete []c;
    return ata;
}
#endif // USE_BLAS_LEVEL3

// ===========================================================================
// AAT(): Compute A A^T matrix x matrix product

template<class MT>
TSymMatrix<MT> AAT (const TDenseMatrix<MT> &A)
{
    int r, c, k, nr = A.nRows(), nc = A.nCols();
    MT *Ar, *Ac;
    TSymMatrix<MT> aat(nr);
    for (r = 0; r < nr; r++) {
	for (c = 0; c <= r; c++) {
	    Ar = A.val + r*nc;
	    Ac = A.val + c*nc;
	    MT &aat_rc = aat(r,c);
	    for (k = 0; k < nc; k++)
		aat_rc += *Ar++ * *Ac++;
	}
    }
    return aat;
}
#ifdef USE_BLAS_LEVEL3 // interface to BLAS level 3 xSYRK functions
template<>
MATHLIB TSymMatrix<double> AAT (const TDenseMatrix<double> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'T'; // 'T' actually indicates ATA, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request AAT
    static double alpha = 1.0;
    static double beta = 0.0;

    int n = A.nRows(), k = A.nCols();
    double *c = new double[n*n];

    dsyrk_(uplo, trans, n, k, alpha, A.val, k, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<double> aat(n);
    double *aat_buf = aat.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    aat_buf[idx++] = c[i*n+j];
    delete []c;
    return aat;
}
template<>
MATHLIB TSymMatrix<float> AAT (const TDenseMatrix<float> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'T'; // 'T' actually indicates ATA, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request AAT
    static float alpha = 1.0;
    static float beta = 0.0;

    int n = A.nRows(), k = A.nCols();
    float *c = new float[n*n];

    ssyrk_(uplo, trans, n, k, alpha, A.val, k, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<float> aat(n);
    float *aat_buf = aat.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    aat_buf[idx++] = c[i*n+j];
    delete []c;
    return aat;
}
template<>
MATHLIB TSymMatrix<toast::complex> AAT (const TDenseMatrix<toast::complex> &A)
{
    static char uplo = 'U';  // return result in upper triangle
    static char trans = 'T'; // 'T' actually indicates ATA, but since fortran
                             // stores in column format, all matrices are
                             // interpreted as transpose, so we actually
                             // request AAT
    static toast::complex alpha = toast::complex (1.0, 0.0);
    static toast::complex beta = toast::complex (0.0, 0.0);

    int n = A.nRows(), k = A.nCols();
    toast::complex *c = new toast::complex[n*n];

    zsyrk_(uplo, trans, n, k, alpha, A.val, k, beta, c, n);

    // now we need to copy the results into a SymMatrix
    TSymMatrix<toast::complex> aat(n);
    toast::complex *aat_buf = aat.data_buffer();
    int i, j, idx;
    for (i = idx = 0; i < n; i++)
        for (j = 0; j <= i; j++)
	    aat_buf[idx++] = c[i*n+j];
    delete []c;
    return aat;
}
#endif // USE_BLAS_LEVEL3

// ===========================================================================

template<class MT>
TDenseMatrix<MT> kron (const TDenseMatrix<MT> &A, const TDenseMatrix<MT> &B)
{
    int ia, ja, ib, jb;
    int na = A.nRows(), ma = A.nCols();
    int nb = B.nRows(), mb = B.nCols();
    int n = na*nb, m = ma*mb;
    MT a_ij;

    TDenseMatrix<MT> C(n,m);
    
    for (ia = 0; ia < na; ia++) {
	for (ja = 0; ja < ma; ja++) {
	    a_ij = A(ia,ja);
	    for (ib = 0; ib < nb; ib++) {
		for (jb = 0; jb < mb; jb++)
		    C(ia*nb+ib, ja*mb+jb) = a_ij * B(ib,jb);
	    }
	}
    }
    return C;
}

// ===========================================================================

template<class MT>
MT TDenseMatrix<MT>::colmult (int c1, int c2)
{
    int dc = c2-c1;
    MT sum = (MT)0;
    MT *v1 = val+c1;
    for (int i = 0; i < this->rows; i++) {
        sum += *v1 * *(v1+dc);
	v1 += this->cols;
    }
    return sum;
}

template<class MT>
TDenseMatrix<MT> &TDenseMatrix<MT>::operator= (const TDenseMatrix<MT> &m)
{
    if (this->rows != m.rows || this->cols != m.cols)
        New_dirty (m.rows, m.cols);     // resize
    memcpy (val, m.val, rc*sizeof(MT)); // copy data vector
    return *this;
}

template<class MT>
TDenseMatrix<MT> &TDenseMatrix<MT>::operator= (MT v)
{
    for (int i = 0; i < this->rows * this->cols; i++) val[i] = v;
    return *this;
}

template<class MT>
TDenseMatrix<MT> TDenseMatrix<MT>::operator+ (const TDenseMatrix<MT> &m) const
{
    dASSERT(this->rows == m.rows && this->cols == m.cols,
	    "Incompatible matrix dimensions");

    TDenseMatrix<MT> tmp(m);
    for (int i = 0; i < rc; i++) tmp.val[i] += val[i];
    return tmp;
}

template<class MT>
TDenseMatrix<MT> TDenseMatrix<MT>::operator- (const TDenseMatrix<MT> &m) const
{
    dASSERT(this->rows == m.rows && this->cols == m.cols,
	    "Incompatible matrix dimensions");

    TDenseMatrix<MT> tmp(*this);
    for (int i = 0; i < rc; i++) tmp.val[i] -= m.val[i];
    return tmp;
}

template<class MT>
MATHLIB TDenseMatrix<MT> transpose (const TDenseMatrix<MT> &A)
{
    int Arows = A.nRows();
    int Acols = A.nCols();

    TDenseMatrix<MT> AT(Acols, Arows);
    for (int r = 0; r < Arows; r++)
        for (int c = 0; c < Acols; c++)
	    AT(c,r) = A(r,c);
    return AT;
}

// ===========================================================================

template<class MT>
MATHLIB TDenseMatrix<MT> cath (const TDenseMatrix<MT> &A,
    const TDenseMatrix<MT> &B)
{
    // Concatenates matrices A and B horizontally

    int i, ac, bc, br, nr, nc;
    ac = A.nCols();
    bc = B.nCols();
    nr = A.nRows();
    br = B.nRows();

    // Check trivial cases of one matrix of zero dimensions
    if (ac == 0 && nr == 0) return B;
    if (bc == 0 && br == 0) return A;

    nc = ac+bc;
    xASSERT (nr == B.nRows(), "Matrix row dimensions do not match");

    TDenseMatrix<MT> C (nr, nc);
    MT *valt  = C.val;
    MT *vals1 = A.val;
    MT *vals2 = B.val;
    
    for (i = 0; i < nr; i++) {
	memcpy (valt, vals1, ac*sizeof(MT));
	valt += ac; vals1 += ac;
	memcpy (valt, vals2, bc*sizeof(MT));
	valt += bc; vals2 += bc;
    }
    return C;
}

// ===========================================================================

template<class MT>
MATHLIB TDenseMatrix<MT> catv (const TDenseMatrix<MT> &A,
    const TDenseMatrix<MT> &B)
{
    // Concatenates matrices A and B vertically
    
    int i, ar, br, bc, nr, nc;
    ar = A.nRows();
    br = B.nRows();
    nc = A.nCols();
    bc = B.nRows();

    // Check trivial cases of one matrix of zero dimensions
    if (nc == 0 && ar == 0) return B;
    if (bc == 0 && br == 0) return A;

    nr = ar+br;
    xASSERT (nc = B.nCols(), "Matrix column dimensions do not match");

    TDenseMatrix<MT> C (nr, nc);
    MT *valt  = C.val;
    MT *vals1 = A.val;
    MT *vals2 = B.val;

    for (i = 0; i < ar; i++) {
	memcpy (valt, vals1, nc*sizeof(MT));
	valt += nc; vals1 += nc;
    }
    for (i = 0; i < br; i++) {
	memcpy (valt, vals2, nc*sizeof(MT));
	valt += nc; vals2 += nc;
    }
    return C;
}

// ===========================================================================

template<class MT>
MATHLIB MT det (const TDenseMatrix<MT> &A, TDenseMatrix<MT> *Ai)
{
    xASSERT(A.nRows() == A.nCols(),
	    "Need square matrix to compute determinant");
    const double detmin = 1e-50;
    int n = A.nCols();
    MT d;
#ifdef FEM_DEBUG
    if (Ai) dASSERT(Ai->rows == n && Ai->cols == n,
		    "Argument 2 matrix has wrong dimension");
#endif

    switch (n) {
    case 0:
        d = (MT)0;
	break;
    case 1:
        d = A(0,0);
	if (Ai) {
	    xASSERT(norm(d) > detmin, "Matrix singular");
	    Ai->Set (0,0,(MT)1/d);
	}
	break;
    case 2:
        d = A(0,0)*A(1,1) - A(0,1)*A(1,0);
	if (Ai) {
	    xASSERT(norm(d) > detmin, "Matrix singular");
	    Ai->Set(0,0, A(1,1)/d);
	    Ai->Set(0,1,-A(0,1)/d);
	    Ai->Set(1,0,-A(1,0)/d);
	    Ai->Set(1,1, A(0,0)/d);
	}
	break;
    case 3:
      d = A(0,0) * (A(1,1)*A(2,2) - A(2,1)*A(1,2)) -
	  A(0,1) * (A(1,0)*A(2,2) - A(2,0)*A(1,2)) +
	  A(0,2) * (A(1,0)*A(2,1) - A(2,0)*A(1,1));
      if (Ai) {
	  xASSERT(norm(d) > detmin, "Matrix singular");
	  Ai->Set(0,0, ( A(1,1)*A(2,2) - A(2,1)*A(1,2))/d);
	  Ai->Set(1,0, (-A(1,0)*A(2,2) + A(2,0)*A(1,2))/d);
	  Ai->Set(2,0, ( A(1,0)*A(2,1) - A(2,0)*A(1,1))/d);
	  Ai->Set(0,1, (-A(0,1)*A(2,2) + A(2,1)*A(0,2))/d);
	  Ai->Set(1,1, ( A(0,0)*A(2,2) - A(2,0)*A(0,2))/d);
	  Ai->Set(2,1, (-A(0,0)*A(2,1) + A(2,0)*A(0,1))/d);
	  Ai->Set(0,2, ( A(0,1)*A(1,2) - A(1,1)*A(0,2))/d);
	  Ai->Set(1,2, (-A(0,0)*A(1,2) + A(1,0)*A(0,2))/d);
	  Ai->Set(2,2, ( A(0,0)*A(1,1) - A(1,0)*A(0,1))/d);
	}
	break;
    default: {
	int i, j, k, jj;
        for (d = 0, k = 0; k < n; k++) {
	    TDenseMatrix<MT> tmp(n-1);       // create submatrix
	    for (j = jj = 0; j < n; j++) {
	        if (j == k) continue;
		for (i = 1; i < n; i++)
		    tmp(i-1,jj) = A(i,j);
		jj++;
	    }
	    if (k%2) d -= A(k,0) * det(tmp); // recursion
	    else     d += A(k,0) * det(tmp);
	}
	xASSERT(Ai == 0,
	   "This routine does not compute inverse for matrix greater than 3x3");
	}
	break;
    }
    return d;
}

template<class MT>
MATHLIB TDenseMatrix<MT> inverse (const TDenseMatrix<MT> &A)
{
    int n = A.nRows();
    xASSERT (n == A.nCols(), "Inverse requires square matrix.");
    TDenseMatrix<MT> Ai (n);

    if (n <= 3) {
	det (A, &Ai);
    } else { // use LU factorisation to invert matrix
        // need to create a copy of A since A is constant
        TDenseMatrix<MT> AA(A);
	IVector indx(n);
	double d;
	int i, j;
	LUFactorize (AA, indx, d);
	for (j = 0; j < n; j++) {
	    TVector<MT> col(n);
  	    col[j] = 1.0;
	    LUSolve (AA, indx, col);
	    for (i = 0; i < n; i++) Ai(i,j) = col[i];
	}
    }
    return Ai;
}

// ==========================================================================
// QR decomposition
// Calculates Householder QR decomposition of (n,m) matrix A
// Mirror vectors v_k, S_k=I-v_k v_k^t are stored in A[k..(n-1)][k]
// Factor R is stored in A[j][(j+1)..(m-1)] and in d[j]=r_jj

template<class MT>
MATHLIB int QRFactorize (TDenseMatrix<MT> &A, TVector<MT> &c, TVector<MT> &d)
{
    int i, j, k;
    int n = A.nRows();
    int m = A.nCols();
    MT sum, b, f;

    for (k = 0; k < m; k++) {
        for (sum = 0, i = k; i < n; i++)
	    sum += A(i,k)*A(i,k);
	d[k] = (A(k,k) < 0 ? -sqrt(sum) : sqrt(sum));
	b = sqrt((MT)2.0*d[k]*(A(k,k) + d[k]));
	A(k,k) = (A(k,k) + d[k])/b;
	for (i = k+1; i < n; i++)
	    A(i,k) /= b;
	for (j = k+1; j < m; j++) {
	    for (sum = 0, i = k; i < n; i++)
	        sum += A(i,k)*A(i,j);
	    f = (MT)2.0*sum;
	    for (i = k; i < n; i++)
	        A(i,j) -= f*A(i,k);
	}
    }
    return 0;
}

template<class MT>
MATHLIB void RSolve (const TDenseMatrix<MT> &A, const TVector<MT> &d, TVector<MT> &b)
{
    int i, j;
    int n = A.nRows();
    int m = A.nCols();
    MT sum;
    b[m-1] /= -d[m-1];
    for (i = m-2; i >= 0; i--) {
        for (sum = 0.0, j = i+1; j < m; j++)
	    sum = sum + A(i,j) * b[j];
	b[i] = (b[i]-sum) / -d[i];
    }
}

template<class MT>
MATHLIB void QRSolve (const TDenseMatrix<MT> &A, const TVector<MT> &c,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x)
{
    int n = A.nRows();
    int m = A.nCols();
    int i, k;
    MT sum, f;

    // Calculates y = Q^T b
    x = b;
    for (k = 0; k < m; k++) {
        for (sum = 0, i = k; i < n; i++)
	    sum += A(i,k)*x[i];
	f = sum*2.0;
	for (i = k; i < n; i++)
	    x[i] -= f*A(i,k);
    }

    // Solves Rx = y
    RSolve (A, d, x);
}

template<class MT>
MATHLIB void LUFactorize (TDenseMatrix<MT> &a, IVector &indx, double &d)
{
    // Implementation from NR. pp 46
    // Modifications: matrix and vector interfaces, zero-based indices

    const double TINY = 1.0e-20;
    double big, sum, dum, temp;
    int i, j, k, imax, n = a.nRows();
    if (n != a.nCols())
        xERROR ("LUFactorize requires square matrix.");
    indx.New(n);
    RVector vv(n);
    d = 1.0;  // no row interchanges yet
    for (i = 0; i < n; i++) {
        big = 0.0;
	for (j = 0; j < n; j++)
	    if ((temp = fabs(a(i,j))) > big) big = temp;
	if (big == 0.0)
	    xERROR ("Singular matrix in LUFactorize");
	vv[i] = 1.0/big;
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
	    sum = a(i,j);
	    for (k = 0; k < i; k++)
	        sum -= a(i,k) * a(k,j);
	    a(i,j) = sum;
	}
	big = 0.0;
	for (i = j; i < n; i++) {
	    sum = a(i,j);
	    for (k = 0; k < j; k++)
	        sum -= a(i,k) * a(k,j);
	    a(i,j) = sum;
	    if ((dum = vv[i]*fabs(sum)) >= big) {
	        big = dum;
		imax = i;
	    }
	}
	if (j != imax) {
	    for (k = 0; k < n; k++) {
	        dum = a(imax,k);
		a(imax,k) = a(j,k);
		a(j,k) = dum;
	    }
	    d = -d;
	    vv[imax] = vv[j];
	}
	indx[j] = imax;
	if (a(j,j) == 0.0) a(j,j) = TINY;
	if (j+1 != n) {
	  dum = 1.0 / a(j,j);
	  for (i = j+1; i < n; i++)
	      a(i,j) *= dum;
	}
    }
}

template<class MT>
MATHLIB void LUSolve (const TDenseMatrix<MT> &a, const IVector &indx, TVector<MT> &b)
{
    int i, j, ip, ii = 0;
    int n = a.nRows();
    MT sum;

    for (i = 0; i < n; i++) {
        ip = indx[i];
	sum = b[ip];
	b[ip] = b[i];
	if (ii)
	    for (j = ii-1; j < i; j++)
	        sum -= a(i,j)*b[j];
	else if (sum) ii = i+1;
	b[i] = sum;
    }
    for (i = n-1; i >= 0; i--) {
        sum = b[i];
	for (j = i+1; j < n; j++)
	    sum -= a(i,j) * b[j];
	b[i] = sum / a(i,i);
    }
}

template<class MT>
MATHLIB istream &operator>> (istream &is, TDenseMatrix<MT> &m)
{
    int i, n;
    char c;
    struct vectag {
	TVector<MT> v;
	vectag *next;
    } *vtfirst = 0, *vtlast = 0;

    do {
	is >> c;
    } while (c != '[');
    for (n = 0;; n++) {
	TVector<MT> v;
	is >> v;
	if (!is.good()) break;
	vectag *vt = new vectag;
	vt->v = v;
	vt->next = 0;
	if (vtlast) vtlast->next = vt;
	else vtfirst = vt;
	vtlast = vt;
    }

    m.New (n, vtfirst->v.Dim());
    for (i = 0; i < n; i++) {
	m.SetRow (i, vtfirst->v);
	vectag *vt = vtfirst;
	vtfirst = vtfirst->next;
	delete vt;
    }
    return is;
}

template<class MT>
MATHLIB ostream &operator<< (ostream &os, const TDenseMatrix<MT> &m)
{
    os << '[';
    for (int i = 0; i < m.rows; i++) {
	if (i) os << '\n';

	os << '[';
	for (int j = 0; j < m.cols; j++) {
	    if (j) os << ' ';
	    os << m(i, j);
	}
	os << ']';
    }
    os << ']';
    return os;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TDenseMatrix<double>;
template class MATHLIB TDenseMatrix<float>;
template class MATHLIB TDenseMatrix<toast::complex>;
template class MATHLIB TDenseMatrix<scomplex>;
template class MATHLIB TDenseMatrix<int>;

template MATHLIB istream &operator>> (istream &is, RDenseMatrix &m);
template MATHLIB istream &operator>> (istream &is, CDenseMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const RDenseMatrix &m);
template MATHLIB ostream &operator<< (ostream &os, const CDenseMatrix &m);

template MATHLIB TDenseMatrix<double> cath (const TDenseMatrix<double> &A,
    const TDenseMatrix<double> &B);
template MATHLIB TDenseMatrix<toast::complex> cath (const TDenseMatrix<toast::complex> &A,
    const TDenseMatrix<toast::complex> &B);

template MATHLIB TDenseMatrix<double> catv (const TDenseMatrix<double> &A,
    const TDenseMatrix<double> &B);
template MATHLIB TDenseMatrix<toast::complex> catv (const TDenseMatrix<toast::complex> &A,
    const TDenseMatrix<toast::complex> &B);

#ifndef USE_BLAS_LEVEL3 // otherwise use BLAS interface specialisations
template MATHLIB TSymMatrix<double> ATA (const TDenseMatrix<double> &A);
template MATHLIB TSymMatrix<float> ATA (const TDenseMatrix<float> &A);
template MATHLIB TSymMatrix<toast::complex> ATA (const TDenseMatrix<toast::complex> &A);
#endif

#ifndef USE_BLAS_LEVEL3 // otherwise use BLAS interface specialisations
template MATHLIB TSymMatrix<double> AAT (const TDenseMatrix<double> &A);
template MATHLIB TSymMatrix<float> AAT (const TDenseMatrix<float> &A);
template MATHLIB TSymMatrix<toast::complex> AAT (const TDenseMatrix<toast::complex> &A);
#endif

template MATHLIB TDenseMatrix<double> kron (const TDenseMatrix<double> &A,
    const TDenseMatrix<double> &B);
template MATHLIB TDenseMatrix<toast::complex> kron (const TDenseMatrix<toast::complex> &A,
    const TDenseMatrix<toast::complex> &B);
template MATHLIB TDenseMatrix<int> kron (const TDenseMatrix<int> &A,
    const TDenseMatrix<int> &B);

template MATHLIB double det (const TDenseMatrix<double> &A, TDenseMatrix<double> *Ai);
template MATHLIB TDenseMatrix<double> inverse (const TDenseMatrix<double> &A);
template MATHLIB TDenseMatrix<double> transpose (const TDenseMatrix<double> &A);
template MATHLIB TDenseMatrix<toast::complex> transpose (const TDenseMatrix<toast::complex> &A);

template MATHLIB int QRFactorize (TDenseMatrix<double> &A, TVector<double> &c,
    TVector<double> &d);
template MATHLIB int QRFactorize (TDenseMatrix<toast::complex> &A, TVector<toast::complex> &c,
    TVector<toast::complex> &d);

template MATHLIB void RSolve (const TDenseMatrix<double> &A, const TVector<double> &d,
    TVector<double> &b);
template MATHLIB void RSolve (const TDenseMatrix<toast::complex> &A, const TVector<toast::complex> &d,
    TVector<toast::complex> &b);

template MATHLIB void QRSolve (const TDenseMatrix<double> &A, const TVector<double> &c,
    const TVector<double> &d, const TVector<double> &b, TVector<double> &x);
template MATHLIB void QRSolve (const TDenseMatrix<toast::complex> &A, const TVector<toast::complex> &c,
    const TVector<toast::complex> &d, const TVector<toast::complex> &b, TVector<toast::complex> &x);


template MATHLIB void LUFactorize (TDenseMatrix<double> &a, IVector &indx, double &d);
template MATHLIB void LUSolve (const TDenseMatrix<double> &a, const IVector &indx,
    TVector<double> &b);

#endif // NEED_EXPLICIT_INSTANTIATION
