// ==========================================================================
// Module mathlib
// File crmatrix_mpi.h
// Declaration of template class TCompRowMatrixMPI
// Distributed (MPI) version of TCompRowMatrix (template sparse matrix in
// compressed row format)
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __CRMATRIXMPI_CC

//#include "toast_mpi.h"
#include "crmatrix_mpi.h"

using namespace std;
using namespace toast;

// ==========================================================================
// Member definitions

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI ()
  : TGenericSparseMatrix<MT> ()
{
    MPIinit();
    rowptr = new int[1];
    rowptr[0] = 0;
}

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI (int rows, int cols)
  : TGenericSparseMatrix<MT> (rows, cols)
{
    MPIinit();
    rowptr = new int[my_nr+1];
    for (int i = 0; i <= my_nr; i++)
        rowptr[i] = 0;
}

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI (int nrows, int ncols,
    const int *_rowptr, const int *_colidx,
    int proc_nrows, const int *proc_rows)
  : TGenericSparseMatrix<MT> (nrows, ncols)
{
    MPIinit();

    int i, j, k, row;

    my_nr = proc_nrows;
    my_r  = new int[proc_nrows];
    for (i = 0; i < proc_nrows; i++)
        my_r[i] = proc_rows[i];

    // Eliminate all rows from the sparsity structure that I'm not responsible
    // for
    int *nentry = new int[nrows];
    for (i = 0; i < nrows; i++) nentry[i] = 0;
    for (i = this->nval = 0; i < my_nr; i++) {
        row = my_r[i];
	this->nval += (nentry[row] = _rowptr[row+1]-_rowptr[row]);
    }
    rowptr = new int[nrows+1];
    colidx = new int[this->nval];
    this->val = new MT[this->nval];
    memset (this->val, 0, this->nval*sizeof(MT));

    rowptr[0] = 0;
    for (i = 0; i < nrows; i++)
        rowptr[i+1] = rowptr[i] + nentry[i];

    for (i = 0; i < my_nr; i++) {
        row = my_r[i];
	k = rowptr[row];
	for (j = _rowptr[row]; j < _rowptr[row+1]; j++)
	    colidx[k++] = _colidx[j];
    }

    delete []nentry;
}

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI (int rows, int cols,
    const int *_rowptr, const int *_colidx, const MT *_data)
  : TGenericSparseMatrix<MT> (rows, cols)
{
    MPIinit();
    this->nval = my_nr = 0;
    Initialise (_rowptr, _colidx, _data);
}

template<class MT>
TCompRowMatrixMPI<MT>::~TCompRowMatrixMPI ()
{
    Unlink();
    delete []rowptr;
    if (my_nr) delete []my_r;
    //if (mpi_r0) delete []mpi_r0;
    //if (mpi_nr) delete []mpi_nr;
}

template<class MT>
void TCompRowMatrixMPI<MT>::Initialise (const int *_rowptr, const int *_colidx,
    const MT *_data)
{
    int i, j;

    int nval_tot = _rowptr[this->rows]; // total number of nonzeros
    int nval_tgt = nval_tot / sze;  // target number of nonzeros per process

    int r;
    mpi_r0[0] = 0;
    mpi_r0[sze] = this->rows;

    for (i = 1; i < sze; i++) {
        int nv_tgt = i*nval_tgt;  // target for upper bound
	for (r = mpi_r0[i-1]; r < this->rows && _rowptr[r] < nv_tgt; r++);
	if (r && (_rowptr[r]-nv_tgt > nv_tgt-_rowptr[r-1])) r--;
	mpi_r0[i] = r;
    }

    for (i = 0; i < sze; i++)
        mpi_nr[i] = mpi_r0[i+1]-mpi_r0[i];

    r0 = mpi_r0[rnk];
    r1 = mpi_r0[rnk+1];

    if (my_nr != r1-r0) {
        if (my_nr) delete []rowptr;
	my_nr = r1-r0;
	if (my_nr) rowptr = new int[my_nr+1];
    }

    int ofs = _rowptr[r0];

    for (i = 0; i <= my_nr; i++)
        rowptr[i] = _rowptr[i+r0] - ofs;

    if (this->nval != rowptr[my_nr]) {
        if (this->nval) {
  	    delete []colidx;
	    delete []this->val;
	}
	this->nval = rowptr[my_nr];
	if (this->nval) {
	    colidx = new int[this->nval];
	    this->val = new MT[this->nval];
	}
    }

    for (i = 0; i < this->nval; i++)
        colidx[i] = _colidx[i+ofs];
    for (i = 0; i < this->nval; i++)
        this->val[i] = (_data ? _data[i+ofs] : (MT)0);

    cout << "P-" << rnk << ": r0=" << r0 << ", r1=" << r1 << ", nval=" << this->nval << endl;
}

template<class MT>
void TCompRowMatrixMPI<MT>::Zero ()
{
    for (int i = 0; i < this->nval; i++)
        this->val[i] = (MT)0;
}

template<class MT>
MT TCompRowMatrixMPI<MT>::Get (int r, int c) const
{
    xERROR("Not implemented yet");
    return (MT)0;
}

template<class MT>
TVector<MT> TCompRowMatrixMPI<MT>::Row (int r) const
{
    int rk, nc = this->cols;
    TVector<MT> rw(nc);

    for (rk = 0; rk < sze; rk++)
        if (r < mpi_r0[rk+1]) break;

    if (rk == rnk) { // it's our row
        int i, i1, i2;
        r -= r0;
	i1 = rowptr[r], i2 = rowptr[r+1];
	for (i = i1; i < i2; i++)
	    rw[colidx[i]] = this->val[i];
    }
    MPI_Bcast (rw.data_buffer(), nc, mpitp, rk, MPI_COMM_WORLD);
    return rw;
}


template<class MT>
int TCompRowMatrixMPI<MT>::SparseRow (int r, int *colidx, MT *val) const
{
    xERROR("Not implemented yet");
    return 0;
}

template<class MT>
TVector<MT> TCompRowMatrixMPI<MT>::Col (int c) const
{
    xERROR("Not implemented yet");
    TVector<MT> col(this->rows);
    return col;
}

template<class MT>
bool TCompRowMatrixMPI<MT>::Exists (int r, int c) const
{
    dASSERT(r < this->rows, "Row index out of range");
    dASSERT(c < this->cols, "Col index out of range");

    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return true;
    return false;
}

template<class MT>
MT &TCompRowMatrixMPI<MT>::operator() (int r, int c)
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
int TCompRowMatrixMPI<MT>::Get_index (int r, int c) const
{
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) return rp;
    return -1;
}

template<class MT>
MT TCompRowMatrixMPI<MT>::GetNext (int &r, int &c) const
{
    xERROR("Not implemented");
    return (MT)0;
}

template<class MT>
void TCompRowMatrixMPI<MT>::RowScale (const TVector<MT> &scale)
{
    xERROR("Not implemented yet");
}

template<class MT>
void TCompRowMatrixMPI<MT>::ColScale (const TVector<MT> &scale)
{
    xERROR("Not implemented yet");
}

template<class MT>
void TCompRowMatrixMPI<MT>::Unlink ()
{
    if (this->nval) {
      //delete []colidx;
      //delete []val;
      //memset (rowptr, 0, (rows+1)*sizeof(int));
	this->nval = 0;
    }
}

template<class MT>
void TCompRowMatrixMPI<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == x.Dim(), Argument 1 wrong size);
    if (b.Dim() != this->rows) b.New(this->rows);

    b.Clear();
    if (!this->rows) return; // nothing to do

    int i, j, j2, r;
    MT sum;

    for (i = 0; i < my_nr; i++) {
        r = my_r[i];
	j2 = rowptr[r+1];
	for (sum = (MT)0, j = rowptr[r]; j < j2; j++)
	    sum += this->val[j] * x[colidx[j]];
	b[r] = sum;
    }

    // Now collect the results and re-distribute to all processes
    MPI_Allreduce (MPI_IN_PLACE, b.data_buffer(), b.Dim(),
		   mpitp, MPI_SUM, MPI_COMM_WORLD);
}

template<class MT>
void TCompRowMatrixMPI<MT>::Ax (const TVector<MT> &x, TVector<MT> &b,
    int i1, int i2) const
{
    xERROR("Not implemented");
}

template<class MT>
void TCompRowMatrixMPI<MT>::Add_proc (int r, int c, MT v)
{
    dASSERT (r >= r0 && r < r1, Row index out of process block range);
    
    r -= r0;
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) {
	    this->val[rp] += v;
	    return;
	}
    xERROR("Specified element not found");
}

template<class MT>
void TCompRowMatrixMPI<MT>::MPIinit ()
{
    rnk = TMPI::Rank();
    sze = TMPI::Size();
    mpitp = TMPI::MPIType<MT>();
    mpi_r0 = new int[sze+1];
    mpi_nr = new int[sze];
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TCompRowMatrixMPI<double>;
template class MATHLIB TCompRowMatrixMPI<float>;
template class MATHLIB TCompRowMatrixMPI<complex>;
template class MATHLIB TCompRowMatrixMPI<scomplex>;
template class MATHLIB TCompRowMatrixMPI<int>;

#endif // NEED_EXPLICIT_INSTANTIATION
