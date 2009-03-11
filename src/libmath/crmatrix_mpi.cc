// ==========================================================================
// Module mathlib
// File crmatrix_mpi.h
// Declaration of template class TCompRowMatrixMPI
// Distributed (MPI) version of TCompRowMatrix (template sparse matrix in
// compressed row format)
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __CRMATRIXMPI_CC

#include "crmatrix_mpi.h"

using namespace std;
using namespace toast;

// ==========================================================================
// Member definitions

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI ()
  : TMatrix<MT> ()
{
    MPIinit();
    rowptr = new int[1];
    rowptr[0] = 0;
}

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI (int rows, int cols)
  : TMatrix<MT> (rows, cols)
{
    MPIinit();
    rowptr = new int[nr+1];
    for (int i = 0; i <= nr; i++)
        rowptr[i] = 0;
}

template<class MT>
TCompRowMatrixMPI<MT>::TCompRowMatrixMPI (int rows, int cols,
    const int *_rowptr, const int *_colidx, const MT *_data)
  : TMatrix<MT> (rows, cols)
{
    MPIinit();
    nval = nr = 0;
    Initialise (_rowptr, _colidx, _data);
}

template<class MT>
TCompRowMatrixMPI<MT>::~TCompRowMatrixMPI ()
{
    Unlink();
    delete []rowptr;
    if (mpi_r0) delete []mpi_r0;
    if (mpi_nr) delete []mpi_nr;
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

    if (nr != r1-r0) {
        if (nr) delete []rowptr;
	nr = r1-r0;
	if (nr) rowptr = new int[nr+1];
    }

    int ofs = _rowptr[r0];

    for (i = 0; i <= nr; i++)
        rowptr[i] = _rowptr[i+r0] - ofs;

    if (nval != rowptr[nr]) {
        if (nval) {
  	    delete []colidx;
	    delete []data;
	}
	nval = rowptr[nr];
	if (nval) {
	    colidx = new int[nval];
	    data = new MT[nval];
	}
    }

    for (i = 0; i < nval; i++)
        colidx[i] = _colidx[i+ofs];
    for (i = 0; i < nval; i++)
        data[i] = (_data ? _data[i+ofs] : (MT)0);

    cout << "P-" << rnk << ": r0=" << r0 << ", r1=" << r1 << ", nval=" << nval << endl;
}

template<class MT>
void TCompRowMatrixMPI<MT>::Zero ()
{
    for (int i = 0; i < nval; i++)
        data[i] = (MT)0;
}

template<class MT>
MT TCompRowMatrixMPI<MT>::Get (int r, int c) const
{
    xERROR(Not implemented yet);
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
	    rw[colidx[i]] = data[i];
    }
    MPI_Bcast (rw.data_buffer(), nc, mpitp, rk, MPI_COMM_WORLD);
    return rw;
}


template<class MT>
int TCompRowMatrixMPI<MT>::SparseRow (int r, int *colidx, MT *val) const
{
    xERROR(Not implemented yet);
    return 0;
}

template<class MT>
TVector<MT> TCompRowMatrixMPI<MT>::Col (int c) const
{
    xERROR(Not implemented yet);
    TVector<MT> col(this->rows);
    return col;
}

template<class MT>
void TCompRowMatrixMPI<MT>::RowScale (const TVector<MT> &scale)
{
    xERROR(Not implemented yet);
}

template<class MT>
void TCompRowMatrixMPI<MT>::ColScale (const TVector<MT> &scale)
{
    xERROR(Not implemented yet);
}

template<class MT>
void TCompRowMatrixMPI<MT>::Unlink ()
{
    if (nval) {
        delete []colidx;
	delete []data;
	memset (rowptr, 0, (nr+1)*sizeof(int));
	nval = 0;
    }
}

template<class MT>
void TCompRowMatrixMPI<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == x.Dim(), Argument 1 wrong size);
    if (b.Dim() != this->rows) b.New(this->rows);

    if (!this->rows) return; // nothing to do

    MT *bbuf = b.data_buffer();

    // Individual solution blocks
    int r, i, i2;
    MT sum;
    for (r = i = 0; r < nr; r++) {
        i2 = rowptr[r+1];
	for (sum = (MT)0; i < i2; i++)
	    sum += data[i] * x[colidx[i]];
	bbuf[r+r0] = sum;
    }

    // Now collect the results and re-distribute to all processes
    MPI_Allgatherv (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bbuf,
    		    mpi_nr, mpi_r0, mpitp, MPI_COMM_WORLD);
}

template<class MT>
void TCompRowMatrixMPI<MT>::Add_proc (int r, int c, MT val)
{
    dASSERT (r >= r0 && r < r1, Row index out of process block range);
    
    r -= r0;
    for (int rp = rowptr[r]; rp < rowptr[r+1]; rp++)
        if (colidx[rp] == c) {
	    data[rp] += val;
	    return;
	}
    xERROR(Specified element not found);
}

template<class MT>
void TCompRowMatrixMPI<MT>::MPIinit ()
{
    rnk = TMPI<MT>::Rank();
    sze = TMPI<MT>::Size();
    mpitp = TMPI<MT>::MPIType();
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
