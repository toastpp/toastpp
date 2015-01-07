// ==========================================================================
// General dense matrix class
// Distributed (MPI) version
// Each MPI process operates on a block of matrix rows (r0 <= r < r1)
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "dnsmatrix_mpi.h"

using namespace std;
using namespace toast;

template<class MT>
TDenseMatrixMPI<MT>::TDenseMatrixMPI (): TMatrix<MT> ()
{
    r0 = 0;
    r1 = 0;
    nr = 0;
    val = 0;
}

template<class MT>
TDenseMatrixMPI<MT>::TDenseMatrixMPI (int r, int c): TMatrix<MT> (r, c)
{
    Alloc (r, c); // use default process subdivision
}

template<class MT>
TDenseMatrixMPI<MT>::~TDenseMatrixMPI<MT> ()
{
    Unlink();
}

template<class MT>
void TDenseMatrixMPI<MT>::New (int r, int c)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (r != this->rows || c != this->cols) {
	TMatrix<MT>::New (r,c);
	Unlink();
	Alloc (r, c);
    } else {
	if (nr)
	    memset (val, 0, nr*this->cols*sizeof(MT));
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

template<class MT>
MT TDenseMatrixMPI<MT>::Get (int r, int c) const
{
    // not implemented yet
    return (MT)0;
}

template<class MT>
MT &TDenseMatrixMPI<MT>::operator() (int r, int c)
{
    dASSERT (r >= 0 && r < this->rows && c >= 0 && c < this->cols, Index out of range);
    dASSERT (r >= r0 && r < r1, Row index not within process support);
    return val[(r-r0)*this->cols + c];
}

template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::Row (int r) const
{
    MPI_Datatype mpitp = MPIType();
    int nc = this->cols;
    TVector<MT> rw(nc);
    if (r >= r0 && r < r1) {
	memcpy(rw.data_buffer(), val+(r-r0)*nc, nc*sizeof(MT));
	MPI_Bcast (rw.data_buffer(), nc, mpitp, rnk, MPI_COMM_WORLD);
	// ERROR: MPI_Bcast needs to be called by ALL processes.
	// Therefore each process needs to know which rank is sending the data
    }
    return rw;
}

template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::Col (int c) const
{
    // not implemented yet
    return TVector<MT>();
}

template<class MT>
int TDenseMatrixMPI<MT>::SparseRow (int r, int *ci, MT *rv) const
{
    // not implemented yet
    return 0;
}

template<class MT>
void TDenseMatrixMPI<MT>::SetRow (int r, const TVector<MT> &rval)
{
    if (r >= this->r0 && r < this->r1) {
	dASSERT(rval.Dim() == this->cols, Argument 2: wrong size);
	memcpy (val + ((r-r0)*this->cols), rval.data_buffer(),
		this->cols*sizeof(MT));
    }
}

template<class MT>
void TDenseMatrixMPI<MT>::ColScale (const TVector<MT> &scale)
{
    dASSERT (scale.Dim() == this->cols, Argument 1: wrong size);
    int r, c;
    for (r = r0; r < r1; r++)
	for (c = 0; c < this->cols; c++)
	    val[(r-r0)*this->cols+c] *= scale[c];
}

template<class MT>
void TDenseMatrixMPI<MT>::RowScale (const TVector<MT> &scale)
{
    int r, c;
    dASSERT (scale.Dim() == this->rows, Argument 1: wrong size);
    for (r = r0; r < r1; r++)
	for (c = 0; c < this->cols; c++)
	    val[(r-r0)*this->cols+c] *= scale[r];
}

template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::RowSum () const
{
    MPI_Datatype mpitp = MPIType();
    MT *buf = new MT[sze*nr_max];

    if (nr) {
	int r, c;
	int nc = this->cols;
	MT sum, *idx = val;
	for (r = 0; r < nr; r++) {
	    sum = (MT)0;
	    for (c = 0; c < nc; c++)
		sum += *idx++;
	    buf[r+r0] = sum;
	}
    }

    MPI_Allgather (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buf,
		   nr_max, mpitp, MPI_COMM_WORLD);
    TVector<MT>res (this->rows);
    memcpy (res.data_buffer(), buf, this->rows*sizeof(MT));
    delete []buf;
    return res;
}

template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::RowSumSq () const
{
    MPI_Datatype mpitp = MPIType();
    MT *buf = new MT[sze*nr_max];

    if (nr) {
	int r, c;
	int nc = this->cols;
	MT sumsq, *idx = val;
	for (r = 0; r < nr; r++) {
	    sumsq = (MT)0;
	    for (c = 0; c < nc; c++) {
		sumsq += (*idx)*(*idx);
		idx++;
	    }
	    buf[r+r0] = sumsq;
	}
    }

    MPI_Allgather (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buf,
		   nr_max, mpitp, MPI_COMM_WORLD);
    TVector<MT> res (this->rows);
    memcpy (res.data_buffer(), buf, this->rows*sizeof(MT));
    delete []buf;
    return res;
}

template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::ColSum () const
{
    MPI_Datatype mpitp = MPIType();
    int nc = this->cols;
    MT *buf = new MT[nc];
    memset (buf, 0, nc*sizeof(MT));
    
    if (nr) {
	int r, c;
	MT *idx = val;
	for (r = 0; r < nr; r++)
	    for (c = 0; c < nc; c++)
		buf[c] += *idx++;
    }

    TVector<MT> res(nc);
    MPI_Allreduce (buf, res.data_buffer(), nc, mpitp, MPI_SUM, MPI_COMM_WORLD);
    delete []buf;
    return res;
}


template<class MT>
TVector<MT> TDenseMatrixMPI<MT>::ColSumSq () const
{
    MPI_Datatype mpitp = MPIType();
    int nc = this->cols;
    MT *buf = new MT[nc];
    memset (buf, 0, nc*sizeof(MT));
    
    if (nr) {
	int r, c;
	MT *idx = val;
	for (r = 0; r < nr; r++)
	    for (c = 0; c < nc; c++) {
		buf[c] += (*idx)*(*idx);
		idx++;
	    }
    }

    TVector<MT> res(nc);
    MPI_Allreduce (buf, res.data_buffer(), nc, mpitp, MPI_SUM, MPI_COMM_WORLD);
    delete []buf;
    return res;
}


template<class MT>
void TDenseMatrixMPI<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == x.Dim(), Argument 1 wrong size);
    dASSERT(this->rows == b.Dim(), Argument 2 wrong size);

    MPI_Datatype mpitp = MPIType();
    MT *bbuf = new MT[sze*nr_max];
    // pad the buffer so we can use MPI_Allgather with equal fragment sizes

    // Individual solution blocks
    if (nr) {
	int r, c;
	int nc = this->cols;
	MT sum, *idx = val;
	for (r = 0; r < nr; r++) {
	    sum = (MT)0;
	    for (c = 0; c < nc; c++)
		sum += *idx++ * x[c];
	    bbuf[r+r0] = sum;
	}
    }

    // Now collect the results and re-distribute to all processes
    MPI_Allgather (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bbuf,
		   nr_max, mpitp, MPI_COMM_WORLD);
    memcpy(b.data_buffer(), bbuf, b.Dim()*sizeof(MT));

    delete []bbuf;
}

template<class MT>
void TDenseMatrixMPI<MT>::ATx (const TVector<MT> &x, TVector<MT> &b) const
{
    dASSERT(this->cols == b.Dim(), Argument 1 wrong size);
    dASSERT(this->rows == x.Dim(), Argument 2 wrong size);

    MPI_Datatype mpitp = MPIType();
    int nc = this->cols;
    MT *bbuf = new MT[nc];
    memset (bbuf, 0, nc*sizeof(MT));

    if (nr) {
	int r, c;
	MT xr, *idx = val;
	for (r = 0; r < nr; r++) {
	    xr = x[r+r0];
	    for (c = 0; c < nc; c++)
		bbuf[c] += *idx++ * xr;
	}
    }

    // Sum up contributions from all processes and re-distribute
    MPI_Allreduce (bbuf, b.data_buffer(), nc, mpitp, MPI_SUM, MPI_COMM_WORLD);

    delete []bbuf;
}

template<class MT>
void TDenseMatrixMPI<MT>::Unlink ()
{
    if (val) delete []val;
}


template<class MT>
void TDenseMatrixMPI<MT>::Alloc (int r, int c)
{
    // Allocate default sub-range for this process, given matrix size r x c.

    // Find range of rows to be managed by this process
    MPI_Comm_rank (MPI_COMM_WORLD, &rnk);
    MPI_Comm_size (MPI_COMM_WORLD, &sze);

    nr_max = max (1, (r+sze-1)/sze);   // rows per processor
    r0 = min (r, rnk*nr_max);          // base index
    r1 = min (r, (rnk+1)*nr_max);      // top index+1
    nr = r1-r0;                               // actual rows for this process

    if (nr) {
	val = new MT[nr*c];
	memset (val, 0, nr*c*sizeof(MT));
    } else {
	val = 0;
    }
}

template<class MT>
MPI_Datatype TDenseMatrixMPI<MT>::MPIType() const
{
    return MPI_INT; // default
}

template<>
MPI_Datatype TDenseMatrixMPI<double>::MPIType() const
{
    return MPI_DOUBLE;
}

template<>
MPI_Datatype TDenseMatrixMPI<complex>::MPIType() const
{
    return MPI_DOUBLE_COMPLEX;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TDenseMatrixMPI<double>;
template class MATHLIB TDenseMatrixMPI<float>;
template class MATHLIB TDenseMatrixMPI<complex>;
template class MATHLIB TDenseMatrixMPI<scomplex>;
template class MATHLIB TDenseMatrixMPI<int>;

#endif // NEED_EXPLICIT_INSTANTIATION
