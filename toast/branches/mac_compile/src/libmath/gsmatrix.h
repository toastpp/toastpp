// -*-C++-*-
// ==========================================================================
// Module mathlib
// File gsmatrix.h
// Declaration of template class TGenericSparseMatrix
// Base class for sparse matrix types
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RGenericSparseMatrix = TGenericSparseMatrix<double>   ('real')
//	FGenericSparseMatrix = TGenericSparseMatrix<float>    ('float')
//	CGenericSparseMatrix = TGenericSparseMatrix<complex>  ('complex')
//	IGenericSparseMatrix = TGenericSparseMatrix<int>      ('integer')
// ==========================================================================

#ifndef __GSMATRIX_H
#define __GSMATRIX_H

#ifdef TOAST_PARALLEL
//#define CG_PARALLEL
#endif

#ifdef DBG_TIMING
extern double cgtime;
#endif

#include "matrix.h"

const int BUFFER_CHUNK_SIZE = 256;
// allocation chunk size for dynamically growing matrices

typedef enum {
    ITMETHOD_CG,
    ITMETHOD_BICG,
    ITMETHOD_BICGSTAB,
    ITMETHOD_GMRES,
    ITMETHOD_GAUSSSEIDEL
} IterativeMethod;

// ==========================================================================
// Nonmember declarations

extern MATHLIB IterativeMethod itmethod_general;
extern MATHLIB IterativeMethod itmethod_complex;
extern MATHLIB int IterCount;

template<class MT>class TGenericSparseMatrix;
template<class MT>class TPreconditioner;

template<class MT>
int QRFactorize (TGenericSparseMatrix<MT> &A, TVector<MT> &c, TVector<MT> &d);

template<class MT>
void RSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &d,
    TVector<MT> &b);

template<class MT>
void QRSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &c,
    const TVector<MT> &d, const TVector<MT> &b, TVector<MT> &x);


template<class MT>
int IterativeSolve (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
void IterativeSolve (const TGenericSparseMatrix<MT> &A,
    const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit = 0,
    TPreconditioner<MT> *precon = 0, IterativeSolverResult *res = 0);

template<class MT>
int CG (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
int BiCG (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon = 0,
    int maxit = 0);

template<class MT>
int ComplexBiCGSolve (const TGenericSparseMatrix<MT> &Are,
    const TGenericSparseMatrix<MT> &Aim, const TVector<MT> &bre,
    const TVector<MT> &bim, TVector<MT> &xre, TVector<MT> &xim,
    double &tol, int maxit = 0);

template<class MT>
int GaussSeidel (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, int maxit = 0);

// ==========================================================================
// class TGenericSparseMatrix
// ==========================================================================
/**
 * \brief Virtual base class for sparse matrix types
 *
 * The following template types are instantiated by default:
 * - TGenericSparseMatrix<double> (RGenericSparseMatrix)
 * - TGenericSparseMatrix<float> (FGenericSparseMatrix)
 * - TGenericSparseMatrix<complex> (CGenericSparseMatrix)
 * - TGenericSparseMatrix<scomplex> (SCGenericSparseMatrix)
 * - TGenericSparseMatrix<int> (IGenericSparseMatrix)
 */
template<class MT> class TGenericSparseMatrix: public TMatrix<MT> {

public:
    /**
     * \brief Create a sparse matrix of size 0 x 0.
     */
    TGenericSparseMatrix ();

    /**
     * \brief Create a sparse matrix of logical size rows x cols.
     * \param rows number of rows (>= 0)
     * \param cols number of columns (>= 0)
     * \param nv number of nonzero entries (0 <= nv <= rows.cols)
     * \param data data array of size >= nv
     * \param cmode if set to SHALLOW_COPY, the provided data array is
     *   used directly instead of copied to a local buffer. The buffer
     *   must remain valid during the lifetime of the matrix.
     * \note If nv>0 and data=0, all allocated elements are initialised to
     *   zero.
     */
    TGenericSparseMatrix (int rows, int cols, int nv=0, const MT *data=0);

    TGenericSparseMatrix (int rows, int cols, int nv, MT *data, CopyMode cmode=DEEP_COPY);

    /**
     * \brief Constructs a matrix as a copy of 'm'
     */
    TGenericSparseMatrix (const TGenericSparseMatrix<MT> &m);

    /**
     * \brief Destructor
     */
    virtual ~TGenericSparseMatrix ();

    /**
     * \brief Reset the matrix dimensions.
     * \param nrows number of matrix rows
     * \param ncols number of matrix columns
     * \note This method unlinks the matrix from its current data block by
     *   calling Unlink(), then resets its logical dimensions by calling
     *   TMatrix::New().
     * \note To allocate a new data block for the matrix after New(), use
     *    Initialise().
     * \sa Unlink, Initialise, TMatrix::New
     */
    virtual void New (int nrows, int ncols);

    virtual MT Get (int r, int c) const = 0;
    // derived classes (also sparse versions) return the value of matrix
    // element (r,c). (Not reference because element might not be stored
    // for some matrix types).

    /**
     * \brief Checks allocation of a matrix element
     * \param r row index (>= 0)
     * \param c column index (>= 0)
     * \return \e true if memory space is allocated for the element, \e false
     *   if not.
     */
    virtual bool Exists (int r, int c) const = 0;

    virtual void Unlink ();
    // deallocate data block

    void Initialise (int nv, const MT *data);
    void Initialise (int nv, MT *data, CopyMode cmode);
    // reallocates data vector of length nv and initialises it with 'data',
    // if given, or zero otherwise

    void Zero ();
    // Reset data to zero, but preserve sparse structure (no deallocation,
    // keep index lists)

    inline int nVal() const { return nval; }

    inline MT &Val (int i)
    { dASSERT(i >= 0 && i < nval, "Index out of range");
      return val[i];
    }
    // return i-th value of data vector

    inline MT *ValPtr () { return val; }
    inline const MT *ValPtr () const { return val; }

    virtual MT &operator() (int r, int c) = 0;
    // element access

    virtual void Add (int r, int c, const MT &val)
    { (*this)(r,c) += val; }

    virtual int Get_index (int r, int c) const = 0;
    // returns offset into data array of entry for row r and column c
    // returns -1 if entry does not exist

    virtual MT GetNext (int &r, int &c) const = 0;
    // Iterator; Returns the next nonzero element, together with its row and
    // column index. In the first call set r < 0 to retrieve first nonzero.
    // In the following calls set r and c to the result of the previous call
    // r < 0 on return indicates end of matrix.
    // Note that elements are not necessarily returned in a particular order

    TGenericSparseMatrix &operator*= (const MT &sc);
    // Multiply *this with scalar sc and return the result

    inline TVector<MT> operator* (const TVector<MT> &x) const
    { TVector<MT> b(this->rows); Ax(x,b); return b; }
    // Multiplies *this with x and returns the result

    virtual void Ax (const TVector<MT> &x, TVector<MT> &b) const = 0;
    virtual void Ax (const TVector<MT> &x, TVector<MT> &b, int i1, int i2)
        const = 0;
    // Returns result of (*this) * x in b
    // The second version processes only rows i1 <= i < i2 of A

    virtual void ATx (const TVector<MT> &x, TVector<MT> &b) const = 0;
    // Returns result of transpose(*this) * x in b

    inline TVector<MT> ATx (const TVector<MT> &x) const
    { return TMatrix<MT>::ATx (x); }

    inline double FillFraction () const
    { return (nval ?
      (double)nval/(double)(this->rows * this->cols) : 0.0); }
    // Returns fractional matrix allocation size

    virtual void Display (std::ostream &os) const;
    // pretty-print matrix in full format. Only advisable for small matrices

    virtual void PrintFillinGraph (const char *fname, int maxdim = 600,
        bool binary = true, bool antialias = true);
    // Generates a binary or ASCII PGM bitmap file 'fname' with the image of
    // the fill-in graph of the matrix (NOT nozeros!)
    // If the matrix dimensions are > maxdim and antialiasing is on
    // then pixel values are interpolated

    /**
     * \brief Write sparse matrix to ASCII output stream.
     *
     * This writes the nonzero elements of the matrix in a row-column-value
     * format
     * \code
     * <r_i> <c_i> <val_i>
     * \endcode
     * containing integer row and column index, and one (double) or two
     * (complex) floating point values, all white-space separated.
     * \note The row and column indices are 1-based.
     * \note This format can be read into MATLAB by using the 'load' command
     * and converted into a MATLAB sparse matrix by using the 'spconvert'
     * command.
     */
    virtual void ExportRCV (std::ostream &os) {}

    static void GlobalSelectIterativeSolver (IterativeMethod method);
    static void GlobalSelectIterativeSolver_complex (IterativeMethod method);
    static IterativeMethod GetGlobalIterativeSolver ();
    static IterativeMethod GetGlobalIterativeSolver_complex ();

    friend int QRFactorize<> (TGenericSparseMatrix<MT> &A, TVector<MT> &c,
        TVector<MT> &d);
    // QR decomposition. Return value != 0 indicates singularity

    friend void RSolve<> (const TGenericSparseMatrix<MT> &A,
	const TVector<MT> &d, TVector<MT> &b);
    // Solves set of linear equations Rx=b where R is upper triangular matrix
    // stored in A and d. On return b is overwritten by the result x.

    friend void QRSolve<> (const TGenericSparseMatrix<MT> &A,
        const TVector<MT> &c, const TVector<MT> &d, const TVector<MT> &b,
        TVector<MT> &x);
    // Solves set of linear equations Ax=b, where A, c and d are the result of
    // a preceeding call to QRFactorize

    friend int IterativeSolve<> (const TGenericSparseMatrix<MT> &A,
        const TVector<MT> &b, TVector<MT> &x, double &tol,
	TPreconditioner<MT> *precon, int maxit);
    // Solves Ax=b using the previously selected iterative solver
    // x contains initial guess on start, solution on exit
    // tol is tolerance limit on start, and final error on exit
    // maxit is iteration limit (0 = no limit)
    // return value is iteration count

    friend void IterativeSolve<> (const TGenericSparseMatrix<MT> &A,
        const TVector<MT> *b, TVector<MT> *x, int nrhs, double tol, int maxit,
        TPreconditioner<MT> *precon, IterativeSolverResult *res);

    friend int CG<> (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
        TVector<MT> &x, double &tol, TPreconditioner<MT> *precon,
        int maxit);
    // preconditioned conjugate gradient method.

    friend int BiCG<> (const TGenericSparseMatrix<MT> &A, const TVector<MT> &b,
        TVector<MT> &x, double &tol, TPreconditioner<MT> *precon,
        int maxit);
    // biconjugate gradient method.

    friend int ComplexBiCGSolve<> (const TGenericSparseMatrix<MT> &Are,
        const TGenericSparseMatrix<MT> &Aim, const TVector<MT> &bre,
	const TVector<MT> &bim, TVector<MT> &xre, TVector<MT> &xim,
	double &tol, int maxit);
    // Solves Ax = b for complex case using biconjugate gradient solver
    // obsolete

    friend int GaussSeidel<> (const TGenericSparseMatrix<MT> &A,
        const TVector<MT> &b, TVector<MT> &x, double &tol, int maxit);
    // Solves Ax = b using Gauss-Seidel iteration

protected:
    void Append (MT v = 0);
    // Adds an additional entry to the end of the data array, reallocating
    // buffers as required. Used by derived classes which allow dynamic growth

    MT *val;   // data block
    int nbuf;  // data block allocation size
    int nval;  // data block entry count (<= nbuf)

#ifdef CG_PARALLEL
    static void cg_loop1(void*,int,int);
    static void cg_loop2(void*,int,int);
    static void cg_loop3(void*,int,int);
#endif

//private:
//    static IterativeMethod itmethod_general;
//    static IterativeMethod itmethod_complex;
};

// ==========================================================================
// typedefs for specific instances of `TGenericSparseMatrix'

typedef TGenericSparseMatrix<double>	RGenericSparseMatrix;	// 'real'
typedef TGenericSparseMatrix<float>	FGenericSparseMatrix;	// 'float'
typedef TGenericSparseMatrix<std::complex<double> > CGenericSparseMatrix; // 'complex'
typedef TGenericSparseMatrix<std::complex<float> >  SCGenericSparseMatrix;	// 's. complex'
typedef TGenericSparseMatrix<int>	IGenericSparseMatrix;	// 'integer'


// ==========================================================================
// ==========================================================================
// Member definitions

//template<class MT>
//IterativeMethod TGenericSparseMatrix<MT>::itmethod_general = ITMETHOD_CG;

//template<class MT>
//IterativeMethod TGenericSparseMatrix<MT>::itmethod_complex = ITMETHOD_BICGSTAB;

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT>::TGenericSparseMatrix ()
: TMatrix<MT> ()
{
    nbuf = nval = 0;
}

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT>::TGenericSparseMatrix (int rows, int cols,
    int nv, const MT *data)
: TMatrix<MT> (rows, cols)
{
    nbuf = nval = 0;
    Initialise (nv, data);
}

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT>::TGenericSparseMatrix (int rows, int cols,
    int nv, MT *data, CopyMode cmode)
: TMatrix<MT> (rows, cols)
{
    nbuf = nval = 0;
    Initialise (nv, data, cmode);
}

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT>::TGenericSparseMatrix (const TGenericSparseMatrix<MT>
    &m)
: TMatrix<MT> (m)
{
    nval = m.nval;
    nbuf = m.nbuf;
    int nsize = (nbuf > nval ? nbuf : nval);
    if (nsize) {
        val = new MT[nsize];
	for (int i = 0; i < nval; i++) val[i] = m.val[i];
    }
}

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT>::~TGenericSparseMatrix ()
{
    TGenericSparseMatrix<MT>::Unlink ();
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::New (int nrows, int ncols)
{
    TGenericSparseMatrix<MT>::Unlink ();
    TMatrix<MT>::New (nrows, ncols);
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Unlink ()
{
    if (nbuf) delete []val;
    nbuf = nval = 0;
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Initialise (int nv, const MT *data)
{
    int i;
    if (nv != nval) {
	if (nbuf) delete []val;
	if ((nbuf = (nval = nv))) val = new MT[nv];
    }
    if (data) for (i = 0; i < nv; i++) val[i] = data[i];
    else      for (i = 0; i < nv; i++) val[i] = (MT)0;
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Initialise (int nv, MT *data, CopyMode cmode)
{
    if (cmode == SHALLOW_COPY) {
	dASSERT (data, "Nonzero data buffer reference required");
	if (nbuf) {
	    delete []val;
	    nbuf = 0;
	}
	val = data;
	nval = nv;
    } else Initialise (nv, data);
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Zero ()
{
    for (int i = 0; i < nval; i++)
        val[i] = (MT)0;
}

// --------------------------------------------------------------------------

template<class MT>
TGenericSparseMatrix<MT> &TGenericSparseMatrix<MT>::operator*= (const MT &sc)
{
    for (int i = 0; i < nval; i++)
        val[i] *= sc;
    return *this;
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Append (MT v)
{
    if (nval == nbuf) { // no slack buffer - reallocation required
        MT *tmp_val = new MT[nbuf+BUFFER_CHUNK_SIZE];
	for (int i = 0; i < nval; i++) tmp_val[i] = val[i];
	if (nbuf) delete []val;
	val = tmp_val;
	nbuf += BUFFER_CHUNK_SIZE;
    }
    val[nval++] = v;
}

// --------------------------------------------------------------------------

template<class MT>
void TGenericSparseMatrix<MT>::Display (std::ostream &os) const
{
    for (int r = 0; r < this->rows; r++)
        for (int c = 0; c < this->cols; c++)
	    os << Get(r,c) << (c < this->cols-1 ? ' ' : '\n');
}

// --------------------------------------------------------------------------
inline int scale (double x)
{ return (int)(x*2.56e7+0.5); }

template<class MT>
void TGenericSparseMatrix<MT>::PrintFillinGraph (const char *fname, int maxdim,
    bool binary, bool antialias)
{
    int i, v, nx, ny, c, r = -1, *img;

    if (this->cols <= maxdim && this->rows <= maxdim) {
        nx = this->cols, ny = this->rows;
	img = new int[nx*ny];
	memset (img, 0, nx*ny*sizeof(int));
	for (;;) {
	    GetNext (r, c);
	    if (r < 0) break; // finished
	    img[r*nx+c] = 255;
	}
    } else {
        int ix, iy;
	double d, x, y, rx, ry, wx0, wx1, wx2, wy0, wy1, wy2, dummy;
        if (this->cols > this->rows)
	    nx = maxdim, ny = (maxdim*this->rows)/this->cols;
	else
	    ny = maxdim, nx = (maxdim*this->cols)/this->rows;
	d = 0.5*(double)nx/(double)this->cols;
	img = new int[nx*ny];
	memset (img, 0, nx*ny*sizeof(int));
	for (;;) {
	    GetNext (r, c);
	    if (r < 0) break; // finished
	    x = (c+0.5)/(double)this->cols * (double)nx;
	    rx = modf (x, &dummy);
	    ix = (int)(dummy+0.5);
	    y = (r+0.5)/(double)this->rows * (double)ny;
	    ry = modf (y, &dummy);
	    iy = (int)(dummy+0.5);
	    if (!antialias) {
	        img[iy*nx+ix] = 255;
		continue;
	    }
	    if (rx < d) {
	        wx0 = d-rx, wx1 = d+rx, wx2 = 0.0;
	    } else if (rx > 1.0-d) {
	        wx0 = 0.0, wx1 = d+1.0-rx, wx2 = d-1.0+rx;
	    } else {
	        wx0 = wx2 = 0.0, wx1 = d+d;
	    }
	    if (ry < d) {
	        wy0 = d-ry, wy1 = d+ry, wy2 = 0.0;
	    } else if (ry > 1.0-d) {
	        wy0 = 0.0, wy1 = d+1.0-ry, wy2 = d-1.0+ry;
	    } else {
	        wy0 = wy2 = 0.0, wy1 = d+d;
	    }
	    img[iy*nx+ix] += scale(wx1*wy1);
	    if (wx0 && ix) {
	        img[iy*nx+ix-1] += scale(wx0*wy1);
		if (wy0 && iy) img[(iy-1)*nx+ix-1] += scale(wx0*wy0);
		else if (wy2 && iy<ny-1) img[(iy+1)*nx+ix-1] += scale(wx0*wy2);
	    } else if (wx2 && ix<nx-1) {
	        img[iy*nx+ix+1] += scale(wx2*wy1);
		if (wy0 && iy) img[(iy-1)*nx+ix+1] += scale(wx2*wy0);
		else if (wy2 && iy<ny-1) img[(iy+1)*nx+ix+1] += scale(wx2*wy2);
	    }
	    if (wy0 && iy) {
	        img[(iy-1)*nx+ix] += scale(wx1*wy0);
	    } else if (wy2 && iy<ny-1) {
	        img[(iy+1)*nx+ix] += scale(wx1*wy2);
	    }
	}
	if (antialias)
	    for (i = 0; i < nx*ny; i++) // scale back down
	        if ((img[i] /= 100000) > 255) img[i] = 255;
    }
    std::ofstream ofs (fname);
    ofs << (binary ? "P5" : "P2") << std::endl;
    ofs << "# CREATOR: TOAST TGenericSparseMatrix::PrintFillinGraph" << std::endl;
    ofs << nx << ' ' << ny << std::endl;
    ofs << 255 << std::endl;
    for (r = i = 0; r < ny; r++) {
        for (c = 0; c < nx; c++) {
	    v = 255-img[r*nx+c];
	    if (binary)  ofs << (unsigned char)v;
	    else         ofs << v << (++i % 35 ? ' ' : '\n');
	}
    }
    delete []img;
}

// --------------------------------------------------------------------------
// select/return default iterative solver

template<class MT>
void TGenericSparseMatrix<MT>::GlobalSelectIterativeSolver
    (IterativeMethod method)
{ itmethod_general = method; }

template<class MT>
void TGenericSparseMatrix<MT>::GlobalSelectIterativeSolver_complex
    (IterativeMethod method)
{ itmethod_complex = method; }

template<class MT>
IterativeMethod TGenericSparseMatrix<MT>::GetGlobalIterativeSolver ()
{ return itmethod_general; }

template<class MT>
IterativeMethod TGenericSparseMatrix<MT>::GetGlobalIterativeSolver_complex ()
{ return itmethod_complex; }

// --------------------------------------------------------------------------

#ifdef CMAT2RMAT
/* add specific function for CGenericSparseMatrix */
RGenericSparseMatrix Re(const CGenericSparseMatrix &);
RGenericSparseMatrix Im(const CGenericSparseMatrix &);
#endif //!CMAT2RMAT

#endif // !__GSMATRIX_H
