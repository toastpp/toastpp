// ==========================================================================
// Module mathlib
// File precon.cc
// Definition of template classes
//     TPreconditioner
//     TPrecon_Null
//     TPrecon_Diag
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#include "mathlib.h"

#ifdef HAVE_ILU
#include "ilutoast.h"
#include <ilupack.h>
#endif // HAVE_ILU

using namespace std;

// ==========================================================================
// class TPreconditioner

template<class MT>
TPreconditioner<MT> *TPreconditioner<MT>::Create (PreconType type)
{
    switch (type) {
    case PRECON_NULL:	      // no preconditioning
        return new TPrecon_Null<MT>;
	break;
    case PRECON_DIAG:	      // diagonal of system matrix
        return new TPrecon_Diag<MT>;
	break;
    case PRECON_ICH:	      // incomplete Cholesky decomposition
        return new TPrecon_IC<MT>;
	break;
    case PRECON_DILU:         // diagonal incomplete LU decomposition
        return new TPrecon_DILU<MT>;
	break;
#ifdef HAVE_ILU
     case PRECON_ILU:         // diagonal incomplete LU decomposition
        return new TPrecon_ILU<MT>;
	break;
#endif // HAVE_ILU
    //case PRECON_CG_MULTIGRID: // CG multigrid
        //return new TPrecon_CG_Multigrid<MT>;
	//break;
    default:
        xERROR ("Unknown preconditioner type");
	return 0;
    }
}

// ==========================================================================
// class TPrecon_Diag

template<class MT>
void TPrecon_Diag<MT>::Reset (const TMatrix<MT> *A)
{
    int dim = (A->nRows() < A->nCols() ? A->nRows() : A->nCols());
    if (idiag.Dim() != dim) idiag.New (dim);
    TVector<MT> diag = A->Diag();
    for (int i = 0; i < dim; i++)
	idiag[i] = (diag[i] != (MT)0 ? (MT)1/diag[i] : (MT)1);
}

template<class MT>
void TPrecon_Diag<MT>::ResetFromDiagonal (const TVector<MT> &diag)
{
    int dim = diag.Dim();
    if (idiag.Dim() != dim) idiag.New (dim);
    idiag = inv (diag);
}

template<class MT>
void TPrecon_Diag<MT>::Apply (const TVector<MT> &r, TVector<MT> &s) 
{
    dASSERT(r.Dim() == idiag.Dim(), "Dimension mismatch");
    dASSERT(s.Dim() == idiag.Dim(), "Dimension mismatch");
    s = r * idiag;
}

template<class MT>
void TPrecon_Diag<MT>::Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const
 {
    dASSERT(r.nRows() == idiag.Dim(), "Dimension mismatch");
    dASSERT(s.nRows() == r.nRows(), "Dimension mismatch");
    dASSERT(s.nCols() == r.nCols(), "Dimension mismatch");
    int nr = r.nRows();
    int nc = r.nCols();
    for(int i=0; i < nr; i++)
	for(int j=0; j < nc; j++)
	    s.Set(i, j, idiag[i]*r.Get(i, j));
}


// ==========================================================================
// class TPrecon_IC

template<class MT>
void TPrecon_IC<MT>::Reset (const TMatrix<MT> *A)
{
    // *** This assumes A is of type TCompRowMatrix ***

    xASSERT(A->isSparse(), "Sparse matrix type required for IC preconditioner");
    const TCompRowMatrix<MT> *spA = (const TCompRowMatrix<MT>*)A;

    int nrows = A->nRows(), ncols = A->nCols();
    if (L.nRows() != nrows || L.nCols() != ncols) {
        L.New (nrows, ncols);
	d.New (nrows);
	int *rowptr, *colidx;
	spA->CalculateIncompleteCholeskyFillin (rowptr, colidx);
	L.Initialise (rowptr, colidx);
	delete []rowptr;
	delete []colidx;
    }
    IncompleteCholeskyFactorize (*spA, L, d, true);
    L.SetColAccess();
}

template<class MT>
void TPrecon_IC<MT>::Apply (const TVector<MT> &r, TVector<MT> &s)
{
    CholeskySolve (L, d, r, s);
}

// ==========================================================================
// class TPrecon_DILU

template<class MT>
void TPrecon_DILU<MT>::Reset (const TMatrix<MT> *_A)
{
    // assumes symmetric A
    A = _A;

    const MT unit = (MT)1;
    int i, r, c, nz;
    dim = (A->nRows() < A->nCols() ? A->nRows() : A->nCols());
    TVector<MT> Ar(A->nCols());
    int *ci = new int[dim];
    MT *rv = new MT[dim];
    if (ipivot.Dim() != dim) ipivot.New (dim);
    ipivot = A->Diag();
    for (r = 0; r < dim; r++) {
        ipivot[r] = unit/ipivot[r];
	nz = A->SparseRow (r, ci, rv);
	for (i = 0; i < nz; i++)
	    if ((c = ci[i]) > r) ipivot[c] -= rv[i]*rv[i]*ipivot[r];
    }
    delete []ci;
    delete []rv;
}

template<class MT>
void TPrecon_DILU<MT>::Apply (const TVector<MT> &r, TVector<MT> &s) 
{
    int i, j, k, nz;
    MT sum;
    TVector<MT> row(dim);
    TVector<MT> z(dim);
    int *ci = new int[dim];
    MT  *rv = new MT[dim];
    for (i = 0; i < dim; i++) {
        nz = A->SparseRow (i, ci, rv);
	sum = (MT)0;
	for (k = 0; k < nz; k++)
	    if ((j = ci[k]) < i) sum += rv[k]*z[j];
	z[i] = ipivot[i]*(r[i]-sum);
    }
    for (i = dim-1; i >= 0; i--) {
        nz = A->SparseRow (i, ci, rv);
	sum = (MT)0;
	for (k = 0; k < nz; k++)
	    if ((j = ci[k]) > i) sum += rv[k]*s[j];
	s[i] = z[i] - ipivot[i]*sum;
    }
    delete []ci;
    delete []rv;
}

// ==========================================================================
// class TPrecon_CG_Multigrid

template<class MT>
void TPrecon_CG_Multigrid<MT>::Reset (const TMatrix<MT> *_A)
{
    xASSERT(_A->StorageType() == MATRIX_COMPROW,
	    "Only implemented for CompRow matrix type");
    A = (TCompRowMatrix<MT>*)_A;
}

template<class MT>
void TPrecon_CG_Multigrid<MT>::Apply (const TVector<MT> &r, TVector<MT> &s)
{
    MGM (r, s, 0);
}

template<class MT>
double TPrecon_CG_Multigrid<MT>::MGM (const TVector<MT> &r, TVector<MT> &x,
    int level) const
{
    double tol = 1e-8, err = 1;

    int j, nn, n = x.Dim();
    cout << "MGM level " << level << " multigrid cycles " << nmg
	 << " smoothing steps " << ngs << endl;

    if (level == maxlevel-1) {
        CG (A[level], r, x, tol, precon[level], n);
      	return tol;
    } else {
        nn = P[level]->nCols(); // dimension at next coarser level
        Smoother (A[level], r, x,ngs);
	TVector<MT> r2(n);
	TVector<MT> x2(n);
	TVector<MT> r1(nn);
	TVector<MT> x1(nn);
	A[level].Ax (x, r2);
	r2 = r -r2;                    // residual at this level 
	r1 = (*R[level])* r2;            // restriction
	for (j = 0; j < nmg && err > tol; j++) 
	  err = MGM (r1, x1, level+1);
	x += (*P[level])*x1;           // prolongation
	r2 = r - A[level]*x;          //    r2 = b2  - A2.f2;    (* residual *)
        Smoother (A[level], r2, x2,ngs);      // smooth
	x += x2;
	r2 = r - A[level]*x;          //    r2 = b2  - A2.f2;    (* residual *)
	err = l2normsq (r2);
	cout << "\terror at MGM level " << level << ':'<< err << endl;
      	return err;
    }
}

template<class MT>
void TPrecon_CG_Multigrid<MT>::Smoother (const TCompRowMatrix<MT> &A,
    const TVector<MT> &r, TVector<MT> &x, int itermax) const
{
    double tol = 1e-10;
    GaussSeidel (A, r, x, tol, itermax);
}

// ==========================================================================
// class SCPreconditionerMixed and derived classes

inline SCPreconditionerMixed *SCPreconditionerMixed::Create (PreconType type)
{
    switch (type) {
    case PRECON_NULL:	      // no preconditioning
        return new SCPreconMixed_Null;
	break;
    case PRECON_DIAG:	      // diagonal of system matrix
        return new SCPreconMixed_Diag;
	break;
	//case PRECON_ICH:	      // incomplete Cholesky decomposition
        //return new TPrecon_IC<MT>;
	//break;
    case PRECON_DILU:         // diagonal incomplete LU decomposition
        return new SCPreconMixed_DILU;
	break;
    //case PRECON_CG_MULTIGRID: // CG multigrid
        //return new TPrecon_CG_Multigrid<MT>;
	//break;
    default:
        xERROR ("Unknown preconditioner type");
	return 0;
    }
}

inline void SCPreconMixed_Diag::Reset (const SCCompRowMatrixMixed *A)
{
    int dim = (A->nRows() < A->nCols() ? A->nRows() : A->nCols());
    if (idiag.Dim() != dim) idiag.New (dim);
    idiag = MakeCVector(inv (A->Diag()));
}

inline void SCPreconMixed_Diag::Apply (const CVector &r, CVector &s)
{
    dASSERT(r.Dim() == idiag.Dim(), "Dimension mismatch");
    dASSERT(s.Dim() == idiag.Dim(), "Dimension mismatch");
    s = r * idiag;
}

inline void SCPreconMixed_DILU::Reset (const SCCompRowMatrixMixed *_A)
{
    // assumes symmetric A
    A = _A;

    const std::complex<double> unit = std::complex<double>(1,0);
    int i, r, c, nz;
    dim = (A->nRows() < A->nCols() ? A->nRows() : A->nCols());
    CVector Ar(A->nCols());
    int *ci = new int[dim];
    std::complex<double> *rv = new std::complex<double>[dim];
    if (ipivot.Dim() != dim) ipivot.New (dim);
    ipivot = MakeCVector (A->Diag());
    for (r = 0; r < dim; r++) {
        ipivot[r] = unit/ipivot[r];
	nz = A->SparseRow (r, ci, rv);
	for (i = 0; i < nz; i++)
	    if ((c = ci[i]) > r) ipivot[c] -= rv[i]*rv[i]*ipivot[r];
    }
    delete []ci;
    delete []rv;
}

inline void SCPreconMixed_DILU::Apply (const CVector &r, CVector &s)
{
    int i, j, k, nz;
    std::complex<double> sum;
    CVector row(dim);
    CVector z(dim);
    int *ci = new int[dim];
    std::complex<double> *rv = new std::complex<double>[dim];
    for (i = 0; i < dim; i++) {
        nz = A->SparseRow (i, ci, rv);
	sum = std::complex<double>(0,0);
	for (k = 0; k < nz; k++)
	    if ((j = ci[k]) < i) sum += rv[k]*z[j];
	z[i] = ipivot[i]*(r[i]-sum);
    }
    for (i = dim-1; i >= 0; i--) {
        nz = A->SparseRow (i, ci, rv);
	sum = std::complex<double>(0,0);
	for (k = 0; k < nz; k++)
	    if ((j = ci[k]) > i) sum += rv[k]*s[j];
	s[i] = z[i] - ipivot[i]*sum;
    }
    delete []ci;
    delete []rv;
}

// ==========================================================================
// class TPrecon_ILU

#ifdef HAVE_ILU

template<class MT>
void TPrecon_ILU<MT>::Reset (TCompRowMatrix<MT> &, int matching, char *ordering, double droptol, int condest, int elbow)
{
    ERROR_UNDEF;
}

template<class MT>
TPrecon_ILU<MT>::~TPrecon_ILU ()
{
	ZGNLAMGdelete(&A,&PRE,&param);
	delete []rhs;
	delete []sol;
	delete []A.ia;
	delete []A.ja;
	delete []A.a;
}

template<class MT>
void TPrecon_ILU<MT>::Apply (const TVector<MT> &rh, TVector<MT> &s)
{
    ERROR_UNDEF;
}

#endif // HAVE_ILU

// ==========================================================================
// class and friend instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TPreconditioner<float>;
template class MATHLIB TPrecon_Null<float>;
template class MATHLIB TPrecon_Diag<float>;
template class MATHLIB TPrecon_IC<float>;
template class MATHLIB TPrecon_DILU<float>;
template class MATHLIB TPrecon_CG_Multigrid<float>;

template class MATHLIB TPreconditioner<double>;
template class MATHLIB TPrecon_Null<double>;
template class MATHLIB TPrecon_Diag<double>;
template class MATHLIB TPrecon_IC<double>;
template class MATHLIB TPrecon_DILU<double>;
template class MATHLIB TPrecon_CG_Multigrid<double>;

template class MATHLIB TPreconditioner<toast::complex>;
template class MATHLIB TPrecon_Null<toast::complex>;
template class MATHLIB TPrecon_Diag<toast::complex>;
template class MATHLIB TPrecon_IC<toast::complex>;
template class MATHLIB TPrecon_DILU<toast::complex>;
template class MATHLIB TPrecon_CG_Multigrid<toast::complex>;

template class MATHLIB TPreconditioner<scomplex>;
template class MATHLIB TPrecon_Null<scomplex>;
template class MATHLIB TPrecon_Diag<scomplex>;
template class MATHLIB TPrecon_IC<scomplex>;
template class MATHLIB TPrecon_DILU<scomplex>;
template class MATHLIB TPrecon_CG_Multigrid<scomplex>;

#ifdef HAVE_ILU
template class MATHLIB TPrecon_ILU<float>;
template class MATHLIB TPrecon_ILU<double>;
template class MATHLIB TPrecon_ILU<toast::complex>;
template class MATHLIB TPrecon_ILU<scomplex>;
#endif // HAVE_ILU

#endif // NEED_EXPLICIT_INSTANTIATION
