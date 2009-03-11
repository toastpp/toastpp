// converts matrix formats between toast and ilupack
// Both store their matrices in compressed row format,
// but toast has zero-based index lists, ilupack has
// 1-based.

#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"

#ifdef HAVE_ILU

#define _DOUBLE_REAL_
#include "ilupack.hpp"
//#include "ilupackmacros.h"
#include "ilutoast.h"

static int permtype = ILUPERM_NULL;

void ILUSetPermtype (int ptype)
{
    permtype = ptype;
}

// Create an ILU Dmat from a TOAST RCompRowMatrix
// note that this allocates index lists, but uses a
// reference to the original data array, since that
// remains valid

#ifdef UNDEF
void CreateDmat (const RCompRowMatrix &A, Dmat *mat)
{
    mat->nr = A.nRows();
    mat->nc = A.nCols();
    mat->ia = new int[A.nRows()+1];
    mat->ja = new int[A.nVal()];

    mat->a = new double[A.nVal()];
    memcpy (mat->a, A.ValPtr(), A.nVal()*sizeof(double));
    // need to copy the value array because it is corrupted by
    // the solver

    for (int i = 0; i <= A.nRows(); i++)
	mat->ia[i] = A.rowptr[i]+1;
    for (int j = 0; j < A.nVal(); j++)
	mat->ja[j] = A.colidx[j]+1;
}

// Delete a Dmat by deallocating the index lists

void DeleteDmat (Dmat *mat)
{
    delete []mat->ia;
    delete []mat->ja;
    delete []mat->a;
}

void ILUSolveDGNL (const RCompRowMatrix &A, const RVector &b, RVector &x)
{
    Dmat Ailu;
    DAMGlevelmat PRE;
    DILUPACKparam param;
    int flags, elbow, max_it, ierr, nrestart;
    int n = A.nRows();
    int nlev = 0;
    doubleprecision condest, restol, droptols[2];
    integer (*perm0)(Dmat, double*,double*, int*, int*, int*,
		     DILUPACKparam*);
    integer (*perm)(Dmat, double*,double*, int*, int*, int*,
		    DILUPACKparam*);
    integer (*permf)(Dmat, double*,double*, int*, int*, int*,
		     DILUPACKparam*);

    // convert matrices and vectors
    CreateDmat (A, &Ailu);
    double *rhs = (double*)b.data_buffer();
    double *sol = x.data_buffer();

    size_t ndbuff = 3*n;
    doubleprecision *dbuff = (doubleprecision*)malloc
	(ndbuff*sizeof(doubleprecision));

    DGNLAMGinit (Ailu, &param);
    param.dbuff = dbuff;
    param.ndbuff = ndbuff;

    perm0 = DGNLperm_mmd;
    perm  = DGNLperm_mmd;
    permf = DGNLperm_mmd;

    DGNLAMGgetparams (param, &flags, &elbow, droptols, &condest, &restol,
		      &max_it, &nrestart);
    droptols[1] = 1e-2; elbow = 5;
    DGNLAMGsetparams (Ailu, &param, flags, elbow, droptols, condest, restol,
		      max_it, nrestart);
    ierr = DGNLAMGfactor (&Ailu, &PRE, &nlev, &param, perm0, perm, permf);
    ierr = DGNLAMGsolver (Ailu, PRE, nlev, &param, sol, rhs);
    //delete []dbuff;
}
#endif

void CreateZmat (const CCompRowMatrix &A, Zmat *mat)
{
    int i, j, nv = A.nVal();

    mat->nr = A.nRows();
    mat->nc = A.nCols();
    mat->ia = new int[A.nRows()+1];
    mat->ja = new int[nv];
    mat->a  = new doublecomplex[nv];

    const toast::complex *vptr = A.ValPtr();
    for (i = 0; i < nv; i++) {
	mat->a[i].r = vptr[i].re;
	mat->a[i].i = vptr[i].im;
    }

    for (i = 0; i <= A.nRows(); i++)
	mat->ia[i] = A.rowptr[i]+1;
    for (j = 0; j < A.nVal(); j++)
	mat->ja[j] = A.colidx[j]+1;
}

void CreateZSYMmat (const CCompRowMatrix &A, Zmat *mat)
{
    // creates a symmetric ILU Zmat from a general TOAST matrix
    // simply takes upper triangle without testing for symmetry

    int i, j, r, c, nz, nv = A.nVal();
    const toast::complex *vptr = A.ValPtr();

    mat->nr = A.nRows();
    mat->nc = A.nCols();
    mat->ia = new int[A.nRows()+1];
    mat->ia[0] = 1;

    // pass 1: scan for nonzeros in upper triangle
    for (r = 0; r < A.nRows(); r++) {
	nz = 0;
	for (i = A.rowptr[r]; i < A.rowptr[r+1]; i++) {
	    c = A.colidx[i];
	    if (c >= r) nz++; // upper triangle
	}
	mat->ia[r+1] = mat->ia[r]+nz;
    }
    nz = mat->ia[A.nRows()]-1;
    mat->ja = new int[nz];
    mat->a  = new doublecomplex[nz];

    // pass 2: fill column index and value arrays
    for (r = nz = 0; r < A.nRows(); r++) {
	for (i = A.rowptr[r]; i < A.rowptr[r+1]; i++) {
	    c = A.colidx[i];
	    if (c >= r) {
		mat->ja[nz] = c+1;
		mat->a[nz].r = vptr[i].re;
		mat->a[nz].i = vptr[i].im;
		nz++;
	    }
	}
    }
}

void DeleteZmat (Zmat *mat)
{
    delete []mat->ia;
    delete []mat->ja;
    delete []mat->a;
}

int ILUSolveZGNL (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol, double droptol, int maxit)
{
    Zmat Ailu;
    ZAMGlevelmat PRE;
    ZILUPACKparam param;
    int flags, elbow, max_it, ierr, nrestart;
    int n = A.nRows();
    int nlev = 0;
    doubleprecision condest, restol, droptols[2];
    integer (*perm0)(Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);
    integer (*perm) (Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);
    integer (*permf)(Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);

    // convert matrices and vectors
    CreateZmat (A, &Ailu);
    doublecomplex *rhs = (doublecomplex*)b.data_buffer();
    doublecomplex *sol = (doublecomplex*)x.data_buffer();

    size_t ndbuff = 3*n;
    doublecomplex *dbuff = (doublecomplex*)malloc
	(ndbuff*sizeof(doublecomplex));

    ZGNLAMGinit (Ailu, &param);
    param.dbuff = dbuff;
    param.ndbuff = ndbuff;

    switch (permtype) {
    case ILUPERM_MMD:
	perm0 = ZGNLperm_mmd;
	perm  = ZGNLperm_mmd;
	permf = ZGNLperm_mmd;
	break;
    case ILUPERM_RCM:
	perm0 = ZGNLperm_rcm;
	perm  = ZGNLperm_rcm;
	permf = ZGNLperm_rcm;
	break;
    case ILUPERM_ND:
	perm0 = ZGNLperm_nd;
	perm  = ZGNLperm_nd;
	permf = ZGNLperm_nd;
	break;
    case ILUPERM_INDSET:
	perm0 = ZGNLperm_indset;
	perm  = ZGNLperm_indset;
	permf = ZGNLperm_indset;
	break;
    case ILUPERM_AMF:
	perm0 = ZGNLperm_amf;
	perm  = ZGNLperm_amf;
	permf = ZGNLperm_amf;
	break;
    case ILUPERM_AMD:
	perm0 = ZGNLperm_amd;
	perm  = ZGNLperm_amd;
	permf = ZGNLperm_amd;
	break;
    case ILUPERM_PQ:
	perm0 = ZGNLperm_pq;
	perm  = ZGNLperm_pq;
	permf = ZGNLperm_pq;
	break;
    default:
	perm0 = ZGNLperm_null;
	perm  = ZGNLperm_null;
	permf = ZGNLperm_null;
	break;
    }

    ZGNLAMGgetparams (param, &flags, &elbow, droptols, &condest, &restol,
		      &max_it, &nrestart);
    droptols[1] = droptol; elbow = 5;
    restol = tol;
    max_it = maxit;
    ZGNLAMGsetparams (Ailu, &param, flags, elbow, droptols, condest, restol,
		      max_it, nrestart);
    ierr = ZGNLAMGfactor (&Ailu, &PRE, &nlev, &param, perm0, perm, permf);
    ierr = ZGNLAMGsolver (Ailu, PRE, nlev, &param, sol, rhs);

    return ierr;
    //delete []dbuff;
    
}

int ILUSolveZSYM (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol, double droptol, int maxit)
{
    Zmat Ailu;
    ZAMGlevelmat PRE;
    ZILUPACKparam param;
    int flags, elbow, max_it, ierr, nrestart;
    int n = A.nRows();
    int nlev = 0;
    doubleprecision condest, restol, droptols[2];
    integer (*perm0)(Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);
    integer (*perm) (Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);
    integer (*permf)(Zmat, doublecomplex*, doublecomplex*, int*, int*, int*,
		     ZILUPACKparam*);

    // convert matrices and vectors
    CreateZSYMmat (A, &Ailu);
    doublecomplex *rhs = (doublecomplex*)b.data_buffer();
    doublecomplex *sol = (doublecomplex*)x.data_buffer();

    size_t ndbuff = 3*n;
    doublecomplex *dbuff = (doublecomplex*)malloc
	(ndbuff*sizeof(doublecomplex));

    ZSYMAMGinit (Ailu, &param);
    param.dbuff = dbuff;
    param.ndbuff = ndbuff;

    switch (permtype) {
    case ILUPERM_MMD:
	perm0 = ZSYMperm_mmd;
	perm  = ZSYMperm_mmd;
	permf = ZSYMperm_mmd;
	break;
    case ILUPERM_RCM:
	perm0 = ZSYMperm_rcm;
	perm  = ZSYMperm_rcm;
	permf = ZSYMperm_rcm;
	break;
    case ILUPERM_ND:
	perm0 = ZSYMperm_nd;
	perm  = ZSYMperm_nd;
	permf = ZSYMperm_nd;
	break;
    case ILUPERM_INDSET:
	perm0 = ZSYMperm_indset;
	perm  = ZSYMperm_indset;
	permf = ZSYMperm_indset;
	break;
    case ILUPERM_AMF:
	perm0 = ZSYMperm_amf;
	perm  = ZSYMperm_amf;
	permf = ZSYMperm_amf;
	break;
    case ILUPERM_AMD:
	perm0 = ZSYMperm_amd;
	perm  = ZSYMperm_amd;
	permf = ZSYMperm_amd;
	break;
    default:
	perm0 = ZSYMperm_null;
	perm  = ZSYMperm_null;
	permf = ZSYMperm_null;
	break;
    }

    ZSYMAMGgetparams (param, &flags, &elbow, droptols, &condest, &restol,
		      &max_it);
    droptols[1] = droptol; elbow = 3;
    restol = tol;
    max_it = maxit;
    ZSYMAMGsetparams (Ailu, &param, flags, elbow, droptols, condest, restol,
		      max_it);
    ierr = ZSYMAMGfactor (&Ailu, &PRE, &nlev, &param, perm0, perm, permf);
    ierr = ZSYMAMGsolver (Ailu, PRE, nlev, &param, sol, rhs);

    return ierr;
    //delete []dbuff;
    

}


#else // !HAVE_ILU

void ILUSetPermtype (int ptype)
{
    xERROR (No ILUPack support provided);
}

int ILUSolveZGNL (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol, double droptol, int maxit)
{
    xERROR (No ILUPack support provided);
    return 0;
}

int ILUSolveZSYM (const CCompRowMatrix &A, const CVector &b, CVector &x,
    double tol, double droptol, int maxit)
{
    xERROR (No ILUPack support provided);
    return 0;
}

#endif // HAVE_ILU
