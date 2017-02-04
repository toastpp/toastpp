// Single-complex part of fwdsolver
// Required in a separate file because SuperLU header files for different
// types cannot be combined

#define STOASTLIB_IMPLEMENTATION
#include "slu_zdefs.h"
#include "supermatrix.h"
#include "stoastlib.h"
#include "fwdsolver_zslu.h"

// =========================================================================
// =========================================================================

ZSuperLU::ZSuperLU (int n)
{
    A         = 0;
    perm_c    = 0;
    perm_r    = 0;
    etree     = 0;
    R         = 0;
    C         = 0;
    allocated = false;

    set_default_options(&options);
    options.Fact = DOFACT;
    options.Equil = NO;
    options.ColPerm = NATURAL;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
}

// =========================================================================

ZSuperLU::~ZSuperLU ()
{
    Deallocate();
}

// =========================================================================

void ZSuperLU::AllocMatrix (const CCompRowMatrix *F)
{
    int m = F->nRows();
    int n = F->nCols();
    int nz = F->nVal();
    doublecomplex *cdat = (doublecomplex*)F->ValPtr();
    zCreate_CompCol_Matrix (&sA, m, n, nz, cdat, 
        F->colidx, F->rowptr, SLU_NR, SLU_Z, SLU_GE);    
}

// =========================================================================

void ZSuperLU::Deallocate()
{
    if (allocated) {
	delete []perm_c;
	delete []perm_r;
	delete []etree;
	delete []R;
	delete []C;
	Destroy_SuperMatrix_Store (&sA);
	allocated = false;
	if (options.Fact != DOFACT) {
	    Destroy_SuperNode_Matrix (&sL);
	    Destroy_CompCol_Matrix (&sU);
	}
    }
}

// =========================================================================

void ZSuperLU::Reset (const CCompRowMatrix *F)
{
    A = F;
    Deallocate();
    AllocMatrix (F);
    int m = F->nRows();
    int n = F->nCols();
    perm_r = new int[m];
    perm_c = new int[n];
    etree  = new int[n];
    R      = new double[m];
    C      = new double[n];
    get_perm_c (0, &sA, perm_c);
    options.Fact = DOFACT;
    allocated = true;
}

// =========================================================================

void ZSuperLU::Solve (SuperMatrix *sB, SuperMatrix *sX) const
{
    int info;
    mem_usage_t mem_usage;
    char equed = 'N';
    double R = 0.0;
    double C = 0.0;
    GlobalLU_t Glu;
    int nrhs = sB->ncol;
	double *ferr = new double[nrhs];
	double *berr = new double[nrhs];
    double recip_pivot_growth, rcond;

    SuperLUStat_t stat;
    StatInit (&stat);
    
    zgssvx (&options, &sA, perm_c, perm_r, etree, &equed, &R, &C,
	    &sL, &sU, 0, 0, sB, sX, &recip_pivot_growth, &rcond,
	    ferr, berr, &Glu, &mem_usage, &stat, &info);

    options.Fact = FACTORED;
    StatFree (&stat);
	delete []ferr;
	delete []berr;
}

// =========================================================================

void ZSuperLU::CalcField (const CVector &qvec, CVector &phi,
    IterativeSolverResult *res, int en) const
{
    SuperMatrix sB, sX;
    int m = qvec.Dim();
    int n = phi.Dim();

    dASSERT (A, "Matrix not defined");
    dASSERT (A->nRows() == m && A->nCols() == n, "Invalid vector sizes");

    doublecomplex *rhsbuf = (doublecomplex*)qvec.data_buffer();
    doublecomplex *xbuf   = (doublecomplex*)phi.data_buffer();
    zCreate_Dense_Matrix (&sB, m, 1, rhsbuf, m, SLU_DN, SLU_Z,SLU_GE);
    zCreate_Dense_Matrix (&sX, n, 1, xbuf, n, SLU_DN, SLU_Z, SLU_GE);

    Solve (&sB, &sX);
    
    Destroy_SuperMatrix_Store (&sB);
    Destroy_SuperMatrix_Store (&sX);
}

// =========================================================================

void ZSuperLU::CalcFields (const CCompRowMatrix &qvec, CVector *phi,
    IterativeSolverResult *res, int en) const
{
    SuperMatrix sB, sX;
    int i;

    int m = qvec.nCols();
    int n = phi[0].Dim();
    int nrhs = qvec.nRows();

    dASSERT (A, "Matrix not defined");
    dASSERT (A->nRows() == m && A->nCols() == n, "Invalid vector sizes");

    // write sparse source vector array into dense matrix
    CDenseMatrix qd(qvec);
    doublecomplex *qbuf = (doublecomplex*)qd.ValPtr();
    zCreate_Dense_Matrix (&sB, m, nrhs, qbuf, m, SLU_DN, SLU_Z, SLU_GE);

    std::complex<double> *x = new std::complex<double>[n*nrhs];
    for (i = 0; i < nrhs; i++)
        memcpy (x+(i*n), phi[i].data_buffer(), n*sizeof(std::complex<double>));
    doublecomplex *xbuf = (doublecomplex*)x;
    zCreate_Dense_Matrix (&sX, n, nrhs, xbuf, n, SLU_DN, SLU_Z, SLU_GE);

    Solve (&sB, &sX);
    
    for (i = 0; i < nrhs; i++)
        memcpy (phi[i].data_buffer(), x+(i*n), n*sizeof(std::complex<double>));

    Destroy_SuperMatrix_Store (&sB);
    Destroy_SuperMatrix_Store (&sX);
    
    delete []x;
}

