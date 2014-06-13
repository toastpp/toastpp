// Single-complex part of fwdsolver
// Required in a separate file because SuperLU header files for different
// types cannot be combined

#define STOASTLIB_IMPLEMENTATION
#include "slu_zdefs.h"
#include "supermatrix.h"
#include "stoastlib.h"
#include "fwdsolver_zslu.h"

// =========================================================================

class ZSuperLU_engine {
public:
    ZSuperLU_engine();
    ~ZSuperLU_engine();
    void Reset (const CCompRowMatrix *F);
    void Solve (SuperMatrix *B, SuperMatrix *X);
    SuperMatrix A, L, U;
    int *perm_c;
    int *perm_r;
    int *etree;
    double *R, *C;
    superlu_options_t options;
    
private:
    void AllocMatrix (const CCompRowMatrix *F);
    void Deallocate ();
    bool allocated;
};

// =========================================================================

ZSuperLU_engine::ZSuperLU_engine ()
{
    perm_c = 0;
    perm_r = 0;
    etree  = 0;
    R      = 0;
    C      = 0;
    allocated = false;
    
    set_default_options(&options);
    options.Fact = DOFACT;
    options.Equil = NO;
    options.ColPerm = NATURAL;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
}

// =========================================================================

ZSuperLU_engine::~ZSuperLU_engine ()
{
    Deallocate();
}

// =========================================================================

void ZSuperLU_engine::Deallocate ()
{
    if (allocated) {
	delete []perm_c;
	delete []perm_r;
	delete []etree;
	delete []R;
	delete []C;
	Destroy_SuperMatrix_Store (&A);
	allocated = false;
	if (options.Fact != DOFACT) {
	    Destroy_SuperNode_Matrix (&L);
	    Destroy_CompCol_Matrix (&U);
	}
    }
}

// =========================================================================

void ZSuperLU_engine::Reset (const CCompRowMatrix *F)
{
    Deallocate();
    AllocMatrix (F);
    int m = F->nRows();
    int n = F->nCols();
    perm_r = new int[m];
    perm_c = new int[n];
    etree  = new int[n];
    R      = new double[m];
    C      = new double[n];
    get_perm_c (0, &A, perm_c);
    options.Fact = DOFACT;
    allocated = true;
}

// =========================================================================

void ZSuperLU_engine::AllocMatrix (const CCompRowMatrix *F)
{
    int m = F->nRows();
    int n = F->nCols();
    int nz = F->nVal();
    doublecomplex *cdat = (doublecomplex*)F->ValPtr();
    zCreate_CompCol_Matrix (&A, m, n, nz, cdat, 
        F->colidx, F->rowptr, SLU_NR, SLU_Z, SLU_GE);
}

// =========================================================================

void ZSuperLU_engine::Solve (SuperMatrix *B, SuperMatrix *X)
{
    int info;
    mem_usage_t mem_usage;
    char equed = 'N';
    double R = 0.0;
    double C = 0.0;
    static int nrhs_max = 1;
    static double *ferr = new double[nrhs_max];
    static double *berr = new double[nrhs_max];
    double recip_pivot_growth, rcond;
    int nrhs = B->ncol;
    if (nrhs > nrhs_max) {
        nrhs_max = nrhs;
	delete []ferr;   ferr = new double[nrhs_max];
	delete []berr;   berr = new double[nrhs_max];
    }


    SuperLUStat_t stat;
    StatInit (&stat);
    
    zgssvx (&options, &A, perm_c, perm_r, etree, &equed, &R, &C,
    		  &L, &U, 0, 0, B, X, &recip_pivot_growth, &rcond,
    		  ferr, berr, &mem_usage, &stat, &info);

    options.Fact = FACTORED;
    StatFree (&stat);
}

// =========================================================================
// =========================================================================

ZSuperLU::ZSuperLU ()
{
    A = 0;
    engine = new ZSuperLU_engine;
}

// =========================================================================

ZSuperLU::~ZSuperLU ()
{
    delete engine;
}

// =========================================================================

void ZSuperLU::Reset (const CCompRowMatrix *F)
{
    A = F;
    engine->Reset (F);
}

// =========================================================================

void ZSuperLU::CalcField (const CVector &qvec, CVector &phi,
    IterativeSolverResult *res) const
{
    SuperMatrix B, X;
    int m = qvec.Dim();
    int n = phi.Dim();

    dASSERT (A, "Matrix not defined");
    dASSERT (A->nRows() == m && A->nCols() == n, "Invalid vector sizes");

    doublecomplex *rhsbuf = (doublecomplex*)qvec.data_buffer();
    doublecomplex *xbuf   = (doublecomplex*)phi.data_buffer();
    zCreate_Dense_Matrix (&B, m, 1, rhsbuf, m, SLU_DN, SLU_Z,SLU_GE);
    zCreate_Dense_Matrix (&X, n, 1, xbuf, n, SLU_DN, SLU_Z, SLU_GE);

    engine->Solve (&B, &X);
    
    Destroy_SuperMatrix_Store (&B);
    Destroy_SuperMatrix_Store (&X);
}

// =========================================================================

void ZSuperLU::CalcFields (const CCompRowMatrix &qvec, CVector *phi,
    IterativeSolverResult *res) const
{
    SuperMatrix B, X;
    int i, r;
    const idxtype *rptr, *cidx;

    int m = qvec.nCols();
    int n = phi[0].Dim();
    int nrhs = qvec.nRows();

    dASSERT (A, "Matrix not defined");
    dASSERT (A->nRows() == m && A->nCols() == n, "Invalid vector sizes");

    // write sparse source vector array into dense matrix
    CDenseMatrix qd(qvec);
    doublecomplex *qbuf = (doublecomplex*)qd.ValPtr();
    zCreate_Dense_Matrix (&B, m, nrhs, qbuf, m, SLU_DN, SLU_Z, SLU_GE);

    std::complex<double> *x = new std::complex<double>[n*nrhs];
    for (i = 0; i < nrhs; i++)
        memcpy (x+(i*n), phi[i].data_buffer(), n*sizeof(std::complex<double>));
    doublecomplex *xbuf = (doublecomplex*)x;
    zCreate_Dense_Matrix (&X, n, nrhs, xbuf, n, SLU_DN, SLU_Z, SLU_GE);

    engine->Solve (&B, &X);
    
    for (i = 0; i < nrhs; i++)
        memcpy (phi[i].data_buffer(), x+(i*n), n*sizeof(std::complex<double>));

    Destroy_SuperMatrix_Store (&B);
    Destroy_SuperMatrix_Store (&X);
    
    delete []x;
}

