// Single-complex part of fwdsolver
// Required in a separate file because SuperLU header files for different
// types cannot be combined

#define STOASTLIB_IMPLEMENTATION
#include "slu_cdefs.h"
#include "supermatrix.h"
#include "stoastlib.h"
#include "fwdsolver_cslu.h"

// =========================================================================

class CSuperLU_engine {
public:
    CSuperLU_engine();
    ~CSuperLU_engine();
    void Reset (const SCCompRowMatrix *F);
    void Solve (SuperMatrix *B, SuperMatrix *X);
    SuperMatrix A, L, U;
    int *perm_c;
    int *perm_r;
    int *etree;
    float *R, *C;
    superlu_options_t options;
    
private:
    void AllocMatrix (const SCCompRowMatrix *F);
    void Deallocate ();
    bool allocated;
};

// =========================================================================

CSuperLU_engine::CSuperLU_engine ()
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

CSuperLU_engine::~CSuperLU_engine ()
{
    Deallocate();
}

// =========================================================================

void CSuperLU_engine::Deallocate ()
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

void CSuperLU_engine::Reset (const SCCompRowMatrix *F)
{
    Deallocate();
    AllocMatrix (F);
    int m = F->nRows();
    int n = F->nCols();
    perm_r = new int[m];
    perm_c = new int[n];
    etree  = new int[n];
    R      = new float[m];
    C      = new float[n];
    get_perm_c (0, &A, perm_c);
    options.Fact = DOFACT;
    allocated = true;
}

// =========================================================================

void CSuperLU_engine::AllocMatrix (const SCCompRowMatrix *F)
{
    int m = F->nRows();
    int n = F->nCols();
    int nz = F->nVal();
    ::complex *cdat = (::complex*)F->ValPtr();
    cCreate_CompCol_Matrix (&A, m, n, nz, cdat, 
        F->colidx, F->rowptr, SLU_NR, SLU_C, SLU_GE);
}

// =========================================================================

void CSuperLU_engine::Solve (SuperMatrix *B, SuperMatrix *X)
{
    int info;
    mem_usage_t mem_usage;
    char equed = 'N';
    float R = 0.0;
    float C = 0.0;
    GlobalLU_t Glu;
    float ferr, berr;
    float recip_pivot_growth, rcond;
    SuperLUStat_t stat;
    StatInit (&stat);

    cgssvx (&options, &A, perm_c, perm_r, etree, &equed, &R, &C,
	    &L, &U, 0, 0, B, X, &recip_pivot_growth, &rcond,
	    &ferr, &berr, &Glu, &mem_usage, &stat, &info);
    options.Fact = FACTORED;
    StatFree (&stat);
}

// =========================================================================
// =========================================================================

CSuperLU::CSuperLU ()
{
    A = 0;
    engine = new CSuperLU_engine;
}

// =========================================================================

CSuperLU::~CSuperLU ()
{
    delete engine;
}

// =========================================================================

void CSuperLU::Reset (const SCCompRowMatrix *F)
{
    A = F;
    engine->Reset (F);
}

// =========================================================================

void CSuperLU::CalcField (const SCVector &qvec, SCVector &phi,
    IterativeSolverResult *res) const
{
    SuperMatrix B, X;
    int m = qvec.Dim();
    int n = phi.Dim();

    dASSERT (A, "Matrix not defined");
    dASSERT (A->nRows() == m && A->nCols() == n, "Invalid vector sizes");

    ::complex *rhsbuf = (::complex*)qvec.data_buffer();
    ::complex *xbuf   = (::complex*)phi.data_buffer();
    cCreate_Dense_Matrix (&B, m, 1, rhsbuf, m, SLU_DN, SLU_C,SLU_GE);
    cCreate_Dense_Matrix (&X, n, 1, xbuf, n, SLU_DN, SLU_C, SLU_GE);

    engine->Solve (&B, &X);
    
    Destroy_SuperMatrix_Store (&B);
    Destroy_SuperMatrix_Store (&X);
}
