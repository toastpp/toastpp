// ==========================================================================
// Module mathlib
// File crmatrix.cc
// Definition of template class TCompRowMatrix ('template compressed-row
//  matrix')
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#define __CRMATRIX_CC

#include "mathlib.h"
#include "slu_zdefs.h"  // SuperLU 4
#include "ilutoast.h"

#ifdef USE_CUDA_FLOAT
#include "toastcuda.h"
#include "toastspmv.h"
#endif

#ifdef ML_INTERFACE
#include "ml_defs.h"
#include "ml_operator.h"
#endif // ML_INTERFACE

using namespace std;

template<>
void LU (TCompRowMatrix<std::complex<double> > &A,
    const TVector<std::complex<double> > &b,
    TVector<std::complex<double> > &x)
{
    int n = b.Dim();
    char equed  = 'N';
    int *perm_c = new int[n];
    int *perm_r = new int[n];
    int *etree  = new int[n];
    double R = 0.0;
    double C = 0.0;
    double ferr, berr;
    double recip_pivot_growth, rcond;
    int info;
    mem_usage_t mem_usage;
    SuperMatrix smA, L, U, B, X;
    superlu_options_t options;
    SuperLUStat_t stat;
    GlobalLU_t Glu;
    doublecomplex *cdat = (doublecomplex*)A.ValPtr();
    zCreate_CompCol_Matrix (&smA, A.nRows(), A.nCols(), A.nVal(),
			    cdat, A.colidx, A.rowptr, SLU_NR, SLU_Z, SLU_GE);

    doublecomplex *bdat = (doublecomplex*)b.data_buffer();
    doublecomplex *xdat = (doublecomplex*)x.data_buffer();
    zCreate_Dense_Matrix (&B, n, 1, bdat, n, SLU_DN, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix (&X, n, 1, xdat, n, SLU_DN, SLU_Z, SLU_GE);

    get_perm_c (0, &smA, perm_c);

    set_default_options(&options);
    options.Fact = DOFACT;
    options.Equil = NO;
    options.ColPerm = NATURAL;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;

    StatInit (&stat);

    zgssvx (&options, &smA, perm_c, perm_r, etree, &equed, &R, &C,
	    &L, &U, 0, 0, &B, &X, &recip_pivot_growth, &rcond,
	    &ferr, &berr, &Glu, &mem_usage, &stat, &info);

    StatFree (&stat);

    //toast_zgssvx (&fact, &trans, &refact, &smA, 0, perm_c, perm_r, etree,
    //    &equed, &R, &C, &L, &U, 0, 0, &B, &X,
    //    &recip_pivot_growth, &rcond, &ferr, &berr, &mem_usage, &info);

    Destroy_SuperMatrix_Store (&B);
    Destroy_SuperMatrix_Store (&X);
    delete []perm_c;
    delete []perm_r;
    delete []etree;
    Destroy_SuperMatrix_Store (&smA);
    Destroy_SuperNode_Matrix (&L);
    Destroy_CompCol_Matrix (&U);
}

