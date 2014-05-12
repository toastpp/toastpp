/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "dsp_defs.h"
#include "util.h"
#include "Cnames.h"

int
c_bridge_dgssv_(int *n, int *nnz, int *nrhs, double *values, int *rowind,
		int *colptr, double *b, int *ldb, int *info)

{
    SuperMatrix A, B, L, U;
    SCformat *Lstore;
    NCformat *Ustore;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      panel_size, permc_spec, i;
    mem_usage_t   mem_usage;

    /* Adjust to 0-based indexing */
    for (i = 0; i < *nnz; ++i) --rowind[i];
    for (i = 0; i <= *n; ++i) --colptr[i];

    dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr, NC, _D, GE);
    dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, DN, _D, GE);

    if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     */    	
    permc_spec = 0;
    get_perm_c(permc_spec, &A, perm_c);

    panel_size = sp_ienv(1);
    
    dgssv(&A, perm_c, perm_r, &L, &U, &B, info);
    
    if ( *info == 0 ) {

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	
	dQuerySpace(&L, &U, panel_size, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
	
    } else {
	printf("dgssv() error returns INFO= %d\n", *info);
	if ( *info <= *n ) { /* factorization completes */
	    dQuerySpace(&L, &U, panel_size, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		   mem_usage.expansions);
	}
    }

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    /* Restore to 1-based indexing */
    for (i = 0; i < *nnz; ++i) ++rowind[i];
    for (i = 0; i <= *n; ++i) ++colptr[i];

}

