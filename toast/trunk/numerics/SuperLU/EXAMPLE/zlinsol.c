

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include "zsp_defs.h"
#include "util.h"

main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    doublecomplex   *a;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, panel_size, m, n, nnz, permc_spec;
    char     trans[1];
    doublecomplex   *xact, *rhs;
    mem_usage_t   mem_usage;

    nrhs   = 1;
    *trans = 'N';
    
    zreadhb(&m, &n, &nnz, &a, &asub, &xa);

    zCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, _Z, GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    if ( !(rhs = doublecomplexMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    zCreate_Dense_Matrix(&B, m, nrhs, rhs, m, DN, _Z, GE);
    xact = doublecomplexMalloc(n * nrhs);
    ldx = n;
    zGenXtrue(n, nrhs, xact, ldx);
    zFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree on structure of A'*A
     *   permc_spec = 2: minimum degree on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec = 0;
    get_perm_c(permc_spec, &A, perm_c);

    panel_size = sp_ienv(1);
    
    zgssv(&A, perm_c, perm_r, &L, &U, &B, &info);
    
    if ( info == 0 ) {

	zinf_norm_error(nrhs, &B, xact); /* Inf. norm of the error */

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	
	zQuerySpace(&L, &U, panel_size, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
	
    } else {
	printf("zgssv() error returns INFO= %d\n", info);
	if ( info <= n ) { /* factorization completes */
	    zQuerySpace(&L, &U, panel_size, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		   mem_usage.expansions);
	}
    }

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
}

