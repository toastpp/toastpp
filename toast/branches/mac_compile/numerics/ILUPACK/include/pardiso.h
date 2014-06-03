/* ============================================================================================ *
 *												*
 *  smat.h -- Sparse matrix routines for csr matrices.						*
 *												*
 *  Author: 	Michael Hagemann								*
 *  Date:	2003-06-21 20:21								*
 *												*
 *  Copyright (c) 2002, 2003, 2004, Michael Hagemann, Olaf Schenk.						*
 *												*
 * -------------------------------------------------------------------------------------------- *
 *												*
 * ============================================================================================ */

		/* $Id: pardiso.h,v 1.2 2005/06/30 07:25:11 cchcak Exp $ */


/* @@LICENCE@@ */


/* -------------------------------------------------------------------------------------------- */

#ifndef __SMAT_H
#define __SMAT_H

#ifdef __cplusplus
extern "C" {
#endif


/* -------------------------------------------------------------------------------------------- */

#include "long_integer.h"
#include "namespardiso.h"

#define pint integer

// not really correct
#define perm_t integer

	/**
	 *  Compressed sparse row matrix structure.
	 */
	struct smat_struct {
		unsigned integer     m;		/*!< Dimension rows. */
		unsigned integer     n;		/*!< Dimension columns. */
		unsigned integer     nnz;	/*!< Number of nonzero elements.  Equal to
					 * <tt>ia[m+1]</tt>. */

		integer	     sym;	/*!< symmetric? */

		pint	    *ia;	/*!< Index to rows in a and ja.  Has n+1 entries (length
					 * of last row). */
		pint	    *ja;	/*!< Column of each value in a. */
		double	    *a;		/*!< Values. */
	};

	/**
	 *  Compressed sparse row matrix type.
	 *
	 *  @see smat_struct
	 */
	typedef struct smat_struct smat_t;


	/**
	 *  Structure to allow column access of matrix.
	 */
	typedef struct {
		smat_t	    *mat;	/* Corresponding matrix. */
		unsigned integer    *col_start;	/* Index of first pointer to i-th column in col_ptrs. */
		unsigned integer    *row_ind;	/* Row of each element pointed to. */
		unsigned integer    *row_ptr;	/* Pointer to ia, ja in matrix. */
	} col_index_t;


	/* ------------------------------------------------------------------------------------ */

	/**
	 *  Contructs new matrix.
	 *
	 *  Memory is allocated for the row indices only!
	 *
	 *  @param	m		Row dimension
	 *  @param	n		Column dimension
	 *  @param	is_symmetric	Is the matrix symmetric?
	 *
	 *  @returns	A new matrix structure.
	 */
	smat_t	* smat_new		(unsigned integer	 m,
					 unsigned integer	 n,
					 integer		 is_symmetric);


	/**
	 *  Contructs new matrix and allocates all necessary memory.
	 *
	 *  @param	m		Row dimension
	 *  @param	n		Column dimension
	 *  @param	nnz		Designated number of non-zero elements.  See 
	 *				{@link smat_realloc smat_realloc} for later reallocation.
	 *  @param	is_symmetric	Is the matrix symmetric?
	 *
	 *  @returns	A new matrix structure.
	 */
	smat_t	* smat_new_nnz		(unsigned integer	 m,
					 unsigned integer	 n,
					 unsigned integer	 nnz,
					 integer		 is_symmetric);


	/**
	 *  Reallocates memory for the values and column indices according to <tt>A->nnz</tt>.
	 */
	void	  smat_realloc		(smat_t		*A,
					 unsigned integer	new_nnz);


	/**
	 *  Contructs a matrix from an existing csr structure.  No copies will
	 *  be made.
	 *
	 *  <P> <B>Warning:</B> If the base has to be converted, the original
	 *  structure will be altered!
	 *
	 *  @param	m		Row dimension
	 *  @param	n		Column dimension
	 *  @param	ia		Indices to Rows
	 *  @param	ja		Column indices
	 *  @param	a		Values
	 *  @param	is_symmetric	Is the matrix symmetric?
	 *  @param	base_index	Is the matrix 0-based (C) or 1-based (Fortran)?
	 *
	 *  @returns	A new matrix structure.
	 */
	smat_t	*smat_new_from		(unsigned integer	 m,
					 unsigned integer	 n,
					 pint		*ia,
					 pint		*ja,
					 double		*a,
					 integer		 is_symmetric,
					 integer		 base_index);

	/**
	 *  Deallocates the matrix structure and all arrays.
	 *
	 *  @param	A	Matrix to liberate.
	 */
	void	 smat_free		(smat_t		*A);


	/**
	 *  Copies a matrix.
	 *
	 *  @param	src_mat		Source Matrix.
	 *
	 *  @returns	A deep copy of src_mat.  Everything is copied.
	 */
	smat_t	*smat_copy		(smat_t		*src_mat);


	/**
	 *  Copies the structure of a matrix.
	 *
	 *  @param	src		Source Matrix.
	 *
	 *  @returns	Copy of structure of matrix src.  Values are not allocated.
	 */
	smat_t	*smat_copy_structure	(smat_t		*src);


	/**
	 *  Convert row and column indices to a 1-based scheme.  This is handy if you call
	 *  Fortran code.  NOTE: 1-based indices are not allowed in sagg functions.  Convert the
	 *  matrix back to C-indexing after the call to Fortran.
	 *
	 *  @param	A		Source Matrix.
	 *
	 *  @returns	Nothing, matrix is changed in-place.
	 */
	void	 smat_to_fortran_indexing (smat_t	*mat);


	/**
	 *  Converts a matrix from Fortran 1-based indexing to C 0-based indexing.
	 *
	 *  @param	A		Source Matrix.
	 *
	 *  @returns	Nothing, matrix is changed in-place.
	 */
	void	 smat_to_c_indexing	(smat_t		*mat);


/* -------------------------------------------------------------------------------------------- */

	/**
	 *  Various consistency checks for the matrix structure.
	 *
	 *  @param	A		Matrix.
	 *
	 *  @returns	0, if one of the tests failed, 1 otherwise.  Reports
	 *		error messages to console.
	 */
	integer	 smat_is_consistent	(smat_t		*A);

	/**
	 *  Check that row indices have a legal pattern (start from zero, increase steadily).
	 */
	integer	 smat_check_ia		(smat_t		*A);

	/**
	 *  Check that column indices have a legal pattern (i.e. they are bounded by zero and
	 *  n).
	 */
	integer	 smat_check_ja		(smat_t		*A);

	/**
	 *  Check that structurally symmetric matrices are really symmetric.
	 */
	integer	 smat_check_symmetry	(smat_t		*A);

	/**
	 *  Check that diagonal values exist.
	 */
	integer	 smat_check_diagonal	(smat_t		*A);


	/**
	 *  Print internals of matrix (a, ia, ja, ...).  Long lists are abbreviated by ellipses.
	 */
	void	 smat_print_raw		(smat_t		*A);


/* -------------------------------------------------------------------------------------------- */
//#define __PARDISO_2__
#ifdef __PARDISO_2__
        int mps_pardiso(int n, int * ia, int *ja, double *a, int *p, double *rowscal, double *colsal, int i);
#else
        int mps_pardiso(int contraint1, int contraint2, int n, int *ia, int *ja, double *a, int *p, double *rowscal, double *colsal, int i);
#endif
	smat_t	*smat_copy_permute		(smat_t		*src,
						 perm_t		*perm,
						 integer		 copy_vals);


	void	 smat_reordering_metis_ddist	(smat_t		*A,
						 perm_t		*perm,
						 integer		 nproc,
						 perm_t		*ddist);

	void	 smat_reordering_amd		(smat_t		*A,
						 perm_t		*perm);

	void	 smat_reordering_gepmmd		(smat_t		*A,
						 perm_t		*perm);

	void	 smat_reordering_gepcmd		(smat_t		*A,
						 perm_t		*perm,
						 pint		 nparts);


/* -------------------------------------------------------------------------------------------- */

	/**
	 *  Create a column index of the matrix.  
	 *
	 *  You can traverse a column index along the lines of:
	 *
	 *  <pre>
	 *  for (j = ci->col_start[col];
	 *       j < ci->col_start[col+1]; j++)
	 *  {
	 *  	integer   row = ci->row_ind[j];
	 *	real  val = A->a[ ci->row_ptr[j] ];
	 *	printf ("row: %d, col: %d, value: %e \n", row, col, val);
	 *  }
	 *  </pre>
	 *
	 *  @see col_index_t
	 *
	 *  @returns	Column index.
	 */
	col_index_t *smat_col_index_new	(smat_t		*A);


	/**
	 *  Free a column index.
	 */
	void	 smat_col_index_free	(col_index_t	*ci);


/* -------------------------------------------------------------------------------------------- */

#ifdef __cplusplus
}
#endif

#endif  /* __SMAT_H */
