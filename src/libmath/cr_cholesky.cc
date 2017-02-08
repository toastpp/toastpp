// ============================================================================
// Cholesky factorisation routines for SCR (sparse compressed row) matrices
// ============================================================================

#define MATHLIB_IMPLEMENTATION
#include <stdio.h>
#include <string.h>
#include "mathlib.h"

// ============================================================================
// col_parent
// Given a column `col' with `len' entries and row index list `rowidx',
// return the column parent, defined as the smallest row index > col
// Assumes sorted rowidx

int col_parent (int col, int *rowidx, int len)
{
    for (int r = 0; r < len; r++)
        if (rowidx[r] > col) return rowidx[r];
    return -1; // no parent found
}

// ============================================================================
// merge_cols
// Merge row indices `rowidx1' (of length `len1') and `rowidx2' (of length
// `len2') and write the resulting list into buffer `buf' (of length `buflen').
// The first entry in rowidx1 (corresponding to the diagonal index of its
// parent) is ignored. Return value is the length of the merged list or -1 if
// buflen is exceeded. Assumes sorted row indices, and returns a sorted row
// index

int merge_cols (int *rowidx1, int len1, int *rowidx2, int len2,
		int *buf, int buflen)
{
    int r1, r2, i1 = 1, i2 = 0, n = 0;
    if (i1 < len1) r1 = rowidx1[i1];
    if (i2 < len2) r2 = rowidx2[i2];
    while (i1 < len1 && i2 < len2) {
        if (r1 < r2) {
	    if (n == buflen) return -1;
	    buf[n++] = r1;
	    r1 = rowidx1[++i1];
	} else if (r1 > r2) {
	    if (n == buflen) return -1;
	    buf[n++] = r2;
	    r2 = rowidx2[++i2];
	} else {
	    if (n == buflen) return -1;
	    buf[n++] = r1;
	    r1 = rowidx1[++i1];
	    r2 = rowidx2[++i2];
	}
    }
    while (i1 < len1) {
        if (n == buflen) return -1;
	buf[n++] = r1;
	r1 = rowidx1[++i1];
    }
    while (i2 < len2) {
        if (n == buflen) return -1;
	buf[n++] = r2;
	r2 = rowidx2[++i2];
    }
    return n;
}

// ============================================================================
// symbolic_cholesky_factor
// Given a matrix of dimension dim x dim, with structure rowptr and colidx
// generate the structure of the Cholesky factor and return it in
// frowptr and fcolidx. Return value is the number of entries in the factor.

MATHLIB int symbolic_cholesky_factor (int dim, idxtype *rowptr, idxtype *colidx,
    idxtype *&frowptr, idxtype *&fcolidx)
{
    const int fill_factor = 10; // expected fill-in
    int nval = rowptr[dim];
    int i, c, r, ci, ri, n, parent, newlen;

    const char *swapname = "/tmp/chswap.tmp";
    FILE *swap = 0; // column swap file
    int *colbuf, swap_next = 0;

    // create work buffer
    int buflen = nval*fill_factor;
    int *buf = new int[buflen];
    int bufalloc = 0;

    // create column entry count list
    int *nval_col = new int[dim];
    for (c = 0; c < dim; c++) nval_col[c] = 0;
    for (r = 0; r < dim; r++)
        for (ri = rowptr[r]; ri < rowptr[r+1]; ri++) {
	    c = colidx[ri];
	    if (r > c) nval_col[c]++; // only store lower triangle
	}

    // create column pointer array
    int **colref = new int*[dim];
    for (c = bufalloc = 0; c < dim; c++) {
        colref[c] = buf+bufalloc;
	bufalloc += nval_col[c];
    }

    // init row index array
    for (c = 0; c < dim; c++) nval_col[c] = 0;
    for (r = 0; r < dim; r++) {
        for (ri = rowptr[r]; ri < rowptr[r+1]; ri++) {
	    c = colidx[ri];
	    if (r > c) {
	        *(colref[c]+nval_col[c]) = r;
		nval_col[c]++;
	    }
	}
    }

    // sort row index array
    for (c = 0; c < dim; c++) {
        int *cptr = colref[c];
	int i, j, r, n = nval_col[c], inc = 1;
	do { inc *= 3; inc++; } while (inc < n);
	do {
	    inc /= 3;
	    for (i = inc; i < n; i++) {
	        r = cptr[i];
		j = i;
		while (cptr[j-inc] > r) {
		    cptr[j] = cptr[j-inc];
		    j -= inc;
		    if (j < inc) break;
		}
		cptr[j] = r;
	    }
	} while (inc > 1);
    }

    for (c = 0; c < dim; c++) {
        if ((parent = col_parent (c, colref[c], nval_col[c])) < 0)
	    continue;
	do {
	    newlen = merge_cols (colref[c], nval_col[c], colref[parent],
			     nval_col[parent], buf+bufalloc, buflen-bufalloc);
	    if (newlen < 0) { // insufficent space -> swap
	        if (!swap) swap = fopen (swapname, "wb");
		// swap all columns before current
		for (; swap_next < c; swap_next++) {
		    fwrite (colref[swap_next], sizeof(int),
			    nval_col[swap_next], swap);
		    colref[swap_next] = 0;
		}
		// defrag the rest. This should really be done in-place
		int *nbuf = new int[buflen];
		int nbufalloc = 0;
		for (i = c; i < dim; i++) {
		   memcpy (nbuf+nbufalloc, colref[i], nval_col[i]*sizeof(int));
		   colref[i] = nbuf+nbufalloc;
		   nbufalloc += nval_col[i];
		}
		delete []buf;
		buf = nbuf;
		bufalloc = nbufalloc;
	    }
	} while (newlen < 0);
	colref[parent] = buf+bufalloc;
        nval_col[parent] = newlen;
	bufalloc += newlen;
    }

    if (swap) {
        fclose (swap);
	swap = fopen (swapname, "rb");
	colbuf = new int[dim];
    }

    // write results back into row storage
    nval = 0;
    for (c = 0; c < dim; c++) nval += nval_col[c];

    // row entry count
    int *nval_row = new int[dim];
    for (r = 0; r < dim; r++) nval_row[r] = 0;
    for (c = 0; c < dim; c++) {
        int *colp;
        if (c < swap_next) {
	    n = fread (colbuf, sizeof(int), nval_col[c], swap);
	    xASSERT(n == nval_col[c],
		    "Unexpected number of elements read from stream");
	    colp = colbuf;
	} else {
	    colp = colref[c];
	}
        for (ci = 0; ci < nval_col[c]; ci++)
	    nval_row[colp[ci]]++;
    }

    // init row pointer array
    frowptr = new idxtype[dim+1];
    for (r = 0, frowptr[0] = 0; r < dim; r++)
        frowptr[r+1] = frowptr[r] + nval_row[r];

    if (swap) rewind (swap);
    // init column index array
    fcolidx = new idxtype[nval];
    for (r = 0; r < dim; r++) nval_row[r] = 0;
    for (c = 0; c < dim; c++) {
        int *colp;
	if (c < swap_next) {
	    n = fread (colbuf, sizeof(int), nval_col[c], swap);
	    xASSERT(n == nval_col[c],
		    "Unexpected number of elements read from stream");
	    colp = colbuf;
	} else {
	    colp = colref[c];
	}
        for (ci = 0; ci < nval_col[c]; ci++) {
	    r = colp[ci];
	    fcolidx[frowptr[r]+nval_row[r]] = c;
	    nval_row[r]++;
	}
    }

    // cleanup
    delete []nval_row;
    delete []nval_col;
    delete []colref;
    delete []buf;
    if (swap) {
	delete []colbuf;
	fclose (swap);
	remove (swapname);
    }
    return nval;
}
