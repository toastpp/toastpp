// ============================================================================
// meshopt
// Mesh optimisation routines
// ============================================================================

#define FELIB_IMPLEMENTATION

#include "felib.h"
#include "ordmmd.h"
#include "rcm.h"

using namespace std;

// ============================================================================
// Local prototypes

int tinney_ (int *l, int *nl, int *fl, int *il, int *jl, int *eltop,
	     int ieltop, int jeltop, int *totels, int *nodsort, int *newnod,
	     int *totnod, int tsz);

// ============================================================================
// Sort boundary nodes to the end of the node list
// This only reorders the permutation vector `perm', not the node list itself
// Returns the number of boundary nodes found

FELIB int SortBndToEnd (Mesh &mesh, int *perm)
{
    int src, dst, tmp, nbnd = 0;
    NodeList &nlist = mesh.nlist;

    for (src = dst = mesh.nlen()-1; src >= 0; src--)
        if (nlist[perm[src]].isBnd()) {
	    tmp = perm[src], perm[src] = perm[dst], perm[dst] = tmp;
	    dst--;
	    nbnd++;
	}
    return nbnd;
}

// ============================================================================
// Minimum bandwidth optimisation
// return value is error flag
//   0 = success

FELIB int Optimise_MinBandwidth (Mesh &mesh, int *perm, int ofs, int len)
{
    int i, nzero, nfirst, maskval, nds = mesh.nlen();
	idxtype *rowptr, *colidx;
    int iperm0 = 0;
    int *mask = new int[nds];
    int *riord = new int[nds];
    int levels[512], nlev, init;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    
    // adjust lists to 1-based, required by fortran
    for (i = 0; i <= nds; i++) {
        rowptr[i]++;
	mask[i] = 1;
    }
    for (i = 0; i < nzero; i++)
        colidx[i]++;
    nfirst = 1;
    riord[0] = 1;
    //for (i = 0; i < nds; i++) {
    //   if (mesh.nlist[i].isBnd()) {
    //	    riord[1] = i+1;
    //	    break;
    //	}
    //}
    maskval = 1;
    init = 1;
    perphn_ (&nds, colidx, rowptr, &init, &iperm0, mask, &maskval, &nlev,
    	     riord, levels);
    //bfs_ (&nds, colidx, rowptr, &nfirst, &iperm0, mask, &maskval, riord,
    //	  levels, &nlev);
    reversp_ (&nds, riord); // not sure whether this is necessary
    for (i = 0; i < nds; i++) perm[i] = riord[i]-1;
    delete []mask;
    // delete []riord;
    // Deleting riord causes a segmentation fault for a reason I can't
    // imagine. Probably perphn does something nasty to it. For the moment
    // I don't delete riord, but this is sub-optimal!
    return 0;
}

// ============================================================================
// Minimum degree optimisation
// return value is error flag:
//   0 = success
//  -1 = insufficient working space (should never happen)

FELIB int Optimise_MMD (Mesh &mesh, int *perm, int ofs, int len)
{
    int i, j, k, idx, iflag, nofsub, nds = mesh.nlen();
    int iwsiz = 4*nds, *iwork = new int[iwsiz]; // work space
    idxtype *rowptr, *colidx;
	int nzero;

    mesh.SparseRowStructure (rowptr, colidx, nzero);

    // eliminate diagonal (note: this may not be necessary)
    for (i = idx = 0; i < nds; i++)
	for (j = rowptr[i]; j < rowptr[i+1]; j++)
	    if ((k = colidx[j]) != i) colidx[idx++] = k;
    for (i = 1; i <= nds; i++) rowptr[i] -= i;

    for (i = 0; i < nzero; i++) colidx[i]++; // convert to 1-based
    for (i = 0; i <= nds; i++)  rowptr[i]++; // convert to 1-based

    int *iperm = new int[nds];
    for (i = 0; i < nds; i++) {
        iperm[perm[i]] = i+1;
        perm[i]++; // ordmmd expects 1-based
    }
    ordmmd_ (&nds, rowptr, colidx, iperm, perm, &iwsiz, iwork, &nofsub,
	     &iflag);
    for (i = 0; i < nds; i++) perm[i]--; // back to 0-based

    delete []iperm;
    delete []iwork;
    delete []rowptr;
    delete []colidx;

    return iflag;
}

// ============================================================================
// Tinney 2 optimisation scheme
// return value is error flag:
//   0 = success
//   1 = mesh contains elements with different number of nodes
//   2 = fatal error in tinney subroutine

FELIB int Optimise_Tinney2 (Mesh &mesh, int *perm, int ofs, int len)
{
    int nds = mesh.nlen();
    int els = mesh.elen();
    int i, j, res, nnd = mesh.elist[0]->nNode();
    for (i = 1; i < els; i++)
        if (mesh.elist[i]->nNode() != nnd) return 1;
    // we only allow meshes with elements that have the same number of nodes
    int *l, tsz, il, jl = nds;
    int *nl = new int[jl];
    int *fl = new int[jl];
    int ieltop = els;
    int jeltop = nnd;
    int *eltop = new int[ieltop*jeltop];

    for (i = 0; i < els; i++)
        for (j = 0; j < nnd; j++)
	    eltop[i+ieltop*j] = mesh.elist[i]->Node[j];

    int *iperm = new int[nds];
    for (i = 0; i < nds; i++) iperm[perm[i]] = i;

    il = 131072;
    tsz = 1024;
    do {
        l = new int[il];
        res = tinney_ (l, nl, fl, &il, &jl, eltop, ieltop, jeltop, &els,
		       perm, iperm, &nds, tsz);
	delete []l;
	switch (res) {
	case 1:
	    il *= 2;
	    cerr << "Insufficient l-buffer. Increasing to " << il
		 << " and retrying ..." << endl;
	    break;
	case 2:
	    tsz *= 2;
	    cerr << "Insufficient t-buffer. Increasing to " << tsz
		 << " and retrying ..." << endl;
	    break;
	case 3:
	    cerr << "Panic! Unrecoverable error in Tinney. Terminating.\n";
	    return 2;
	}
    } while (res != 0);

    delete []nl;
    delete []fl;
    delete []eltop;
    delete []iperm;
    return 0;
}


// ============================================================================
// ============================================================================
// Tinney2 optimisation routine. f2c'd from Fortran original

/* Builtin functions */
inline double d_int(double *x) {return *x;};

/* Subroutine */ 
int symadd_(int *nl, int *l, int *il, int nnode, bool *new__)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;


/*     INSERT SYMBOL IN COLUMN NNODE OF THIS ROW */


    /* Parameter adjustments */
    --l;

    /* Function Body */
    *new__ = true;
    i__1 = *nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (l[i__] == nnode) {
	    *new__ = false;
	}
/* L100: */
    }
    if (*new__) {
	++(*nl);
	l[*nl] = nnode;
    }

    return 0;
} /* symadd_ */


// ============================================================================
// WARNING: Not completely stable yet. Fails on a larger TET4 cubic mesh
// ============================================================================

/* Subroutine */ 
int tinney_ (int *l, int *nl, int *fl, int *il, int *jl, int *eltop,
	     int ieltop, int jeltop, int *totels, int *nodsort, int *newnod,
	     int *totnod, int tsz)
{
    // System generated locals
    int i__1, i__2, i__3;

    // Local variables
    static int iele, lcol;
    static double rcol;
    static int temp, i__, j, k, inode, jnode, knode, mnode, found,
        mnval, mxval, ti, nt;
    static double sumcol;
    static int fli;
    static bool new__;

    int result = 0;

    // static int t[200];
    // `200' was in the original but is not sufficient for large meshes
    int *t = new int[tsz];

    // static int valence[5000];
    // `5000' was in the original but must be >= number of nodes!
    int *valence = new int[*totnod];  // M.S. 24.9.99

    static int sum, lcolold, nxtsift;

    // Fortran I/O blocks
/*
    static cilist io___162 = { 0, 6, 0, 0, 0 };
    static cilist io___164 = { 0, 6, 0, 0, 0 };
    static cilist io___166 = { 0, 6, 0, 0, 0 };
    static cilist io___167 = { 0, 6, 0, 0, 0 };
    static cilist io___168 = { 0, 6, 0, 0, 0 };
    static cilist io___179 = { 0, 6, 0, 0, 0 };
    static cilist io___180 = { 0, 6, 0, 0, 0 };
*/

    // CALCULATE THE TINNEY SCHEME 2 (WEBSTER PAGE 141) NODE RENUMBERING 
    // FOR MINIMUM FILL-IN DURING CHOLESKI FACTORISATION

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    //   NOW THE REAL WORK, CALCULATE THE TINNEY SCHEME 2 RENUMBERING
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    //     INITIALISE

    // Parameter adjustments
    //--l;
    //--fl;
    //--nl;
    //--newnod;
    //--nodsort;

    // Function Body
    rcol = (double) (*il) / (double) (*totnod);
    lcol = (int) d_int(&rcol);
    sumcol = 0.0;
    fli = 0;
    i__1 = *totnod;
    for (inode = 0; inode < i__1; ++inode) {
	fl[inode] = fli;
	sumcol += rcol;
	fli = (int) d_int(&sumcol);
	nl[inode] = 0;
	valence[inode] = 0;
    }
    // BUILD SYMBOLIC SYSTEM STIFFNESS MATRIX
    cerr << "Building symbolic system matrix.\n-> rcol " << rcol << ", lcol "
	 << lcol  << endl;
    i__1 = *totels;
    for (iele = 0; iele < i__1; ++iele) {
	for (i__ = 0; i__ < jeltop; ++i__) {
	    inode = newnod[eltop[iele + i__ * ieltop]];
	    for (j = i__; j < jeltop; ++j) {
		jnode = newnod[eltop[iele + j * ieltop]];

		if (jnode > inode) {
		    symadd_(&nl[inode], &l[fl[inode]], &lcol, jnode+1,
			    &new__);
		} else {
		    symadd_(&nl[jnode], &l[fl[jnode]], &lcol, inode+1,
			    &new__);
		}
		if (new__) {
		    if (inode == jnode) {
			++valence[inode];
		    } else {
			++valence[inode];
			++valence[jnode];
		    }
		}
/* L280: */
	    }
/* L290: */
	}
/* L300: */
    }
    // CALCULATE THE SPARSITY

    sum = 0;
    i__1 = *totnod;
    for (i__ = 0; i__ < i__1; ++i__)
	sum += valence[i__];
    cerr << "Non-zeros before factorisation " << sum  << endl;
/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/*     FOR EACH ROW(COLUMN) OF THE MATRIX: */
/*       RE-ORDER SO THAT THE ROW(COLUMN) WITH THE SMALLEST VALENCE */
/*         IS THE NEXT ROW(COLUMN) */
/*       THEN PERFORM ONE STEP IN THE SYMBOLIC FACTORISATION */

    nxtsift = *totnod / 4;
    i__1 = *totnod;
    for (inode = 0; inode < i__1; ++inode) {

/*       IF IT IS TIME PERFORM A SIFT OF THE MATRIX */

	if (inode+1 == nxtsift) {
/*
	    s_wsle(&io___164);
	    do_lio(&c__9, &c__1, "SIFTING AT ", 11L);
	    do_lio(&c__3, &c__1, (char *)&inode, (ftnlen)sizeof(int));
	    e_wsle();
*/
	    lcolold = lcol;
	    rcol = (double) (*il) / (double) (*totnod - nxtsift + 1);
	    lcol = (int) d_int(&rcol);
	    sumcol = (float)0.;
	    fli = (int) sumcol;
	    i__2 = *totnod;
	    for (i__ = inode; i__ < i__2; ++i__) {
		if (nl[i__] >= lcolold) {
		    result = 1; // il too small. retry
		    goto done;
		}
		if (fli > fl[i__]) {
		    result = 3; // Panic. abort
		    goto done;
		}
		i__3 = nl[i__];
		for (j = 0; j < i__3; ++j) {
		    l[fli+j] = l[fl[i__]+j];
/* L380: */
		}
		fl[i__] = fli;
		sumcol += rcol;
		fli = (int) d_int(&sumcol);
/* L390: */
	    }
	    if (lcol < *totnod - (inode+1)) {
		nxtsift += (*totnod - nxtsift) / 4;
	    }
	}

	// FIND THE ROW(COLUMN) WITH SMALLEST VALENCE

	mnode = inode;
	mnval = *totnod;
	i__2 = *totnod;
	for (i__ = inode; i__ < i__2; ++i__) {
	    if (valence[i__] < mnval) {
		mnode = i__;
		mnval = valence[i__];
	    }
/* L400: */
	}

	// SWAP NODE NUMBERS INODE AND MNODE

	if (inode != mnode) {

	    temp = nodsort[inode];
	    nodsort[inode] = nodsort[mnode];
	    nodsort[mnode] = temp;
/*         ... SWAP IN NEW NODE NUMBER ARRAY */

	    temp = valence[inode];
	    valence[inode] = valence[mnode];
	    valence[mnode] = temp;
/*         ... SWAP COLUMN VALENCIES */

/*         THIS IS A VERY MESSY OPERATION! WE NEED TO DO A COLUMN AND ROW */
/*         SWAP ON THE MATRIX STORED AS A SPARSE, LOWER TRIANGULER PART OF */
/*         A SYMETRIC MATRIX.  IT IS TRICKY BUT FAST.  THE ONLY HELP I CAN */
/*         GIVE IS TO URGE UTMOST CONCENTRATION. */

	    nt = nl[inode];
	    if (nt > tsz) {
		result = 2; // t-buffer too small. retry
		goto done;
	    }
	    nl[inode] = 0;
	    i__2 = nt;
	    for (i__ = 0; i__ < i__2; ++i__) {
		ti = l[i__ + fl[inode]];
		t[i__] = ti;
/* L600: */
	    }
/*         ...  MOVE COLUMN INODE INTO TEMPORARY STORAGE */

	    i__2 = mnode;
	    for (i__ = inode+1; i__ < i__2; ++i__) {
		found = 0;
		i__3 = nl[i__];
		for (j = 0; j < i__3; ++j) {
		    if (l[j + fl[i__]] == mnode+1) {
			++nl[inode];
			l[nl[inode] + fl[inode]-1] = i__+1;
			found = j+1;
		    }
/* L680: */
		}
		if (found > 0) {
		    i__3 = nl[i__] - 1;
		    for (j = found; j <= i__3; ++j) {
			l[j + fl[i__]-1] = l[j + fl[i__]];
/* L690: */
		    }
		    --nl[i__];
		}
/* L700: */
	    }
/*         ... FIX THE COLUMNS BETWEEN INODE AND MNODE */

	    i__2 = nl[mnode];
	    for (i__ = 0; i__ < i__2; ++i__) {
		++nl[inode];
		if (l[i__ + fl[mnode]] == mnode+1) {
		    l[nl[inode] + fl[inode]-1] = inode+1;
		} else {
		    l[nl[inode] + fl[inode]-1] = l[i__ + fl[mnode]];
		}
/* L800: */
	    }
	    nl[mnode] = 0;
/*         ... FIX COLUMN MNODE */

	    i__2 = nt;
	    for (i__ = 0; i__ < i__2; ++i__) {
		ti = t[i__];
		if (ti == inode+1) {
		    ++nl[mnode];
		    l[nl[mnode] + fl[mnode]-1] = mnode+1;
		} else if (ti < mnode+1) {
		    ++nl[ti-1];
		    l[nl[ti-1] + fl[ti-1]-1] = mnode+1;
		} else if (ti > mnode+1) {
		    ++nl[mnode];
		    l[nl[mnode] + fl[mnode]-1] = ti;
		} else {
		    ++nl[inode];
		    l[nl[inode] + fl[inode]-1] = mnode+1;
		}
/* L900: */
	    }
/*         ... FIX THE OLD COLUMN INODE */

	}

/*       PERFORM THE NEXT STEP IN THE SYMBOLIC FACTORISATION */

	i__2 = nl[inode];
	for (j = 0; j < i__2; ++j) {
	    jnode = l[j + fl[inode]]-1;
	    if (jnode > inode) {
		i__3 = nl[inode];
		for (k = 0; k < i__3; ++k) {
		    knode = l[k + fl[inode]]-1;
		    if (knode > jnode) {

			symadd_(&nl[jnode], &l[fl[jnode]], &lcol, knode+1, 
				&new__);
			if (new__) {
			    if (jnode == knode) {
				++valence[jnode];
			    } else {
				++valence[jnode];
				++valence[knode];
			    }
			}
		    }
/* L990: */
		}
	    }
/* L1000: */
	}
/* L2000: */
    }

/*     CALCULATE THE FILL-IN */

    mxval = 0;
    sum = 0;
    i__1 = *totnod;
    for (i__ = 0; i__ < i__1; ++i__) {
	sum += valence[i__];
	if (nl[i__] > mxval)
	    mxval = nl[i__];
/* L3000: */
    }
    cerr <<  "Non-zeros after  factorisation " << sum << endl;
    cerr <<  "Maximum number of non-zeros in factor is " << mxval << endl;

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

    // Cleanup

 done:
    delete []t;
    delete []valence; // M.S. 24.9.99
    return result;

} /* tinney_ */

