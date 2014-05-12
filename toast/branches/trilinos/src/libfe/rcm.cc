/* rcm.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <iostream>
#include <math.h>
#include "rcm.h"
#define abs(x) fabs(x)
#define min(x,y) ((x)<(y) ? (x):(y))
#define max(x,y) ((x)>(y) ? (x):(y))

/* Table of constant values */

static int c__1 = 1;

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*               ROERDERING ROUTINES -- REORD MODULE                    c */
/* ----------------------------------------------------------------------c */
/* BSF     : Breadth-First Search traversal (Cuthill Mc Kee ordering)   c */
/* dblstr  : two-way dissection partitioning -- with equal size domains c */
/* stripes : routine used by dblstr to assign points                    c */
/* perphn  : finds a peripheral node and does a BFS search from it.     c */
/* add_lvst: routine for adding a new level set in BFS algorithm        c */
/* reversp : routine to reverse a given permuation (e.g., for RCMK)     c */
/* maskdeg : int function to compute the `masked' of a node         c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int bfs_(int *n, int *ja, int *ia, int *
	nfirst, int *iperm, int *mask, int *maskval, int *
	riord, int *levels, int *nlev)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int iend;
    extern /* Subroutine */ int add_lvst__(int *, int *, int *, 
	    int *, int *, int *, int *, int *);
    static int j, ii, istart;
    static bool permut;
    static int nod;

/* -----------------------------------------------------------------------
 */
/* finds the level-structure (breadth-first-search or CMK) ordering for a 
*/
/* given sparse matrix. Uses add_lvst. Allows an set of nodes to be */
/* the initial level (instead of just one node). Allows masked nodes. */
/* -------------------------parameters------------------------------------
 */
/* on entry: */
/* ---------- */
/* n      = number of nodes in the graph */
/* ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 
*/
/*          structure) */
/* nfirst = number of nodes in the first level that is input in riord */
/* iperm  = int array indicating in which order to  traverse the graph
 */
/*          in order to generate all connected components. */
/*          The nodes will be traversed in order iperm(1),....,iperm(n) */
/*          Convention: */
/*          if iperm(1) .eq. 0 on entry then BFS will traverse the */
/*          nodes in the  order 1,2,...,n. */

/* riord  = (also an ouput argument). on entry riord contains the labels 
*/
/*          of the nfirst nodes that constitute the first level. */

/* mask   = array used to indicate whether or not a node should be */
/*          condidered in the graph. see maskval. */
/*          mask is also used as a marker of  visited nodes. */

/* maskval= consider node i only when:  mask(i) .eq. maskval */
/*          maskval must be .gt. 0. */
/*          thus, to consider all nodes, take mask(1:n) = 1. */
/*          maskval=1 (for example) */

/* on return */
/* --------- */
/* mask   = on return mask is restored to its initial state. */
/* riord  = `reverse permutation array'. Contains the labels of the nodes 
*/
/*          constituting all the levels found, from the first level to */
/*          the last. */
/* levels = pointer array for the level structure. If lev is a level */
/*          number, and k1=levels(lev),k2=levels(lev+1)-1, then */
/*          all the nodes of level number lev are: */
/*          riord(k1),riord(k1+1),...,riord(k2) */
/* nlev   = number of levels found */
/* -----------------------------------------------------------------------
 */
/* Notes on possible usage */
/* ------------------------- */
/* 1. if you want a CMK ordering from a known node, say node init then */
/*    call BFS with nfirst=1,iperm(1) =0, mask(1:n) =1, maskval =1, */
/*    riord(1) = init. */
/* 2. if you want the RCMK ordering and you have a preferred initial node 
*/
/*     then use above call followed by reversp(n,riord) */
/* 3. Similarly to 1, and 2, but you know a good LEVEL SET to start from 
*/
/*    (nfirst = number if nodes in the level, riord(1:nfirst) contains */
/*    the nodes. */
/* 4. If you do not know how to select a good initial node in 1 and 2, */
/*    then you should use perphn instead. */

/* -----------------------------------------------------------------------
 */
/*     local variables -- */
    /* Parameter adjustments */
    --mask;
    --iperm;
    --ja;
    --ia;
    --riord;
    --levels;

    /* Function Body */
    permut = iperm[1] != 0;

/*     start pointer structure to levels */

    *nlev = 0;

/*     previous end */

    istart = 0;
    ii = 0;

/*     current end */

    iend = *nfirst;

/*     intialize masks to zero -- except nodes of first level -- */

    i__1 = *nfirst;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = 0;
/* L12: */
    }
/* -----------------------------------------------------------------------
 */
/* L13: */

L1:
    ++(*nlev);
    levels[*nlev] = istart + 1;
    add_lvst__(&istart, &iend, nlev, &riord[1], &ja[1], &ia[1], &mask[1], 
	    maskval);
    if (istart < iend) {
	goto L1;
    }
L2:
    ++ii;
    if (ii <= *n) {
	nod = ii;
	if (permut) {
	    nod = iperm[nod];
	}
	if (mask[nod] == *maskval) {

/*     start a new level */

	    istart = iend;
	    ++iend;
	    riord[iend] = nod;
	    mask[nod] = 0;
	    goto L1;
	} else {
	    goto L2;
	}
    }
/* -----------------------------------------------------------------------
 */
/* L3: */
    levels[*nlev + 1] = iend + 1;
    i__1 = iend;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = *maskval;
    }
    return 0;
/* -----------------------------------------------------------------------
 */
/* -----end-of-BFS--------------------------------------------------------
 */
} /* bfs_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int dblstr_(int *n, int *ja, int *ia, int *
	ip1, int *ip2, int *nfirst, int *riord, int *ndom, 
	int *map, int *mapptr, int *mask, int *levels, 
	int *iwk)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int idom, jdom, kdom, init, nlev, j, k;
    extern /* Subroutine */ int perphn_(int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *);
    static int numnod;
    extern /* Subroutine */ int bfs_(int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *);
    static int maskval, nextdom, ndp1;
    extern /* Subroutine */ int stripes_(int *, int *, int *, 
	    int *, int *, int *, int *);

/* -----------------------------------------------------------------------
 */
/* this routine performs a two-way partitioning of a graph using */
/* level sets recursively. First a coarse set is found by a */
/* simple Cuthill-Mc Kee type algorithm. Them each of the large */
/* domains is further partitioned into subsets using the same */
/* technique. The ip1 and ip2 parameters indicate the desired number */
/* number of partitions 'in each direction'. So the total number of */
/* partitions on return ought to be equal (or close) to ip1*ip2 */
/*----------------------parameters---------------------------------------
-*/
/* on entry: */
/* --------- */
/* n      = row dimension of matrix == number of vertices in graph */
/* ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data 
*/
/*          structure) */
/* ip1    = int indicating the number of large partitions ('number of 
*/
/*          paritions in first direction') */
/* ip2    = int indicating the number of smaller partitions, per */
/*          large partition, ('number of partitions in second direction') 
*/
/* nfirst = number of nodes in the first level that is input in riord */
/* riord  = (also an ouput argument). on entry riord contains the labels 
*/
/*          of the nfirst nodes that constitute the first level. */
/* on return: */
/* ----------- */
/* ndom   = total number of partitions found */
/* map    = list of nodes listed partition by pertition from partition 1 
*/
/*          to paritition ndom. */
/* mapptr = pointer array for map. All nodes from position */
/*          k1=mapptr(idom),to position k2=mapptr(idom+1)-1 in map belong 
*/
/*          to partition idom. */
/* work arrays: */
/* ------------- */
/* mask   = array of length n, used to hold the partition number of each 
*/
/*          node for the first (large) partitioning. */
/*          mask is also used as a marker of  visited nodes. */
/* levels = int array of length .le. n used to hold the pointer */
/*          arrays for the various level structures obtained from BFS. */
/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --iwk;
    --levels;
    --mask;
    --mapptr;
    --map;
    --riord;
    --ia;
    --ja;

    /* Function Body */
    maskval = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	mask[j] = maskval;
    }
    iwk[1] = 0;
    bfs_(n, &ja[1], &ia[1], nfirst, &iwk[1], &mask[1], &maskval, &riord[1], &
	    levels[1], &nlev);
/*      init = riord(1) */
/*      call perphn (ja,ia,mask,maskval,init,nlev,riord,levels) */
    stripes_(&nlev, &riord[1], &levels[1], ip1, &map[1], &mapptr[1], ndom);
/* -----------------------------------------------------------------------
 */
    if (*ip2 == 1) {
	return 0;
    }
    ndp1 = *ndom + 1;

/*     pack info into array iwk */

    i__1 = *ndom + 1;
    for (j = 1; j <= i__1; ++j) {
	iwk[j] = ndp1 + mapptr[j];
    }
    i__1 = mapptr[*ndom + 1] - 1;
    for (j = 1; j <= i__1; ++j) {
	iwk[ndp1 + j] = map[j];
    }
/* -----------------------------------------------------------------------
 */
    i__1 = *ndom;
    for (idom = 1; idom <= i__1; ++idom) {
	i__2 = mapptr[idom + 1] - 1;
	for (k = mapptr[idom]; k <= i__2; ++k) {
	    mask[map[k]] = idom;
	}
    }
    nextdom = 1;

/*     jdom = counter for total number of (small) subdomains */

    jdom = 1;
    mapptr[jdom] = 1;
/* -----------------------------------------------------------------------
 */
    i__1 = *ndom;
    for (idom = 1; idom <= i__1; ++idom) {
	maskval = idom;
	*nfirst = 1;
	numnod = iwk[idom + 1] - iwk[idom];
	j = iwk[idom];
	init = iwk[j];
	nextdom = mapptr[jdom];
	perphn_(&numnod, &ja[1], &ia[1], &init, &iwk[j], &mask[1], &maskval, &
		nlev, &riord[1], &levels[1]);
	stripes_(&nlev, &riord[1], &levels[1], ip2, &map[nextdom], &mapptr[
		jdom], &kdom);
	mapptr[jdom] = nextdom;
	i__2 = jdom + kdom - 1;
	for (j = jdom; j <= i__2; ++j) {
	    mapptr[j + 1] = nextdom + mapptr[j + 1] - 1;
	}
	jdom += kdom;
    }

    *ndom = jdom - 1;
    return 0;
} /* dblstr_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int perphn_(int *n, int *ja, int *ia, int *
	init, int *iperm, int *mask, int *maskval, int *nlev, 
	int *riord, int *levels)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j, nlevp, mindeg, nfirst, deg;
    extern /* Subroutine */ int bfs_(int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *);
    static int nod;
    extern int maskdeg_(int *, int *, int *, int *, 
	    int *);

/* -----------------------------------------------------------------------
 */
/*     finds a pseudo-peripheral node and does a BFS search from it. */
/* -----------------------------------------------------------------------
 */
/* see routine  dblstr for description of parameters */
/* input: */
/* ------- */
/* ja, ia  = list pointer array for the adjacency graph */
/* mask    = array used for masking nodes -- see maskval */
/* maskval = value to be checked against for determing whether or */
/*           not a node is masked. If mask(k) .ne. maskval then */
/*           node k is not considered. */
/* init    = init node in the pseudo-peripheral node algorithm. */

/* output: */
/* ------- */
/* init    = actual pseudo-peripherial node found. */
/* nlev    = number of levels in the final BFS traversal. */
/* riord   = */
/* levels  = */
/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --levels;
    --riord;
    --mask;
    --iperm;
    --ia;
    --ja;

    /* Function Body */
    nlevp = 0;
L1:
    riord[1] = *init;
    nfirst = 1;
    bfs_(n, &ja[1], &ia[1], &nfirst, &iperm[1], &mask[1], maskval, &riord[1], 
	    &levels[1], nlev);
    if (*nlev > nlevp) {
	mindeg = levels[*nlev + 1] - 1;
	i__1 = levels[*nlev + 1] - 1;
	for (j = levels[*nlev]; j <= i__1; ++j) {
	    nod = riord[j];
	    deg = maskdeg_(&ja[1], &ia[1], &nod, &mask[1], maskval);
	    if (deg < mindeg) {
		*init = nod;
		mindeg = deg;
	    }
	}
	nlevp = *nlev;
	goto L1;
    }
    return 0;
} /* perphn_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int add_lvst__(int *istart, int *iend, int *nlev,
	 int *riord, int *ja, int *ia, int *mask, int *
	maskval)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, j, k, ir, nod;

/* ---------------------------------------------------------------------- 
*/
/* adds one level set to the previous sets. span all nodes of previous */
/* set. Uses Mask to mark those already visited. */
/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --mask;
    --ia;
    --ja;
    --riord;

    /* Function Body */
    nod = *iend;
    i__1 = *iend;
    for (ir = *istart + 1; ir <= i__1; ++ir) {
	i__ = riord[ir];
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    if (mask[j] == *maskval) {
		++nod;
		mask[j] = 0;
		riord[nod] = j;
	    }
/* L24: */
	}
/* L25: */
    }
    *istart = *iend;
    *iend = nod;
    return 0;
/* -----------------------------------------------------------------------
 */
} /* add_lvst__ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int stripes_(int *nlev, int *riord, int *levels, 
	int *ip, int *map, int *mapptr, int *ndom)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    static int ilev, nsiz, psiz, k, ib, ktr;

/* -----------------------------------------------------------------------
 */
/*    this is a post processor to BFS. stripes uses the output of BFS to 
*/
/*    find a decomposition of the adjacency graph by stripes. It fills */
/*    the stripes level by level until a number of nodes .gt. ip is */
/*    is reached. */
/*---------------------------parameters----------------------------------
-*/
/* on entry: */
/* -------- */
/* nlev   = number of levels as found by BFS */
/* riord  = reverse permutation array produced by BFS -- */
/* levels = pointer array for the level structure as computed by BFS. If 
*/
/*          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1, 
*/
/*          then all the nodes of level number lev are: */
/*                      riord(k1),riord(k1+1),...,riord(k2) */
/*  ip    = number of desired partitions (subdomains) of about equal size.
 */

/* on return */
/* --------- */
/* ndom     = number of subgraphs (subdomains) found */
/* map      = node per processor list. The nodes are listed contiguously 
*/
/*            from proc 1 to nproc = mpx*mpy. */
/* mapptr   = pointer array for array map. list for proc. i starts at */
/*            mapptr(i) and ends at mapptr(i+1)-1 in array map. */
/* -----------------------------------------------------------------------
 */
/* local variables. */

    /* Parameter adjustments */
    --levels;
    --riord;
    --map;
    --mapptr;

    /* Function Body */
    *ndom = 1;
    ib = 1;
/* to add: if (ip .le. 1) then ... */
    nsiz = levels[*nlev + 1] - levels[1];
/* Computing MAX */
    i__1 = 1, i__2 = *ip - *ndom + 1;
    psiz = (nsiz - ib) / max(i__1,i__2) + 1;
    mapptr[*ndom] = ib;
    ktr = 0;
    i__1 = *nlev;
    for (ilev = 1; ilev <= i__1; ++ilev) {

/*     add all nodes of this level to domain */

	i__2 = levels[ilev + 1] - 1;
	for (k = levels[ilev]; k <= i__2; ++k) {
	    map[ib] = riord[k];
	    ++ib;
	    ++ktr;
	    if (ktr >= psiz || k >= nsiz) {
		++(*ndom);
		mapptr[*ndom] = ib;
/* Computing MAX */
		i__3 = 1, i__4 = *ip - *ndom + 1;
		psiz = (nsiz - ib) / max(i__3,i__4) + 1;
		ktr = 0;
	    }

/* L3: */
	}
/* L10: */
    }
    --(*ndom);
    return 0;
/* -----------------------------------------------------------------------
 */
/* -----end-of-stripes----------------------------------------------------
 */
} /* stripes_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int reversp_(int *n, int *riord)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j, k;

/* -----------------------------------------------------------------------
 */
/*     this routine does an in-place reversing of the permutation array */
/*     riord -- */
/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --riord;

    /* Function Body */
    i__1 = *n / 2;
    for (j = 1; j <= i__1; ++j) {
	k = riord[j];
	riord[j] = riord[*n - j + 1];
	riord[*n - j + 1] = k;
/* L26: */
    }
    return 0;
} /* reversp_ */

/* ----------------------------------------------------------------------- */
int maskdeg_(int *ja, int *ia, int *nod, int *mask, 
	int *maskval)
{
    /* System generated locals */
    int ret_val, i__1;

    /* Local variables */
    static int k, deg;

/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --mask;
    --ia;
    --ja;

    /* Function Body */
    deg = 0;
    i__1 = ia[*nod + 1] - 1;
    for (k = ia[*nod]; k <= i__1; ++k) {
	if (mask[ja[k]] == *maskval) {
	    ++deg;
	}
    }
    ret_val = deg;
    return ret_val;
} /* maskdeg_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int rperm_(int *nrow, doublereal *a, int *ja, 
	int *ia, doublereal *ao, int *jao, int *iao, int *
	perm, int *job)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, j, k, ii, ko;
    static bool values;

/* -----------------------------------------------------------------------
 */
/* this subroutine permutes the rows of a matrix in CSR format. */
/* rperm  computes B = P A  where P is a permutation matrix. */
/* the permutation P is defined through the array perm: for each j, */
/* perm(j) represents the destination row number of row number j. */
/* Youcef Saad -- recoded Jan 28, 1991. */
/* -----------------------------------------------------------------------
 */
/* on entry: */
/* ---------- */
/* n 	= dimension of the matrix */
/* a, ja, ia = input matrix in csr format */
/* perm 	= int array of length nrow containing the permutation arrays 
*/
/* 	  for the rows: perm(i) is the destination of row i in the */
/*         permuted matrix. */
/*         ---> a(i,j) in the original matrix becomes a(perm(i),j) */
/*         in the output  matrix. */

/* job	= int indicating the work to be done: */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/*                       (including the copying of real values ao and */
/*                       the array iao). */
/* 		job .ne. 1 :  ignore real values. */
/*                     (in which case arrays a and ao are not needed nor 
*/
/*                      used). */

/* ------------ */
/* on return: */
/* ------------ */
/* ao, jao, iao = input matrix in a, ja, ia format */
/* note : */
/*        if (job.ne.1)  then the arrays a and ao are not used. */
/* ----------------------------------------------------------------------c
 */
/*           Y. Saad, May  2, 1990                                      c 
*/
/* ----------------------------------------------------------------------c
 */
    /* Parameter adjustments */
    --perm;
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;

    /* Function Body */
    values = *job == 1;

/*     determine pointers for output matix. */

    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	i__ = perm[j];
	iao[i__ + 1] = ia[j + 1] - ia[j];
/* L50: */
    }

/* get pointers from lengths */

    iao[1] = 1;
    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	iao[j + 1] += iao[j];
/* L51: */
    }

/* copying */

    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/* old row = ii  -- new row = iperm(ii) -- ko = new pointer */

	ko = iao[perm[ii]];
	i__2 = ia[ii + 1] - 1;
	for (k = ia[ii]; k <= i__2; ++k) {
	    jao[ko] = ja[k];
	    if (values) {
		ao[ko] = a[k];
	    }
	    ++ko;
/* L60: */
	}
/* L100: */
    }

    return 0;
/* ---------end-of-rperm -------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* rperm_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int cperm_(int *nrow, doublereal *a, int *ja, 
	int *ia, doublereal *ao, int *jao, int *iao, int *
	perm, int *job)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, k, nnz;

/* -----------------------------------------------------------------------
 */
/* this subroutine permutes the columns of a matrix a, ja, ia. */
/* the result is written in the output matrix  ao, jao, iao. */
/* cperm computes B = A P, where  P is a permutation matrix */
/* that maps column j into column perm(j), i.e., on return */
/*      a(i,j) becomes a(i,perm(j)) in new matrix */
/* Y. Saad, May 2, 1990 / modified Jan. 28, 1991. */
/* -----------------------------------------------------------------------
 */
/* on entry: */
/* ---------- */
/* nrow 	= row dimension of the matrix */

/* a, ja, ia = input matrix in csr format. */

/* perm	= int array of length ncol (number of columns of A */
/*         containing the permutation array  the columns: */
/*         a(i,j) in the original matrix becomes a(i,perm(j)) */
/*         in the output matrix. */

/* job	= int indicating the work to be done: */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/*                       (including the copying of real values ao and */
/*                       the array iao). */
/* 		job .ne. 1 :  ignore real values ao and ignore iao. */

/* ------------ */
/* on return: */
/* ------------ */
/* ao, jao, iao = input matrix in a, ja, ia format (array ao not needed) 
*/

/* Notes: */
/* ------- */
/* 1. if job=1 then ao, iao are not used. */
/* 2. This routine is in place: ja, jao can be the same. */
/* 3. If the matrix is initially sorted (by increasing column number) */
/*    then ao,jao,iao  may not be on return. */

/* ----------------------------------------------------------------------c
 */
/* local parameters: */

    /* Parameter adjustments */
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;
    --perm;

    /* Function Body */
    nnz = ia[*nrow + 1] - 1;
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	jao[k] = perm[ja[k]];
/* L100: */
    }

/*     done with ja array. return if no need to touch values. */

    if (*job != 1) {
	return 0;
    }

/* else get new pointers -- and copy values too. */

    i__1 = *nrow + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iao[i__] = ia[i__];
/* L1: */
    }

    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
/* L2: */
    }

    return 0;
/* ---------end-of-cperm--------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* cperm_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int dperm_(int *nrow, doublereal *a, int *ja, 
	int *ia, doublereal *ao, int *jao, int *iao, int *
	perm, int *qperm, int *job)
{
    extern /* Subroutine */ int cperm_(int *, doublereal *, int *, 
	    int *, doublereal *, int *, int *, int *, int 
	    *), rperm_(int *, doublereal *, int *, int *, 
	    doublereal *, int *, int *, int *, int *);
    static int locjob;

/* -----------------------------------------------------------------------
 */
/* This routine permutes the rows and columns of a matrix stored in CSR */
/* format. i.e., it computes P A Q, where P, Q are permutation matrices. 
*/
/* P maps row i into row perm(i) and Q maps column j into column qperm(j):
 */
/*      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix */
/* In the particular case where Q is the transpose of P (symmetric */
/* permutation of A) then qperm is not needed. */
/* note that qperm should be of length ncol (number of columns) but this 
*/
/* is not checked. */
/* -----------------------------------------------------------------------
 */
/* Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. */
/* -----------------------------------------------------------------------
 */
/* on entry: */
/* ---------- */
/* n 	= dimension of the matrix */
/* a, ja, */
/*    ia = input matrix in a, ja, ia format */
/* perm 	= int array of length n containing the permutation arrays */
/* 	  for the rows: perm(i) is the destination of row i in the */
/*         permuted matrix -- also the destination of column i in case */
/*         permutation is symmetric (job .le. 2) */

/* qperm	= same thing for the columns. This should be provided only */
/*         if job=3 or job=4, i.e., only in the case of a nonsymmetric */
/* 	  permutation of rows and columns. Otherwise qperm is a dummy */

/* job	= int indicating the work to be done: */
/* * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P) */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/* 		job = 2 permute matrix ignoring real values. */
/* * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q */
/* 		job = 3	permute a, ja, ia into ao, jao, iao */
/* 		job = 4 permute matrix ignoring real values. */

/* on return: */
/* ----------- */
/* ao, jao, iao = input matrix in a, ja, ia format */

/* in case job .eq. 2 or job .eq. 4, a and ao are never referred to */
/* and can be dummy arguments. */
/* Notes: */
/* ------- */
/*  1) algorithm is in place */
/*  2) column indices may not be sorted on return even  though they may be
 */
/*     on entry. */
/* ----------------------------------------------------------------------c
 */
/* local variables */

/*     locjob indicates whether or not real values must be copied. */

    /* Parameter adjustments */
    --perm;
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;
    --qperm;

    /* Function Body */
    locjob = *job % 2;

/* permute rows first */

    rperm_(nrow, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1], &perm[1], &
	    locjob);

/* then permute columns */

    locjob = 0;

    if (*job <= 2) {
	cperm_(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		perm[1], &locjob);
    } else {
	cperm_(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		qperm[1], &locjob);
    }

    return 0;
/* -------end-of-dperm----------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* dperm_ */

/* Subroutine */ int csrcoo_(int *nrow, int *job, int *nzmax, 
	doublereal *a, int *ja, int *ia, int *nnz, doublereal *ao,
	 int *ir, int *jc, int *ierr)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, k, k1, k2;

/* -----------------------------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
/*  Compressed Sparse Row      to      Coordinate */
/* -----------------------------------------------------------------------
 */
/* converts a matrix that is stored in coordinate format */
/*  a, ir, jc into a row general sparse ao, jao, iao format. */

/* on entry: */
/* --------- */
/* nrow	= dimension of the matrix. */
/* job   = int serving as a job indicator. */
/*         if job = 1 fill in only the array ir, ignore jc, and ao. */
/*         if job = 2 fill in ir, and jc but not ao */
/*         if job = 3 fill in everything. */
/*         The reason why these options are provided is that on return */
/*         ao and jc are the same as a, ja. So when job = 3, a and ja are 
*/
/*         simply copied into ao, jc.  When job=2, only jc and ir are */
/*         returned. With job=1 only the array ir is returned. Moreover, 
*/
/*         the algorithm is in place: */
/* 	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) */
/*         will write the output matrix in coordinate format on a, ja,ia. 
*/
/*         (Important: note the order in the output arrays a, ja, ia. ) */
/*         i.e., ao can be the same as a, ir can be the same as ia */
/*         and jc can be the same as ja. */

/* a, */
/* ja, */
/* ia    = matrix in compressed sparse row format. */
/* nzmax = length of space available in ao, ir, jc. */
/*         the code will stop immediatly if the number of */
/*         nonzero elements found in input matrix exceeds nzmax. */

/* on return: */
/* ----------- */
/* ao, ir, jc = matrix in coordinate format. */

/* nnz        = number of nonzero elements in matrix. */
/* ierr       = int error indicator. */
/*         ierr .eq. 0 means normal retur */
/*         ierr .eq. 1 means that the the code stopped */
/*         because there was no space in ao, ir, jc */
/*         (according to the value of  nzmax). */

/*-----------------------------------------------------------------------
-*/
    /* Parameter adjustments */
    --jc;
    --ir;
    --ao;
    --ia;
    --ja;
    --a;

    /* Function Body */
    *ierr = 0;
    *nnz = ia[*nrow + 1] - 1;
    if (*nnz > *nzmax) {
	*ierr = 1;
	return 0;
    }
/*-----------------------------------------------------------------------
-*/
    switch (*job) {
	case 1:  goto L3;
	case 2:  goto L2;
	case 3:  goto L1;
    }
L1:
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
/* L10: */
    }
L2:
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	jc[k] = ja[k];
/* L11: */
    }
/* copy backward to allow */
L3:
    for (i__ = *nrow; i__ >= 1; --i__) {
	k1 = ia[i__ + 1] - 1;
	k2 = ia[i__];
	i__1 = k2;
	for (k = k1; k >= i__1; --k) {
	    ir[k] = i__;
/* L12: */
	}
/* L13: */
    }
    return 0;
/* ------------- end of csrcoo -------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* csrcoo_ */

/* Subroutine */ int amudia_(int *nrow, int *job, doublereal *a, 
	int *ja, int *ia, doublereal *diag, doublereal *b, int *
	jb, int *ib)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int k, k1, k2, ii;

/* -----------------------------------------------------------------------
 */
/* performs the matrix by matrix product B = A * Diag  (in place) */
/* -----------------------------------------------------------------------
 */
/* on entry: */
/* --------- */
/* nrow	= int. The row dimension of A */

/* job   = int. job indicator. Job=0 means get array b only */
/*         job = 1 means get b, and the int arrays ib, jb. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* diag = diagonal matrix stored as a vector dig(1:n) */

/* on return: */
/* ---------- */

/* b, */
/* jb, */
/* ib	= resulting matrix B in compressed sparse row sparse format. */

/* Notes: */
/* ------- */
/* 1)        The column dimension of A is not needed. */
/* 2)        algorithm in place (B can take the place of A). */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --ib;
    --diag;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/*     scale each element */

	k1 = ia[ii];
	k2 = ia[ii + 1] - 1;
	i__2 = k2;
	for (k = k1; k <= i__2; ++k) {
	    b[k] = a[k] * diag[ja[k]];
/* L2: */
	}
/* L1: */
    }

    if (*job == 0) {
	return 0;
    }

    ib[1] = ia[1];
    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {
	ib[ii] = ia[ii];
	i__2 = ia[ii + 1] - 1;
	for (k = ia[ii]; k <= i__2; ++k) {
	    jb[k] = ja[k];
/* L31: */
	}
/* L3: */
    }
    return 0;
/* -----------------------------------------------------------------------
 */
/* -----------end-of-amudiag----------------------------------------------
 */
} /* amudia_ */

/* Subroutine */ int ilutp_(int *n, doublereal *a, int *ja, int *
	ia, int *lfil, doublereal *droptol, doublereal *permtol, int *
	mbloc, doublereal *alu, int *jlu, int *ju, int *iwk, 
	doublereal *wu, doublereal *wl, int *jr, int *jwl, int *
	jwu, int *iperm, int *ierr)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal fact;
    static int lenl, imax, lenu, icut, jpos;
    static doublereal xmax;
    static int jrow, lenl0, lenu0;
    static doublereal xmax0;
    static int i__, j, k;
    static doublereal s, t;
    static int j1, j2;
    static doublereal tnorm;
    static int ii, jj, nl;
    extern /* Subroutine */ int qsplit_(doublereal *, int *, int *, 
	    int *);
    static int ju0, len;
    static doublereal tmp;

/* -----------------------------------------------------------------------
 */
/*     SPARSKIT ROUTINE   --  ILUT with PIVOTING -- */
/* -----------------------------------------------------------------------
 */
/*      implicit none */
/* ----------------------------------------------------------------------*
 */
/*                      *** ILUT preconditioner ***                     * 
*/
/*                      ---------------------------                     * 
*/
/*      incomplete LU factorization with dual truncation mechanism      * 
*/
/*      VERSION 2 : sorting  done for both L and U.                     * 
*/
/*                                                                      * 
*/
/* ----------------------------------------------------------------------*
 */
/* ---- coded by Youcef Saad Sep 8, 1993 . ------------------------------*
 */
/* ---- Dual drop-off strategy works as follows.                         *
 */
/*                                                                      * 
*/
/* 1) Theresholding in L and U as set by droptol. Any element whose size* 
*/
/*    is less than some tolerance (relative to the norm of current      * 
*/
/*    row in u) is dropped.                                             * 
*/
/*                                                                      * 
*/
/* 2) Keeping only the largest lfil+il(i) elements in the i-th row      * 
*/
/*    of L and the largest lfil+iu(i) elements in the i-th row of       * 
*/
/*    U where il(i), iu(i) are the original number of nonzero           * 
*/
/*    elements of the L-part and the U-part of the i-th row of A        * 
*/
/*                                                                      * 
*/
/* column permuting is used --                                          * 
*/
/*  see also comments in ilut                                           * 
*/
/* ----------------------------------------------------------------------*
 */
/* PARAMETERS */
/* ----------- */

/* on entry: */
/* ========== */
/* n       = int. The dimension of the matrix A. */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */
/*           ONE RETURN THE COLUMNS OF A ARE PERMUTED. */

/* lfil    = int. The fill-in parameter. Each row of L and */
/*           each row of U will have a maximum of lfil elements */
/*           in addition to their original number of nonzero elements. */
/*           Thus storage can be determined beforehand. */
/*           lfil must be .ge. 0. */

/* droptol = tolerance used for dropping elements in L and U. */
/*           elements are dropped if they are .lt. norm(row) x droptol */
/*           row = row being eliminated */

/* permtol = tolerance ratio used for determning whether to permute */
/*           two columns. We will permute two columns only if */
/*           a(i,j)*permtol .gt. a(i,i) [good values 0.1 to 0.01] */

/* mbloc   = if desired, permuting can be done only within the diagonal */
/*           blocks of size mbloc. Useful for PDE problems with several */
/*           degrees of freedom.. If feature not wanted take mbloc=n. */

/* iwk     = int. The minimum length of arrays alu and jlu */
/*           to work properly, the code requires that iwk be */

/*                      .ge. nnz + 2*lfil*n + 2 */

/*           where nnz = original number of nonzero elements in A. */
/*           if iwk is not large enough the code will stop prematurely */
/*           with ierr = -2 or ierr = -3 (see below). */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing 
*/
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
*/
/*           contains the i-th row of L (excluding the diagonal entry=1) 
*/
/*           followed by the i-th row of U. */

/* ju      = int array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */
/* iperm   = contains the permutation arrays .. */
/*           iperm(1:n) = old numbers of unknowns */
/*           iperm(n+1:2*n) = reverse permutation = new unknowns. */

/* ierr    = int. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. 
*/
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered. */

/* work arrays: */
/* ============= */
/* jr,jwu,jwl 	  = int work arrays of length n. */
/* wu, wl          = real work arrays of length n+1, and n resp. */

/* Notes: */
/* ------ */
/* A must have all nonzero diagonal elements. */
/* U -- MATRIX   STORED IN UNPERMUTED FORMAT TO AVOID PERMUTATION ARRAYS 
*/
/* code working -- CODED BY Y. SAAD, SEPT 9, 1993. */
/* -----------------------------------------------------------------------
 */
/*     local variables */


    /* Parameter adjustments */
    --iperm;
    --jwu;
    --jwl;
    --jr;
    --wl;
    --wu;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;

    /* Function Body */
    if (*lfil < 0) {
	goto L998;
    }
/* ------------------------------- */
/* initialize ju0 (points to next element to be added to alu,jlu) */
/* and pointer. */
/* -----------------------------------------------------------------------
 */
    ju0 = *n + 2;
    jlu[1] = ju0;

/*  int double pointer array. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jr[j] = 0;
	iperm[j] = j;
	iperm[*n + j] = j;
/* L1: */
    }
/* -----------------------------------------------------------------------
 */
/*  beginning of main loop. */
/* -----------------------------------------------------------------------
 */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	j1 = ia[ii];
	j2 = ia[ii + 1] - 1;
	tnorm = 0.;
	i__2 = j2;
	for (k = j1; k <= i__2; ++k) {
	    tnorm += (d__1 = a[k], abs(d__1));
/* L501: */
	}
	if (tnorm == (float)0.) {
	    goto L999;
	}
	tnorm /= j2 - j1 + 1;

/* --- unpack L-part and U-part of row of A in arrays wl, wu -- */

	lenu = 1;
	lenl = 0;
	jwu[1] = ii;
	wu[1] = (float)0.;
	jr[ii] = 1;
/* ------------------------------------------------------------------
----- */
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    k = iperm[*n + ja[j]];
	    t = a[j];
	    if (abs(t) < *droptol * tnorm && k != ii) {
		goto L170;
	    }
	    if (k < ii) {
		++lenl;
		jwl[lenl] = k;
		wl[lenl] = t;
		jr[k] = lenl;
	    } else if (k == ii) {
		wu[1] = t;
	    } else {
		++lenu;
		jwu[lenu] = k;
		wu[lenu] = t;
		jr[k] = lenu;
	    }
L170:
	    ;
	}
/* ------------------------------------------------------------------
----- */
	tnorm /= j2 - j1 + 1;
	lenl0 = lenl;
	lenu0 = lenu;
	jj = 0;
	nl = 0;
/* -------------------------------------------------------------------
 */
/* ---------------------- eliminate previous rows --------------------
 */
/* -------------------------------------------------------------------
 */
L150:
	++jj;
	if (jj > lenl) {
	    goto L160;
	}
/* -------------------------------------------------------------------
 */
/* in order to do the elimination in the correct order we need to */
/* exchange the current row number with the one that has */
/* smallest column number, among jj,jj+1,...,lenl. */
/* -------------------------------------------------------------------
 */
	jrow = jwl[jj];
	k = jj;

/* determine smallest column index */

	i__2 = lenl;
	for (j = jj + 1; j <= i__2; ++j) {
	    if (jwl[j] < jrow) {
		jrow = jwl[j];
		k = j;
	    }
/* L151: */
	}

/*     exchange in jwl */

	if (k != jj) {
	    j = jwl[jj];
	    jwl[jj] = jwl[k];
	    jwl[k] = j;

/*     exchange in jr */

	    jr[jrow] = jj;
	    jr[j] = k;

/*     exchange in wl */

	    s = wl[jj];
	    wl[jj] = wl[k];
	    wl[k] = s;
	}
/* ------------------------------------------------------------------
----- */
	if (jrow >= ii) {
	    goto L160;
	}

/*     get the multiplier for row to be eliminated: jrow */

	fact = wl[jj] * alu[jrow];
/*     zero out element in row by setting jr(jrow) = 0 */

	jr[jrow] = 0;

	if (abs(fact) * wu[*n + 2 - jrow] <= *droptol * tnorm) {
	    goto L150;
	}
/* -------------------------------------------------------------------
 */
/* ------------combine current row and row jrow --------------------- 
*/
/* -------------------------------------------------------------------
 */
	i__2 = jlu[jrow + 1] - 1;
	for (k = ju[jrow]; k <= i__2; ++k) {
	    s = fact * alu[k];
/*     new column number */
	    j = iperm[*n + jlu[k]];
	    jpos = jr[j];

/*     if fill-in element is small then disregard: */

	    if (abs(s) < *droptol * tnorm && jpos == 0) {
		goto L203;
	    }
	    if (j >= ii) {

/*     dealing with upper part. */

		if (jpos == 0) {
/*     this is a fill-in element */
		    ++lenu;
		    if (lenu > *n) {
			goto L995;
		    }
		    jwu[lenu] = j;
		    jr[j] = lenu;
		    wu[lenu] = -s;
		} else {
/*     no fill-in element -- */
		    wu[jpos] -= s;
		}
	    } else {

/*     dealing with lower part. */

		if (jpos == 0) {
/*     this is a fill-in element */
		    ++lenl;
		    if (lenl > *n) {
			goto L995;
		    }
		    jwl[lenl] = j;
		    jr[j] = lenl;
		    wl[lenl] = -s;
		} else {
/*     no fill-in element -- */
		    wl[jpos] -= s;
		}
	    }
L203:
	    ;
	}
	++nl;
	wl[nl] = fact;
	jwl[nl] = jrow;
	goto L150;
/* ---------------------------------------------------------- */
/* ------------ update l-matrix ----------------------------- */
/* ---------------------------------------------------------- */
L160:
/* Computing MIN */
	i__2 = nl, i__3 = lenl0 + *lfil;
	len = min(i__2,i__3);
/*     160    len = min0(nl,lfil) */
	qsplit_(&wl[1], &jwl[1], &nl, &len);

/*     store L-part -- in original coordinates .. */

	i__2 = len;
	for (k = 1; k <= i__2; ++k) {
	    if (ju0 > *iwk) {
		goto L996;
	    }
	    alu[ju0] = wl[k];
	    jlu[ju0] = iperm[jwl[k]];
/*     jlu(ju0) = jwl(k) */
	    ++ju0;
/* L204: */
	}

/*     save pointer to beginning of row ii of U */

	ju[ii] = ju0;

/*     reset double-pointer jr to zero (L-part - except first */
/*     jj-1 elements which have already been reset) */

	i__2 = lenl;
	for (k = jj; k <= i__2; ++k) {
	    jr[jwl[k]] = 0;
/* L306: */
	}
/* ---------------------------------------------------------- */
/* ------------update u-matrix ----------------------------- */
/* ---------------------------------------------------------- */
/* Computing MIN */
	i__2 = lenu, i__3 = lenu0 + *lfil;
	len = min(i__2,i__3);
/*     len = min0(lenu,lfil) */
	i__2 = lenu - 1;
	qsplit_(&wu[2], &jwu[2], &i__2, &len);

	imax = 1;
	xmax = (d__1 = wu[imax], abs(d__1));
	xmax0 = xmax;

	icut = ii - 1 + *mbloc - (ii - 1) % *mbloc;

	i__2 = len;
	for (k = 2; k <= i__2; ++k) {
	    t = (d__1 = wu[k], abs(d__1));
	    if (t > xmax && t * *permtol > xmax0 && jwu[k] <= icut) {
		imax = k;
		xmax = t;
	    }
	}

/*     exchange wu's */

	tmp = wu[1];
	wu[1] = wu[imax];
	wu[imax] = tmp;

/*     update iperm and reverse iperm */

	j = jwu[imax];
	i__ = iperm[ii];
	iperm[ii] = iperm[j];
	iperm[j] = i__;
/*     reverse iperm */
	iperm[*n + iperm[ii]] = ii;
	iperm[*n + iperm[j]] = j;

	t = (d__1 = wu[k], abs(d__1));
	if (len + ju0 > *iwk) {
	    goto L997;
	}

/*     store U-part in original coordinates */

	i__2 = len;
	for (k = 2; k <= i__2; ++k) {
	    jlu[ju0] = iperm[jwu[k]];
	    alu[ju0] = wu[k];
	    t += (d__1 = wu[k], abs(d__1));
	    ++ju0;
/* L302: */
	}

/*     save norm in wu (backwards). Norm is in fact average abs value 
*/

	wu[*n + 2 - ii] = t / (len + 1);

/*     store inverse of diagonal element of u */

	if (wu[1] == (float)0.) {
	    wu[1] = (*droptol + 1e-4) * tnorm;
	}

	alu[ii] = 1. / wu[1];

/*     update pointer to beginning of next row of U. */

	jlu[ii + 1] = ju0;

/*     reset double-pointer jr to zero (U-part) */

	i__2 = lenu;
	for (k = 1; k <= i__2; ++k) {
	    jr[jwu[k]] = 0;
/* L308: */
	}
/* ------------------------------------------------------------------
----- */
/*     end main loop */
/* ------------------------------------------------------------------
----- */
/* L500: */
    }

/*     permute all column indices of LU ... */


/*     call dvperm(n,alu,iperm(n+1)) */

/*     do ii =1, n */
/*     do k = ju(ii), jlu(ii+1)-1 */
/*     jlu(k) = iperm(jlu(k)) */
/*     enddo */
/*     enddo */
/* -----------------------------------------------------------------------
 */
    i__1 = jlu[*n + 1] - 1;
    for (k = jlu[1]; k <= i__1; ++k) {
	jlu[k] = iperm[*n + jlu[k]];
    }

/*     ...and A */

    i__1 = ia[*n + 1] - 1;
    for (k = 1; k <= i__1; ++k) {
	ja[k] = iperm[*n + ja[k]];
    }

    *ierr = 0;
    return 0;

/*     zero pivot : */

/*     900    ierr = ii */
/*     return */

/*     incomprehensible error. Matrix must be wrong. */

L995:
    *ierr = -1;
    return 0;

/*     insufficient storage in L. */

L996:
    *ierr = -2;
    return 0;

/*     insufficient storage in U. */

L997:
    *ierr = -3;
    return 0;

/*     illegal lfil entered. */

L998:
    *ierr = -4;
    return 0;

/*     zero row encountered */

L999:
    *ierr = -5;
    return 0;
/* ----------------end-of-ilutp-------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* ilutp_ */

/* ---------------------------------------------------------------------- */
/* Subroutine */ int qsplit_(doublereal *a, int *ind, int *n, int 
	*ncut)
{
    /* System generated locals */
    int i__1;
    doublereal d__1;

    /* Local variables */
    static int last, itmp, j, first;
    static doublereal abskey;
    static int mid;
    static doublereal tmp;

/* -----------------------------------------------------------------------
 */
/*     does a quick-sort split of a real array. */
/*     on input a(1:n). is a real array */
/*     on output a(1:n) is permuted such that its elements satisfy: */
/*     a(i) .le. a(ncut) for i .le. ncut and */
/*     a(i) .ge. a(ncut) for i .ge. ncut */
/*    ind(1:n) is an int array which permuted in the same way as a(*).
*/
/* -----------------------------------------------------------------------
 */
/* ----- */
    /* Parameter adjustments */
    --ind;
    --a;

    /* Function Body */
    first = 1;
    last = *n;
    if (*ncut < first || *ncut > last) {
	return 0;
    }

/*     outer loop -- while mid .ne. ncut do */

L1:
    mid = first;
    abskey = (d__1 = a[mid], abs(d__1));
    i__1 = last;
    for (j = first + 1; j <= i__1; ++j) {
	if ((d__1 = a[j], abs(d__1)) > abskey) {
	    ++mid;
/*     interchange */
	    tmp = a[mid];
	    itmp = ind[mid];
	    a[mid] = a[j];
	    ind[mid] = ind[j];
	    a[j] = tmp;
	    ind[j] = itmp;
	}
/* L2: */
    }

/*     interchange */

    tmp = a[mid];
    a[mid] = a[first];
    a[first] = tmp;

    itmp = ind[mid];
    ind[mid] = ind[first];
    ind[first] = itmp;

/*     test for while loop */

    if (mid == *ncut) {
	return 0;
    }
    if (mid > *ncut) {
	last = mid - 1;
    } else {
	first = mid + 1;
    }
    goto L1;
/* ----------------end-of-qsplit------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* qsplit_ */

/* Subroutine */ int atob_(int *n, doublereal *a, int *ja, int *
	ia, doublereal *b, int *jb, int *ib)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

/* ... Copy matrix a,ja,ia to b,jb,ib.  Both matrices are in */
/*     compressed sparse row format. */
/* ... Input arguments: */
/* ... Output arguments: */
/* ... Local variable: */
    /* Parameter adjustments */
    --ib;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    i__1 = ia[*n + 1] - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = a[i__];
	jb[i__] = ja[i__];
    }
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib[i__] = ia[i__];
    }
    return 0;
} /* atob_ */

/*  end of atob */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int djreord_(int *neq, int *neqp1, int *nnzmx, 
	char *premeth, doublereal *jac, int *ja, int *ia, doublereal *
	awk, int *jwk, int *iwk, int *perm, int *qperm, 
	int *levels, int *mask, int *ireorder)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    extern /* Subroutine */ int atob_(int *, doublereal *, int *, 
	    int *, doublereal *, int *, int *);
    static int nlev, i__;
    extern /* Subroutine */ int dperm_(int *, doublereal *, int *, 
	    int *, doublereal *, int *, int *, int *, int 
	    *, int *);
    static int nfirst;
    extern /* Subroutine */ int bfs_(int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *);
    static int maskval;
    extern /* Subroutine */ int reversp_(int *, int *);

/* ... Version of 10-6-95 */
/* ... If desired, reorder the Jacobian matrix. */
/* ... Input arguments: */
/* total number of equations */
/* NEQ + 1 */
/* maximum number of nonzeroes in Jacobian */
/* nonzero Jacobian elements */
/* column indices of nonzero Jacobian elements */
/* indices of 1st nonzero element in each row */
/* ... Work-array arguments: */
/* used in reordering the rows and columns of */
/* the Jacobian matrix. */
/* Integer array containing the permutation */
/* permutation in array perm. */
/* Integer array holding the inverse of the */
/* subroutine.   See subroutine BFS for */
/* more details. */
/* Work array used by the bfs reordering */
/* subroutine.  See BFS subroutine. */
/* Work array used by the BFS reordering */
/* of the Jacobian matrix is desired. */
/* = 1 means a reverse Cuthill-McKee */
/*     reordering of the rows and columns */
/*     of the Jacobian is done. */
/* = 0 means no reordering. */
/* ... Local variables: */
/* Flag used to determine if a reordering */
/* See subroutine BFS for more details. */
/* Number of levels in levels array. */
/* Scalar used with MASK. */
    /* Parameter adjustments */
    --mask;
    --levels;
    --qperm;
    --perm;
    --iwk;
    --ia;
    --jwk;
    --awk;
    --ja;
    --jac;

    /* Function Body */
    if (*ireorder == 1) {
/* ... Copy JAC, JA, and IA to AWK, JWK, and IWK. */
	atob_(neq, &jac[1], &ja[1], &ia[1], &awk[1], &jwk[1], &iwk[1]);
/* ... Perform a Cuthill-McKee reordering of the Jacobian. */
	nfirst = 1;
	perm[1] = 0;
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mask[i__] = 1;
	}
	maskval = 1;
	qperm[1] = 1;
	bfs_(neq, &jwk[1], &iwk[1], &nfirst, &perm[1], &mask[1], &maskval, &
		qperm[1], &levels[1], &nlev);
/* ... Reverse the permutation to obtain the reverse Cuthill-McKee */
/*     reordering. */
	reversp_(neq, &qperm[1]);
/* ... Calculate the inverse of QPERM and put it in PERM. */
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    perm[qperm[i__]] = i__;
	}
/* ... Permute rows and columns of Jacobian using PERM. */
	dperm_(neq, &awk[1], &jwk[1], &iwk[1], &jac[1], &ja[1], &ia[1], &perm[
		1], &perm[1], &c__1);
/* ... End of If block */
    }
    return 0;
/* ------------  End of Subroutine DJREORD  ------------------------------
 */
} /* djreord_ */
