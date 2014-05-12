/* mc21.f -- translated by f2c (version 20060506).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* ######DATE 8 Oct 1992 COPYRIGHT Rutherford Appleton Laboratory */
/* ######8/10/92 Toolpack tool decs employed. */
/* ######8/10/92 D version created by name change only. */
/* Subroutine */ int mc21ad_(integer *n, integer *icn, integer *licn, integer 
	*ip, integer *lenr, integer *iperm, integer *numnz, integer *iw)
{
    /* System generated locals */
    integer iw_dim1, iw_offset;

    /* Local variables */
    extern /* Subroutine */ int mc21bd_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    iw_dim1 = *n;
    iw_offset = 1 + iw_dim1;
    iw -= iw_offset;
    --iperm;
    --lenr;
    --ip;
    --icn;

    /* Function Body */
    mc21bd_(n, &icn[1], licn, &ip[1], &lenr[1], &iperm[1], numnz, &iw[iw_dim1 
	    + 1], &iw[(iw_dim1 << 1) + 1], &iw[iw_dim1 * 3 + 1], &iw[(iw_dim1 
	    << 2) + 1]);
    return 0;

} /* mc21ad_ */

/* Subroutine */ int mc21bd_(integer *n, integer *icn, integer *licn, integer 
	*ip, integer *lenr, integer *iperm, integer *numnz, integer *pr, 
	integer *arp, integer *cv, integer *out)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, j1, ii, kk, in1, in2, jord, ioutk;

/*   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH. */
/* IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM. */
/*   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE */
/* ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE */
/* (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES. */
/*   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I */
/* WAS VISITED. */
/*   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I */
/* WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT. */
/*   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I */
/* WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP. */

/*   INITIALIZATION OF ARRAYS. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --out;
    --cv;
    --arp;
    --pr;
    --iperm;
    --lenr;
    --ip;
    --icn;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arp[i__] = lenr[i__] - 1;
	cv[i__] = 0;
	iperm[i__] = 0;
/* L10: */
    }
    *numnz = 0;


/*   MAIN LOOP. */
/*   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT */
/* OR GIVES A ROW WITH NO ASSIGNMENT. */
    i__1 = *n;
    for (jord = 1; jord <= i__1; ++jord) {
	j = jord;
	pr[j] = -1;
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
/* LOOK FOR A CHEAP ASSIGNMENT */
	    in1 = arp[j];
	    if (in1 < 0) {
		goto L30;
	    }
	    in2 = ip[j] + lenr[j] - 1;
	    in1 = in2 - in1;
	    i__3 = in2;
	    for (ii = in1; ii <= i__3; ++ii) {
		i__ = icn[ii];
		if (iperm[i__] == 0) {
		    goto L80;
		}
/* L20: */
	    }
/*   NO CHEAP ASSIGNMENT IN ROW. */
	    arp[j] = -1;
/*   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J. */
L30:
	    out[j] = lenr[j] - 1;
/* INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS. */
	    i__3 = jord;
	    for (kk = 1; kk <= i__3; ++kk) {
		in1 = out[j];
		if (in1 < 0) {
		    goto L50;
		}
		in2 = ip[j] + lenr[j] - 1;
		in1 = in2 - in1;
/* FORWARD SCAN. */
		i__4 = in2;
		for (ii = in1; ii <= i__4; ++ii) {
		    i__ = icn[ii];
		    if (cv[i__] == jord) {
			goto L40;
		    }
/*   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS. */
		    j1 = j;
		    j = iperm[i__];
		    cv[i__] = jord;
		    pr[j] = j1;
		    out[j1] = in2 - ii - 1;
		    goto L70;

L40:
		    ;
		}

/*   BACKTRACKING STEP. */
L50:
		j = pr[j];
		if (j == -1) {
		    goto L100;
		}
/* L60: */
	    }

L70:
	    ;
	}

/*   NEW ASSIGNMENT IS MADE. */
L80:
	iperm[i__] = j;
	arp[j] = in2 - ii - 1;
	++(*numnz);
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
	    j = pr[j];
	    if (j == -1) {
		goto L100;
	    }
	    ii = ip[j] + lenr[j] - out[j] - 2;
	    i__ = icn[ii];
	    iperm[i__] = j;
/* L90: */
	}

L100:
	;
    }

/*   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE */
/* PERMUTATION IPERM. */
    if (*numnz == *n) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arp[i__] = 0;
/* L110: */
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iperm[i__] != 0) {
	    goto L120;
	}
	++k;
	out[k] = i__;
	goto L130;

L120:
	j = iperm[i__];
	arp[j] = i__;
L130:
	;
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (arp[i__] != 0) {
	    goto L140;
	}
	++k;
	ioutk = out[k];
	iperm[ioutk] = i__;
L140:
	;
    }
    return 0;

} /* mc21bd_ */

