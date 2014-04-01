// =========================================================================
// Generic inexact line search function
// Uses a callback function to evaluate objective function at trial steps
// =========================================================================

#define STOASTLIB_IMPLEMENTATION
#include "lsearch.h"

STOASTLIB int LineSearch (const RVector &x0, const RVector &dx, double s0, double of0,
		OF_CLBK of, double *smin, double *ofmin, void *context)
{
    const int MAXIT = 16;

    double h0 = 0.0; // lower step bound
    double h2 = s0;  // upper step bound
    double hm;       // intermediate step
    double of2;      // objective function at upper bound
    double ofm;      // objective function at midpoint
    double of_sub[2], ofm_sub[2];
    bool do_interpolation = true;

    // phase 1: bracket the minimum
    RVector x = x0 + dx*h2;
    of2 = of(x, of_sub, context);
    while (of2 < 0.0) { // error flag
	LOGOUT ("Parameters out of range in trial step");
	h2 *= 0.5;
	x = x0 + dx*h2;
	of2 = of (x, of_sub, context);
    }
    LOGOUT("Lsearch: STEP %g OF %g [LH %g PR %g]", 
	   h2,of2,of_sub[0],of_sub[1]);
    if (of2 < of0) { // increase interval
	hm = h2;  ofm = of2;
	ofm_sub[0] = of_sub[0]; ofm_sub[1] = of_sub[1];
	h2 *= 2.0;
	x = x0 + dx*h2;
	of2 = of (x, of_sub, context);
	if (of2 >= 0.0) {
	    LOGOUT("Lsearch: STEP %g OF %g [LH %g PR %g]",
		   h2,of2,of_sub[0],of_sub[1]);
	} else {
	    LOGOUT("Parameters out of range in trial step");
	    of2 = ofm*4.0; // stop growing interval
	    do_interpolation = false;
	}
	while (of2 < ofm) {
	    h0 = hm;  of0 = ofm;
	    hm = h2;  ofm = of2;
	    ofm_sub[0] = of_sub[0]; ofm_sub[1] = of_sub[1];
	    h2 *= 2.0;
	    x = x0 + dx*h2;
	    of2 = of (x, of_sub, context);
	    if (of2 >= 0.0) {
		LOGOUT("Lsearch: STEP %g OF %g [LH %g PR %g]",
		       h2,of2,of_sub[0],of_sub[1]);
	    } else {
		LOGOUT("Parameters out of range in trial step");
		of2 = ofm*4.0; // stop growing interval
		do_interpolation = false;
	    }
	}
    } else { // decrease interval
	hm = 0.5*h2;
	x = x0 + dx*hm;
	ofm = of (x, of_sub, context);
	LOGOUT("Lsearch: STEP %g OF %g [LH %g PR %g]",
	       hm,ofm,of_sub[0],of_sub[1]);
	int itcount = 0;
	while (ofm > of0) {
	    if (++itcount > MAXIT) return 1;
	    h2 = hm;  of2 = ofm;
	    hm = 0.5*h2;
	    x = x0 + dx*hm;
	    ofm = of (x, of_sub, context);
	    LOGOUT("Lsearch: STEP %g OF %g [LH %g PR %g]",
		   hm,ofm,of_sub[0],of_sub[1]);
	}
    }

    // phase 2: quadratic interpolation
    if (do_interpolation) {
        double a = ((of0-of2)/(h0-h2) - (of0-ofm)/(h0-hm)) / (h2-hm);
	double b = (of0-of2)/(h0-h2) - a*(h0+h2);
	*smin = -b/(2.0*a);
	x = x0 + dx*(*smin);
	ofm_sub[0] = of_sub[0], ofm_sub[1] = of_sub[1];
	*ofmin = of (x, of_sub, context);
	if (*ofmin > ofm) { // interpolation didn't give improvement
	    *smin = hm;
	    *ofmin = ofm;
	    of_sub[0] = ofm_sub[0], of_sub[1] = ofm_sub[1];
	}
    } else {
        *smin = hm;
	*ofmin = ofm;
	of_sub[0] = ofm_sub[0]; of_sub[1] = ofm_sub[1];
    }

    LOGOUT("Lsearch final: STEP %g OF %g [LH %g PR %g]",
	   *smin,*ofmin,of_sub[0],of_sub[1]);
    return 0;
}
