// =========================================================================
// Generic inexact line search function
// Uses a callback function to evaluate objective function at trial steps
// =========================================================================

#define STOASTLIB_IMPLEMENTATION
#include "lsearch.h"

int LineSearch (const RVector &x0, const RVector &dx, double s0, double of0,
		OF_CLBK of, double *smin, double *ofmin, void *context)
{
    double h0 = 0.0; // lower step bound
    double h2 = s0;  // upper step bound
    double hm;       // intermediate step
    double of2;      // objective function at upper bound
    double ofm;      // objective function at midpoint

    // phase 1: bracket the minimum
    RVector x = x0 + dx*h2;
    of2 = of(x, context);
    while (of2 < 0.0) { // error flag
	h2 *= 0.5;
	x = x0 + dx*h2;
	of2 = of (x, context);
    }
    if (of2 < of0) { // increase interval
	hm = h2;  ofm = of2;
	h2 *= 2.0;
	x = x0 + dx*h2;
	of2 = of (x, context);
	while (of2 < ofm) {
	    h0 = hm;  of0 = ofm;
	    hm = h2;  ofm = of2;
	    h2 *= 2.0;
	    x = x0 + dx*h2;
	    of2 = of (x, context);
	}
    } else { // decrease interval
	hm = 0.5*h2;
	x = x0 + dx*hm;
	ofm = of (x, context);
	while (ofm > of0) {
	    h2 = hm;  of2 = ofm;
	    hm = 0.5*h2;
	    x = x0 + dx*hm;
	    ofm = of (x, context);
	}
    }

    // phase 2: quadratic interpolation
    double a = ((of0-of2)/(h0-h2) - (of0-ofm)/(h0-hm)) / (h2-hm);
    double b = (of0-of2)/(h0-h2) - a*(h0+h2);
    *smin = -b/(2.0*a);
    x = x0 + dx*(*smin);
    *ofmin = of (x, context);
    if (*ofmin > ofm) { // interpolation didn't give improvement
	*smin = hm;
	*ofmin = ofm;
    }
    return 0;
}
