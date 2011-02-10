// -*-C++-*-
// =========================================================================
// Generic inexact line search function
// Uses a callback function to evaluate objective function at trial steps
// =========================================================================

#ifndef __LSEARCH_H
#define __LSEARCH_H

#include "stoastlib.h"

typedef double (*OF_CLBK)(const RVector &x, void *context);

int LineSearch (const RVector &x0, const RVector &dx, double s0, double of0,
		OF_CLBK of, double *smin, double *ofmin, void *context = 0);

#endif // !__LSEARCH_H
