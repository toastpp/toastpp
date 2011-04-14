// -*-C++-*-
// =========================================================================
// Generic inexact line search function
// Uses a callback function to evaluate objective function at trial steps
// =========================================================================

#ifndef __LSEARCH_H
#define __LSEARCH_H

#include "stoastlib.h"

typedef double (*OF_CLBK)(const RVector &x, double *ofparts, void *context);


/**
 * \brief Generic linesearch function
 * \param x0 current state vector
 * \param dx search direction
 * \param s0 initial step length
 * \param of0 value of objective function at current state
 * \param of callback function for evaluation of objective function
 * \param smin pointer to variable receiving the final step size
 * \param ofmin pointer to variable receiving the final objective value
 * \param context optional pointer to variable passed to callback function
 * \note The callback function must have the interface
 *     double of(const RVector &x, void *context)
 *  where x is the trial state vector, context is the context variable
 * defined in the call to LineSearch, and the return value is the value
 * of the objective function
 */
STOASTLIB int LineSearch (const RVector &x0, const RVector &dx, double s0, double of0,
		OF_CLBK of, double *smin, double *ofmin, void *context = 0);

#endif // !__LSEARCH_H
