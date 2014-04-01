// -*-C++-*-
// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe
// File: lin_qr
//
// Gaussian quadrature rules over the 1D interval [0,1]
// ==========================================================================

// All quadrature rules return the weights in array 'wght', the integration
// abscissae in array 'absc', and the number of points as the function
// return value

// Functions are of the form QRule_lin_X where X is the degree of the
// rule (and the number of abscissae)

#ifndef __LIN_QR_H
#define __LIN_QR_H

#include <math.h>

int QRule_lin_1 (const double **wght, const double **absc);
int QRule_lin_2 (const double **wght, const double **absc);
int QRule_lin_3 (const double **wght, const double **absc);
int QRule_lin_4 (const double **wght, const double **absc);
int QRule_lin_5 (const double **wght, const double **absc);
int QRule_lin_6 (const double **wght, const double **absc);

#endif // !__LIN_QR_H
