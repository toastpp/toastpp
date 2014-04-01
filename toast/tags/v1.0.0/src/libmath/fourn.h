#ifndef __MATH_FOURN_H
#define __MATH_FOURN_H

#ifdef CVERSION_FFT
extern "C" void fourn(float *data, int *nn, int ndim, int isign);
#else
MATHLIB void fourn(float *data, int *nn, int ndim, int isign);
MATHLIB void fourn(double *data, int *nn, int ndim, int isign);
#endif
#endif
