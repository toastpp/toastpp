/* Fortran to C header file */

#ifndef __F2C__
#define __F2C__

#include "long_integer.h"

#define logical                     integer
#define ftnlen                      integer
#define real                        float
#define doubleprecision             double
#define character                   char
typedef struct { real r, i; } complex;

// original
//typedef struct { doubleprecision r, i; } doublecomplex;

// modified
typedef struct { double r, i; } ilu_doublecomplex;

#endif /* __F2C__ */
