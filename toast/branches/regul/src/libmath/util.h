// -*-C++-*-
// ============================================================================
// TOAST                                              (c) Martin Schweiger 2000
// Library: libmath
// File:    util
//
// General-purpose utility routines
// ============================================================================

#ifndef __UTIL_H
#define __UTIL_H

#include "vector.h"

// string and file parsing auxiliary routines
char *ToupperString (char *str);
void SplitEqString (char *str, char **cat, char **val);

// ============================================================================
// Write an image (stored in 1-D data array 'img') of dimension xdim,ydim
// to a PPM pixmap file. Image is scaled to range scalemin-scalemax.
// If scalemin and/or scalemax is NULL, then autscaling is used.

void WritePPM (const RVector &img, int xdim, int ydim,
    double *scalemin, double *scalemax, char *fname);

#endif // !__UTIL_H
