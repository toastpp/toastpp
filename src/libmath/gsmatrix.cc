// ==========================================================================
// Module mathlib
// File gsmatrix.cc
// Definition of template class TGenericSparseMatrix
// ==========================================================================

#define __GSMATRIX_CC
#define MATHLIB_IMPLEMENTATION
#define DBG_TIMING

#include "mathlib.h"
#include "timing.h"

#ifdef DBG_TIMING
double cgtime = 0.0; // global CG timer
#endif

MATHLIB IterativeMethod itmethod_general = ITMETHOD_CG;
MATHLIB IterativeMethod itmethod_complex = ITMETHOD_BICGSTAB;
MATHLIB int IterCount = 0;