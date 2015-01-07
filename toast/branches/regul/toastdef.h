// -*-C++-*-
// ==========================================================================
// A list of DEFINEs which fine tune TOAST compilation options
// ==========================================================================

#ifndef __TOASTDEF_H
#define __TOASTDEF_H

#include "arch.h"

#if defined(WIN32)||defined(WIN64)
#define NEED_EXPLICIT_INSTANTIATION   // JK in Windows (VS2005) explicit instantiation is required
#endif

#if defined(WIN32)||defined(WIN64)
#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)
#define DLLEXTERN extern
#else
#define DLLEXPORT
#define DLLIMPORT
#define DLLEXTERN
#endif

#ifdef TOASTLIB
#define DLL_LIB DLLEXPORT
#else
#define DLL_LIB DLLIMPORT
#endif

//#define NEED_FRIEND_PT


#define EXPIRY_YEAR 2012
#define EXPIRY_MONTH 12
// This defines the date after which the CHECK_EXPIRY macro will
// terminate a program.

#define VERBOSE_LEVEL 1
// Log file verbosity level


#ifdef TOAST_THREAD

#define THREAD_LEVEL 2 // 0=none, 1=fine-grain, 2=coarse-grain
#define TOAST_THREAD_MATLAB_GRADIENT  // parallelise Matlab toastGradient
#define TOAST_THREAD_MATLAB_QMVEC     // parallelise Matlab Mvec
#define TOAST_THREAD_ASSEMBLE         // parallelise system matrix assembly

#else

#define THREAD_LEVEL 0

#endif


#include "blasnames.h"

#if (!defined(WIN32))&&(!defined(WIN64))
//#define USE_BLAS_LEVEL1
// Use external BLAS level 1 routines (vector-vector operations)
// instead of local C++ implementations where possible.
// This is not guaranteed to improve performance, since the function
// call overhead may outweigh any savings from optimised BLAS.

#endif

//#define USE_BLAS_LEVEL2
// Use external BLAS level 2 routines (matrix-vector operations)
// instead of local C++ implementations where possible

//#define USE_BLAS_LEVEL3
// Use external BLAS level 3 routines (matrix-matrix operations)
// instead of local C++ implementations where possible

#if !defined(WIN32)&&!defined(WIN64) // for now
#define USE_SPBLAS
// Use external sparse BLAS routines from the NIST SPARSE BLAS toolkit
#endif

#define TRI6IP_STORE_COORDS
#define TRI10IP_STORE_COORDS
#define TET10IP_STORE_COORDS

#define TRI10_STORE_COORDS

#define NSUBSAMPLE 50
// The number of sub-samples (in 1-D) used for numerical integration
// over elements. The total number of samples is n(n+1)/2 in a triangle, and
// ((n+3)*n+2)*n/6 in a tetrahedron, where n=NSUBSAMPLE

// GCC version number
#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
		     + __GNUC_MINOR__ * 100 \
		     + __GNUC_PATCHLEVEL__)
#else
#define GCC_VERSION 0
#endif

#ifdef TOAST_MPI

// Uncomment to enable MPI bindings in TFwdSolver class
//#define MPI_FWDSOLVER

#endif

#endif // !__TOASTDEF_H