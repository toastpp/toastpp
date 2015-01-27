// -*-C++-*-
// ==========================================================================
// Module mathlib
// File mathdef.h
// General maths types and macros
// ==========================================================================

#ifndef __MATHDEF_H
#define __MATHDEF_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <complex>

#ifdef FEM_DEBUG
#define MATH_DEBUG  		// switch on error checks
#define NO_TEMPLATE_INLINE	// switch off inlines in template classes
#endif

//#define COMPUTE_FLOPS         // record floating point operations for
                                // selected calculations

// matrix/vector index types
typedef int idxtype;
//typedef long idxtype;

// some constants ===========================================================

const double Pi   = 3.1415926535897932384626433832795;
const double Pi2  = Pi*2.0;
const double Pi05 = Pi*0.5;
const double Rad  = Pi/180.0;
const double Deg  = 180.0/Pi;

// type definitions =========================================================

typedef unsigned char BYTE;

// macros ===================================================================

//#ifndef __MINMAX_DEFINED
//namespace toast {
//    template<class T>
//    inline T min (const T x, const T y) { return (x < y) ? x : y; }
//    template<class T>
//    inline T max (const T x, const T y) { return (x > y) ? x : y; }
//}
//#endif

typedef enum {
    DEEP_COPY,
    SHALLOW_COPY
} CopyMode;

template<class T>
inline T sqr (const T a) { return (a == 0.0 ? (T)0 : a*a); }

inline double sign (const double a, const double b)
{ return (b >= 0 ? fabs(a) : -fabs(a)); }
// return a with the sign of b

inline double norm (double f)
{ return fabs (f); }

inline double norm2 (double f)
{ return f*f; }

inline int norm (int i)
{ return abs (i); }

inline bool pow (bool base, bool e)
{ return base; }

inline int pow (int base, int e)
{ return (int)pow ((double)base, (double)e); }

// inline float pow (float base, float e)
// { return (float)pow ((double)base, (double)e); }




//inline float conj(float a)
//{ return a; }

//inline double conj(double a)
//{ return a; }



// SP 24.01.15: OS X LIBC++ provides definitions of conj
// for real types which return a std::complex. This is incompatible
// with the gmres defined in gmres_imp.hpp. The following function
// replicate original TOAST (LIBSTDC++) behaviour.

namespace toast {
    
    inline float conj(float a)
    { return a; }
    
    inline double conj(double a)
    { return a; }
    
    template<typename T> inline std::complex<T> conj(std::complex<T> a)
    { return std::conj(a);}
    
}



inline bool iszero (double d)
{ return d == 0.0; }

inline bool iszero (std::complex<double> d)
{
    return d.real() == 0 && d.imag() == 0;
}

inline bool operator> (const std::complex<double> &z1,
		       const std::complex<double> &z2)
{
    return std::norm(z1) > std::norm(z2);
}

inline bool operator< (const std::complex<double> &z1,
		       const std::complex<double> &z2)
{
    return std::norm(z1) < std::norm(z2);
}

inline bool operator! (const std::complex<double> &z)
{
    return z.real() == 0 && z.imag() == 0;
}

inline std::complex<double> hadamard (const std::complex<double>& a,
    const std::complex<double>& b)
{
    return std::complex<double> (a.real()*b.real(), a.imag()*b.imag());
}

inline int binomial_coeff (int n, int r)
{
    int i, sum = 1;
	for (i = n; i > std::max (r, (n-r)); i--) sum *= i;
	for (i = std::min (r, (n-r)); i > 1; i--) sum /= i;
    return sum;
}

inline double rad (double deg) { return deg*Rad; }
inline double deg (double rad) { return rad*Deg; }

#ifdef __BORLANDC__			// BorlandC doesn't know about atanh
inline double atanh (double x)
{ return 0.5 * log ((1.0+x) / (1.0-x)); }	// valid for |x| < 1
#endif

#endif
