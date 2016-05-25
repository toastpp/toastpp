/**************************************************************************
 * 
 * Copyright 2008 Tungsten Graphics, Inc., Cedar Park, Texas.
 * All Rights Reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sub license, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial portions
 * of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT.
 * IN NO EVENT SHALL TUNGSTEN GRAPHICS AND/OR ITS SUPPLIERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 **************************************************************************/


/**
 * Math utilities and approximations for common math functions.
 * Reduced precision is usually acceptable in shaders...
 *
 * "fast" is used in the names of functions which are low-precision,
 * or at least lower-precision than the normal C lib functions.
 */


#ifndef U_MATH_H
#define U_MATH_H


#include "pipe/p_compiler.h"
#include "util/u_debug.h"


#ifdef __cplusplus
extern "C" {
#endif


#if defined(PIPE_SUBSYSTEM_WINDOWS_MINIPORT)
__inline double ceil(double val)
{
   double ceil_val;

   if((val - (long) val) == 0) {
      ceil_val = val;
   }
   else {
      if(val > 0) {
         ceil_val = (long) val + 1;
      }
      else {
         ceil_val = (long) val;
      }
   }

   return ceil_val;
}

#ifndef PIPE_SUBSYSTEM_WINDOWS_CE_OGL
__inline double floor(double val)
{
   double floor_val;

   if((val - (long) val) == 0) {
      floor_val = val;
   }
   else {
      if(val > 0) {
         floor_val = (long) val;
      }
      else {
         floor_val = (long) val - 1;
      }
   }

   return floor_val;
}
#endif

#pragma function(pow)
__inline double __cdecl pow(double val, double exponent)
{
   /* XXX */
   assert(0);
   return 0;
}

#pragma function(log)
__inline double __cdecl log(double val)
{
   /* XXX */
   assert(0);
   return 0;
}

#pragma function(atan2)
__inline double __cdecl atan2(double val)
{
   /* XXX */
   assert(0);
   return 0;
}
#else
#include <math.h>
#include <stdarg.h>
#endif


#if defined(_MSC_VER) 

#if _MSC_VER < 1400 && !defined(__cplusplus) || defined(PIPE_SUBSYSTEM_WINDOWS_CE)
 
static INLINE float cosf( float f ) 
{
   return (float) cos( (double) f );
}

static INLINE float sinf( float f ) 
{
   return (float) sin( (double) f );
}

static INLINE float ceilf( float f ) 
{
   return (float) ceil( (double) f );
}

static INLINE float floorf( float f ) 
{
   return (float) floor( (double) f );
}

static INLINE float powf( float f, float g ) 
{
   return (float) pow( (double) f, (double) g );
}

static INLINE float sqrtf( float f ) 
{
   return (float) sqrt( (double) f );
}

static INLINE float fabsf( float f ) 
{
   return (float) fabs( (double) f );
}

static INLINE float logf( float f ) 
{
   return (float) log( (double) f );
}

#else
/* Work-around an extra semi-colon in VS 2005 logf definition */
#ifdef logf
#undef logf
#define logf(x) ((float)log((double)(x)))
#endif /* logf */
#endif

static INLINE double log2( double x )
{
   const double invln2 = 1.442695041;
   return log( x ) * invln2;
}

#endif /* _MSC_VER */





#define POW2_TABLE_SIZE_LOG2 9
#define POW2_TABLE_SIZE (1 << POW2_TABLE_SIZE_LOG2)
#define POW2_TABLE_OFFSET (POW2_TABLE_SIZE/2)
#define POW2_TABLE_SCALE ((float)(POW2_TABLE_SIZE/2))
extern float pow2_table[POW2_TABLE_SIZE];



extern void
util_init_math(void);


union fi {
   float f;
   int32_t i;
   uint32_t ui;
};


/**
 * Fast version of 2^x
 * Identity: exp2(a + b) = exp2(a) * exp2(b)
 * Let ipart = int(x)
 * Let fpart = x - ipart;
 * So, exp2(x) = exp2(ipart) * exp2(fpart)
 * Compute exp2(ipart) with i << ipart
 * Compute exp2(fpart) with lookup table.
 */
static INLINE float
util_fast_exp2(float x)
{
   int32_t ipart;
   float fpart, mpart;
   union fi epart;
   
   if(x > 129.00000f)
      return 3.402823466e+38f;
   
   if(x < -126.99999f)
      return 0.0f;

   ipart = (int32_t) x;
   fpart = x - (float) ipart;
   
   /* same as
    *   epart.f = (float) (1 << ipart)
    * but faster and without integer overflow for ipart > 31 */
   epart.i = (ipart + 127 ) << 23;
   
   mpart = pow2_table[POW2_TABLE_OFFSET + (int)(fpart * POW2_TABLE_SCALE)];
   
   return epart.f * mpart;
}


/**
 * Fast approximation to exp(x).
 */
static INLINE float
util_fast_exp(float x)
{
   const float k = 1.44269f; /* = log2(e) */
   return util_fast_exp2(k * x);
}


#define LOG2_TABLE_SIZE_LOG2 16
#define LOG2_TABLE_SCALE (1 << LOG2_TABLE_SIZE_LOG2)
#define LOG2_TABLE_SIZE (LOG2_TABLE_SCALE + 1)
extern float log2_table[LOG2_TABLE_SIZE];


static INLINE float
util_fast_log2(float x)
{
   union fi num;
   float epart, mpart;
   num.f = x;
   epart = (float)(((num.i & 0x7f800000) >> 23) - 127);
   /* mpart = log2_table[mantissa*LOG2_TABLE_SCALE + 0.5] */
   mpart = log2_table[((num.i & 0x007fffff) + (1 << (22 - LOG2_TABLE_SIZE_LOG2))) >> (23 - LOG2_TABLE_SIZE_LOG2)];
   return epart + mpart;
}


static INLINE float
util_fast_pow(float x, float y)
{
   return util_fast_exp2(util_fast_log2(x) * y);
}



/**
 * Floor(x), returned as int.
 */
static INLINE int
util_ifloor(float f)
{
   int ai, bi;
   double af, bf;
   union fi u;
   af = (3 << 22) + 0.5 + (double)f;
   bf = (3 << 22) + 0.5 - (double)f;
   u.f = (float) af;  ai = u.i;
   u.f = (float) bf;  bi = u.i;
   return (ai - bi) >> 1;
}


/**
 * Round float to nearest int.
 */
static INLINE int
util_iround(float f)
{
#if defined(PIPE_CC_GCC) && defined(PIPE_ARCH_X86) 
   int r;
   __asm__ ("fistpl %0" : "=m" (r) : "t" (f) : "st");
   return r;
#elif defined(PIPE_CC_MSVC) && defined(PIPE_ARCH_X86)
   int r;
   _asm {
	 fld f
	 fistp r
	}
   return r;
#else
   if (f >= 0.0f)
      return (int) (f + 0.5f);
   else
      return (int) (f - 0.5f);
#endif
}



/**
 * Test if x is NaN or +/- infinity.
 */
static INLINE boolean
util_is_inf_or_nan(float x)
{
   union fi tmp;
   tmp.f = x;
   return !(int)((unsigned int)((tmp.i & 0x7fffffff)-0x7f800000) >> 31);
}


/**
 * Find first bit set in word.  Least significant bit is 1.
 * Return 0 if no bits set.
 */
#if defined(_MSC_VER) && _MSC_VER >= 1300
static INLINE
unsigned long ffs( unsigned long u )
{
   unsigned long i;
   if(_BitScanForward(&i, u))
      return i + 1;
   else
      return 0;
}
#elif defined(PIPE_CC_MSVC) && defined(PIPE_ARCH_X86)
static INLINE
unsigned ffs( unsigned u )
{
   unsigned i;

   if( u == 0 ) {
      return 0;
   }

   __asm bsf eax, [u]
   __asm inc eax
   __asm mov [i], eax

   return i;
}
#elif defined(__MINGW32__)
#define ffs __builtin_ffs
#endif


/**
 * Return float bits.
 */
static INLINE unsigned
fui( float f )
{
   union fi fi;
   fi.f = f;
   return fi.ui;
}



static INLINE float
ubyte_to_float(ubyte ub)
{
   return (float) ub * (1.0f / 255.0f);
}


/**
 * Convert float in [0,1] to ubyte in [0,255] with clamping.
 */
static INLINE ubyte
float_to_ubyte(float f)
{
   const int ieee_0996 = 0x3f7f0000;   /* 0.996 or so */
   union fi tmp;

   tmp.f = f;
   if (tmp.i < 0) {
      return (ubyte) 0;
   }
   else if (tmp.i >= ieee_0996) {
      return (ubyte) 255;
   }
   else {
      tmp.f = tmp.f * (255.0f/256.0f) + 32768.0f;
      return (ubyte) tmp.i;
   }
}



#define CLAMP( X, MIN, MAX )  ( (X)<(MIN) ? (MIN) : ((X)>(MAX) ? (MAX) : (X)) )

#define MIN2( A, B )   ( (A)<(B) ? (A) : (B) )
#define MAX2( A, B )   ( (A)>(B) ? (A) : (B) )


static INLINE int
align(int value, int alignment)
{
   return (value + alignment - 1) & ~(alignment - 1);
}


#ifndef COPY_4V
#define COPY_4V( DST, SRC )         \
do {                                \
   (DST)[0] = (SRC)[0];             \
   (DST)[1] = (SRC)[1];             \
   (DST)[2] = (SRC)[2];             \
   (DST)[3] = (SRC)[3];             \
} while (0)
#endif


#ifndef COPY_4FV
#define COPY_4FV( DST, SRC )  COPY_4V(DST, SRC)
#endif


#ifndef ASSIGN_4V
#define ASSIGN_4V( DST, V0, V1, V2, V3 ) \
do {                                     \
   (DST)[0] = (V0);                      \
   (DST)[1] = (V1);                      \
   (DST)[2] = (V2);                      \
   (DST)[3] = (V3);                      \
} while (0)
#endif


#ifdef __cplusplus
}
#endif

#endif /* U_MATH_H */
