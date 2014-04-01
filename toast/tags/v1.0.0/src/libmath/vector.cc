// ==========================================================================
// Module mathlib
// File vector.cc
// Definition of template class TVector ('template vector')
// ==========================================================================

#define __VECTOR_CC
#define MATHLIB_IMPLEMENTATION

#include <math.h>
#include <cstdlib>
#include "mathlib.h"

using namespace std;
using namespace toast;

/* Explicit complex conversions */
/* These ought to be friends really, except that I can't see how to do that
when using template */

MATHLIB TVector<float> Re (const TVector<scomplex> &vec)
{
    TVector<float> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

MATHLIB TVector<double> Re (const TVector<toast::complex> &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

MATHLIB TVector<float> Im (const TVector<scomplex> &vec)
{
    TVector<float> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

MATHLIB TVector<double> Im (const TVector<toast::complex> &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

MATHLIB TVector<double> Mod (const TVector<toast::complex> &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = mod(vec[i]);
    return tmp;
}

MATHLIB TVector<double> LogMod (const TVector<toast::complex> &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = log(mod(vec[i]));
    return tmp;
}

MATHLIB TVector<double> Arg (const TVector<toast::complex> &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = arg(vec[i]);
    return tmp;
}

MATHLIB TVector<toast::complex> Conj (const TVector<toast::complex> &vec)
{
	TVector<toast::complex> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = conj(vec[i]);
    return tmp;
}

MATHLIB void SelfConj (const TVector<toast::complex> &vec)
{
  /* version converts this vector to conjugate */
    for (int i = 0; i < vec.Dim(); i++) vec[i] = conj(vec[i]);
}

MATHLIB TVector<toast::complex> Hadamard (const TVector<toast::complex> &a, const TVector<toast::complex> &b)
{
    dASSERT(a.Dim() == b.Dim(), "Dimension mismatch");
	TVector<toast::complex> tmp(a.Dim());
    for (int i = 0; i < a.Dim(); i++) tmp[i] = hadamard(a[i],b[i]);
    return tmp;
}

MATHLIB void SetReal (TVector<toast::complex> &z, const TVector<double> &zre)
{
    dASSERT(z.Dim() == zre.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++) z[i].re = zre[i];
}

MATHLIB void SetImag (TVector<toast::complex> &z, const TVector<double> &zim)
{
    dASSERT(z.Dim() == zim.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++) z[i].im = zim[i];
}

MATHLIB TVector<toast::complex> MakeCVector (const TVector<scomplex> &v)
{
	TVector<toast::complex> c(v.Dim());
    for (int i = 0; i < v.Dim(); i++)
	c[i] = (toast::complex)v[i];
    return c;
}
