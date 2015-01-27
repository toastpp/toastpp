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

/* Explicit complex conversions */
/* These ought to be friends really, except that I can't see how to do that
when using template */

MATHLIB TVector<float> Re (const TVector<std::complex<float> > &vec)
{
    TVector<float> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].real();
    return tmp;
}

MATHLIB TVector<double> Re (const TVector<std::complex<double> > &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].real();
    return tmp;
}

MATHLIB TVector<float> Im (const TVector<std::complex<float> > &vec)
{
    TVector<float> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].imag();
    return tmp;
}

MATHLIB TVector<double> Im (const TVector<std::complex<double> > &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].imag();
    return tmp;
}

MATHLIB TVector<double> Mod (const TVector<std::complex<double> > &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = std::abs(vec[i]);
    return tmp;
}

MATHLIB TVector<double> LogMod (const TVector<std::complex<double> > &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = log(std::abs(vec[i]));
    return tmp;
}

MATHLIB TVector<double> Arg (const TVector<std::complex<double> > &vec)
{
    TVector<double> tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = arg(vec[i]);
    return tmp;
}

MATHLIB TVector<std::complex<double> > Conj (
    const TVector<std::complex<double> > &vec)
{
    TVector<std::complex<double> > tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = toast::conj(vec[i]);
    return tmp;
}

MATHLIB void SelfConj (const TVector<std::complex<double> > &vec)
{
  /* version converts this vector to conjugate */
    for (int i = 0; i < vec.Dim(); i++) vec[i] = toast::conj(vec[i]);
}

MATHLIB TVector<std::complex<double> > Hadamard (
    const TVector<std::complex<double> > &a,
    const TVector<std::complex<double> > &b)
{
    dASSERT(a.Dim() == b.Dim(), "Dimension mismatch");
    TVector<std::complex<double> > tmp(a.Dim());
    for (int i = 0; i < a.Dim(); i++) tmp[i] = hadamard(a[i],b[i]);
    return tmp;
}

MATHLIB void SetReal (TVector<std::complex<double> > &z,
    const TVector<double> &zre)
{
    dASSERT(z.Dim() == zre.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++)
        z[i] = zre[i];
}

MATHLIB void SetImag (TVector<std::complex<double> > &z,
    const TVector<double> &zim)
{
    dASSERT(z.Dim() == zim.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++)
      z[i] = std::complex<double> (z[i].real(), zim[i]);
}

MATHLIB TVector<std::complex<double> > MakeCVector (
    const TVector<std::complex<float> > &v)
{
    TVector<std::complex<double> > c(v.Dim());
    for (int i = 0; i < v.Dim(); i++)
        c[i] = (std::complex<double>)v[i];
    return c;
}
