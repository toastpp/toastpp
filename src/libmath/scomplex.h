// ==========================================================================
// Module mathlib 				   Martin Schweiger - 21.5.96
// File scomplex.h                                                   
// Declaration and inline definitions of class scomplex.
// Out-of-line definitions are in scomplex.cc
// single precision version of complex class	   Simon Arridge - 8.5.05
// ==========================================================================

#ifndef __SCOMPLEX_H
#define __SCOMPLEX_H

#include <iostream>
#include "complex.h"
#include <math.h>

// ==========================================================================
// declaration
// ==========================================================================

class scomplex {

public:

    float re, im;                          // real and imaginary components

    scomplex ();                             // default constructor
    scomplex (const scomplex& z);            // copy constructor
    scomplex (const toast::complex& z);      // copy constructor
    scomplex (double r, double i=0.0);       // elementwise initialisation

    operator toast::complex()
	{ return toast::complex((double)re, (double)im); }
    // explicit cast operator scomplex->complex

    scomplex& operator= (const scomplex &z);    // assignment
    //scomplex& operator= (const complex &z);    // assignment

    // arithmetic operators
    scomplex operator+ (const scomplex &z) const;
    scomplex operator- (const scomplex &z) const;
    scomplex operator- () const;		     // unary -
    scomplex operator* (const scomplex &z) const;
    scomplex operator/ (const scomplex &z) const;
    scomplex operator* (const double &x) const;
    scomplex operator/ (const double &x) const;

    scomplex operator+= (const scomplex &z);   // scomplex addition/assignment
    scomplex operator-= (const scomplex &z);   // scomplex subtraction/ass.
    scomplex operator*= (const scomplex &z);   // scomplex multiplication/ass.
    scomplex operator/= (const scomplex &z);   // scomplex division/ass.
    scomplex operator*= (const double &x);   // scomplex multiplication/ass.
    scomplex operator/= (const double &x);   // scomplex division/ass.

    unsigned char operator== (scomplex z) const { return re==z.re && im==z.im; }
    unsigned char operator!= (scomplex z) const { return re!=z.re || im!=z.im; }
    unsigned char operator! () const { return re==0 && im==0; }
    /* some hacked comparision operators DANGER !!! */
    unsigned char operator> (scomplex z) const { return re>z.re && im>z.im; }
    unsigned char operator< (scomplex z) const { return re<z.re && im<z.im; }
    // arithmetic operators double -> single
    scomplex operator+ (const toast::complex &z) const;
    scomplex operator- (const toast::complex &z) const;
    scomplex operator* (const toast::complex &z) const;
    scomplex operator/ (const toast::complex &z) const;

    scomplex operator+= (const toast::complex &z);
    // scomplex addition/assignment
    scomplex operator-= (const toast::complex &z);
    // scomplex subtraction/ass.
    scomplex operator*= (const toast::complex &z);
    // scomplex multiplication/ass.
    scomplex operator/= (const toast::complex &z);
    // scomplex division/ass.

    unsigned char operator== (toast::complex z) const
	{ return re==z.re && im==z.im; }
    unsigned char operator!= (toast::complex z) const
	{ return re!=z.re || im!=z.im; }
    /* some hacked comparision operators DANGER !!! */
    unsigned char operator> (toast::complex z) const
	{ return re>z.re && im>z.im; }
    // friends
    friend double re (const scomplex &z);     // return real part
    friend double im (const scomplex &z);     // return imaginary part.
    friend scomplex conj (const scomplex &z);  // scomplex conjugate
    friend double norm2 (const scomplex &z);  // square of norm
    friend double norm (const scomplex &z);
    friend scomplex sqrt (const scomplex &z);
    friend double arg (const scomplex &z);
    friend double mod (const scomplex &z);
    friend scomplex log (const scomplex &z);
    friend scomplex pow (const scomplex &zb, const scomplex &ze);
    friend scomplex pow (const scomplex &zb, double e);
    friend scomplex exp (scomplex &z);
    friend bool iszero (scomplex &z) { return z.re == 0.0 && z.im == 0.0; }
    friend std::istream &operator>> (std::istream &is, scomplex &z);
    friend std::ostream &operator<< (std::ostream &os, const scomplex &z);
};


// ==========================================================================
// member definitions
// ==========================================================================

inline scomplex::scomplex () { re=im=0.0; }
inline scomplex::scomplex (const scomplex& z): re(z.re), im(z.im) {}
inline scomplex::scomplex (const toast::complex& z): re((float)z.re), im((float)z.im) {}
inline scomplex::scomplex (double r, double i): re((float)r), im((float)i) {}

inline scomplex& scomplex::operator= (const scomplex& z)
{ re=z.re, im=z.im; return *this; }

inline scomplex scomplex::operator+ (const scomplex& z) const
{ return scomplex (re+z.re, im+z.im); }

inline scomplex scomplex::operator- (const scomplex& z) const
{ return scomplex (re-z.re, im-z.im); }

inline scomplex scomplex::operator- () const
{ return scomplex (-re, -im); }

inline scomplex scomplex::operator* (const scomplex &z) const
{ return scomplex (re*z.re - im*z.im, re*z.im + z.re*im); }

inline scomplex scomplex::operator* (const double &x) const
{ return scomplex (re*x, im*x); }

inline scomplex scomplex::operator/ (const double &x) const
{ return scomplex (re/x, im/x); }

inline scomplex scomplex::operator+= (const scomplex& z)
{ re += z.re, im += z.im; return *this; }

inline scomplex scomplex::operator-= (const scomplex& z)
{ re -= z.re, im -= z.im; return *this; }

inline scomplex scomplex::operator*= (const scomplex& z)
{ float r = re*z.re-im*z.im;
  im = re*z.im+z.re*im, re = r; return *this; }

inline scomplex scomplex::operator*= (const double& x)
{  re *= (float)x; im *= (float)x; return *this; }

inline scomplex scomplex::operator/= (const double& x)
{  re /= (float)x; im /= (float)x; return *this; }

//--------------------- duble -> single versions 
//inline scomplex& scomplex::operator= (const complex& z)
//{ re=z.re, im=z.im; return *this; }

inline scomplex scomplex::operator+ (const toast::complex& z) const
{ return scomplex (re+z.re, im+z.im); }

inline scomplex scomplex::operator- (const toast::complex& z) const
{ return scomplex (re-z.re, im-z.im); }

inline scomplex scomplex::operator* (const toast::complex &z) const
{ return scomplex (re*z.re - im*z.im, re*z.im + z.re*im); }

inline scomplex scomplex::operator+= (const toast::complex& z)
{ re += (float)z.re, im += (float)z.im; return *this; }

inline scomplex scomplex::operator-= (const toast::complex& z)
{ re -= (float)z.re, im -= (float)z.im; return *this; }

inline scomplex scomplex::operator*= (const toast::complex& z)
{ double r = re*z.re-im*z.im;
  im = (float)(re*z.im+z.re*im), re = (float)r; return *this; }


// ==========================================================================
// friend definitions
// ==========================================================================

inline double re (const scomplex &z)
{ return z.re; }

inline double im (const scomplex &z)
{ return z.im; }


inline scomplex conj (const scomplex& z)
{ return scomplex (z.re, -z.im); }
/* --- these are defined in complex.h 
inline double re (const double r)
{ return r; }

inline double im (const double r)
{ return 0.0; }

inline double conj (const double r)
{ return r; }
*/

inline double arg (const scomplex& z)
{ return atan2 (z.im, z.re); }

inline double mod (const scomplex& z)
{ return hypot (z.re, z.im); }

inline scomplex log (const scomplex &z)
{ return scomplex (log (mod (z)), arg (z)); }

inline double norm2 (const scomplex& z)
{ return z.re*z.re + z.im*z.im; }

inline double norm (const scomplex& z)
{ return hypot (z.re, z.im); }

inline scomplex exp (scomplex& z)
{
  double r = ::exp (z.re);
  return scomplex (r * cos (z.im), r * sin (z.im));
}
inline scomplex hadamard (const scomplex& a, const scomplex & b)
{ return scomplex (a.re * b.re, a.im * b.im); }

#endif // !__SCOMPLEX_H
