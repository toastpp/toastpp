// ==========================================================================
// Module mathlib 				   Martin Schweiger - 21.5.96
// File complex.h                                                   
// Declaration and inline definitions of class complex.
// Out-of-line definitions are in complex.cc
// ==========================================================================

#ifndef __COMPLEX_H
#define __COMPLEX_H

#include <iostream>
#include <math.h>
#include "arch.h"

// ==========================================================================
// class toast::complex
// ==========================================================================

namespace toast {

/**
 * A class representing double-precision complex numbers.
 * See also \ref complex_func.
 */
class MATHLIB complex {

public:

    double re;  ///< real component
    double im;  ///< imaginary component

    /**
     * \brief Default constructor: Creates a complex number with zero real and
     *   imaginary components.
     */
    inline complex ();

    /**
     * \brief Copy constructor: Creates a complex number as a copy of another.
     * \param z source number
     */
    inline complex (const complex& z);

    /**
     * \brief Explicit constructor: Creates a complex number from the provided
     *   real and imaginary components.
     * \param r real part
     * \param i imaginary part
     */
    inline complex (double r, double i=0.0);

    /**
     * \brief Assignment operator. Implements a=b for complex numbers.
     * \param z Right-hand side operand.
     * \return Assignment result (to allow assignment chains, e.g. a=b=c)
     */
    inline complex& operator= (const complex &z); // assignment

    /**
     * \brief Complex addition.
     * \param z Right-hand side operand.
     * \return Summation result.
     * \note This operator returns *this + z, but does not modify *this.
     */
    inline complex operator+ (const complex &z) const;

    /**
     * \brief Complex subtraction.
     * \param z Right-hand side operand.
     * \return Subtraction result.
     * \note This operator returns *this - z, but does not modify *this.
     */
    inline complex operator- (const complex &z) const;

    /**
     * \brief Unary minus.
     * \return The negative of *this.
     * \note Returns (-this>-re,-this->im) but does not modify *this.
     */
    inline complex operator- () const;

    /**
     * \brief Complex multiplication.
     * \param z Right-hand side operand.
     * \return Multiplication result.
     * \note This operator returns the product *this * z, but does not modify
     *   *this.
     */
    inline complex operator* (const complex &z) const;

    /**
     * \brief Complex division.
     * \param z Right-hand side operand.
     * \return Division result.
     * \note This operator returns the quotient *this / z, but does not
     *   modify *this.
     */
    complex operator/ (const complex &z) const;

    /**
     * \brief Multiplication with real number.
     * \param x Right-hand side operand (real)
     * \return Multiplication result.
     * \note This operator returns the product *this * x, but does not modify
     *   *this.
     */
    inline complex operator* (const double &x) const;

    /**
     * \brief Division by real number.
     * \param x Right-hand side operand (real)
     * \return Division result
     * \note This operator returns the quotient *this / z, but does not
     *   modify *this.
     */
    inline complex operator/ (const double &x) const;

    /**
     * \brief Addition/assignment operator.
     * \param z Right-hand side operand
     * \return Addition result (to allow chain constructs like a = b += c)
     * \note This operator performs *this = *this+z and returns the modified
     *   *this.
     */
    inline complex operator+= (const complex &z);

    /**
     * \brief Subtraction/assignment operator.
     * \param z Right-hand side operand
     * \return Subtraction result (to allow chain constructs like a = b-= c)
     * \note This operator performs *this = *this-z and returns the modified
     *   *this.
     */
    inline complex operator-= (const complex &z);

    /**
     * \brief Multiplication/assignment operator.
     * \param z Right-hand side operand
     * \return Multiplication result (to allow chain constructs like
     *   a = b *= c)
     * \note This operator performs *this = *this * z and returns the modified
     *   *this.
     */
    inline complex operator*= (const complex &z);

    /**
     * \brief Division/assignment operator.
     * \param z Right-hand side operand
     * \return Division result (to allow chain constructs like (a = b /= c)
     * \note This operator performs *this = *this / z and returns the modified
     *   *this.
     */
    complex operator/= (const complex &z);

    /**
     * \brief Multiplication/assignment with real operand.
     * \param x Right-hand side operand (real)
     * \return Multiplication result (to allow chain constructs like
     *   a = b *= x)
     * \note This operator performs *this = *this * x and returns the modified
     *   *this.
     */
    inline complex operator*= (const double &x);

    /**
     * \brief Division/assignment with real operand.
     * \param x Right-hand side operand (real)
     * \return Division result (to allow constructs like a = b /= x)
     * \note This operator performs *this = *this / x and returns the modified
     *   *this.
     */
    inline complex operator/= (const double &x);

    /**
     * \brief Equality operator.
     * \param z Right-hand side operand
     * \return \e true if both operands are equal, \e false otherwise.
     * \note Two complex numbers are equal if their real and imaginary
     *   components are equal.
     */
    inline bool operator== (const complex z) const;

    /**
     * \brief Inequality operator.
     * \param z Right-hand side operand
     * \return \e false if both operands are equal, \e true otherwise.
     * \note Two complex numbers are unequal if either their real or
     *   imaginary parts are unequal.
     */
    inline bool operator!= (const complex z) const;

    /**
     * \brief Logical negation operator.
     * \return \e true if real and imaginary components are zero, \e false
     *   otherwise.
     */
    inline bool operator! () const;

    /**
     * \brief Relational comparison operator: greater than.
     * \param z Right-hand side operand.
     * \return \e true if left-hand operator is greater than right-hand
     *   operator, \e false otherwise.
     * \note A complex number is considered greater than another if its
     *   norm is greater.
     * \note This definition differs from Matlab, which only appears to use
     *   the real part for the comparison.
     */
    inline bool operator> (const complex z) const;

    /**
     * \brief Relational comparison operator: less than.
     * \param z Right-hand side operand.
     * \return \e true if left-hand operator is less than right-hand
     *   operator, \e false otherwise.
     * \note A complex number is considered less than another if its
     *   norm is smaller.
     * \note This definition differs from Matlab, which only appears to use
     *   the real part for the comparison.
     */
    inline bool operator< (const complex z) const;
};

} // end namespace toast


// ==========================================================================
// nonmember declarations
// ==========================================================================
/**
 * \defgroup complex_func Functions operating on complex numbers
 *
 * This section contains global-namespace functions that operate on
 * complex numbers. Where applicable, equivalent functions operating on
 * real-valued arguments are provided as well for consistency (unless
 * such functions are already part of C++, e.g. \e log).
 */
// ==========================================================================
/// @{

/**
 * \brief Returns the real part of a complex number.
 * \param z complex number
 * \return real part (z.re)
 */
inline double re (const toast::complex &z);

/**
 * \brief Returns the imaginary part of a complex number.
 * \param z complex number
 * \return imaginary part (z.im)
 */
inline double im (const toast::complex &z);

/**
 * \brief Returns the complex conjugate of a complex number.
 * \param z complex number
 * \return complex conjugate (z.re - i * z.im)
 */
inline toast::complex conj (const toast::complex &z);

/**
 * \brief Returns the argument of a complex number.
 * \param z complex number
 * \return arc tan (z.im / z.re)
 */
inline double arg (const toast::complex &z);

/**
 * \brief Returns the modulus of a complex number.
 * \param z complex number
 * \return |z| = sqrt(z.re^2 + z.im^2)
 * \sa norm
 */
inline double mod (const toast::complex &z);

/**
 * \brief Returns the norm of a complex number.
 * \param z complex number
 * \return |z| = sqrt(z.re^2 + z.im^2)
 * \sa mod
 */
inline double norm (const toast::complex &z);

/**
 * \brief Returns the square of the norm of a complex number.
 * \param z complex number
 * \return |z|^2 = z.re^2 + z.im^2
 * \sa norm
 */
inline double norm2 (const toast::complex &z);

/**
 * \brief Returns the natural logarithm of a complex number.
 * \param z complex number
 * \return log(z) = log(mod(z)) + i * arg(z)
 */
inline toast::complex log (const toast::complex &z);

/**
 * \brief Returns the square root of a complex number.
 * \param z complex number
 * \return sqrt(z)
 */
MATHLIB toast::complex sqrt (const toast::complex &z);

/**
 * \brief Returns the complex power of a complex number.
 * \param zb complex basis
 * \param ze complex exponent
 * \return zb^ze
 * \bug This function is currently not implemented and returns 0.
 * \sa pow(const toast::complex&,double)
 */
toast::complex pow (const toast::complex &zb, const toast::complex &ze);

/**
 * \brief Returns the real power of a complex number.
 * \param zb complex basis
 * \param e real exponent
 * \return zb^e
 * \bug This function is currently not implemented and returns 0.
 * \sa pow(const toast::complex&,const toast::complex&)
 */
toast::complex pow (const toast::complex &zb, double e);

/**
 * \brief Returns the exponential of a complex number.
 * \param z complex number
 * \return e^z = e^z.re (cos(z.im) + i * sin(z.im))
 */
inline toast::complex exp (toast::complex &z);

/**
 * \brief Checks a complex number for zero.
 * \param z complex number
 * \return \e true if z==0, \e false otherwise.
 */
bool iszero (toast::complex &z);

/**
 * \brief Stream input operator.
 * \param is stream object
 * \param z complex number to be assigned from the stream.
 * \return stream object
 * \note The expected format for representing a complex number in a stream
 *   is "<re im>", i.e. real and imaginary part separated by a space and
 *   enclosed in less-than and greater-than characters. If the tokens read
 *   from the stream do not match this format, the read operation fails.
 */
MATHLIB std::istream &operator>> (std::istream &is, toast::complex &z);

/**
 * \brief Stream output operator.
 * \param os stream object
 * \param z complex number to be written to the stream.
 * \return stream object
 * \note The complex number is written to the stream in the format "<re im>",
 *   i.e. real and imaginary part separated by a space and enclosed in
 *   less-than and greater-than characters.
 */
MATHLIB std::ostream &operator<< (std::ostream &os, const toast::complex &z);

///@}

// ==========================================================================
// member definitions
// ==========================================================================

inline toast::complex::complex () { re=im=0.0; }
inline toast::complex::complex (const toast::complex& z): re(z.re), im(z.im) {}
inline toast::complex::complex (double r, double i): re(r), im(i) {}

inline toast::complex& toast::complex::operator= (const toast::complex& z)
{ re=z.re, im=z.im; return *this; }

inline toast::complex toast::complex::operator+ (const toast::complex& z) const
{ return toast::complex (re+z.re, im+z.im); }

inline toast::complex toast::complex::operator- (const toast::complex& z) const
{ return toast::complex (re-z.re, im-z.im); }

inline toast::complex toast::complex::operator- () const
{ return toast::complex (-re, -im); }

inline toast::complex toast::complex::operator* (const toast::complex &z) const
{ return toast::complex (re*z.re - im*z.im, re*z.im + z.re*im); }

inline toast::complex toast::complex::operator* (const double &x) const
{ return toast::complex (re*x, im*x); }

inline toast::complex toast::complex::operator/ (const double &x) const
{ return toast::complex (re/x, im/x); }

inline toast::complex toast::complex::operator+= (const toast::complex& z)
{ re += z.re, im += z.im; return *this; }

inline toast::complex toast::complex::operator-= (const toast::complex& z)
{ re -= z.re, im -= z.im; return *this; }

inline toast::complex toast::complex::operator*= (const toast::complex& z)
{ double r = re*z.re-im*z.im;
  im = re*z.im+z.re*im, re = r; return *this; }

inline toast::complex toast::complex::operator*= (const double& x)
{  re *=x; im *= x; return *this; }

inline toast::complex toast::complex::operator/= (const double& x)
{  re /=x; im /= x; return *this; }

inline bool toast::complex::operator== (const complex z) const
{  return re == z.re && im == z.im; }

inline bool toast::complex::operator!= (const complex z) const
{  return re != z.re || im != z.im; }

inline bool toast::complex::operator! () const
{  return re == 0 && im == 0; }

inline bool toast::complex::operator> (const complex z) const
{  return norm(*this) > norm(z); }

inline bool toast::complex::operator< (const complex z) const
{  return norm(*this) < norm(z); }


// ==========================================================================
// nonmember definitions
// ==========================================================================

inline double re (const toast::complex &z)
{  return z.re; }

inline double im (const toast::complex &z)
{  return z.im; }

inline double re (const double r)
{  return r; }

inline double im (const double r)
{  return 0.0; }

inline toast::complex conj (const toast::complex& z)
{  return toast::complex (z.re, -z.im); }

inline double conj (const double r)
{  return r; }

inline double arg (const toast::complex& z)
{  return atan2 (z.im, z.re); }

inline double mod (const toast::complex& z)
{  return hypot (z.re, z.im); }

inline toast::complex log (const toast::complex &z)
{ return toast::complex (::log (mod (z)), arg (z)); }

inline double norm (const toast::complex& z)
{ return mod(z); }

inline double norm2 (const toast::complex& z)
{ return z.re*z.re + z.im*z.im; }

inline toast::complex exp (toast::complex& z)
{
  double r = ::exp (z.re);
  return toast::complex (r * cos (z.im), r * sin (z.im));
}

inline bool iszero (toast::complex &z)
{ return z.re == 0.0 && z.im == 0.0; }

inline toast::complex hadamard (const toast::complex& a,
    const toast::complex & b)
{ return toast::complex (a.re * b.re, a.im * b.im); }

#endif // !__COMPLEX_H
