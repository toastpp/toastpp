#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include "mathlib.h"

using namespace std;
using namespace toast;

// division and sqrt algorithms from NR p.176 and p.948

scomplex scomplex::operator/ (const scomplex &z) const
{
    double r, den;
    if (fabs (z.re) >= fabs (z.im)) {
	r = z.im / z.re;
	den = z.re + r*z.im;
	return scomplex ((re+r*im)/den, (im-r*re)/den);
    } else {
	r = z.re / z.im;
	den = z.im + r*z.re;
	return scomplex ((re*r+im)/den, (im*r-re)/den);
    }
}

scomplex scomplex::operator/= (const scomplex& z)
{
    float r, den, nr;
    if (fabs (z.re) >= fabs (z.im)) {
	r = z.im / z.re;
	den = z.re + r*z.im;
	nr = (re+r*im)/den;
	im = (im-r*re)/den;
	re = nr;
    } else {
	r = z.re / z.im;
	den = z.im + r*z.re;
	nr = (re*r+im)/den;
	im = (im*r-re)/den;
	re = nr;
    }
    return *this;
}

scomplex sqrt (const scomplex& z)
{
    if (!z.re && !z.im) return scomplex(0.0, 0.0);
    double r, w, x = fabs (z.re), y = fabs (z.im);
    if (x >= y) {
	r = y/x;
	w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0+r*r)));
    } else {
	r = x/y;
	w = sqrt (y) * sqrt (0.5 * (r + sqrt (1.0+r*r)));
    }
    if (z.re >= 0.0) return scomplex (w, z.im/(2.0*w));
    if (z.im >= 0.0) return scomplex (z.im/(2.0*w), w);
    return scomplex (-z.im/(2.0*w), -w);
}

scomplex pow (const scomplex& /*zb*/, const scomplex& /*ze*/)
{
    xERROR(Function not implemented!);
    return (scomplex(0,0));
}

scomplex pow (const scomplex& /*zb*/, double /*e*/)
{
    xERROR(Function not implemented!);
    return (scomplex(0,0));
}

ostream &operator<< (ostream &os, const scomplex &z)
{
    os << '<' << z.re << ' ' << z.im << '>';
    return os;
}

istream &operator>> (istream& is, scomplex& z)
{
    char c;
    for (;;) {
	is >> c;
	if (c == '<') break;
	if (c == ' ' || c == '\t' || c == '\n') continue;
	is.unget();
	is.setstate (ios::failbit);
	return is;
    }

    is >> z.re >> z.im;
    do { is >> c; } while (is.good() && c != '>');
    return is;
}
/*--------------- double -> single versions */
scomplex scomplex::operator/ (const complex &z) const
{
    double r, den;
    if (fabs (z.re) >= fabs (z.im)) {
	r = z.im / z.re;
	den = z.re + r*z.im;
	return scomplex ((re+r*im)/den, (im-r*re)/den);
    } else {
	r = z.re / z.im;
	den = z.im + r*z.re;
	return scomplex ((re*r+im)/den, (im*r-re)/den);
    }
}

scomplex scomplex::operator/= (const complex& z)
{
    double r, den, nr;
    if (fabs (z.re) >= fabs (z.im)) {
	r = z.im / z.re;
	den = z.re + r*z.im;
	nr = (re+r*im)/den;
	im = (float)((im-r*re)/den);
	re = (float)nr;
    } else {
	r = z.re / z.im;
	den = z.im + r*z.re;
	nr = (re*r+im)/den;
	im = (float)((im*r-re)/den);
	re = (float)nr;
    }
    return *this;
}
