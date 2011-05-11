#ifndef __SPHERICALHARMONIC_ALGEBRA_H
#define __SPHERICALHARMONIC_ALGEBRA_H
#include <mathlib.h>
#include <felib.h>
#include "PN_incl.h"
#include "PN_angint.h"
#define Alm(l, m) (sqrt(toast::complex((((l)+(m))*((l)-(m)))/((double)(2*(l)+1)*(2*(l)-1)), 0)))
#define Blm(l, m) (sqrt(toast::complex((((l)+(m))*((l)+(m)-1))/((double)(2*(l)+1)*(2*(l)-1)), 0)))

double factorial(int n);

double doublefactorial(int n);

/** Computes C^{n} where C is a 'complex' number and 'n' is a positive integer
**/
toast::complex cpowi(toast::complex &c, const int m);

/** Computes (-1)^{m}
**/
int sign(int m);

/** Returns 0, 1, and 2 respectively based on whether m=0, m<0 and m>0.
**/
int signum(int m);

/* evaluates \int_{S^{2}} Y_{\ell_{1}, m_{1}} Y_{\ell_{2}, m_{2}} where both Y's are complex spherical harmonics
*/
toast::complex kronD(const IVector &a, const IVector &b);

/** Computes associated Legendre polynomials of a given order on a set of points 
*	l -> Maximum order
*	numpts -> number of points
*	pt -> the three dimensional point set
*	LT[] -> Legendre polynomials in a table form
**/
void LegendreTable(const int l, const int numpts, const RDenseMatrix &pt, RDenseMatrix* &LT);

/**Real spherical harmonics computed on a set of points
* order -> Maximum order of spherical harmonics
* numpts ->  number of points on which to evaluate
* pt -> 3D point set
* sphHarm[l] -> spherical harmonics over point set of order 'l'
**/
void sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt, RDenseMatrix* &sphHarm);

/** Y_{l, m}^{R} = a_{m}Y_{l, m} + b_{m}Y_{l, -m} where the supercript 'R' denotes the real-valued spherical harmonics.
* This function gives the value of a_{m}
**/
toast::complex am(int m);

/** Y_{l, m}^{R} = a_{m}Y_{l, m} + b_{m}Y_{l, -m} where the supercript 'R' denotes the real-valued spherical harmonics.
* This function gives the value of b_{m}
**/
toast::complex bm(int m);

/* Sin(\theta)Cos(\phi)Y_{l, m}^{R} = Sin(\theta)Cos(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}.
*/
void sincosY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c);

/* Sin(\theta)Sin(\phi)Y_{l, m}^{R} = Sin(\theta)Sin(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}. 
*/
void sinsinY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c);

/* Cos(\theta)Y_{l, m}^{R} = Cos(\theta)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m}, Y_{l+1, m}, Y_{l-1, -m}, Y_{l+1, -m}. 
*/
void cosY(const int l, const int m, CVector& e, CVector& f, IDenseMatrix& e1c, IDenseMatrix& f1c);

/* Y_{l, m}^{R} = (a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l, m}, Y_{l, -m}. 
*/
void sphY(const int l, const int m, CVector& p, IDenseMatrix& p1c);

/*
Computes integral of the form
\int (a_{0}Y_{l1, m1} + a_{1}Y_{l2, m2})(b_{0}Y_{l3, m3} + b_{1}Y_{l4, m4})
*/
double Integrate2(CVector &a, CVector &b, IDenseMatrix &a1c, IDenseMatrix &b1c);

#endif






