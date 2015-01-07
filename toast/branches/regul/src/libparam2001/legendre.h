#ifndef __LEGENDRE_H
#define __LEGENDRE_H


//***************************************************************************
//******** Legendre Accosiated Polynomials and Spherical Harmonic ***********
//***************************************************************************
//***************************************************************************

int factorial(int i);

double plgndr(int l, int m, double x);

std::complex<double> SphericalHarmonic(int l,int m, double thi, double fi); // calculates the spherical harmonic
                                                               //of degree l order m
                                                               // at the (thi,fi) point
#endif
