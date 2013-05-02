// ========================================================================= //
// TOAST v.15                                      (c) Martin Schweiger 1999 //
// Library: libfe     File: tri10_ip.cc                                      //
//                                                                           //
// 10-noded isoparametric triangle for cubic shape functions and curved      //
// boundaries.                                                               //
// ========================================================================= //

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "lin_qr.h"
#include "tri_qr.h"

Triangle10_ip::Triangle10_ip (const Triangle10_ip &el)
: Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Triangle10_ip::Copy ()
{
    return new Triangle10_ip(*this);
}

void Triangle10_ip::Initialise (const NodeList &nlist)
{
#ifdef TRI10IP_STORE_COORDS
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
    n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1];
    n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1];
    n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1];
    n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1];
    n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1];
    n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1];
    n9x = nlist[Node[9]][0], n9y = nlist[Node[9]][1];
#endif

    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1];

    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    Element_Unstructured_2D::Initialise (nlist);

    intff.Zero (10);
    intff = ComputeIntFF (nlist);
}

int Triangle10_ip::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Side index out of range");
    dASSERT(node >= 0 && node < 4, "Node index out of range");
    static int SN[3][4] = {{0,1,3,4},{1,2,5,6},{2,0,7,8}};
    return SN[side][node];
}

Point Triangle10_ip::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    static const double EPS = 1e-4;
    static const int MAXIT = 10;

    Point loc(2), lglob(2);
    RVector dloc(2), fun(10);
    RDenseMatrix IJ(2,2);

    // First do linear inverse mapping
    double scale = 1.0/(a0+a1+a2);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1]) * scale;

    for (int i = 0; i < MAXIT; i++) {
        // map forward again
        fun = LocalShapeF (loc);
	IJacobian (loc, IJ);
	lglob[0] = lglob[1] = 0.0;
	for (int j = 0; j < 10; j++)
	    lglob += nlist[Node[j]] * fun[j];

	// apply correction
	dloc = transpose(IJ) * (glob-lglob);
	loc += dloc;
	if (length(dloc) < EPS) return loc;
    }
    loc[0] = loc[1] = -1.0;
    // if no convergence, assume point is outside element - HACK!
    return loc;
}

Point Triangle10_ip::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 10, "Node index out of range");
    static Point2D nloc[10] = {
        Point2D(0,0), Point2D(1,0), Point2D(0,1),
	Point2D(1.0/3.0,0),       Point2D(2.0/3.0,0),
	Point2D(2.0/3.0,1.0/3.0), Point2D(1.0/3.0,2.0/3.0),
	Point2D(0,2.0/3.0),       Point2D(0,1.0/3.0),
	Point2D(1.0/3.0,1.0/3.0)
    };
    return nloc[node];
}

RVector Triangle10_ip::DirectionCosine (int side, RDenseMatrix &jacin)
{
    // INVALID
    RVector cosin(2);
    switch(side) {
    case 0: cosin[0] = -b2, cosin[1] = -c2; break;
    case 1: cosin[0] = -b0, cosin[1] = -c0; break;
    case 2: cosin[0] = -b1, cosin[1] = -c1; break;
    default: xERROR("Side index out of range");
    }
    return cosin/length(cosin);
}

const RVector &Triangle10_ip::LNormal (int side) const
{
    static const RVector lnm0 = RVector (2, "0 -1");
    static const RVector lnm1 = RVector (2, "0.7071067814 0.7071067814");
    static const RVector lnm2 = RVector (2, "-1 0");
    static const RVector *lnm[3] = {
      &lnm0, &lnm1, &lnm2
    };
    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
    return *lnm[side];
}

bool Triangle10_ip::LContains (const Point& loc, bool pad) const
{
    dASSERT(loc.Dim() == 2, "Local point must be 2D.");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 &&
	    loc[0]+loc[1]-EPS <= 1.0);
    } else {
 	return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[0]+loc[1] <= 1.0);
   }
}

bool Triangle10_ip::GContains (const Point& glob, const NodeList& nlist) const
{
    return LContains (Local (nlist, glob));
}

RVector Triangle10_ip::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RVector fun(10);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    double f0 = L0*(3.0*L0-1.0);
    double f1 = L1*(3.0*L1-1.0);
    double f2 = L2*(3.0*L2-1.0);
    fun[0] = 0.5 * f0 * (3.0*L0-2.0);
    fun[1] = 0.5 * f1 * (3.0*L1-2.0);
    fun[2] = 0.5 * f2 * (3.0*L2-2.0);
    fun[3] = 4.5 * L1 * f0;
    fun[4] = 4.5 * L0 * f1;
    fun[5] = 4.5 * L2 * f1;
    fun[6] = 4.5 * L1 * f2;
    fun[7] = 4.5 * L0 * f2;
    fun[8] = 4.5 * L2 * f0;
    fun[9] = 27.0 * L0 * L1 * L2;
    return fun;
}

RDenseMatrix Triangle10_ip::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der (2,10);
    double xi = loc[0], eta = loc[1];
    double s = 1.0-xi-eta;
    double f1 = 3.0*s-1.0,   g1 = 3.0*s-2.0,    h1 = 13.5*xi*eta;
    double f2 = 3.0*xi-1.0,  g2 = 3.0*xi-2.0,   h2 = 13.5*s*xi;
    double f3 = 3.0*eta-1.0, g3 = 3.0*eta-2.0,  h3 = 13.5*s*eta;

    der(0,0) = der(1,0) = -0.5*g1*f1 - 1.5*s*(g1+f1);
    der(0,1) = 0.5*f2*g2 + 1.5*xi*(f2+g2);
    // der(1,1) = 0.0;
    // der(0,2) = 0.0;
    der(1,2) = 0.5*g3*f3 + 1.5*eta*(g3+f3);

    f1 *= 4.5, f2 *= 4.5, f3 *= 4.5;

    der(0,3) =  f1*(s-xi) - h2;
    der(1,3) = -f1*xi - h2;
    der(0,4) =  f2*(s-xi) + h2;
    der(1,4) = -f2*xi;
    der(0,5) =  f2*eta + h1;
    der(1,5) =  f2*xi;
    der(0,6) =  f3*eta;
    der(1,6) =  f3*xi + h1;
    der(0,7) = -f3*eta;
    der(1,7) =  f3*(s-eta) + h3;
    der(0,8) = -f1*eta - h3;
    der(1,8) =  f1*(s-eta) - h3;
    der(0,9) =  27.0*eta*(s-xi);
    der(1,9) =  27.0*xi*(s-eta);
    return der;
}

RSymMatrix Triangle10_ip::IntFF() const
{
    return intff;
}

double Triangle10_ip::IntFF (int i, int j) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
    return intff(i,j);
}

double Triangle10_ip::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10 && k >= 0 && k < 10);
    int p, np;
    double ff = 0.0;
    RVector fun(10);
    const double *wght;
    const Point  *absc;
    np = QRule_tri_9_19 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	ff += wght[p] * fun[i] * fun[j] * fun[k] * DetJ (absc[p]);
    }
    return ff;
}

RSymMatrix Triangle10_ip::IntPFF (const RVector &P) const
{
    RSymMatrix pff(10);
    RVector fun(10);
    double djac, v;
    int i, j, k, p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_9_19 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	djac = DetJ (absc[p]) * wght[p];
	for (k = 0; k < 10; k++) {
	    v = P[Node[k]] * fun[k] * djac;
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    pff(i,j) += v * fun[i] * fun[j];
	}
    }
    return pff;
}

double Triangle10_ip::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
    RVector fun(10);
    int k, p, np;
    double v, pff = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_9_19 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = fun[i] * fun[j] * DetJ (absc[p]) * wght[p];
	for (k = 0; k < 10; k++)
	    pff += P[Node[k]] * fun[k] * v;
    }
    return pff;
}

double Triangle10_ip::IntFDD (int i, int j, int k) const
{
    RVector fun(10);
    RDenseMatrix IJ(2,2), der(2,10);
    int p, np;
    double djac, fdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_7_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	der = LocalShapeD (absc[p]);
        djac = IJacobian (absc[p], IJ) * wght[p];
	fdd += fun[i] * djac *
	  ((IJ(0,0)*der(0,j) + IJ(0,1)*der(1,j)) *
	   (IJ(0,0)*der(0,k) + IJ(0,1)*der(1,k)) +
	   (IJ(1,0)*der(0,j) + IJ(1,1)*der(1,j)) *
	   (IJ(1,0)*der(0,k) + IJ(1,1)*der(1,k)));
    }
    return fdd;
}

double Triangle10_ip::IntPDD (int i, int j, const RVector &P) const
{
    RVector fun(10);
    RDenseMatrix IJ(2,2), der(2,10);
    int k, p, np;
    double djac, dd, fdd, pdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_7_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	der = LocalShapeD (absc[p]);
        djac = IJacobian (absc[p], IJ) * wght[p];
	dd = djac *
	  ((IJ(0,0)*der(0,i) + IJ(0,1)*der(1,i)) *
	   (IJ(0,0)*der(0,j) + IJ(0,1)*der(1,j)) +
	   (IJ(1,0)*der(0,i) + IJ(1,1)*der(1,i)) *
	   (IJ(1,0)*der(0,j) + IJ(1,1)*der(1,j)));
	for (k = 0, fdd = 0.0; k < 10; k++)
	    fdd += fun[k] * P[Node[k]];
	pdd += fdd * dd;
    }
    return pdd;
}

double Triangle10_ip::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double v, xi, eta, dFx, dFy, res = 0.0;
    Point loc(2);
    RVector fun(10);
    int p, np;
    const double *wght;
    const double *absc;
    np = QRule_lin_6 (&wght, &absc);

#ifndef TRI10IP_STORE_COORDS
    xERROR("Requires definition of TRI10IP_STORE_COORDS");
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y,
           n5x, n5y, n6x, n6y, n7x, n7y, n8x, n8y, n9x, n9y;
#endif

    if (bndside[0]) {
	loc[1] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t1, t3, t4, f = 13.5*xi*xi;
	    dFx = (t0 = (-5.5 + 18.0*xi -     f)) * n0x +
	          (t1 = ( 1.0 -  9.0*xi +     f)) * n1x +
	          (t3 = ( 9.0 - 45.0*xi + 3.0*f)) * n3x +
	          (t4 = (-4.5 + 36.0*xi - 3.0*f)) * n4x;
	    dFy = t0*n0y + t1*n1y + t3*n3y + t4*n4y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[0]]*fun[0] + P[Node[1]]*fun[1] +
		    P[Node[3]]*fun[3] + P[Node[4]]*fun[4])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[1]) {
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p]; loc[1] = eta = 1.0-loc[0];
	    fun = LocalShapeF (loc);
	    double f = 13.5*xi*xi;
	    dFx = ( 1.0 -  9.0*xi +     f) * n1x +
	          (-5.5 + 18.0*xi -     f) * n2x +
	          (-4.5 + 36.0*xi - 3.0*f) * n5x +
	          ( 9.0 - 45.0*xi + 3.0*f) * n6x;
	    f = 13.5*eta*eta;
	    dFy = (-5.5 + 18.0*eta -     f) * n1y +
	          ( 1.0 -  9.0*eta +     f) * n2y +
		  ( 9.0 - 45.0*eta + 3.0*f) * n5y +
		  (-4.5 + 36.0*eta - 3.0*f) * n6y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[2]]*fun[2] + P[Node[1]]*fun[1] + 
		    P[Node[6]]*fun[6] + P[Node[5]]*fun[5])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[2]) {
	loc[0] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t2, t7, t8, f = 13.5*eta*eta;
	    dFx = (t0 = (-5.5 + 18.0*eta -     f)) * n0x +
	          (t2 = ( 1.0 -  9.0*eta +     f)) * n2x +
	          (t7 = (-4.5 + 36.0*eta - 3.0*f)) * n7x +
	          (t8 = ( 9.0 - 45.0*eta + 3.0*f)) * n8x;
	    dFy = t0*n0y + t2*n2y + t7*n7y + t8*n8y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[0]]*fun[0] + P[Node[2]]*fun[2] +
		    P[Node[7]]*fun[7] + P[Node[8]]*fun[8])
	      * fun[i] * fun[j] * v;
	}
    }
    return res;
}

double Triangle10_ip::ComputeSize (const NodeList &nlist) const
{
    double sz = 0.0;
    int p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_1_1 (&wght, &absc);
    for (p = 0; p < np; p++) {
        sz += wght[p] * DetJ (absc[p], &nlist);
    }
    return sz;
}

RSymMatrix Triangle10_ip::ComputeIntFF (const NodeList &nlist) const
{
    int i, j, p, np;
    RSymMatrix ff(10);
    RVector fun(10);
    double v;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_6_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = wght[p] * DetJ (absc[p], &nlist);
        for (i = 0; i < 10; i++)
	    for (j = 0; j <= i; j++)
	        ff(i,j) += v * fun[i] * fun[j];
    }
    return ff;
}

RSymMatrix Triangle10_ip::ComputeIntDD (const NodeList &nlist) const
{
    RSymMatrix dd(10);
    RDenseMatrix IJ(2,2);
    double djac;
    int p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_4_6 (&wght, &absc);
    for (p = 0; p < np; p++) {
        djac = IJacobian (absc[p], IJ);
	dd += ATA (IJ*LocalShapeD (absc[p])) * (djac*wght[p]);
    }
    return dd;
}

RSymMatrix Triangle10_ip::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(10);
    if (!bndel) return bff;  // not a boundary element -> return zero matrix
    Point loc(2);
    RVector fun(10);
    int i, j, p, np;
    double v, xi, eta, dFx, dFy;
    const double *wght;
    const double *absc;
    np = QRule_lin_6 (&wght, &absc);

    if (bndside[0]) {
#ifndef TRI10IP_STORE_COORDS
        double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
	double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
	double n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1];
	double n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1];
#endif
	loc[1] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t1, t3, t4, f = 13.5*xi*xi;
	    dFx = (t0 = (-5.5 + 18.0*xi -     f)) * n0x +
	          (t1 = ( 1.0 -  9.0*xi +     f)) * n1x +
	          (t3 = ( 9.0 - 45.0*xi + 3.0*f)) * n3x +
	          (t4 = (-4.5 + 36.0*xi - 3.0*f)) * n4x;
	    dFy = t0*n0y + t1*n1y + t3*n3y + t4*n4y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[1]) {
#ifndef TRI10IP_STORE_COORDS
        double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
	double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
	double n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1];
	double n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1];
#endif
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p]; loc[1] = eta = 1.0-loc[0];
	    fun = LocalShapeF (loc);
	    double f = 13.5*xi*xi;
	    dFx = ( 1.0 -  9.0*xi +     f) * n1x +
	          (-5.5 + 18.0*xi -     f) * n2x +
	          (-4.5 + 36.0*xi - 3.0*f) * n5x +
	          ( 9.0 - 45.0*xi + 3.0*f) * n6x;
	    f = 13.5*eta*eta;
	    dFy = (-5.5 + 18.0*eta -     f) * n1y +
	          ( 1.0 -  9.0*eta +     f) * n2y +
		  ( 9.0 - 45.0*eta + 3.0*f) * n5y +
		  (-4.5 + 36.0*eta - 3.0*f) * n6y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[2]) {
#ifndef TRI10IP_STORE_COORDS
        double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
	double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
	double n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1];
	double n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1];
#endif
	loc[0] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t2, t7, t8, f = 13.5*eta*eta;
	    dFx = (t0 = (-5.5 + 18.0*eta -     f)) * n0x +
	          (t2 = ( 1.0 -  9.0*eta +     f)) * n2x +
	          (t7 = (-4.5 + 36.0*eta - 3.0*f)) * n7x +
	          (t8 = ( 9.0 - 45.0*eta + 3.0*f)) * n8x;
	    dFy = t0*n0y + t2*n2y + t7*n7y + t8*n8y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    return bff;
}

double Triangle10_ip::Jacobian (const Point &loc, RDenseMatrix &J) const
{
#ifndef TRI10IP_STORE_COORDS
    xERROR("Requires definition of TRI10IP_STORE_COORDS");
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y,
           n5x, n5y, n6x, n6y, n7x, n7y, n8x, n8y, n9x, n9y;
#endif
    dASSERT (loc.Dim() == 2, "Parameter 1 wrong dimension");
    dASSERT(J.nRows() == 2 && J.nCols() == 2, "Parameter 3 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);

    for (int i = 0; i < 2; i++) {
        J(i,0) = der(i,0)*n0x + der(i,1)*n1x + der(i,2)*n2x +
	         der(i,3)*n3x + der(i,4)*n4x + der(i,5)*n5x +
	         der(i,6)*n6x + der(i,7)*n7x + der(i,8)*n8x +
                 der(i,9)*n9x;
	J(i,1) = der(i,0)*n0y + der(i,1)*n1y + der(i,2)*n2y +
                 der(i,3)*n3y + der(i,4)*n4y + der(i,5)*n5y +
                 der(i,6)*n6y + der(i,7)*n7y + der(i,8)*n8y +
                 der(i,9)*n9y;
    }
    return J(0,0)*J(1,1) - J(1,0)*J(0,1);
}

double Triangle10_ip::IJacobian (const Point &loc, RDenseMatrix &IJ) const
{
    double d = Jacobian (loc, IJ);
    double id = 1.0/d;
    double tmp = IJ(0,0);
    IJ(0,0) = IJ(1,1)*id;
    IJ(0,1) *= -id;
    IJ(1,0) *= -id;
    IJ(1,1) = tmp*id;
    return d;
}

double Triangle10_ip::DetJ (const Point &loc, const NodeList *nlist) const
{
#ifndef TRI10IP_STORE_COORDS
    dASSERT (nlist != 0, "Node list required");
    double n0x = (*nlist)[Node[0]][0], n0y = (*nlist)[Node[0]][1];
    double n1x = (*nlist)[Node[1]][0], n1y = (*nlist)[Node[1]][1];
    double n2x = (*nlist)[Node[2]][0], n2y = (*nlist)[Node[2]][1];
    double n3x = (*nlist)[Node[3]][0], n3y = (*nlist)[Node[3]][1];
    double n4x = (*nlist)[Node[4]][0], n4y = (*nlist)[Node[4]][1];
    double n5x = (*nlist)[Node[5]][0], n5y = (*nlist)[Node[5]][1];
    double n6x = (*nlist)[Node[6]][0], n6y = (*nlist)[Node[6]][1];
    double n7x = (*nlist)[Node[7]][0], n7y = (*nlist)[Node[7]][1];
    double n8x = (*nlist)[Node[8]][0], n8y = (*nlist)[Node[8]][1];
    double n9x = (*nlist)[Node[9]][0], n9y = (*nlist)[Node[9]][1];
#endif
    dASSERT (loc.Dim() == 2, "Parameter 1 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);
    double j00, j01, j10, j11;

    j00 = der(0,0)*n0x + der(0,1)*n1x + der(0,2)*n2x +
          der(0,3)*n3x + der(0,4)*n4x + der(0,5)*n5x +
          der(0,6)*n6x + der(0,7)*n7x + der(0,8)*n8x +
          der(0,9)*n9x;
    j01 = der(0,0)*n0y + der(0,1)*n1y + der(0,2)*n2y +
          der(0,3)*n3y + der(0,4)*n4y + der(0,5)*n5y +
          der(0,6)*n6y + der(0,7)*n7y + der(0,8)*n8y +
          der(0,9)*n9y;
    j10 = der(1,0)*n0x + der(1,1)*n1x + der(1,2)*n2x +
          der(1,3)*n3x + der(1,4)*n4x + der(1,5)*n5x +
          der(1,6)*n6x + der(1,7)*n7x + der(1,8)*n8x +
          der(1,9)*n9x;
    j11 = der(1,0)*n0y + der(1,1)*n1y + der(1,2)*n2y +
          der(1,3)*n3y + der(1,4)*n4y + der(1,5)*n5y +
          der(1,6)*n6y + der(1,7)*n7y + der(1,8)*n8y +
          der(1,9)*n9y;

    return j00*j11 - j10*j01;
}

