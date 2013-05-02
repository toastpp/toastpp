// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe      File: tri6_ip.cc
//
// Definition of class Triangle6_ip
// 2nd order isoparametric triangle, using numerical integration
// ==========================================================================

// A note on numerical integration:
// Each integration of shape functions or shape function derivatives
// requires the appropriate degree integration rule
// Shape functions are 2nd order, therefore IntFF needs 4th order rule,
// IntFFF needs 6th order rule, etc.
// Shape function derivatives are 1st order, therefore IntDD needs 2nd order
// rule, IntFDD needs 4th order rule, etc.

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "lin_qr.h"
#include "tri_qr.h"

using namespace std;

Triangle6_ip::Triangle6_ip (const Triangle6_ip &el)
: Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Triangle6_ip::Copy ()
{
    return new Triangle6_ip(*this);
}

void Triangle6_ip::Initialise (const NodeList &nlist)
{
#ifdef TRI6IP_STORE_COORDS
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
    n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1];
    n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1];
    n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1];
#endif

    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1];

    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    Element_Unstructured_2D::Initialise (nlist);

    intff.Zero(6);
    intff = ComputeIntFF (nlist);
}

int Triangle6_ip::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Side index out of range");
    dASSERT(node >= 0 && node < 3, "Node index out of range");
    static int SN[3][3] = {{0,1,3},{1,2,4},{2,0,5}};
    return SN[side][node];
}

Point Triangle6_ip::Local (const NodeList &nlist, const Point& glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    static const double EPS = 1e-6;
    static const int MAXIT = 10;

    Point loc(2), lglob(2);
    RVector dloc(2), fun(6);
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
	for (int j = 0; j < 6; j++)
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

Point Triangle6_ip::NodeLocal (int node) const
{
    Point nloc(2);
    switch (node) {
        case 0:  nloc[0] = nloc[1] = 0.0; break;
        case 1:  nloc[0] = 1.0; nloc[1] = 0.0; break;
        case 2:  nloc[0] = 0.0; nloc[1] = 1.0; break;
        case 3:  nloc[0] = 0.5; nloc[1] = 0.0; break;
        case 4:  nloc[0] = nloc[1] = 0.5; break;
        case 5:  nloc[0] = 0.0; nloc[1] = 0.5; break;
        default: xERROR("Node index out of range"); break;
    }
    return nloc;
}

RVector Triangle6_ip::DirectionCosine (int side, RDenseMatrix &/*jacin*/)
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

const RVector &Triangle6_ip::LNormal (int side) const
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

bool Triangle6_ip::LContains (const Point& loc, bool pad) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 &&
	    loc[0]+loc[1]-EPS <= 1.0);
    } else {
	return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[0]+loc[1] <= 1.0);
    }
}

bool Triangle6_ip::GContains (const Point& glob, const NodeList& nlist) const
{
    return LContains (Local (nlist, glob));
}

RVector Triangle6_ip::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RVector fun(6);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = 4.0*L0*L1;
    fun[4] = 4.0*L1*L2;
    fun[5] = 4.0*L2*L0;
    return fun;
}

RDenseMatrix Triangle6_ip::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der(2,6);
    double lx = loc[0], ly = loc[1];
    der(0,0) = der(1,0) = 4.0*(lx+ly) - 3.0;
    der(0,1) = 4.0*lx - 1.0;
    //der(1,1) = 0.0;
    //der(0,2) = 0.0;
    der(1,2) = 4.0*ly - 1.0;
    der(0,3) = 4.0*(1.0-2.0*lx-ly);
    der(1,3) = -4.0*lx;
    der(0,4) = 4.0*ly;
    der(1,4) = 4.0*lx;
    der(0,5) = -4.0*ly;
    der(1,5) = 4.0*(1.0-lx-2.0*ly);
    return der;
}

RSymMatrix Triangle6_ip::IntFF() const
{
    return intff;
}

double Triangle6_ip::IntFF (int i, int j) const
{
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6);
    return intff(i,j);
}

double Triangle6_ip::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6 && k >= 0 && k < 6);
    int p, np;
    double ff = 0.0;
    RVector fun(6);
    const double *wght;
    const Point  *absc;
    np = QRule_tri_6_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	ff += wght[p] * fun[i] * fun[j] * fun[k] * DetJ (absc[p]);
    }
    return ff;
}

RSymMatrix Triangle6_ip::IntPFF (const RVector &P) const
{
    RSymMatrix pff(6);
    RVector fun(6);
    double djac, v;
    int i, j, k, p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_6_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	djac = DetJ (absc[p]) * wght[p];
	for (k = 0; k < 6; k++) {
	    v = P[Node[k]] * fun[k] * djac;
	    for (i = 0; i < 6; i++)
	        for (j = 0; j <= i; j++)
		    pff(i,j) += v * fun[i] * fun[j];
	}
    }
    return pff;
}

double Triangle6_ip::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6);
    RVector fun(6);
    int k, p, np;
    double v, pff = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_6_12 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = fun[i] * fun[j] * DetJ (absc[p]) * wght[p];
	for (k = 0; k < 6; k++)
	    pff += P[Node[k]] * fun[k] * v;
    }
    return pff;
}

double Triangle6_ip::IntFDD (int i, int j, int k) const
{
    RVector fun(6);
    RDenseMatrix IJ(2,2), der(2,6);
    int p, np;
    double djac, fdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_4_6 (&wght, &absc);
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

double Triangle6_ip::IntPDD (int i, int j, const RVector &P) const
{
    RVector fun(6);
    RDenseMatrix IJ(2,2), der(2,6);
    int k, p, np;
    double djac, dd, fdd, pdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_4_6 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	der = LocalShapeD (absc[p]);
        djac = IJacobian (absc[p], IJ) * wght[p];
	dd = djac *
	  ((IJ(0,0)*der(0,i) + IJ(0,1)*der(1,i)) *
	   (IJ(0,0)*der(0,j) + IJ(0,1)*der(1,j)) +
	   (IJ(1,0)*der(0,i) + IJ(1,1)*der(1,i)) *
	   (IJ(1,0)*der(0,j) + IJ(1,1)*der(1,j)));
	for (k = 0, fdd = 0.0; k < 6; k++)
	    fdd += fun[k] * P[Node[k]];
	pdd += fdd * dd;
    }
    return pdd;
}

double Triangle6_ip::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double v, xi, eta, dFx, dFy, res = 0.0;
    Point loc(2);
    RVector fun(6);
    int p, np;
    const double *wght;
    const double *absc;
    np = QRule_lin_6 (&wght, &absc);

#ifndef TRI6IP_STORE_COORDS
    xERROR("Requires definition of TRI6IP_STORE_COORDS");
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y, n5x, n5y;
#endif

    if (bndside[0]) {
	loc[1] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t1, t3;
	    dFx = (t0 = ( 4.0*xi-3.0)) * n0x +
	          (t1 = ( 4.0*xi-1.0)) * n1x +
	          (t3 = (-8.0*xi+4.0)) * n3x;
	    dFy = t0*n0y + t1*n1y + t3*n3y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[0]]*fun[0] + P[Node[1]]*fun[1] + P[Node[3]]*fun[3])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[1]) {
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p]; loc[1] = eta = 1.0-loc[0];
	    fun = LocalShapeF (loc);
	    dFx = (-1.0+4.0* xi) * n1x +
	          (-3.0+4.0* xi) * n2x +
	          ( 4.0-8.0* xi) * n4x;
	    dFy = (-3.0+4.0*eta) * n1y +
	          (-1.0+4.0*eta) * n2y +
	          ( 4.0-8.0*eta) * n4y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[2]]*fun[2] + P[Node[1]]*fun[1] + P[Node[4]]*fun[4])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[2]) {
	loc[0] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t2, t5;
	    dFx = (t0 = (-3.0+4.0*eta)) * n0x +
	          (t2 = (-1.0+4.0*eta)) * n2x +
	          (t5 = ( 4.0-8.0*eta)) * n5x;
	    dFy = t0*n0y + t2*n2y + t5*n5y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    res += (P[Node[0]]*fun[0] + P[Node[2]]*fun[2] + P[Node[5]]*fun[5])
	      * fun[i] * fun[j] * v;
	}
    }
    return res;
}

double Triangle6_ip::ComputeSize (const NodeList &nlist) const
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

RSymMatrix Triangle6_ip::ComputeIntFF (const NodeList &nlist) const
{
    int i, j, p, np;
    RSymMatrix ff(6);
    RVector fun(6);
    double v;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_4_6 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = wght[p] * DetJ (absc[p], &nlist);
        for (i = 0; i < 6; i++)
	    for (j = 0; j <= i; j++)
	        ff(i,j) += v * fun[i] * fun[j];
    }
    return ff;
}

RSymMatrix Triangle6_ip::ComputeIntDD (const NodeList &nlist) const
{
    RSymMatrix dd(6);
    RDenseMatrix IJ(2,2);
    double djac;
    int p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_2_3 (&wght, &absc);
    for (p = 0; p < np; p++) {
        djac = IJacobian (absc[p], IJ);
	dd += ATA (IJ*LocalShapeD (absc[p])) * (djac*wght[p]);
    }
    return dd;
}

RSymMatrix Triangle6_ip::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(6);
    if (!bndel) return bff;  // not a boundary element -> return zero matrix
    Point loc(2);
    RVector fun(6);
    int i, j, p, np;
    double v, xi, eta, dFx, dFy;
    const double *wght;
    const double *absc;
    np = QRule_lin_4 (&wght, &absc);

    if (bndside[0]) {
#ifndef TRI6IP_STORE_COORDS
        double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
	double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
	double n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1];
#endif
	loc[1] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t1, t3;
	    dFx = (t0 = ( 4.0*xi-3.0)) * n0x +
	          (t1 = ( 4.0*xi-1.0)) * n1x +
	          (t3 = (-8.0*xi+4.0)) * n3x;
	    dFy = t0*n0y + t1*n1y + t3*n3y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 6; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[1]) {
#ifndef TRI6IP_STORE_COORDS
        double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
	double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
	double n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1];
#endif
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p]; loc[1] = eta = 1.0-loc[0];
	    fun = LocalShapeF (loc);
	    dFx = (-1.0+4.0* xi) * n1x +
	          (-3.0+4.0* xi) * n2x +
	          ( 4.0-8.0* xi) * n4x;
	    dFy = (-3.0+4.0*eta) * n1y +
	          (-1.0+4.0*eta) * n2y +
	          ( 4.0-8.0*eta) * n4y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 6; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[2]) {
#ifndef TRI6IP_STORE_COORDS
        double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
	double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
	double n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1];
#endif
	loc[0] = 0.0;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p];
	    fun = LocalShapeF (loc);
	    double t0, t2, t5;
	    dFx = (t0 = (-3.0+4.0*eta)) * n0x +
	          (t2 = (-1.0+4.0*eta)) * n2x +
	          (t5 = ( 4.0-8.0*eta)) * n5x;
	    dFy = t0*n0y + t2*n2y + t5*n5y;
	    v = sqrt(dFx*dFx+dFy*dFy) * wght[p];
	    for (i = 0; i < 6; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    return bff;
}

double Triangle6_ip::Jacobian (const Point &loc, RDenseMatrix &J) const
{
#ifndef TRI6IP_STORE_COORDS
    xERROR("Requires definition of TRI6IP_STORE_COORDS");
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y, n5x, n5y;
#endif
    dASSERT (loc.Dim() == 2, "Parameter 1 wrong dimension");
    dASSERT(J.nRows() == 2 && J.nCols() == 2, "Parameter 2 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);

    for (int i = 0; i < 2; i++) {
        J(i,0) = der(i,0)*n0x + der(i,1)*n1x + der(i,2)*n2x +
	         der(i,3)*n3x + der(i,4)*n4x + der(i,5)*n5x;
        J(i,1) = der(i,0)*n0y + der(i,1)*n1y + der(i,2)*n2y +
	         der(i,3)*n3y + der(i,4)*n4y + der(i,5)*n5y;
    }
    return J(0,0)*J(1,1) - J(1,0)*J(0,1);
}

double Triangle6_ip::IJacobian (const Point &loc, RDenseMatrix &IJ) const
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

double Triangle6_ip::DetJ (const Point &loc, const NodeList *nlist) const
{
#ifndef TRI6IP_STORE_COORDS
    dASSERT (nlist != 0, "Node list required");
    double n0x = (*nlist)[Node[0]][0], n0y = (*nlist)[Node[0]][1];
    double n1x = (*nlist)[Node[1]][0], n1y = (*nlist)[Node[1]][1];
    double n2x = (*nlist)[Node[2]][0], n2y = (*nlist)[Node[2]][1];
    double n3x = (*nlist)[Node[3]][0], n3y = (*nlist)[Node[3]][1];
    double n4x = (*nlist)[Node[4]][0], n4y = (*nlist)[Node[4]][1];
    double n5x = (*nlist)[Node[5]][0], n5y = (*nlist)[Node[5]][1];
#endif
    dASSERT (loc.Dim() == 2, "Parameter 1 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);
    double j00, j01, j10, j11;

    j00 = der(0,0)*n0x + der(0,1)*n1x + der(0,2)*n2x +
          der(0,3)*n3x + der(0,4)*n4x + der(0,5)*n5x;
    j01 = der(0,0)*n0y + der(0,1)*n1y + der(0,2)*n2y +
          der(0,3)*n3y + der(0,4)*n4y + der(0,5)*n5y;
    j10 = der(1,0)*n0x + der(1,1)*n1x + der(1,2)*n2x +
          der(1,3)*n3x + der(1,4)*n4x + der(1,5)*n5x;
    j11 = der(1,0)*n0y + der(1,1)*n1y + der(1,2)*n2y +
          der(1,3)*n3y + der(1,4)*n4y + der(1,5)*n5y;

    return j00*j11 - j10*j01;
}
