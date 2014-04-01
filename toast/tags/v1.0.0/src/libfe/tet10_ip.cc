// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe      File: tet10_ip.cc
//
// Definition of class Tetrahedron10_ip
// 2nd order isoparametric tetrahedron, using numerical integration
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "tet_qr.h"
#include "tri_qr.h"

using namespace std;

Tetrahedron10_ip::Tetrahedron10_ip (const Tetrahedron10_ip &el)
: Element_Unstructured_3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Tetrahedron10_ip::Copy ()
{
    return new Tetrahedron10_ip(*this);
}

void Tetrahedron10_ip::Initialise (const NodeList &nlist)
{
    extern double TriangleArea (const Point &p1,
				const Point &p2, const Point &p3);

#ifdef TET10IP_STORE_COORDS
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1], n0z = nlist[Node[0]][2];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1], n1z = nlist[Node[1]][2];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1], n2z = nlist[Node[2]][2];
    n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1], n3z = nlist[Node[3]][2];
    n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1], n4z = nlist[Node[4]][2];
    n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1], n5z = nlist[Node[5]][2];
    n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1], n6z = nlist[Node[6]][2];
    n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1], n7z = nlist[Node[7]][2];
    n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1], n8z = nlist[Node[8]][2];
    n9x = nlist[Node[9]][0], n9y = nlist[Node[9]][1], n9z = nlist[Node[9]][2];
#endif

    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1],
           z0 = nlist[Node[0]][2];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1],
           z1 = nlist[Node[1]][2];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1],
           z2 = nlist[Node[2]][2];
    double x3 = nlist[Node[3]][0], y3 = nlist[Node[3]][1],
           z3 = nlist[Node[3]][2];

    a0 = x1*(y2*z3-y3*z2) - x2*(y1*z3-y3*z1) + x3*(y1*z2-y2*z1);
    a1 = x0*(y3*z2-y2*z3) - x2*(y3*z0-y0*z3) + x3*(y2*z0-y0*z2);
    a2 = x0*(y1*z3-y3*z1) - x1*(y0*z3-y3*z0) + x3*(y0*z1-y1*z0);
    a3 = x0*(y2*z1-y1*z2) - x1*(y2*z0-y0*z2) + x2*(y1*z0-y0*z1);
    b0 = y1*(z3-z2) + y2*(z1-z3) + y3*(z2-z1);
    b1 = y0*(z2-z3) + y2*(z3-z0) + y3*(z0-z2);
    b2 = y0*(z3-z1) + y1*(z0-z3) + y3*(z1-z0);
    b3 = y0*(z1-z2) + y1*(z2-z0) + y2*(z0-z1);
    c0 = x1*(z2-z3) + x2*(z3-z1) + x3*(z1-z2);
    c1 = x0*(z3-z2) + x2*(z0-z3) + x3*(z2-z0);
    c2 = x0*(z1-z3) + x1*(z3-z0) + x3*(z0-z1);
    c3 = x0*(z2-z1) + x1*(z0-z2) + x2*(z1-z0);
    d0 = x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1);
    d1 = x0*(y2-y3) + x2*(y3-y0) + x3*(y0-y2);
    d2 = x0*(y3-y1) + x1*(y0-y3) + x3*(y1-y0);
    d3 = x0*(y1-y2) + x1*(y2-y0) + x2*(y0-y1);

    // calculate surface triangle areas
    for (int sd = 0; sd < 4; sd++)
        side_size[sd] = TriangleArea (nlist[Node[SideNode(sd,0)]],
			 	      nlist[Node[SideNode(sd,1)]],
				      nlist[Node[SideNode(sd,2)]]);

    Element_Unstructured_3D::Initialise (nlist);

    intff.Zero (10);
    intff = ComputeIntFF (nlist);
}

int Tetrahedron10_ip::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 4, "Side index out of range");
    dASSERT(node >= 0 && node < 6, "Node index out of range");
    static int SN[4][6] = {{0,1,2,4,7,5},{0,3,1,6,8,4},
			   {0,2,3,5,9,6},{1,3,2,8,9,7}};
    return SN[side][node];
}

Point Tetrahedron10_ip::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");
    static const double EPS = 1e-4;
    static const int MAXIT = 10;

    Point loc(3), lglob(3);
    RVector dloc(3), fun(10);
    RDenseMatrix IJ(3,3);

    // First do linear inverse mapping
    double scale = 1.0/(a0+a1+a2+a3);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]) * scale;
    loc[2] = (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]) * scale;

    for (int i = 0; i < MAXIT; i++) {
        // map forward again
        fun = LocalShapeF (loc);
	IJacobian (loc, IJ);
	lglob[0] = lglob[1] = lglob[2] = 0.0;
	for (int j = 0; j < 10; j++)
	    lglob += nlist[Node[j]] * fun[j];

	// apply correction
	dloc = transpose(IJ) * (glob-lglob);
	loc += dloc;
	if (length(dloc) < EPS) return loc;
    }
    loc[0] = loc[1] = loc[2] = -1.0;
    // if no convergence, assume point is outside element - HACK!
    return loc;
}

Point Tetrahedron10_ip::NodeLocal (int node) const
{
    Point nloc(3);
    switch (node) {
    case 0: nloc[0] = nloc[1] = nloc[2] = 0.0; break;
    case 1: nloc[0] = 1.0; nloc[1] = nloc[2] = 0.0; break;
    case 2: nloc[0] = nloc[2] = 0.0; nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 0.0; nloc[2] = 1.0; break;
    case 4: nloc[0] = 0.5; nloc[1] = nloc[2] = 0.0; break;
    case 5: nloc[0] = nloc[2] = 0.0; nloc[1] = 0.5; break;
    case 6: nloc[0] = nloc[1] = 0.0; nloc[2] = 0.5; break;
    case 7: nloc[0] = nloc[1] = 0.5; nloc[2] = 0.0; break;
    case 8: nloc[0] = nloc[2] = 0.5; nloc[1] = 0.0; break;
    case 9: nloc[0] = 0.0; nloc[1] = nloc[2] = 0.5; break;
    default: xERROR("Node index out of range");
    }
    return nloc;
}

RVector Tetrahedron10_ip::DirectionCosine (int side, RDenseMatrix &jacin)
{
  // INVALID
    RVector cosin(3);
    switch (side) {
    case 0: cosin[0] = -b3, cosin[1] = -c3, cosin[2] = -d3; break;
    case 1: cosin[0] = -b2, cosin[1] = -c2, cosin[2] = -d2; break;
    case 2: cosin[0] = -b1, cosin[1] = -c1, cosin[2] = -d1; break;
    case 3: cosin[0] = -b0, cosin[1] = -c0, cosin[2] = -d0; break;
    default: xERROR("Side index out of range");
    }
    return cosin/length(cosin);
}

const RVector &Tetrahedron10_ip::LNormal (int side) const
{
    static const RVector lnm0 = RVector (3, "0 0 -1");
    static const RVector lnm1 = RVector (3, "0 -1 0");
    static const RVector lnm2 = RVector (3, "-1 0 0");
    static const RVector lnm3 = RVector (3,
         "0.577350269 0.577350269 0.577350269");
    static const RVector *lnm[4] = {
      &lnm0, &lnm1, &lnm2, &lnm3
    };
    dASSERT(side >= 0 && side < 4, "Argument 1 index out of range");
    return *lnm[side];
}

bool Tetrahedron10_ip::LContains (const Point &loc, bool pad) const
{
    dASSERT(loc.Dim() == 3, "Local point must be 3D");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 && loc[2]+EPS >= 0.0 &&
	    loc[0]+loc[1]+loc[2]-EPS <= 1.0);
    } else {
	return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[2] >= 0.0 &&
	    loc[0]+loc[1]+loc[2] <= 1.0);
    }
}

RVector Tetrahedron10_ip::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    RVector fun(10);
    double L0 = 1.0-loc[0]-loc[1]-loc[2];
    double L1 = loc[0];
    double L2 = loc[1];
    double L3 = loc[2];
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = L3 * (2.0*L3 - 1.0);
    fun[4] = 4.0*L0*L1;
    fun[5] = 4.0*L0*L2;
    fun[6] = 4.0*L0*L3;
    fun[7] = 4.0*L1*L2;
    fun[8] = 4.0*L1*L3;
    fun[9] = 4.0*L2*L3;
    return fun;
}

RDenseMatrix Tetrahedron10_ip::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    RDenseMatrix der(3,10);
    double lx = loc[0], ly = loc[1], lz = loc[2];

    // elements not set are zero
    der(0,0) = der(1,0) = der(2,0) = 4.0*(lx+ly+lz)-3.0;
    der(0,1) = 4.0*lx - 1.0;
    der(1,2) = 4.0*ly - 1.0;
    der(2,3) = 4.0*lz - 1.0;
    der(0,4) = 4.0*(1.0-2.0*lx-ly-lz);
    der(1,4) = der(2,4) = -4.0*lx;
    der(0,5) = der(2,5) = -4.0*ly;
    der(1,5) = 4.0*(1.0-lx-2.0*ly-lz);
    der(0,6) = der(1,6) = -4.0*lz;
    der(2,6) = 4.0*(1.0-lx-ly-2.0*lz);
    der(0,7) = der(2,9) = 4.0*ly;
    der(1,7) = der(2,8) = 4.0*lx;
    der(0,8) = der(1,9) = 4.0*lz;
    return der;
}

RSymMatrix Tetrahedron10_ip::IntFF () const {
    return intff;
}

double Tetrahedron10_ip::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
    return intff(i,j);
}

double Tetrahedron10_ip::IntFFF (int i, int j, int k) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10 && k >= 0 && k < 10);
    int p, np;
    double ff = 0.0;
    RVector fun(10);
    const double *wght;
    const Point  *absc;
    np = QRule_tet_6_24 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	ff += wght[p] * fun[i] * fun[j] * fun[k] * DetJ (absc[p]);
    }
    return ff;
}

RSymMatrix Tetrahedron10_ip::IntPFF (const RVector &P) const
{
    RSymMatrix pff(10);
    RVector fun(10);
    double djac, v;
    int i, j, k, p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_6_24 (&wght, &absc);
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

double Tetrahedron10_ip::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
    RVector fun(10);
    int k, p, np;
    double v, pff = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_6_24 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = fun[i] * fun[j] * DetJ (absc[p]) * wght[p];
	for (k = 0; k < 10; k++)
	    pff += P[Node[k]] * fun[k] * v;
    }
    return pff;
}

double Tetrahedron10_ip::IntFDD (int i, int j, int k) const 
{
    RVector fun(10);
    RDenseMatrix IJ(3,3), der(3,10);
    int p, np;
    double djac, fdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_4_14 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	der = LocalShapeD (absc[p]);
        djac = IJacobian (absc[p], IJ) * wght[p];
	fdd += fun[i] * djac *
	  ((IJ(0,0)*der(0,j) + IJ(0,1)*der(1,j) + IJ(0,2)*der(2,j)) *
	   (IJ(0,0)*der(0,k) + IJ(0,1)*der(1,k) + IJ(0,2)*der(2,k)) +
	   (IJ(1,0)*der(0,j) + IJ(1,1)*der(1,j) + IJ(1,2)*der(2,j)) *
	   (IJ(1,0)*der(0,k) + IJ(1,1)*der(1,k) + IJ(1,2)*der(2,k)) +
	   (IJ(2,0)*der(0,j) + IJ(2,1)*der(1,j) + IJ(2,2)*der(2,j)) *
	   (IJ(2,0)*der(0,k) + IJ(2,1)*der(1,k) + IJ(2,2)*der(2,k)));
    }
    return fdd;
}

double Tetrahedron10_ip::IntPDD (int i, int j, const RVector &P) const
{
    RVector fun(10);
    RDenseMatrix IJ(3,3), der(3,10);
    int k, p, np;
    double djac, dd, fdd, pdd = 0.0;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_4_14 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	der = LocalShapeD (absc[p]);
        djac = IJacobian (absc[p], IJ) * wght[p];
	dd = djac *
	  ((IJ(0,0)*der(0,i) + IJ(0,1)*der(1,i) + IJ(0,2)*der(2,i)) *
	   (IJ(0,0)*der(0,j) + IJ(0,1)*der(1,j) + IJ(0,2)*der(2,j)) +
	   (IJ(1,0)*der(0,i) + IJ(1,1)*der(1,i) + IJ(1,2)*der(2,i)) *
	   (IJ(1,0)*der(0,j) + IJ(1,1)*der(1,j) + IJ(1,2)*der(2,j)) +
	   (IJ(2,0)*der(0,i) + IJ(2,1)*der(1,i) + IJ(2,2)*der(2,i)) *
	   (IJ(2,0)*der(0,j) + IJ(2,1)*der(1,j) + IJ(2,2)*der(2,j)));
	for (k = 0, fdd = 0.0; k < 10; k++)
	    fdd += fun[k] * P[Node[k]];
	pdd += fdd * dd;
    }
    return pdd;
}

double Tetrahedron10_ip::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    int p, np;
    Point loc(3);
    RVector fun(10);
    double v, vx, vy, vz, xi, eta, zeta;
    double dFdxi_x, dFdxi_y, dFdxi_z, dFdeta_x, dFdeta_y, dFdeta_z,
      dFdzeta_x, dFdzeta_y, dFdzeta_z;
    const double *wght;
    const Point  *absc;

    np = QRule_tri_6_12 (&wght, &absc);

    if (bndside[0]) { // z = 0
        loc[2] = 0.0;
	double t0, t1, t2, t4, t5, t7;
	for (p = 0; p < np; p++) {
	    loc[0] = xi  = absc[p][0];
	    loc[1] = eta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t0 = (-3.0+4.0*(xi+eta)))   * n0x +
		      (t1 = (-1.0+4.0*xi))         * n1x +
	              (t4 = ( 4.0-8.0*xi-4.0*eta)) * n4x +
	              (t5 = (-4.0*eta))            * n5x +
	              (t7 = ( 4.0*eta))            * n7x;
	    dFdxi_y = t0*n0y + t1*n1y + t4*n4y + t5*n5y + t7*n7y;
	    dFdxi_z = t0*n0z + t1*n1z + t4*n4z + t5*n5z + t7*n7z;

	    dFdeta_x =  t0                          * n0x +
		       (t2 = (-1.0+4.0*eta))        * n2x +
	               (t4 = (-4.0*xi))             * n4x +
	               (t5 = ( 4.0-4.0*xi-8.0*eta)) * n5x +
	               (t7 = ( 4.0*xi))             * n7x;
	    dFdeta_y = t0*n0y + t2*n2y + t4*n4y + t5*n5y + t7*n7y;
	    dFdeta_z = t0*n0z + t2*n2z + t4*n4z + t5*n5z + t7*n7z;

	    vx = dFdxi_y * dFdeta_z - dFdxi_z * dFdeta_y;
	    vy = dFdxi_z * dFdeta_x - dFdxi_x * dFdeta_z;
	    vz = dFdxi_x * dFdeta_y - dFdxi_y * dFdeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];

	    res += (P[Node[0]]*fun[0] + P[Node[1]]*fun[1] + P[Node[2]]*fun[2] +
		    P[Node[4]]*fun[4] + P[Node[5]]*fun[5] + P[Node[7]]*fun[7])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[1]) { // y = 0
        loc[1] = 0.0;
	double t0, t1, t3, t4, t6, t8;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p][0];
	    loc[2] = zeta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t0 = (-3.0+4.0*(xi+zeta)))   * n0x +
	              (t1 = (-1.0+4.0*xi))          * n1x +
	              (t4 = ( 4.0-8.0*xi-4.0*zeta)) * n4x +
	              (t6 = (-4.0*zeta))            * n6x +
	              (t8 = ( 4.0*zeta))            * n8x;
	    dFdxi_y = t0*n0y + t1*n1y + t4*n4y + t6*n6y + t8*n8y;
	    dFdxi_z = t0*n0z + t1*n1z + t4*n4z + t6*n6z + t8*n8z;

	    dFdzeta_x = t0                          * n0x +
	              (t3 = (-1.0+4.0*zeta))        * n3x +
	              (t4 = (-4.0*xi))              * n4x +
	              (t6 = ( 4.0-4.0*xi-8.0*zeta)) * n6x +
	              (t8 = ( 4.0*xi))              * n8x;
	    dFdzeta_y = t0*n0y + t3*n3y + t4*n4y + t6*n6y + t8*n8y;
	    dFdzeta_z = t0*n0z + t3*n3z + t4*n4z + t6*n6z + t8*n8z;

	    vx = dFdxi_y * dFdzeta_z - dFdxi_z * dFdzeta_y;
	    vy = dFdxi_z * dFdzeta_x - dFdxi_x * dFdzeta_z;
	    vz = dFdxi_x * dFdzeta_y - dFdxi_y * dFdzeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];
	    
	    res += (P[Node[0]]*fun[0] + P[Node[1]]*fun[1] + P[Node[3]]*fun[3] +
		    P[Node[4]]*fun[4] + P[Node[6]]*fun[6] + P[Node[8]]*fun[8])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[2]) { // x = 0
        loc[0] = 0.0;
	double t0, t2, t3, t5, t6, t9;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p][0];
	    loc[2] = zeta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdeta_x = (t0 = (-3.0+4.0*(eta+zeta)))   * n0x +
	               (t2 = (-1.0+4.0*eta))          * n2x +
	               (t5 = ( 4.0-8.0*eta-4.0*zeta)) * n5x +
	               (t6 = (-4.0*zeta))             * n6x +
	               (t9 = ( 4.0*zeta))             * n9x;
	    dFdeta_y = t0*n0y + t2*n2y + t5*n5y + t6*n6y + t9*n9y;
	    dFdeta_z = t0*n0z + t2*n2z + t5*n5z + t6*n6z + t9*n9z;

	    dFdzeta_x = t0                            * n0x +
	               (t3 = (-1.0+4.0*zeta))         * n3x +
	               (t5 = (-4.0*eta))              * n5x +
	               (t6 = ( 4.0-4.0*eta-8.0*zeta)) * n6x +
	               (t9 = ( 4.0*eta))              * n9x;
	    dFdzeta_y = t0*n0y + t3*n3y + t5*n5y + t6*n6y + t9*n9y;
	    dFdzeta_z = t0*n0z + t3*n3z + t5*n5z + t6*n6z + t9*n9z;

	    vx = dFdeta_y * dFdzeta_z - dFdeta_z * dFdzeta_y;
	    vy = dFdeta_z * dFdzeta_x - dFdeta_x * dFdzeta_z;
	    vz = dFdeta_x * dFdzeta_y - dFdeta_y * dFdzeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];
	    
	    res += (P[Node[0]]*fun[0] + P[Node[2]]*fun[2] + P[Node[3]]*fun[3] +
		    P[Node[5]]*fun[5] + P[Node[6]]*fun[6] + P[Node[9]]*fun[9])
	      * fun[i] * fun[j] * v;
	}
    }
    if (bndside[3]) { // x+y+z = 1
	double t1, t2, t3, t7, t8, t9;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p][0];
	    loc[1] = eta = absc[p][1];
	    loc[2] = zeta = 1-xi-eta;
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t1 = (-1.0+4.0*xi))         * n1x +
	              (t3 = (-3.0+4.0*(xi+eta)))   * n3x +
	              (t7 = ( 4.0*eta))            * n7x +
	              (t8 = ( 4.0-8.0*xi-4.0*eta)) * n8x +
	              (t9 = (-4.0*eta))            * n9x;
	    dFdxi_y = t1*n1y + t3*n3y + t7*n7y + t8*n8y + t9*n9y;
	    dFdxi_z = t1*n1z + t3*n3z + t7*n7z + t8*n8z + t9*n9z;

	    dFdeta_x = (t2 = (-1.0+4.0*eta))        * n2x +
	               (t3 = (-3.0+4.0*(xi+eta)))   * n3x +
	               (t7 = ( 4.0*xi))             * n7x +
	               (t8 = (-4.0*xi))             * n8x +
	               (t9 = ( 4.0-4.0*xi-8.0*eta)) * n9x;
	    dFdeta_y = t2*n2y + t3*n3y + t7*n7y + t8*n8y + t9*n9y;
	    dFdeta_z = t2*n2z + t3*n3z + t7*n7z + t8*n8z + t9*n9z;

	    vx = dFdxi_y * dFdeta_z - dFdxi_z * dFdeta_y;
	    vy = dFdxi_z * dFdeta_x - dFdxi_x * dFdeta_z;
	    vz = dFdxi_x * dFdeta_y - dFdxi_y * dFdeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];

	    res += (P[Node[1]]*fun[i] + P[Node[2]]*fun[2] + P[Node[3]]*fun[3] +
		    P[Node[7]]*fun[7] + P[Node[8]]*fun[8] + P[Node[9]]*fun[9])
	      * fun[i] * fun[j] * v;
	}
    }
    return res;
}

double Tetrahedron10_ip::ComputeSize (const NodeList &nlist) const
{
    double sz = 0.0;
    int p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_1_1 (&wght, &absc);
    for (p = 0; p < np; p++)
        sz += wght[p] * DetJ (absc[p], &nlist);
    return sz;
}

RSymMatrix Tetrahedron10_ip::ComputeIntFF (const NodeList &nlist) const
{
    int i, j, p, np;
    RSymMatrix ff(10);
    RVector fun(10);
    double v;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_4_14 (&wght, &absc);
    for (p = 0; p < np; p++) {
        fun = LocalShapeF (absc[p]);
	v = wght[p] * DetJ (absc[p], &nlist);
	for (i = 0; i < 10; i++)
	    for (j = 0; j <= i; j++)
	        ff(i,j) += v * fun[i] * fun[j];
    }
    return ff;
}

RSymMatrix Tetrahedron10_ip::ComputeIntDD (const NodeList &nlist) const
{
    RSymMatrix dd(10);
    RDenseMatrix IJ(3,3);
    double djac;
    int p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tet_2_4 (&wght, &absc);
    for (p = 0; p < np; p++) {
        djac = IJacobian (absc[p], IJ);
	dd += ATA (IJ*LocalShapeD (absc[p])) * (djac*wght[p]);
    }
    return dd;
}

RSymMatrix Tetrahedron10_ip::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(10);
    if (!bndel) return bff;

    Point loc(3);
    RVector fun(10);
    double v, vx, vy, vz, xi, eta, zeta;
    double dFdxi_x, dFdxi_y, dFdxi_z, dFdeta_x, dFdeta_y, dFdeta_z,
      dFdzeta_x, dFdzeta_y, dFdzeta_z;
    int i, j, p, np;
    const double *wght;
    const Point  *absc;
    np = QRule_tri_4_6 (&wght, &absc);

    if (bndside[0]) { // z = 0
#ifndef TET10IP_STORE_COORDS
    double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1],
           n0z = nlist[Node[0]][2];
    double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1],
           n1z = nlist[Node[1]][2];
    double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1],
           n2z = nlist[Node[2]][2];
    double n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1],
           n4z = nlist[Node[4]][2];
    double n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1],
           n5z = nlist[Node[5]][2];
    double n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1],
           n7z = nlist[Node[7]][2];
#endif
        loc[2] = 0.0;
	double t0, t1, t2, t4, t5, t7;
	for (p = 0; p < np; p++) {
	    loc[0] = xi  = absc[p][0];
	    loc[1] = eta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t0 = (-3.0+4.0*(xi+eta)))   * n0x +
		      (t1 = (-1.0+4.0*xi))         * n1x +
	              (t4 = ( 4.0-8.0*xi-4.0*eta)) * n4x +
	              (t5 = (-4.0*eta))            * n5x +
	              (t7 = ( 4.0*eta))            * n7x;
	    dFdxi_y = t0*n0y + t1*n1y + t4*n4y + t5*n5y + t7*n7y;
	    dFdxi_z = t0*n0z + t1*n1z + t4*n4z + t5*n5z + t7*n7z;

	    dFdeta_x =  t0                          * n0x +
		       (t2 = (-1.0+4.0*eta))        * n2x +
	               (t4 = (-4.0*xi))             * n4x +
	               (t5 = ( 4.0-4.0*xi-8.0*eta)) * n5x +
	               (t7 = ( 4.0*xi))             * n7x;
	    dFdeta_y = t0*n0y + t2*n2y + t4*n4y + t5*n5y + t7*n7y;
	    dFdeta_z = t0*n0z + t2*n2z + t4*n4z + t5*n5z + t7*n7z;

	    vx = dFdxi_y * dFdeta_z - dFdxi_z * dFdeta_y;
	    vy = dFdxi_z * dFdeta_x - dFdxi_x * dFdeta_z;
	    vz = dFdxi_x * dFdeta_y - dFdxi_y * dFdeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];

	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[1]) { // y = 0
#ifndef TET10IP_STORE_COORDS
    double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1],
           n0z = nlist[Node[0]][2];
    double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1],
           n1z = nlist[Node[1]][2];
    double n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1],
           n3z = nlist[Node[3]][2];
    double n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1],
           n4z = nlist[Node[4]][2];
    double n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1],
           n6z = nlist[Node[6]][2];
    double n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1],
           n8z = nlist[Node[8]][2];
#endif
        loc[1] = 0.0;
	double t0, t1, t3, t4, t6, t8;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p][0];
	    loc[2] = zeta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t0 = (-3.0+4.0*(xi+zeta)))   * n0x +
	              (t1 = (-1.0+4.0*xi))          * n1x +
	              (t4 = ( 4.0-8.0*xi-4.0*zeta)) * n4x +
	              (t6 = (-4.0*zeta))            * n6x +
	              (t8 = ( 4.0*zeta))            * n8x;
	    dFdxi_y = t0*n0y + t1*n1y + t4*n4y + t6*n6y + t8*n8y;
	    dFdxi_z = t0*n0z + t1*n1z + t4*n4z + t6*n6z + t8*n8z;

	    dFdzeta_x = t0                          * n0x +
	              (t3 = (-1.0+4.0*zeta))        * n3x +
	              (t4 = (-4.0*xi))              * n4x +
	              (t6 = ( 4.0-4.0*xi-8.0*zeta)) * n6x +
	              (t8 = ( 4.0*xi))              * n8x;
	    dFdzeta_y = t0*n0y + t3*n3y + t4*n4y + t6*n6y + t8*n8y;
	    dFdzeta_z = t0*n0z + t3*n3z + t4*n4z + t6*n6z + t8*n8z;

	    vx = dFdxi_y * dFdzeta_z - dFdxi_z * dFdzeta_y;
	    vy = dFdxi_z * dFdzeta_x - dFdxi_x * dFdzeta_z;
	    vz = dFdxi_x * dFdzeta_y - dFdxi_y * dFdzeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];
	    
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[2]) { // x = 0
#ifndef TET10IP_STORE_COORDS
    double n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1],
           n0z = nlist[Node[0]][2];
    double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1],
           n2z = nlist[Node[2]][2];
    double n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1],
           n3z = nlist[Node[3]][2];
    double n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1],
           n5z = nlist[Node[5]][2];
    double n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1],
           n6z = nlist[Node[6]][2];
    double n9x = nlist[Node[9]][0], n9y = nlist[Node[9]][1],
           n9z = nlist[Node[9]][2];
#endif
        loc[0] = 0.0;
	double t0, t2, t3, t5, t6, t9;
	for (p = 0; p < np; p++) {
	    loc[1] = eta = absc[p][0];
	    loc[2] = zeta = absc[p][1];
	    fun = LocalShapeF (loc);
	    dFdeta_x = (t0 = (-3.0+4.0*(eta+zeta)))   * n0x +
	               (t2 = (-1.0+4.0*eta))          * n2x +
	               (t5 = ( 4.0-8.0*eta-4.0*zeta)) * n5x +
	               (t6 = (-4.0*zeta))             * n6x +
	               (t9 = ( 4.0*zeta))             * n9x;
	    dFdeta_y = t0*n0y + t2*n2y + t5*n5y + t6*n6y + t9*n9y;
	    dFdeta_z = t0*n0z + t2*n2z + t5*n5z + t6*n6z + t9*n9z;

	    dFdzeta_x = t0                            * n0x +
	               (t3 = (-1.0+4.0*zeta))         * n3x +
	               (t5 = (-4.0*eta))              * n5x +
	               (t6 = ( 4.0-4.0*eta-8.0*zeta)) * n6x +
	               (t9 = ( 4.0*eta))              * n9x;
	    dFdzeta_y = t0*n0y + t3*n3y + t5*n5y + t6*n6y + t9*n9y;
	    dFdzeta_z = t0*n0z + t3*n3z + t5*n5z + t6*n6z + t9*n9z;

	    vx = dFdeta_y * dFdzeta_z - dFdeta_z * dFdzeta_y;
	    vy = dFdeta_z * dFdzeta_x - dFdeta_x * dFdzeta_z;
	    vz = dFdeta_x * dFdzeta_y - dFdeta_y * dFdzeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];
	    
	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }
    if (bndside[3]) { // x+y+z = 1
#ifndef TET10IP_STORE_COORDS
    double n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1],
           n1z = nlist[Node[1]][2];
    double n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1],
           n2z = nlist[Node[2]][2];
    double n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1],
           n3z = nlist[Node[3]][2];
    double n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1],
           n7z = nlist[Node[7]][2];
    double n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1],
           n8z = nlist[Node[8]][2];
    double n9x = nlist[Node[9]][0], n9y = nlist[Node[9]][1],
           n9z = nlist[Node[9]][2];
#endif
	double t1, t2, t3, t7, t8, t9;
	for (p = 0; p < np; p++) {
	    loc[0] = xi = absc[p][0];
	    loc[1] = eta = absc[p][1];
	    loc[2] = zeta = 1-xi-eta;
	    fun = LocalShapeF (loc);
	    dFdxi_x = (t1 = (-1.0+4.0*xi))         * n1x +
	              (t3 = (-3.0+4.0*(xi+eta)))   * n3x +
	              (t7 = ( 4.0*eta))            * n7x +
	              (t8 = ( 4.0-8.0*xi-4.0*eta)) * n8x +
	              (t9 = (-4.0*eta))            * n9x;
	    dFdxi_y = t1*n1y + t3*n3y + t7*n7y + t8*n8y + t9*n9y;
	    dFdxi_z = t1*n1z + t3*n3z + t7*n7z + t8*n8z + t9*n9z;

	    dFdeta_x = (t2 = (-1.0+4.0*eta))        * n2x +
	               (t3 = (-3.0+4.0*(xi+eta)))   * n3x +
	               (t7 = ( 4.0*xi))             * n7x +
	               (t8 = (-4.0*xi))             * n8x +
	               (t9 = ( 4.0-4.0*xi-8.0*eta)) * n9x;
	    dFdeta_y = t2*n2y + t3*n3y + t7*n7y + t8*n8y + t9*n9y;
	    dFdeta_z = t2*n2z + t3*n3z + t7*n7z + t8*n8z + t9*n9z;

	    vx = dFdxi_y * dFdeta_z - dFdxi_z * dFdeta_y;
	    vy = dFdxi_z * dFdeta_x - dFdxi_x * dFdeta_z;
	    vz = dFdxi_x * dFdeta_y - dFdxi_y * dFdeta_x;
	    v = sqrt (vx*vx + vy*vy + vz*vz) * wght[p];

	    for (i = 0; i < 10; i++)
	        for (j = 0; j <= i; j++)
		    bff(i,j) += v * fun[i] * fun[j];
	}
    }

    return bff;
}

double Tetrahedron10_ip::Jacobian (const Point &loc, RDenseMatrix &J) const
{
#ifndef TET10IP_STORE_COORDS
    xERROR("Requires definition of TRI6IP_STORE_COORDS");
    double n0x, n0y, n0z, n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z;
    double n4x, n4y, n4z, n5x, n5y, n5z, n6x, n6y, n6z, n7x, n7y, n7z;
    double n8x, n8y, n8z, n9x, n9y, n9z;
#endif
    dASSERT (loc.Dim() == 3, "Parameter 1 wrong dimension");
    dASSERT(J.nRows() == 3 && J.nCols() == 3, "Parameter 2 wrong dimension");
    double t1 = 4.0*loc[0];
    double t2 = 4.0*loc[1];
    double t3 = 4.0*loc[2];
    double t4 = t1-1.0;
    double t5 = t2-1.0;
    double t6 = t3-1.0;
    double t7 = t4+t5+t6;
    double t8 = 2.0*t4+t5+t6;
    double t9 = t4+2.0*t5+t6;
    double t10 = t4+t5+2.0*t6;

    J(0,0) = t7*n0x + t4*n1x - t8*n4x + t2*(n7x-n5x) + t3*(n8x-n6x);
    J(0,1) = t7*n0y + t4*n1y - t8*n4y + t2*(n7y-n5y) + t3*(n8y-n6y);
    J(0,2) = t7*n0z + t4*n1z - t8*n4z + t2*(n7z-n5z) + t3*(n8z-n6z);

    J(1,0) = t7*n0x + t5*n2x - t9*n5x + t1*(n7x-n4x) + t3*(n9x-n6x);
    J(1,1) = t7*n0y + t5*n2y - t9*n5y + t1*(n7y-n4y) + t3*(n9y-n6y);
    J(1,2) = t7*n0z + t5*n2z - t9*n5z + t1*(n7z-n4z) + t3*(n9z-n6z);

    J(2,0) = t7*n0x + t6*n3x - t10*n6x + t1*(n8x-n4x) + t2*(n9x-n5x);
    J(2,1) = t7*n0y + t6*n3y - t10*n6y + t1*(n8y-n4y) + t2*(n9y-n5y);
    J(2,2) = t7*n0z + t6*n3z - t10*n6z + t1*(n8z-n4z) + t2*(n9z-n5z);

    return J(0,0) * (J(1,1)*J(2,2) - J(2,1)*J(1,2)) -
           J(0,1) * (J(1,0)*J(2,2) - J(2,0)*J(1,2)) +
           J(0,2) * (J(1,0)*J(2,1) - J(2,0)*J(1,1));
}

double Tetrahedron10_ip::IJacobian (const Point &loc, RDenseMatrix &IJ) const
{
    dASSERT(IJ.nRows() == 3 && IJ.nCols() == 3, "Parameter 2 wrong dimension");
    static RDenseMatrix J(3,3);
    double d = Jacobian (loc, J);
    double id = 1.0/d;
    IJ(0,0) = ( J(1,1)*J(2,2) - J(2,1)*J(1,2)) * id;
    IJ(1,0) = (-J(1,0)*J(2,2) + J(2,0)*J(1,2)) * id;
    IJ(2,0) = ( J(1,0)*J(2,1) - J(2,0)*J(1,1)) * id;
    IJ(0,1) = (-J(0,1)*J(2,2) + J(2,1)*J(0,2)) * id;
    IJ(1,1) = ( J(0,0)*J(2,2) - J(2,0)*J(0,2)) * id;
    IJ(2,1) = (-J(0,0)*J(2,1) + J(2,0)*J(0,1)) * id;
    IJ(0,2) = ( J(0,1)*J(1,2) - J(1,1)*J(0,2)) * id;
    IJ(1,2) = (-J(0,0)*J(1,2) + J(1,0)*J(0,2)) * id;
    IJ(2,2) = ( J(0,0)*J(1,1) - J(1,0)*J(0,1)) * id;
    return d;
}

double Tetrahedron10_ip::DetJ (const Point &loc, const NodeList *nlist) const
{
#ifndef TET10IP_STORE_COORDS
    dASSERT (nlist != 0, "Node list required");
    double n0x = (*nlist)[Node[0]][0], n0y = (*nlist)[Node[0]][1],
           n0z = (*nlist)[Node[0]][2];
    double n1x = (*nlist)[Node[1]][0], n1y = (*nlist)[Node[1]][1],
           n1z = (*nlist)[Node[1]][2];
    double n2x = (*nlist)[Node[2]][0], n2y = (*nlist)[Node[2]][1],
           n2z = (*nlist)[Node[2]][2];
    double n3x = (*nlist)[Node[3]][0], n3y = (*nlist)[Node[3]][1],
           n3z = (*nlist)[Node[3]][2];
    double n4x = (*nlist)[Node[4]][0], n4y = (*nlist)[Node[4]][1],
           n4z = (*nlist)[Node[4]][2];
    double n5x = (*nlist)[Node[5]][0], n5y = (*nlist)[Node[5]][1],
           n5z = (*nlist)[Node[5]][2];
    double n6x = (*nlist)[Node[6]][0], n6y = (*nlist)[Node[6]][1],
           n6z = (*nlist)[Node[6]][2];
    double n7x = (*nlist)[Node[7]][0], n7y = (*nlist)[Node[7]][1],
           n7z = (*nlist)[Node[7]][2];
    double n8x = (*nlist)[Node[8]][0], n8y = (*nlist)[Node[8]][1],
           n8z = (*nlist)[Node[8]][2];
    double n9x = (*nlist)[Node[9]][0], n9y = (*nlist)[Node[9]][1],
           n9z = (*nlist)[Node[9]][2];
#endif
    dASSERT (loc.Dim() == 3, "Parameter 1 wrong dimension");
    double t1 = 4.0*loc[0];
    double t2 = 4.0*loc[1];
    double t3 = 4.0*loc[2];
    double t4 = t1-1.0;
    double t5 = t2-1.0;
    double t6 = t3-1.0;
    double t7 = t4+t5+t6;
    double t8 = 2.0*t4+t5+t6;
    double t9 = t4+2.0*t5+t6;
    double t10 = t4+t5+2.0*t6;

    double j00 = t7*n0x + t4*n1x - t8*n4x + t2*(n7x-n5x) + t3*(n8x-n6x);
    double j01 = t7*n0y + t4*n1y - t8*n4y + t2*(n7y-n5y) + t3*(n8y-n6y);
    double j02 = t7*n0z + t4*n1z - t8*n4z + t2*(n7z-n5z) + t3*(n8z-n6z);

    double j10 = t7*n0x + t5*n2x - t9*n5x + t1*(n7x-n4x) + t3*(n9x-n6x);
    double j11 = t7*n0y + t5*n2y - t9*n5y + t1*(n7y-n4y) + t3*(n9y-n6y);
    double j12 = t7*n0z + t5*n2z - t9*n5z + t1*(n7z-n4z) + t3*(n9z-n6z);

    double j20 = t7*n0x + t6*n3x - t10*n6x + t1*(n8x-n4x) + t2*(n9x-n5x);
    double j21 = t7*n0y + t6*n3y - t10*n6y + t1*(n8y-n4y) + t2*(n9y-n5y);
    double j22 = t7*n0z + t6*n3z - t10*n6z + t1*(n8z-n4z) + t2*(n9z-n5z);

    return j00 * (j11*j22 - j21*j12) -
           j01 * (j10*j22 - j20*j12) +
           j02 * (j10*j21 - j20*j11);
}
