// -*-C++-*-
// =========================================================================
// TOAST v.15                                      (c) Martin Schweiger 1999
// Library: libfe     File: tri6_ip.h
//
// 6-noded isoparametric triangle for quadratic shape functions and curved
// boundaries.
//
//            ^                       Side 0: eta    = 0
//            |                       Side 1: xi+eta = 1
// -------                            Side 2: xi     = 0
// Local     N2
// element    +                       Node coordinates:
// -------    | \                     N0 = (0,0)
//            |   \                   N1 = (1,0)
//            |     \                 N2 = (0,1)
//          N5+       +N4             N3 = (1/2,0)
//            |         \             N4 = (1/2,1/2)
//            |           \           N5 = (0,1/2)
//            |             \
//          N0+-------+-------+N1
//                   N3
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
//           ---> Triangle6_ip
// =========================================================================

#ifndef __TRI6_IP_H
#define __TRI6_IP_H

#include "toastdef.h"

class FELIB Triangle6_ip: public Element_Unstructured_2D {
public:

    Triangle6_ip () { Node = new int[nNode()]; }
    Triangle6_ip (const Triangle6_ip &el);
    ~Triangle6_ip () { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList& nlist);

    BYTE Type() const { return ELID_TRI6_IP; }
    unsigned long GetCaps () const { return ELCAPS_CURVED_BOUNDARY; }
    int nNode() const { return 6; }
    int nSide() const { return 3; }
    int nSideNode (int /*side*/) const { return 3; }
    int SideNode (int side, int node) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;
    bool GContains (const Point& glob, const NodeList& nlist) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;

    double IntF (int i) const
    { ERROR_UNDEF; return 0; }
    RSymMatrix IntFF () const;
    double IntFF (int i, int j) const;
    double IntFFF (int i, int j, int k) const;
    RSymMatrix IntPFF (const RVector& P) const;
    double IntPFF (int i, int j, const RVector& P) const;
    double IntFDD (int i, int j, int k) const;
    RSymMatrix IntPDD (const RVector& P) const
    { ERROR_UNDEF; return RSymMatrix(); }
    double IntPDD (int i, int j, const RVector &P) const;
    double SurfIntF (int i, int sd) const
    { ERROR_UNDEF; return 0; }
    double SurfIntFF (int i, int j, int sd) const
    { ERROR_UNDEF; return 0; }
    RSymMatrix BndIntPFF (const RVector &P) const
    { ERROR_UNDEF; return RSymMatrix(); }
    double BndIntPFF (int i, int j, const RVector &P) const;
    int Intersection (const Point &p1, const Point &p2, Point *s,
	bool add_endpoints, bool boundary_only)
    { ERROR_UNDEF; return 0; }


protected:

    double Jacobian (const Point &loc, RDenseMatrix &J) const;
    // Set J to the Jacobian of the element at local point loc (J must be of
    // dimension 2x2 on input). Return value is det J

    double IJacobian (const Point &loc, RDenseMatrix &IJ) const;
    // Set IJ to the inverse of the Jacobian at local point loc (IJ must be
    // of dimension 2x2 on input). Return value ist det J

    double DetJ (const Point &loc, const NodeList *nlist = 0) const;
    // Return value of det J at local point loc
    // nlist is required ifndef TRI6IP_STORE_COORDS

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    double a0, b0, c0, a1, b1, c1, a2, b2, c2;

#ifdef TRI6IP_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y, n5x, n5y;
    // Global node coordinates
#endif

private:

    RSymMatrix ComputeIntFF (const NodeList &nlist) const;

    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
};

#endif // !__TRI6_IP_H
