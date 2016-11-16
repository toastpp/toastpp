// -*-C++-*-
// =========================================================================
// TOAST v.15                                      (c) Martin Schweiger 1999
// Library: libfe     File: tri10_ip.h
//
// 10-noded isoparametric triangle for cubic shape functions and curved
// boundaries.
//
//            ^ xi                Side 0: xi     = 0
//            |                   Side 1: xi+eta = 1
// -------                        Side 2: eta    = 0
// Local     N2
// element    +                   Node coordinates:
// -------    | \                 N0 = (0,0)    N5 = (2/3,1/3)
//            |   \               N1 = (1,0)    N6 = (1/3,2/3)
//          N7+     +N6           N2 = (0,1)    N7 = (0,2/3)
//            |       \           N3 = (1/3,0)  N8 = (0,1/3)
//            |         \         N4 = (2/3,0)  N9 = (1/3,1/3)
//          N8+     +     +N5
//            |     N9      \
//            |               \
//          N0+-----+-----+-----+N1 --> xi
//                 N3     N4
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
//           ---> Triangle10_ip
// ========================================================================= //

#ifndef __TRI10_IP_H
#define __TRI10_IP_H

#include "toastdef.h"

class FELIB Triangle10_ip: public Element_Unstructured_2D {
public:

    Triangle10_ip () { Node = new int [nNode()]; }
    Triangle10_ip (const Triangle10_ip &el);
    ~Triangle10_ip () { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);

    BYTE Type () const { return ELID_TRI10_IP; }
    unsigned long GetCaps () const { return ELCAPS_CURVED_BOUNDARY; }
    int nNode () const { return 10; }
    int nSide () const { return 3; }
    int nSideNode (int /*side*/) const { return 4; }
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
    // nlist is required ifndef TRI10IP_STORE_COORDS

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    double a0, b0, c0, a1, b1, c1, a2, b2, c2;

#ifdef TRI10IP_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y;
    double n5x, n5y, n6x, n6y, n7x, n7y, n8x, n8y, n9x, n9y;
    // Global node coordinates
#endif

private:

    RSymMatrix ComputeIntFF (const NodeList &nlist) const;

    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
};

#endif // !__TRI10_IP_H
