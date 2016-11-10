// -*-C++-*-
// =========================================================================
// TOAST v.15                                      (c) Martin Schweiger 1999
// Library: libfe     File: tet10_ip.h
//
// Declaration of class Tetrahedron10
// 10-noded 2nd order isoparametric tetrahedron
//
// Node arrangement:
//                                          Side          contains nodes
//           ^ z                            0 (z=0)       0,1,2,4,7,5
//           |                              1 (y=0)       0,3,1,6,8,4
//                      ^ y                 2 (x=0)       0,2,3,5,9,6
//         N3+-_  N9   /                    3 (1-x-y-z=0) 1,3,2,8,9,7
//           |\ -+-_
//           |  \   -+N2                Node coords of local element:
//           |   \  . \                     N0 = (0,0,0)
//           |     \   \                    N1 = (1,0,0)
//         N6+    . +N8 \                   N2 = (0,1,0)
//           | N5+   \   +N7                N3 = (0,0,1)
//           |  .      \  \                 N4 = (1/2,0,0)
//           | .        \  \                N5 = (0,1/2,0)
//           |.           \ \               N6 = (0,0,1/2)
//           +-------+-------+  --> x       N7 = (1/2,1/2,0)
//          N0       N4      N1             N8 = (1/2,0,1/2)
//                                          N9 = (0,1/2,1/2)
// Inheritance:
// ------------
// Element 
// ---> Element_Unstructured
//      ---> Element_Unstructured_3D
//           ---> Tetrahedron10_ip
// =========================================================================

#ifndef __TET10_IP_H
#define __TET10_IP_H

#include "toastdef.h"

class FELIB Tetrahedron10_ip: public Element_Unstructured_3D {
public:

    Tetrahedron10_ip () { Node = new int[nNode()]; }
    Tetrahedron10_ip (const Tetrahedron10_ip &el);
    ~Tetrahedron10_ip () { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_TET10_IP; }
    BYTE VtkType() const { return 24; }
    unsigned long GetCaps () const { return ELCAPS_CURVED_BOUNDARY; }
    int nNode() const { return 10; }
    int nSide() const { return 4; }
    int nSideNode (int /*side*/) const { return 6; }
    int SideNode (int side, int node) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;

    double IntF (int i) const
    { ERROR_UNDEF; return 0; }
    RSymMatrix IntFF() const;
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
    // dimension 3x3 on input). Return value is det J

    double IJacobian (const Point &loc, RDenseMatrix &IJ) const;
    // Set IJ to the inverse of the Jacobian at local point loc (IJ must be
    // of dimension 3x3 on input). Return value ist det J

    double DetJ (const Point &loc, const NodeList *nlist = 0) const;
    // Return value of det J at local point loc
    // nlist is required ifndef TRI10IP_STORE_COORDS

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    double a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    double side_size[4]; // surface triangle areas

#ifdef TET10IP_STORE_COORDS
    double n0x, n0y, n0z, n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z;
    double n4x, n4y, n4z, n5x, n5y, n5z, n6x, n6y, n6z, n7x, n7y, n7z;
    double n8x, n8y, n8z, n9x, n9y, n9z;
    // Global node coordinates
#endif

private:

    RSymMatrix ComputeIntFF (const NodeList &nlist) const;

    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise

};

#endif // !__TET10_H
