// -*-C++-*-
// ==========================================================================
// Module libfe
// File quad4.h
// Declaration of class Quad4
// Quadratic element aligned with the coordinate axes
//
// N2 +-------+ N3     y
//    |       |        ^
//    |       |        |
//    |       |        |
// N0 +-------+ N1     +-----> x
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
//           ---> Quad4
// ==========================================================================

#ifndef __QUAD4_H
#define __QUAD4_H

#include "toastdef.h"

class Surface;

class Quad4: public Element_Unstructured_2D {
public:

    Quad4 () { Node = new int[nNode()]; }
    Quad4 (const Quad4& el);
    ~Quad4 () { delete []Node; }
    // constructors, destructor

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_QUAD4; }

    unsigned long GetCaps () const { return 0; }

    int nNode() const { return 4; }
    int nSide() const { return 4; }
    int nSideNode (int) const { return 2; }
    int SideNode (int side, int node) const;

    Point Local (const NodeList &nlist, const Point &glob) const;
    Point NodeLocal (int node) const;
    bool LContains (const Point& loc, bool pad = true) const;
    bool GContains (const Point& glob, const NodeList& nlist) const;
    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;

private:
    double len; // side length
    Point n0;   // position of node N0
};

#endif // !__QUAD4_H
