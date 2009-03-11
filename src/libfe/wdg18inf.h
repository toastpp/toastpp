// -*-C++-*-
// ==========================================================================
// Module libfe
// File wdg18inf.h
// Declaration of class Wedge18inf
//
// 18-noded infinite wedge element
// Element is finite triangle in x,y, and infinite in z, where
// local coordinates -1,0,1 map to z1,z2,inf
//
// Global element geometry is restricted such that the triangular
// cross section is parallel to the xy-plane, and z-axis is infinite
// ==========================================================================

#ifndef __WDG18INF_H
#define __WDG18INF_H

class Wedge18inf: public Element_Unstructured_3D {
public:
    Wedge18inf () { Node = new int[nNode()]; }
    Wedge18inf (const Wedge18inf &el);
    ~Wedge18inf () { delete []Node; }

    void Initialise (const NodeList &nlist);
    BYTE Type() const { return ELID_WDG18INF; }

    int nNode() const { return 6; }
    // only the 6 nodes of the triangular base are stored

    RSymMatrix IntFF () const;

    std::istream &operator>> (std::istream &i);
    std::ostream &operator<< (std::ostream &o) const;

protected:
    double z0, z1;
    // global z-coordinates of the two nodes in the infinite direction
    // the 3rd node is at infinity

    RDenseMatrix LocalShapeD_Z (double loc) const;
    double JacobianZ (double loc, RDenseMatrix &J) const;

};

#endif // !__WDG18INF_H
