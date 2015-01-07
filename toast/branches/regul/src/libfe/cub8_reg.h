// -*-C++-*-
// ==========================================================================
// Module libfe
// File cub8_reg.h
// Declaration of class Cube8_reg
//
// 8-noded regular cube (aligned with global coordinate axes)
//
// Node arrangement:
//
//     z    N6               N7       Side     contains nodes
//     ^    o-----------------o       0 (z=0)  0,1,2,3
//     |   /|                /|       1 (z=a)  4,6,7,5
//     |  / |               / |       2 (y=0)  0,4,5,1
//       /  |              /  |       3 (y=a)  2,3,7,6
//      /   |             /   |       4 (z=0)  0,2,6,4
//  N4 o----+------------o N5 |       5 (z=a)  1,5,7,3
//     |    |            |    |
//     |    |  ^ y       |    |
//     |    | /          |    |
//     |    |            |    |
//     |    o------------+----o
//     |   / N2          |   / N3
//     |  /              |  /
//     | /               | /
//     |/                |/
//     o-----------------o  --> x
//     N0               N1
//
// Inheritance:
// ------------
// Element ----> Element3D ----> Cube8_reg
// ==========================================================================

#ifndef __CUB8_REG_H
#define __CUB8_REG_H

class Cube8_reg: public Element3D {
public:
    Cube8_reg () { Node = new int[nNode()]; }
    Cube8_reg (const Cube8_reg &el);
    ~Cube8_reg () { delete []Node; }

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_CUB8_REG; }
    int nNode() const { return 8; }
    int nSide() const { return 6; }
    int nSideNode (int /*side*/) const { return 4; }
    int SideNode (int side, int node) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;

    double IntF (int i) const;
    RDenseMatrix ElasticityStiffnessMatrix (double modulus, double pratio)
        const;

private:
    double dx, dy, dz; // cube edge lengths

};

#endif // !__CUB8_REG_H
