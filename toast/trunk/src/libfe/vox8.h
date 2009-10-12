// -*-C++-*-
// ==========================================================================
// Module libfe
// File vox8.h
// Declaration of class Voxel8
//
// 8-noded regular voxel element aligned with global coordinate axes,
// to form a regular mesh grid.
//
// Node arrangement:
//
//               N6               N7        Side          contains nodes
//                +---------------+         0 (z=0)       0,1,2,3
//         z ^   /.              /|         1 (z=1)       4,5,6,7
//           |  / .             / |         2 (y=0)       0,1,4,5 
//           | /  .    ^ y     /  |         3 (y=1)       2,3,6,7 
//            /   .   /       /   |         4 (x=0)       0,2,4,6
//        N4 +---------------+ N5 |         5 (x=1)       1,3,5,7
//           |    .          |    |         ---------------------
//           |    .N2        |    | N3      Node coords of local element:
//           |    +..........|....+         N0 = (0,0,0)
//           |   .           |   /          N1 = (1,0,0)
//           |  .            |  /           N2 = (0,1,0)
//           | .             | /            N3 = (1,1,0)
//           |.              |/             N4 = (0,0,1)
//           +---------------+  --> x       N5 = (1,0,1)
//          N0               N1             N6 = (0,1,1)
//                                          N7 = (1,1,1)
// Inheritance:
// ------------
// Element
// ---> Element_Structured
//      ---> Element_Structured_3D
//           ---> Voxel8
// ==========================================================================

#ifndef __VOX8_H
#define __VOX8_H

class FELIB Voxel8: public Element_Structured_3D {
public:

    Voxel8 () { Node = new int[nNode()]; }
    Voxel8 (const Voxel8 &el);
    ~Voxel8 () { delete []Node; }

    void Initialise (const NodeList &nlist);
    void PostInitialisation ();

    inline BYTE Type() const { return ELID_VOX8; }
    inline unsigned long GetCaps() const { return 0; }
    inline int nNode() const { return 8; }
    inline int nSide() const { return 6; }
    inline int nSideNode (int /*side*/) const { return 4; }
    int SideNode (int side, int node) const;

    inline double Size() const { return size; }

    Point Local (const NodeList &nlist, const Point &glob) const
    { return Local (glob); }
    Point Local (const Point &glob) const;
    Point NodeLocal (int node) const;
    const RVector &LNormal (int side) const;
    inline RVector DirectionCosine (int side, RDenseMatrix& jacin)
    { return LNormal (side); }
    bool LContains (const Point &loc, bool pad = true) const;
    bool GContains (const Point &glob, const NodeList&) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const Point &glob) const
    { return LocalShapeF (Local (glob)); }
    RDenseMatrix GlobalShapeD (const Point &glob) const
    { return LocalShapeD (Local (glob)); }
    RVector GlobalShapeF (const NodeList &nlist, const Point &glob) const
    { return GlobalShapeF (glob); }
    RDenseMatrix GlobalShapeD (const NodeList &nlist, const Point &glob) const
    { return GlobalShapeD (glob); }

    inline double IntF (int i) const
    { return intf; }
    inline double IntFF (int i, int j) const
    { return intff(i,j); }
    inline RSymMatrix IntFF() const
    { return intff; }
    inline double IntFFF (int i, int j, int k) const
    { return intfff[i](j,k); }
    double IntPFF (int i, int j, const RVector& P) const;
    RSymMatrix IntPFF (const RVector& P) const;
    inline RSymMatrix IntDD () const { return intdd; }
    inline double IntDD (int i, int j) const { return intdd(i,j); }
    inline double IntFDD (int i, int j, int k) const { return intfdd[i](j,k); }
    RSymMatrix IntPDD (const RVector& P) const;
    double IntPDD (int i, int j, const RVector &P) const;
    double BndIntFF (int i, int j);
    double BndIntFFSide (int i, int j,int sd);
    RSymMatrix BndIntFF () const
    { xERROR(Not implemented); return RSymMatrix(); }
    RSymMatrix BndIntPFF (const RVector &P) const
    { xERROR(Not implemented); return RSymMatrix(); }
    double BndIntPFF (int i, int j, const RVector &P) const;

    RSymMatrix Intdd() const;
    // Int du_j/dx_l du_k/dx_m dr

    double IntFd (int i, int j, int k) const;
    // Int [u_i du_j/dx_k] dr

    double IntFdd (int i, int j, int k, int l, int m) const;
    // Int u_i du_j/dx_l du_k/dx_m dr
    double IntPdd (const RVector &p, int j, int k, int l, int m) const;
    // Int f(r) du_j/dx_l du_k/dx_m dr
    // where f(r) is given as a nodal vector
    
    double IntFfd (int i, int j, int k, int l) const;
    // Int u_i u_j du_k/dx_l dr
    double IntPfd(const RVector &p,int j,int k,int l) const;
    // Int f(r) u_j du_k/du_l dr

    RDenseMatrix StrainDisplacementMatrix (const Point &glob) const;
    RDenseMatrix ElasticityStiffnessMatrix (double modulus, double pratio)
        const;

    int GlobalIntersection (const NodeList &nlist, const Point &p1,
        const Point &p2, Point **list)
    { xERROR(Not implemented); return 0; }
    int Intersection (const Point &p1, const Point &p2, Point** pi)
    { xERROR(Not implemented); return 0; }

protected:
    void ComputeIntFF () const;
    void ComputeIntFFF () const;
    void ComputeIntDD () const;
    void ComputeIntFDD () const;
    void ComputeBndIntFF () const;
    void ComputeBndIntFFF () const;

private:
    double x0, y0, z0; // global coords of node 0

    // shared properties
    static double dx;
    static double dy;
    static double dz; // voxel edge lengths
    static double size;       // voxel volume
    static double intf;
    static RSymMatrix intff;
    static RSymMatrix intfff[8];
    static RSymMatrix intdd;
    static RSymMatrix intfdd[8];
    static RSymMatrix bndintff[6];
    static RDenseMatrix bndintfff[6][8];

    static bool need_setup;
};

#endif // !__VOX8_H