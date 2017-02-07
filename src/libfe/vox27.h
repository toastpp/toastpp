// -*-C++-*-
// ==========================================================================
// Module libfe
// File vox27.h
// Declaration of class Voxel27
//
// 27-noded regular voxel element aligned with global coordinate axes,
// to form a regular mesh grid. Uses tri-quadratic shape functions
//
// Node arrangement:
//
//               N6      N15      N7        Side          contains nodes
//                +-------+-------+         0 (z=0)       0,1,2,3
//               /.              /|         1 (z=1)       4,5,6,7
//           N13+ .          N14+ |         2 (y=0)       0,1,4,5 
//             /  +N18         /  +N19      3 (y=1)       2,3,6,7 
//            /   . N12       /   |         4 (x=0)       0,2,4,6
//         N4+-------+-------+N5  |         5 (x=1)       1,3,5,7
//           |    .          |    |         ---------------------
//           |    .N2    N11 |    |         Node coords of local element:
//           |    +.......+..|....+N3       N0 = (0,0,0)    N16 = (0,0,1/2)
//        N16+   .        N17+   /          N1 = (1,0,0)    N17 = (1,0,1/2)
//           |  +N9          |  +N10        N2 = (0,1,0)    N18 = (0,1,1/2)
//           | .             | /            N3 = (1,1,0)    N19 = (1,1,1/2)
//     z     |.              |/             N4 = (0,0,1)    N20 = (1/2,1/2,0)
//     ^     +-------+-------+  --> x       N5 = (1,0,1)    N21 = (1/2,0,1/2)
//     |  y N0       N8      N1             N6 = (0,1,1)    N22 = (0,1/2,1/2)
//     | /                                  N7 = (1,1,1)    N23 = (1,1/2,1/2)
//     |/           (nodes 20-26            N8 = (1/2,0,0)  N24 = (1/2,1,1/2)
//     +----->x      not shown)             N9 = (0,1/2,0)  N25 = (1/2,1/2,1)
//                                         N10 = (1,1/2,0)  N26 = (1/2,1/2,1/2)
//                                         N11 = (1/2,1,0)
// Inheritance:                            N12 = (1/2,0,1)
// ------------                            N13 = (0,1/2,1)
// Element                                 N14 = (1,1/2,1)
// ---> Element_Structured                 N15 = (1/2,1,1)
//      ---> Element_Structured_3D
//           ---> Voxel27
// ==========================================================================

#ifndef __VOX27_H
#define __VOX27_H

class Voxel27: public Element_Structured_3D {
public:

    Voxel27 () { Node = new int[nNode()]; }
    Voxel27 (const Voxel27 &el);
    ~Voxel27 () { delete []Node; }

    void Initialise (const NodeList &nlist);

    inline BYTE Type() const { return ELID_VOX27; }
    inline unsigned long GetCaps() const { return 0; }
    inline int nNode() const { return 27; }
    inline int nSide() const { return 6; }
    inline int nSideNode (int /*side*/) const { return 9; }
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
    { return intf[i]; }
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
    
    /**
     * \brief %Surface integral of a product of two shape functions over one of
     *   the sides of the element.
     * \param i first node index (range 0 .. 26)
     * \param j second node index (range 0 .. 26)
     * \param sd side index (range 0 .. 5)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa SurfIntF(int,int)const, BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;

    RSymMatrix BndIntFF () const
    { ERROR_UNDEF; return RSymMatrix(); }
    RSymMatrix BndIntPFF (const RVector &P) const
    { ERROR_UNDEF; return RSymMatrix(); }
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

    int Intersection (const Point &p1, const Point &p2, Point *s,
	bool add_endpoints, bool boundary_only)
    { ERROR_UNDEF; return 0; }

protected:
    void ComputeIntF () const;
    void ComputeIntFF () const;
    void ComputeIntFFF () const;
    void ComputeIntDD () const;
    void ComputeIntFDD () const;
    void ComputeBndIntFF () const;
    void ComputeBndIntFFF () const;

private:
    double x0, y0, z0; // global coords of node 0

    // shared properties
    static double dx, dy, dz; // voxel edge lengths
    static double size;       // voxel volume
    static RVector intf;
    static RSymMatrix intff;
    static RSymMatrix intfff[27];
    static RSymMatrix intdd;
    static RSymMatrix intfdd[27];
    static RSymMatrix bndintff[6];
    static RDenseMatrix bndintfff[6][27];
};

#endif // !__VOX27_H
