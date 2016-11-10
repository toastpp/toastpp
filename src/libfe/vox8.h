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

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);
    void PostInitialisation (const NodeList &nlist);

    inline BYTE Type() const { return ELID_VOX8; }
    BYTE VtkType() const { return 11; }
    inline unsigned long GetCaps() const { return 0; }
    inline int nNode() const { return 8; }
    inline int nSide() const { return 6; }
    inline int nSideNode (int /*side*/) const { return 4; }
    int SideNode (int side, int node) const;

    double Size() const;

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

    double IntF (int i) const;

    double IntFF (int i, int j) const;

    RSymMatrix IntFF() const;

    double IntFFF (int i, int j, int k) const;

    double IntPFF (int i, int j, const RVector& P) const;
    RSymMatrix IntPFF (const RVector& P) const;

    RSymMatrix IntDD () const;

    double IntDD (int i, int j) const;

    double IntFDD (int i, int j, int k) const;

    RSymMatrix IntPDD (const RVector& P) const;
    double IntPDD (int i, int j, const RVector &P) const;

    /**
     * \brief Boundary integral of all shape functions over all boundary sides
     *   of the element.
     * \return Vector of size 8 of integrals
     *   \f[ \int_{\partial\Omega} u_i(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The return vector contains nonzero entries only for nodes that
     *   are part of the mesh surface.
     * \note If the element doesn't contain boundary faces, the returned vector
     *   is all zeros.
     * \sa SurfIntF
     */
    RVector BndIntF () const;

    /**
     * \brief Boundary integral of a shape function over all boundary
     *   sides of the element.
     * \param i node index (range 0 .. 7)
     * \return Integral
     *   \f[ \int_{\partial\Omega} u_i(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The return value is nonzero only if node i is a boundary node.
     * \sa SurfIntF, BndIntFF, SurfIntFF
     */
    double BndIntF (int i);

    /**
     * \brief %Surface integral of a shape function over one of the sides of
     *   the element.
     * \param i node index (range 0 .. 7)
     * \param sd side index (range 0 .. 5)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntF, BndIntFF, SurfIntFF
     */
    double SurfIntF (int i, int sd) const;
    
    /**
     * \brief %Surface integral of a product of two shape functions over one of
     *   the sides of the element.
     * \param i first node index (range 0 .. 7)
     * \param j second node index (range 0 .. 7)
     * \param sd side index (range 0 .. 5)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;

    double BndIntFF (int i, int j);
    
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
	bool add_endpoints, bool boundary_only);

protected:
    void ComputeIntFF () const;
    void ComputeIntFFF () const;
    void ComputeIntDD () const;
    void ComputeIntFDD () const;
    void ComputeBndIntF () const;
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
    static RVector bndintf[6];
    static RSymMatrix bndintff[6];
    static RDenseMatrix bndintfff[6][8];

    static bool need_setup;
};

#endif // !__VOX8_H
