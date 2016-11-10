// -*-C++-*-
// ==========================================================================
// Module libfe
// File tet10.h
// Declaration of class Tetrahedron10
//
// 10-noded tetrahedron for implementation of quadratic shape functions
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
//           ---> Tetrahedron10
// ==========================================================================

#ifndef __TET10_H
#define __TET10_H

class FELIB Tetrahedron10: public Element_Unstructured_3D {
public:

    Tetrahedron10() { Node = new int[nNode()]; }
    Tetrahedron10( const Tetrahedron10 &el);
    ~Tetrahedron10() { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_TET10; }
    BYTE VtkType() const { return 24; }
    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    int nNode() const { return 10; }
    int nSide() const { return 4; }
    int nSideNode (int /*side*/) const { return 6; }
    int SideNode (int side, int node) const;
    double SideSize (int sd, const NodeList &nlist) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    Point SurfToLocal (int side, const Point &p) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;

    double IntF (int i) const;
    RSymMatrix IntFF() const;
    double IntFF (int i, int j) const;
    double IntFFF (int i, int j, int k) const;
    RSymMatrix IntPFF (const RVector& P) const;
    double IntPFF (int i, int j, const RVector& P) const;
    double IntFDD (int i, int j, int k) const;
    RSymMatrix IntPDD (const RVector& P) const
    { ERROR_UNDEF; return RSymMatrix(); }
    double IntPDD (int i, int j, const RVector &P) const;

    /**
     * \brief Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. 9)
     * \param sd side index (range 0 .. 3)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntF, BndIntFF, SurfIntFF
     */
    double SurfIntF (int i, int sd) const;
    
    /**
     * \brief $Surface integral of a product of two shape functions over one of
     *   the sides of the element.
     * \param i first node index (range 0 .. 9)
     * \param j second node index (range 0 .. 9)
     * \param sd side index (range 0 .. 3)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;
    
    RSymMatrix BndIntPFF (const RVector &P) const
    { ERROR_UNDEF; return RSymMatrix(); }
    double BndIntPFF (int i, int j, const RVector &P) const;
    RVector BndIntFX (int side, double (*func)(const Point&),
        const NodeList &nlist) const;
    RVector BndIntFCos (int side, const RVector &cntcos, double a,
        const NodeList &nlist) const;
    RDenseMatrix ElasticityStiffnessMatrix (double modulus, double pratio)
        const;
    RVector InitialStrainVector (double E, double nu, const RVector &e);
    RVector ThermalExpansionVector (double E, double nu, double alpha, double dT);
    int GetLocalSubsampleAbsc (const Point *&absc) const;
    int GetBndSubsampleAbsc (int side, const Point *&absc) const;

    /**
     * \brief Return intersection points of a ray with the element surface.
     * \param p1 first ray endpoint (in local element frame)
     * \param p2 second ray endpoint (in local element frame)
     * \param s pointer to list of intersection points
     * \param add_endpoints flag to add ray endpoints to list if they
     *   are located inside the element
     * \param boundary_only flag to look only for intersection points
     *   with boundary sides
     * \return number of intersection points found
     * \note The point buffer \e s must have been assigned to sufficient
     *   length (2 for convex elements) by the caller.
     * \note If no intersections are found, pi is set to NULL.
     * \note If add_enpoints is true and if the ray starts and/or ends
     *   inside the element, the corresponding end points are added to
     *   the list.
     * \note If boundary_only is true, then only intersections with
     *   boundary sides will be returned.
     * \sa GlobalIntersection
     */
    int Intersection (const Point &p1, const Point &p2, Point *s,
	bool add_endpoints, bool boundary_only);

protected:

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    double a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3;
    double side_size[4]; // surface triangle areas

private:

#ifdef TET10_STORE_INTFF
    RSymMatrix intff;
    // stores integral over element of product of 2 shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // valid after call to Initialise
#endif

};

#endif // !__TET10_H
