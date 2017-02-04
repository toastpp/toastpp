// -*-C++-*-
// ==========================================================================
// Module libfe
// File tet4.h
// Declaration of class Tetrahedron4
//
// 4-noded triangle for implementation of linear shape functions
//
// Node arrangement:
//
//           ^ z
//           |
//                      ^ y                 Side          contains nodes
//         N3+-_       /                    0 (z=0)       0,1,2
//           |\ --_                         1 (y=0)       0,3,1
//           |  \  --+N2                    2 (x=0)       0,2,3
//           |   \  . \                     3 (1-x-y-z=0) 1,3,2
//           |     \   \                    -------------------
//           |    . \   \                   Node coords of local element:
//           |   .   \   \                  N0 = (0,0,0)
//           |  .      \  \                 N1 = (1,0,0)
//           | .        \  \                N2 = (0,1,0)
//           |.           \ \               N3 = (0,0,1)
//           +---------------+  --> x
//          N0               N1
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_3D
//           ---> Tetrahedron4
// ==========================================================================

#ifndef __TET4_H
#define __TET4_H

/**
 * \brief A 4-noded 3-dimensional tetrahedron element with straight edges
 *   and linear shape functions.
 *
 * The node arrangement is given by
 * \code
 *    ^ z                  The local element has node coordinates
 *    |                            N0 = (0,0,0)
 *               ^ y               N1 = (1,0,0)
 *  N3+-_       /                  N2 = (0,1,0)
 *    |\ --_                       N3 = (0,0,1)
 *    |  \  --+N2          
 *    |   \  . \           Sides: side 0 (z=0)       contains N0,N1,N2
 *    |     \   \                 side 1 (y=0)       contains N0,N3,N1
 *    |    . \   \                side 2 (x=0)       contains N0,N2,N3
 *    |   .   \   \               side 3 (1-x-y-z=0) contains N1,N3,N2
 *    |  .      \  \
 *    | .        \  \
 *    |.           \ \
 *    +---------------+  --> x
 *   N0               N1
 * \endcode
 */
class FELIB Tetrahedron4: public Element_Unstructured_3D {
public:

    Tetrahedron4() { Node = new int[nNode()]; }
    Tetrahedron4( const Tetrahedron4 &el);
    ~Tetrahedron4() { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_TET4; }
    BYTE VtkType() const { return 10; }
    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    int nNode() const { return 4; }
    int nSide() const { return 4; }
    int nSideNode (int /*side*/) const { return 3; }
    int SideNode (int side, int node) const;
    double SideSize (int sd, const NodeList &nlist) const;

    double DetJ (const Point &loc, const NodeList *nlist = 0) const
    { return Size()*6.0; }

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    Point SurfToLocal (int side, const Point &p) const;
    void MapToSide (int side, Point &loc) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;

    double IntF (int i) const;
    double IntFF (int i, int j) const;
    RSymMatrix IntFF() const;
    double IntFFF (int i, int j, int k) const;
    RSymMatrix IntPFF (const RVector& P) const;
    double IntPFF (int i, int j, const RVector& P) const;
    double IntFDD (int i, int j, int k) const;
    RSymMatrix IntPDD (const RVector& P) const;
    double IntPDD (int i, int j, const RVector &P) const;

    double IntFd (int i, int j, int k) const;
    // Int [u_i du_k/dx_k] dr

    double IntPd (const RVector &P, int j, int k) const;
    // Returns Int [p(r) du_j/dx_k] dr where p(r) is defined by its nodal
    // basis coefficients P.

    double IntFdd (int i, int j, int k, int l, int m) const;
    double IntPdd (const RVector &p, int j, int k, int l, int m) const;

    RSymMatrix Intdd() const;
    // returns matrix of mixed derivatives

    double Intd (int i, int k) const;

    double IntFfd (int i, int j, int k, int l) const;
    // Int u_i u_j du_k/dx_l dr
    double IntPfd(const RVector &p,int j,int k,int l) const;
    // Int f(r) u_j du_k/du_l dr

    /**
     * \brief Boundary integral of all shape functions over all boundary
     *   sides of the element.
     * \return Vector of size 4, containing the integrals
     *   \f[ \int_{\partial\Omega} u_i(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The returned vector contains nonzero entries at index i only if
     *   node i is a boundary node.
     * \note If the element does not contain boundary sides, the returned
     *   vector is zero.
     * \sa SurfIntF, BndIntFF, SurfIntFF
     */
    RVector BndIntF () const;

    /**
     * \brief Boundary integral of a shape function over all boundary sides
     *   of the element.
     * \param i node index (range 0 .. 3)
     * \return Value of the integral
     *   \f[ \int_{\partial\Omega} u_i(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The returned value is nonzero only if node i is a boundary node.
     */
    double BndIntF (int i) const;
    
    /**
     * \brief %Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. 3)
     * \param sd side index (range 0 .. 3)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntF, BndIntFF, SurfIntFF
     */
    double SurfIntF (int i, int sd) const;

    /**
     * \brief %Surface integral of a product of two shape functions over one of
     *   the sides of the element.
     * \param i first node index (range 0 .. 3)
     * \param j second node index (range 0 .. 3)
     * \param sd side index (range 0 .. 3)
     * \return Value of the integral
     *   \f$ \int_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;

    RSymMatrix BndIntPFF (const RVector &P) const
    { ERROR_UNDEF; return RSymMatrix(); }
    double BndIntPFF (int i, int j, const RVector &P) const;

    /**
     * \brief Calculates integral
     *   \f$ \int_s \frac{partial u_i(s)}{\partial x_j} u_k(s) ds \f$
     *   over an element side.
     * \param sd side index (>= 0)
     * \param i node index 1 (>= 0)
     * \param j direction index (0 <= j < 3)
     * \param k node index 2 (>= 0)
     * \note i and k are element-relative node indices. Side sd must contain
     *   both i and k.
     */
    double BndIntFD (int sd, int i, int j, int k);

    /**
     * \brief Calculates integral
     *   \f$ \int_s \frac{partial u_i(s)}{\partial x_j} u_k(s) ds \f$
     *   over an element side.
     * \param sd side index (>= 0)
     * \param el2 element index for neighbour element
     * \param i node index 1 (>= 0)
     * \param j direction index (0 <= j < 3)
     * \param k node index 2 (>= 0) in neighbour element
     * \note In this version, nodes i and k refer to different elements.
     *   The derivative term du_i/dx_j refers to element *this, while
     *   the shape function term u_k refers to the neighbour element (el2).
     *   The neigbour element must be of the same type (Tetrahedron4), and
     *   the two elements must be joined at side sd.
     */
    //double BndIntFD (Mesh &mesh, int sd, int el2, int sd2, int i, int j, int k);

    RVector BndIntFX (int side, double (*func)(const Point&),
        const NodeList &nlist) const;
    RVector BndIntFCos (int side, const RVector &cntcos, double a,
        const NodeList &nlist) const;
    int GetLocalSubsampleAbsc (const Point *&absc) const;
    int GetBndSubsampleAbsc (int side, const Point *&absc) const;
    RDenseMatrix StrainDisplacementMatrix (const Point &glob) const;
    RDenseMatrix ElasticityStiffnessMatrix (const RDenseMatrix &D) const;
    RDenseMatrix ElasticityStiffnessMatrix (double modulus, double pratio)
        const;
    RVector InitialStrainVector (double E, double nu, const RVector &e0);
    RVector ThermalExpansionVector (double E, double nu, double alpha,
        double dT);
    RVector DThermalExpansionVector (double E, double nu);

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

#ifdef TET4_STORE_INTFF
    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
#endif

};

#endif // !__TET4_H
