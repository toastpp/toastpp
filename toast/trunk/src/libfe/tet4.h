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

    void Initialise (const NodeList &nlist);

    BYTE Type() const { return ELID_TET4; }
    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    int nNode() const { return 4; }
    int nSide() const { return 4; }
    int nSideNode (int /*side*/) const { return 3; }
    int SideNode (int side, int node) const;
    double SideSize (int sd, const NodeList &nlist) const;

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

    RSymMatrix IntFF() const;
    double IntF (int i) const;
    double IntFF (int i, int j) const;
    double IntFFF (int i, int j, int k) const;
    RSymMatrix IntPFF (const RVector& P) const;
    double IntPFF (int i, int j, const RVector& P) const;
    double IntFDD (int i, int j, int k) const;
    RSymMatrix IntPDD (const RVector& P) const
    { xERROR(Not implemented); return RSymMatrix(); }
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

    RVector BndIntF () const;

    double BndIntFSide (int i, int sd);

    double BndIntFFSide (int i, int j, int sd);

    RSymMatrix BndIntPFF (const RVector &P) const
    { xERROR(Not implemented); return RSymMatrix(); }
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

    int GlobalIntersection (const NodeList &nlist, const Point &p1,
	const Point &p2, Point **list);

    /**
     * \brief Calculate intersection of a ray with element surfaces.
     * \param p1 First point defining the ray
     * \param p2 Second point defining the ray
     * \param pi On return, points to list of intersection points
     * \return Number of points found (should be 0 or 2)
     * \note The ray is assumed to be of infinite length, not just the
     *  segment between p1 and p2
     * \note On return, pi points to a static list. It should not be
     *  deallocated by the caller, and it will be overwritten by the next
     *  call to Intersection.
     * \note If no intersection points are found, pi is set to NULL.
     */
    int Intersection (const Point &p1, const Point &p2, Point** pi);

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
