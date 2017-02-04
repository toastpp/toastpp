// -*-C++-*-
// ==========================================================================
// Module libfe
// File tri3.h
// Declaration of unstructured 2-D element class Triangle3
// (3-noded triangle, 1st order shape functions)
//
//
// Node arrangement:   N2
// -----------------    +        The local element has node coordinates
//                      |\               N0 = (0,0)
//                      | \              N1 = (1,0)
//                      |  \             N2 = (0,1)
//                      |   \
//                      |    \   Sides:  side 0 contains N0,N1 (y=0)
//                      |     \          side 1 contains N1,N2 (x+y=1)
//                      |      \         side 2 contains N2,N0 (x=0)
//                    N0+-------+N1
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
//           ---> Triangle3
// ==========================================================================

#ifndef __TRI3_H
#define __TRI3_H

#include "toastdef.h"

class Surface;
class Mesh;

/**
 * \brief A 3-noded (linear) 2-dimensional triangle element.
 *
 * The node arrangement is given by
 * \image html tri3.png "Node arrangement"
 */
class FELIB Triangle3: public Element_Unstructured_2D {
public:

    /**
     * \brief Creates a new triangle, but does not assign any node indices.
     */
    Triangle3();

    /**
     * \brief Creates a new triangle as a copy of an existing one.
     * \param el source element
     */
    Triangle3 (const Triangle3& el);

    /**
     * \brief Destroys the triangle.
     */
    ~Triangle3();

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    /**
     * \brief Initialises the triangle.
     * \param nlist mesh node list
     * \note This method should be called one the mesh node indices have been
     *   assigned, to allow setting up element parameters.
     * \note Calculates Natural coordinates, and subsampling abscissae.
     *   Precalculates some integrals.
     */
    void Initialise (const NodeList& nlist);

    BYTE Type () const { return ELID_TRI3; }
    // returns element type id

    BYTE VtkType() const { return 5; }

    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    // Return element capability flags

    int nNode() const { return 3; }
    // returns number of nodes per element

    int nSide () const { return 3; };
    // number of sides per element

    int nSideNode (int /*side*/) const { return 2; };
    // returns number of nodes attached to side 'side'

    int SideNode (int side, int node) const;
    // returns number of 'node'-th node of side 'side'

    double SideSize (int sd, const NodeList &nlist) const;
    // returns length of side 'sd'

    /**
     * \brief Returns determinant of Jacobian
     * \param loc evaluation point (ignored for this element type)
     * \param nlist mesh node list (not required for this element type)
     * \return 2*Size()
     */
    double DetJ (const Point &loc, const NodeList *nlist = 0) const
    { return Size()*2.0; }

    Point Local (const NodeList& nlist, const Point& glob) const;
    // returns the local coordinate corresponding to global coordinate 'glob'

    Point NodeLocal (int node) const;
    // returns local coordinates of node 'node'

    Point SurfToLocal (int side, const Point &p) const;
    // Map 1-D surface point p (0-1) along side 'side' into 2-D local coords

    RVector DirectionCosine (int side, RDenseMatrix &jacin);
    // returns a vector of direction cosines in global coordinates of the
    // normal to side 'side'.

    const RVector &LNormal (int side) const;

    bool LContains (const Point& loc, bool pad = true) const;
    // returns TRUE if point 'loc' (in local coordinates of the element) is
    // within the standard triangle

    /**
     * \brief Return the values of all shape functions for a point inside
     *   the element, given in local coordinates.
     * \param loc point in local coordinates
     * \return Vector of shape function values u_i(loc) (i=0..2)
     */
    RVector LocalShapeF (const Point &loc) const;

    /**
     * \brief Return the values of all shape functions for a point inside
     *   the element, given in local coordinates.
     * \param [in] loc point in local coordinates
     * \param [out] fun pointer to vector receiving the function values
     * \note fun must point to a valid RVector.
     * \note If the length of fun != 3 on entry, fun is resized to 3.
     */
    void LocalShapeF (const Point &loc, RVector *fun) const;

    RDenseMatrix LocalShapeD (const Point &loc) const;
    // 2x3 matrix of shape function derivatives (d u_j)/(d x_i) at point 'loc'
    // given in local element coordinates

    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    // returns the shape functions for each node at a point, given in global
    // coordinates

    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;
    // returns derivatives of shape functions in global coordinates at point
    // `glob', given in global coordinates

    double IntF (int i) const;

    RSymMatrix IntFF () const;
    // Return integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

    int QuadRule (int order, const double **wght, const Point **absc) const;

    double IntFF (int i, int j) const;
    // Return a single element of IntFF

    double IntFFF (int i, int j, int k) const;
    // returns a single element of integral over element of product of three
    // shape functions:
    // IntFFF = Int_el { F_i(r) * F_j(r) * F_k(r) } dr

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

    void IntFFF (double &iii, double &iij, double &ijk) const;
    // returns the values of the FFF tensor for:
    //   all indices equal (iii)
    //   two indices equal (iij)
    //   all indices different (ijk)
#endif

    RSymMatrix IntPFF (const RVector& P) const;
    // Returns integral over element of product of two shape functions and a
    // function P defined on nodes:
    // PFF = Int_el { P(r) F_i(r) F_J(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPFF (int i, int j, const RVector& P) const;
    // Returns a single element of IntPFF

    RVector IntFD (int i, int j) const;
    // Returns single (vector) element of FD integral over element:
    // IntFD = Int_el { F_i(r) D_j(r) } dr

    double IntFDD (int i, int j, int k) const
    { return intdd.Get(j,k) * 0.3333333333; }
    // returns a single element of integral over element:
    // IntFDD = Int_el { F_i(r) * D_j(r) * D_k(r) } dr

    RSymMatrix IntPDD (const RVector& P) const
    { return intdd * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns integral over element of product of two shape derivatives and a
    // function P defined in nodal basis:
    // PDD = Int_el { P(r) D_i(r) D_j(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPDD (int i, int j, const RVector &P) const
    { return intdd.Get(i,j) * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns a single element of IntPDD

    RSymMatrix BndIntPFF (const RVector &P) const
	{ return intbff * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // This is a hack! Really P would have to be taken into the integral,
    // rather than just averaged.

    double BndIntPFF (int i, int j, const RVector &P) const
    { return intbff.Get(i,j) * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns a single element of BndIntPFF

    double IntFd (int i, int j, int k) const;
    // Returns a single element of Int [u_i du_j/dx_k] dr

    double IntPd (const RVector &P, int j, int k) const;
    // Returns Int [p(r) du_j/dx_k] dr where p(r) is defined by its nodal
    // basis coefficients P.

    double IntFdd (int i, int j, int k, int l, int m) const;
    // Int u_i du_j/dx_l du_k/dx_m dr

    double IntPdd (const RVector &p, int j, int k, int l, int m) const;
    // Int f(r) du_j/dx_l du_k/dx_m dr
    // where f(r) is given as a nodal vector

    RSymMatrix Intdd() const;
    // returns matrix of mixed derivatives

    /**
     * \brief %Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. 2)
     * \param sd side index (range 0 .. 2)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntF, BndIntFF, SurfIntFF
     */
    double SurfIntF (int i, int sd) const;

    /**
     * \brief %Surface integral of a product of two shape functions over one of
     *   the sides of the element.
     * \param i first node index (range 0 .. 2)
     * \param j second node index (range 0 .. 2)
     * \param sd side index (range 0 .. 2)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;

    RVector BndIntFX (int side, double (*func)(const Point&),
        const NodeList &nlist) const;
    // See element.h for description

    RVector BndIntFCos (int side, const Surface *surf, const RVector &cntcos,
        double sigma, double sup, const NodeList &nlist) const;
    // See element.h for description

    RVector BndIntFCos (int side, const RVector &cntcos, double a,
        const NodeList &nlist) const;
    // See element.h for description

    RVector BndIntFDelta (int side, const Surface *surf, const RVector &pos,
        const NodeList &nlist) const;
    // See element.h for description

    int GetLocalSubsampleAbsc (const Point *&absc) const;
    // returns abscissae in local coords for numerical integration by uniform
    // subsampling. Return value is the number of samples returned in 'absc'

    int GetBndSubsampleAbsc (int side, const Point *&absc) const;
    // abscissae for numerical integration over boundary side in local
    // coordinates of boundary element (dim-1). Use SurfToLocal to convert
    // into local element coordinates

    int Intersection (const Point &p1, const Point &p2, Point *s,
        bool add_endpoints, bool boundary_only);
    // creates a list of points where the line defined by p1 and p2 intersects
    // the element (in local coordinates) or starts/ends within the element.
    // Returns the length of the list

    void SplitSide (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2);
    // Subdivide the triangle such that 'side' is split into two parts.
    // If newnode >= 0 and nbr1, nbr2 != NULL, then the side neighbour has
    // already been split. In that case, newnode is the node index of the
    // midpoint node, and nbr1, nbr2 are the new neighbours for the two
    // side segments.
    // Otherwise, the SplitSide method is also called for the neighbour
    // SplitSide either performs a bisection if the current subdivision
    // level is even, otherwise it performs a merge and resplit operation.
    // In the latter case, it also calls SplitSide for any additionally
    // subdivided sides.
    // On return, el1 and el2 are pointers to the two elements that now
    // face the split side (i.e. the neighbours of nbr1 and nbr2, if those
    // were provided)


    void Bisect (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2);

    void MergeAndResplit (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2);

protected:

    double ComputeSize (const NodeList &nlist) const;
    // area of triangle in global coordinates

    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr

    void ComputeIntFD (const NodeList &nlist);
    // Calculates the FD integrals over the element:
    // IntFD = Int_el { F_i(r) D_j(r) } dr

    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

  //private:

    double *intfd_0, *intfd_1;
    // IntFD matrix storage

    double a0, b0, c0, a1, b1, c1, a2, b2, c2;
    // triangle geometry parameters

#ifdef TRI3_STORE_INTFF
    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
#endif

};

#endif //!__TRI3_H
