// -*-C++-*-
// ==========================================================================
// Module libfe
// File tri6.h
// Declaration of class Triangle6
//
// 6-noded triangle for implementation of quadratic shape functions
//
// Node arrangement:   N2           N3 is on side 0
// -----------------    +           N4 is on side 1
//                      |\          N5 is on side 2
//                      | \         ---------------
//                      |  \        The local element has node coordinates
//                    N5+   +N4           N0 = (0,0)
//                      |    \            N1 = (1,0)
//                      |     \           N2 = (0,1)
//                      |      \          N3 = (0.5,0)
//                    N0+---+---+N1       N4 = (0.5,0.5)
//                          N3            N5 = (0,0.5)
//
// Inheritance:
// ------------
// Element
//  ---> Element_Unstructured
//       ---> Element_Unstructured_2D
//            ---> Triangle6
// ==========================================================================

#ifndef __TRI6_H
#define __TRI6_H

/**
 * \brief A 6-noded 2-dimensional triangle element with straight sides
 *   and second-order shape functions.
 *
 * The node arrangement is given by
 * \code
 *   N2             The local element has node coordinates
 *    +                     N0 = (0,0)
 *    |\                    N1 = (1,0)
 *    | \                   N2 = (0,1)
 *    |  \                  N3 = (0.5,0)
 *  N5+   +N4               N4 = (0.5,0.5)
 *    |    \                N5 = (0,0.5)
 *    |     \
 *    |      \      Sides:  side 0 contains N0,N1,N3 (y=0)
 *  N0+---+---+N1           side 1 contains N1,N2,N4 (x+y=1)
 *        N3                side 2 contains N2,N0,N5 (x=0)
 * \endcode
 */
class FELIB Triangle6: public Element_Unstructured_2D {
public:

    Triangle6() { Node = new int[nNode()]; }
    Triangle6 (const Triangle6 &el);
    ~Triangle6() { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList& nlist);

    BYTE Type() const { return ELID_TRI6; }
    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    int nNode() const { return 6; }
    int nSide() const { return 3; }
    int nSideNode (int /*side*/) const { return 3; }
    int SideNode (int side, int node) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;

    /**
     * \brief Return the values of all shape functions for a point inside
     *   the element, given in local coordinates.
     * \param loc point in local coordinates
     * \return Vector of shape function values u_i(loc) (i=0..5)
     */
    RVector LocalShapeF (const Point &loc) const;

    /**
     * \brief Return the values of all shape functions for a point inside
     *   the element, given in local coordinates.
     * \param [in] loc point in local coordinates
     * \param [out] fun pointer to vector receiving the function values
     * \note fun must point to a valid RVector.
     * \note If the length of fun != 6 on entry, fun is resized to 6.
     */
    void LocalShapeF (const Point &loc, RVector *fun) const;

    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;

    double IntF (int i) const;

    RSymMatrix IntFF() const;

    double IntFF (int i, int j) const;

    double IntFFF (int i, int j, int k) const;

    RSymMatrix IntPFF (const RVector& P) const;

    double IntPFF (int i, int j, const RVector& P) const;

    RVector IntFD (int i, int j) const;

    double IntFDD (int i, int j, int k) const;

    RSymMatrix IntPDD (const RVector& P) const
    { ERROR_UNDEF; return RSymMatrix(); }

    double IntPDD (int i, int j, const RVector &P) const;

    /**
     * \brief %Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. 5)
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
     * \param i first node index (range 0 .. 5)
     * \param j second node index (range 0 .. 5)
     * \param sd side index (range 0 .. 2)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    double SurfIntFF (int i, int j, int sd) const;

    RSymMatrix BndIntPFF (const RVector &P) const
    { ERROR_UNDEF; return RSymMatrix(); }

    double BndIntPFF (int i, int j, const RVector &P) const;

    double IntFd (int i, int j, int k) const;
    // Returns a single element of Int [u_i du_j/dx_k] dr

    double IntPd (const RVector &P, int j, int k) const;
    // Returns Int [p(r) du_j/dx_k] dr where p(r) is defined by its nodal
    // basis coefficients P.

    RSymMatrix Intdd() const;
    // returns matrix of mixed derivatives

    RVector BndIntFCos (int side, const Surface *surf, const RVector &cntcos,
        double sigma, double sup, const NodeList &nlist) const;
    RVector BndIntFCos (int side, const RVector &cntcos, double a,
        const NodeList &nlist) const;
    RVector BndIntFDelta (int side, const Surface *surf, const RVector &pos,
        const NodeList &nlist) const;
    int GetLocalSubsampleAbsc (const Point *&absc) const;
    int Intersection (const Point &p1, const Point &p2, Point *s,
        bool add_endpoints, bool boundary_only);

protected:

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    void ComputeIntFD (const NodeList &nlist);
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    double a0, b0, c0, a1, b1, c1, a2, b2, c2;
    // triangle geometry parameters

    RDenseMatrix intfd_0, intfd_1; // IntFD storage

private:

#ifdef TRI6_STORE_INTFF
    SymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
#endif

};

#endif // !__TRI6_H
