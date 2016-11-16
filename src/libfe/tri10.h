// -*-C++-*-
// =========================================================================
// TOAST v.15                                      (c) Martin Schweiger 1999
// Library: libfe     File: tri10.h
//
// 10-noded triangle for cubic shape functions and straight boundaries,
// using analytic integration methods.
//
//            ^ eta               Side 0: eta    = 0
//            |                   Side 1: xi+eta = 1
// -------                        Side 2: xi     = 0
// Local     N2
// element    +                   Node coordinates:
// -------    | \                 N0 = (0,0)    N5 = (2/3,1/3)
//            |   \               N1 = (1,0)    N6 = (1/3,2/3)
//          N7+     +N6           N2 = (0,1)    N7 = (0,2/3)
//            |       \           N3 = (1/3,0)  N8 = (0,1/3)
//            |         \         N4 = (2/3,0)  N9 = (1/3,1/3)
//          N8+     +     +N5
//            |     N9      \
//            |               \
//          N0+-----+-----+-----+N1 --> xi
//                 N3     N4
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
//           ---> Triangle10
// =========================================================================

#ifndef __TRI10_H
#define __TRI10_H

#include "toastdef.h"

/**
 * \brief A 10-noded 2-dimensional triangle element with straight edges
 *   and 3rd order shape functions.
 *
 * The node arrangement is given by
 * \code
 *   N2             The local element has node coordinates
 *    +                     N0 = (0,0)    N5 = (2/3,1/3)
 *    |\                    N1 = (1,0)    N6 = (1/3,2/3)
 *    | \                   N2 = (0,1)    N7 = (0,2/3)
 *  N7+  +N6                N3 = (1/3,0)  N8 = (0,1/3)
 *    |   \                 N4 = (2/3,0)  N9 = (1/3,1/3)
 *    |    \        
 *  N8+  +  +N5     Sides:  side 0 contains N0,N1,N3,N4 (y=0)
 *    |  N9  \              side 1 contains N1,N2,N5,N6 (x+y=1)
 *    |       \             side 2 contains N2,N0,N7,N8 (x=0)
 *  N0+--+--+--+N1           
 *      N3  N4              (N9 is internal)
 * \endcode
 */
class FELIB Triangle10: public Element_Unstructured_2D {
public:

    Triangle10 () { Node = new int [nNode()]; }
    Triangle10 (const Triangle10 &el);
    ~Triangle10 () { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList &nlist);

    BYTE Type () const { return ELID_TRI10; }
    unsigned long GetCaps () const { return 0; }
    int nNode () const { return 10; }
    int nSide () const { return 3; }

    int nSideNode (int /*side*/) const { return 4; }
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
     * \return Vector of shape function values u_i(loc) (i=0..9)
     */
    RVector LocalShapeF (const Point &loc) const;

    /**
     * \brief Return the values of all shape functions for a point inside
     *   the element, given in local coordinates.
     * \param [in] loc point in local coordinates
     * \param [out] fun pointer to vector receiving the function values
     * \note fun must point to a valid RVector.
     * \note If the length of fun != 10 on entry, fun is resized to 10.
     */
    void LocalShapeF (const Point &loc, RVector *fun) const;

    RDenseMatrix LocalShapeD (const Point &loc) const;
    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;

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
     * \brief %Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. 9)
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
     * \param i first node index (range 0 .. 9)
     * \param j second node index (range 0 .. 9)
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
    int Intersection (const Point &p1, const Point &p2, Point *s,
	bool add_endpoints, bool boundary_only)
    { ERROR_UNDEF; return 0; }


protected:

  //double Jacobian (const Point &loc, Matrix &J) const;
    // Set J to the Jacobian of the element at local point loc (J must be of
    // dimension 2x2 on input). Return value is det J

  //double IJacobian (const Point &loc, Matrix &IJ) const;
    // Set IJ to the inverse of the Jacobian at local point loc (IJ must be
    // of dimension 2x2 on input). Return value ist det J

  //double DetJ (const Point &loc, const NodeList *nlist = 0) const;
    // Return value of det J at local point loc
    // nlist is required ifndef TRI10IP_STORE_COORDS

    double ComputeSize (const NodeList &nlist) const;
    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    void ComputeIntFD (const NodeList &nlist);
    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;

    int GetLocalSubsampleAbsc (const Point *&absc) const;

    double a0, b0, c0, a1, b1, c1, a2, b2, c2;
    RDenseMatrix intfd_0, intfd_1; // IntFD storage

#ifdef TRI10_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y;
    double n5x, n5y, n6x, n6y, n7x, n7y, n8x, n8y, n9x, n9y;
    // Global node coordinates
#endif

private:

#ifdef TRI10_STORE_INTFF
    SymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
#endif

};

#endif // !__TRI10_H
