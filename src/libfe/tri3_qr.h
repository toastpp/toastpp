// -*-C++-*-
// ==========================================================================
// Module libfe
// File tri3.h
// Declaration of class Triangle3_qr
// An alternative implementation which uses Gaussian quadrature
//
// Inheritance:
// ------------
// Element ----> Element2D ----> Triangle3
// ==========================================================================

#ifndef __TRI3_QR_H
#define __TRI3_QR_H

#define TRI3_STORE_COORDS

class Triangle3_qr: public Element2D {
public:

    Triangle3_qr () { Node = new int[nNode()]; }
    Triangle3_qr (const Triangle3_qr& el);
    ~Triangle3_qr () { delete []Node; }

    void Initialise (const NodeList &nlist);

    BYTE Type () const { return ELID_TRI3QR; }
    // returns element type id

    int nNode() const { return 3; }
    // returns number of nodes per element

    int nSide () const { return 3; };
    // number of sides per element

    int nSideNode (int /*side*/) const { return 2; };
    // returns number of nodes attached to side 'side'

    int SideNode (int side, int node) const;
    // returns number of 'node'-th node of side 'side'

    const Matrix &Jacobian() const { return jac; }
    // Return the element's 2x2 Jacobian matrix required for mapping
    // between local and global coordinates

    double DetJ() const { return djac; }
    // Determinant of the Jacobian

    Point Local (const NodeList& nlist, const Point& glob) const;
    // returns the local coordinate corresponding to global coordinate 'glob'

    Point NodeLocal (int node) const;
    // returns local coordinates of node 'node'

    RVector DirectionCosine (int side, Matrix& jacin);
    // returns a vector of direction cosines in global coordinates of the
    // normal to side 'side'.

// ##########################################################################
// requires update from here!

    bool LContains (const Point& loc) const;
    // returns TRUE if point 'loc' (in local coordinates of the element) is
    // within the standard triangle

    bool GContains (const Point& glob, const NodeList& nlist) const;
    // returns TRUE if point 'glob' (in global coordinates) is within the
    // element. nlist is the mesh node coordinate list

    RVector LocalShapeF (const Point &loc) const;
    // vector of shape functions u_i(loc) at point 'loc' given in local
    // element coordinates

    Matrix LocalShapeD (const Point &loc) const;
    // 2x3 matrix of shape function derivatives (d u_j)/(d x_i) at point 'loc'
    // given in local element coordinates

    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    // returns the shape functions for each node at a point, given in global
    // coordinates

    Matrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;
    // returns derivatives of shape functions in global coordinates at point
    // `glob', given in global coordinates

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

    SymMatrix IntPFF (const RVector& P) const;
    // Returns integral over element of product of two shape functions and a
    // function P defined on nodes:
    // PFF = Int_el { P(r) F_i(r) F_J(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPFF (int i, int j, const RVector& P) const;
    // Returns a single element of IntPFF

    double IntFDD (int i, int j, int k) const
    { return (j >= k ? intdd[j][k] : intdd[k][j]) * 0.3333333333; }
    // returns a single element of integral over element:
    // IntFDD = Int_el { F_i(r) * D_j(r) * D_k(r) } dr

    SymMatrix IntPDD (const RVector& P) const
    { return intdd * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns integral over element of product of two shape derivatives and a
    // function P defined in nodal basis:
    // PDD = Int_el { P(r) D_i(r) D_j(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPDD (int i, int j, const RVector &P) const {
	RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
	return (i >= j ? intdd[i][j] : intdd[j][i]) * 
	    ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0);
    }
    // Returns a single element of IntPDD

    SymMatrix BndIntPFF (const RVector &P) const
	{ return intbff * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // This is a hack! Really P would have to be taken into the integral,
    // rather than just averaged.

    double BndIntPFF (int i, int j, const RVector &P) const {
	RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
	return (i >= j ? intbff[i][j] : intbff[j][i])
	    * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0);
    }
    // Returns a single element of BndIntPFF

  // SymMatrix FTF_bnd (const Matrix &geom, int side);
    // returns line integral of Fi(T) Fj along side `side' (0-2)

  // int QuadRule (const double **wght, const Point **absc) const;
    // calculates abscissae and weights for numerical quadrature rules on
    // standard element, return number of points

  // int BndQuadRule (int side, double** wght, Point** absc, double* coeff);
    // returns weights and abscissae of a line quadrature rule along `side'
    // `coeff' is a scaling factor. Return value is number of points

  // double UnitLength (Point& /*loc*/, Matrix& geom, int side);
    // returns the unit length along the specified element side; used for line
    // integrations

  // void ConvLineQuadrature (Point** absc, double* labsc, int nqp, int side,
  //      double* coeff);
    // forms line integration rule along side from equivalent 1d integration
    // rule

    int Intersection (const Point &p1, const Point &p2, Point** pi);
    // creates a list of points where the line defined by p1 and p2 intersects
    // the element (in local coordinates) or starts/ends within the element.
    // Returns the length of the list

protected:

    double ComputeSize (const NodeList &nlist) const;
    // area of triangle in global coordinates

    SymMatrix ComputeIntFF (const NodeList &nlist) const;
    // Returns integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

    SymMatrix ComputeIntDD (const NodeList &nlist) const;
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr

    SymMatrix ComputeBndIntFF (const NodeList& nlist) const;
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

private:
    SquareMatrix jac;
    // 2x2 matrix of derivatives of global coordinates with respect to local

    double djac;
    // determinant of jac

    double intfff_scale[3];
    // factors for calculating IntFFF(i,j,k)

#ifdef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
#endif
};

#endif // !__TRI3_QR_H
