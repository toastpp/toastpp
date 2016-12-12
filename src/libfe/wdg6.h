// -*-C++-*-
// ==========================================================================
// Module libfe
// File wdg6.h
// Declaration of class Wedge6
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_3D
//           ---> Wedge6
// ==========================================================================

#ifndef __WDG6_H
#define __WDG6_H

class FELIB Wedge6: public Element_Unstructured_3D {
public:

    Wedge6 () { Node = new int[nNode()]; };
    Wedge6 (const Wedge6 &el);
    ~Wedge6 () { delete []Node; }

    /**
     * \brief Create a copy of the element and return a pointer to it
     */
    Element *Copy();

    void Initialise (const NodeList& nlist);

    BYTE Type (void) const { return ELID_WDG6; }
    unsigned long GetCaps () const { return 0; }
    int nNode (void) const { return 6; }
    int nSide (void) const { return 5; }
    int nSideNode (int side) const;
    int SideNode (int side, int node) const;

    Point Local (const NodeList& nlist, const Point& glob) const;
    Point NodeLocal (int node) const;
    RVector DirectionCosine (int side, RDenseMatrix& jacin);
    const RVector &LNormal (int side) const;
    bool LContains (const Point& loc, bool pad = true) const;

    RVector LocalShapeF (const Point &loc) const;
    RDenseMatrix LocalShapeD (const Point &loc) const;

    double IntF (int i) const
    { ERROR_UNDEF; return 0; }

    RSymMatrix IntFF () const;
    // Return integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

    double IntFF (int i, int j) const;
    // Return a single element of IntFF

    double IntFFF (int i, int j, int k) const
    {
        ERROR_UNDEF;
	return 0.0; // dummy
    }
    // returns a single element of integral over element of product of three
    // shape functions:
    // IntFFF = Int_el { F_i(r) * F_j(r) * F_k(r) } dr

    void IntFFF (double &iii, double &iij, double &ijk) const
    { ERROR_UNDEF; }
    // returns the values of the FFF tensor for:
    //   all indices equal (iii)
    //   two indices equal (iij)
    //   all indices different (ijk)

    RSymMatrix IntPFF (const RVector& P) const
    {
        ERROR_UNDEF;
	return RSymMatrix(); // dummy
    }
    // Returns integral over element of product of two shape functions and a
    // function P defined on nodes:
    // PFF = Int_el { P(r) F_i(r) F_J(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPFF (int i, int j, const RVector& P) const
    {
        ERROR_UNDEF;
	return 0.0; // dummy
    }
    // Returns a single element of IntPFF

    double IntFDD (int i, int j, int k) const
    {
        ERROR_UNDEF;
	return 0.0; // dummy
    }
    // returns a single element of integral over element:
    // IntFDD = Int_el { F_i(r) * D_j(r) * D_k(r) } dr
    // IntFDD = IntDD/3 ? CHECK THIS IS REALLY CORRECT!

    RSymMatrix IntPDD (const RVector& P) const
    {
        ERROR_UNDEF;
	return RSymMatrix(); // dummy
    }
    // Returns integral over element of product of two shape derivatives and a
    // function P defined in nodal basis:
    // PDD = Int_el { P(r) D_i(r) D_j(r) } dr
    // Vector P contains the complete nodal solution for the mesh

    double IntPDD (int i, int j, const RVector &P) const
    {
        ERROR_UNDEF;
	return 0.0; // dummy
    }
    // Returns a single element of IntPDD

    double SurfIntF (int i, int sd) const
    { ERROR_UNDEF; return 0; }
    
    double SurfIntFF (int i, int j, int sd) const
    { ERROR_UNDEF; return 0; }

    RSymMatrix BndIntPFF (const RVector &P) const
    {
        ERROR_UNDEF;
	return RSymMatrix(); // dummy
    }

    double BndIntPFF (int i, int j, const RVector &P) const
    {
        ERROR_UNDEF;
	return 0.0; // dummy
    }
    // Returns a single element of BndIntPFF

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
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr

    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const
    {
        ERROR_UNDEF;
	return RSymMatrix(); // dummy
    }
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

private:

    int QuadRule (const double **wght, const Point **absc) const;
    // calculates abscissae and weights for numerical quadrature rules on
    // standard element, return number of points

    int BndQuadRule (int, double**, Point**, double*)
    { return 0; }
    // not implemented

#ifdef UNDEF
    double UnitLength (Point& /*loc*/, Matrix& geom, int side);
    // returns the unit length along the specified element side; used for line
    // integrations

    void ConvLineQuadrature (Point** absc, double* labsc, int nqp, int side,
        double* coeff);
    // forms line integration rule along side from equivalent 1d integration
    // rule
#endif

    RSymMatrix ComputeIntFF (const NodeList &nlist) const;
    // Returns integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
};

#endif // !__WDG6_H
