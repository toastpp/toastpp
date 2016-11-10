// -*-C++-*-
// ==========================================================================
// Module libfe
// File tri3D3.h
// Declaration of unstructured 2-D element class Triangle3D3
// (3-noded triangle, 1st order shape functions, nodes are three dimensional)
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
//               ---> Triangle3D3
// ==========================================================================

#ifndef __TRI3D3_H
#define __TRI3D3_H

#include "toastdef.h"

class Surface;

class FELIB Triangle3D3: public Triangle3 {
public:

    Triangle3D3();
    Triangle3D3 (const Triangle3D3& el);
    ~Triangle3D3();
    // constructors, destructor

    void Initialise (const NodeList& nlist);

    BYTE Type () const { return ELID_TRI3D3; }
    // returns element type id

    int Dimension () const { return 3; }

  //    unsigned long GetCaps () const { return ELCAPS_SUBSAMPLING; }
    // Return element capability flags

  //    int nNode() const { return 3; }
    // returns number of nodes per element

  //    int nSide () const { return 3; };
    // number of sides per element

  //   int nSideNode (int /*side*/) const { return 2; };
    // returns number of nodes attached to side 'side'

  //  int SideNode (int side, int node) const;
    // returns number of 'node'-th node of side 'side'

  //  double SideSize (int sd, const NodeList &nlist) const;
    // returns length of side 'sd'

    Point Local (const NodeList& nlist, const Point& glob) const;
    // returns the local coordinate corresponding to global coordinate 'glob'

  //    Point NodeLocal (int node) const;
    // returns local coordinates of node 'node'

     Point SurfToLocal (int side, const Point &p) const;
    // Map 1-D surface point p (0-1) along side 'side' into 2-D local coords

    RVector DirectionCosine (int side, RDenseMatrix &jacin);
    // returns a vector of direction cosines in global coordinates of the
    // normal to side 'side'.

  //  const RVector &LNormal (int side) const;

  //    bool LContains (const Point& loc, bool pad = true) const;
    // returns TRUE if point 'loc' (in local coordinates of the element) is
    // within the standard triangle

    bool GContains (const Point& glob, const NodeList& nlist) const;
    // returns TRUE if point 'glob' (in global coordinates) is within the
    // element. nlist is the mesh node coordinate list

  //    RVector LocalShapeF (const Point &loc) const;
    // vector of shape functions u_i(loc) at point 'loc' given in local
    // element coordinates

  //    RDenseMatrix LocalShapeD (const Point &loc) const;
    // 2x3 matrix of shape function derivatives (d u_j)/(d x_i) at point 'loc'
    // given in local element coordinates

    RVector GlobalShapeF (const NodeList& nlist, const Point& glob) const;
    // returns the shape functions for each node at a point, given in global
    // coordinates

    RDenseMatrix GlobalShapeD (const NodeList& nlist, const Point& glob) const;
    // returns derivatives of shape functions in global coordinates at point
    // `glob', given in global coordinates

  //  double IntF (int i) const;

  //  RSymMatrix IntFF () const;
    // Return integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

  //    int QuadRule (int order, const double **wght, const Point **absc) const;

  //  double IntFF (int i, int j) const;
    // Return a single element of IntFF

  //  double IntFFF (int i, int j, int k) const;
    // returns a single element of integral over element of product of three
    // shape functions:
    // IntFFF = Int_el { F_i(r) * F_j(r) * F_k(r) } dr

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

  //    void IntFFF (double &iii, double &iij, double &ijk) const;
    // returns the values of the FFF tensor for:
    //   all indices equal (iii)
    //   two indices equal (iij)
    //   all indices different (ijk)
#endif

  //    RSymMatrix IntPFF (const RVector& P) const;
    // Returns integral over element of product of two shape functions and a
    // function P defined on nodes:
    // PFF = Int_el { P(r) F_i(r) F_J(r) } dr
    // Vector P contains the complete nodal solution for the mesh

  //  double IntPFF (int i, int j, const RVector& P) const;
    // Returns a single element of IntPFF

  //    RVector IntFD (int i, int j) const;
    // Returns single (vector) element of FD integral over element:
    // IntFD = Int_el { F_i(r) D_j(r) } dr

  //   double IntFDD (int i, int j, int k) const
  //  { return intdd.Get(j,k) * 0.3333333333; }
    // returns a single element of integral over element:
    // IntFDD = Int_el { F_i(r) * D_j(r) * D_k(r) } dr

  //    RSymMatrix IntPDD (const RVector& P) const
  //  { return intdd * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns integral over element of product of two shape derivatives and a
    // function P defined in nodal basis:
    // PDD = Int_el { P(r) D_i(r) D_j(r) } dr
    // Vector P contains the complete nodal solution for the mesh

  //  double IntPDD (int i, int j, const RVector &P) const
  //  { return intdd.Get(i,j) * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns a single element of IntPDD

  //    RSymMatrix BndIntPFF (const RVector &P) const
  //	{ return intbff * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // This is a hack! Really P would have to be taken into the integral,
    // rather than just averaged.

  //    double BndIntPFF (int i, int j, const RVector &P) const
  //  { return intbff.Get(i,j) * ((P[Node[0]]+P[Node[1]]+P[Node[2]])/3.0); }
    // Returns a single element of BndIntPFF

  //  double IntFd (int i, int j, int k) const;
    // Returns a single element of Int [u_i du_j/dx_k] dr

  //   double IntPd (const RVector &P, int j, int k) const;
    // Returns Int [p(r) du_j/dx_k] dr where p(r) is defined by its nodal
    // basis coefficients P.

  //    double IntFdd (int i, int j, int k, int l, int m) const;
    // Int u_i du_j/dx_l du_k/dx_m dr

  //    double IntPdd (const RVector &p, int j, int k, int l, int m) const;
    // Int f(r) du_j/dx_l du_k/dx_m dr
    // where f(r) is given as a nodal vector

  //    RSymMatrix Intdd() const;
    // returns matrix of mixed derivatives

  //    double SurfIntFF (int i, int j, int sd);
    // Int [u_i u_j] dr along side sd

  //   RVector BndIntFX (int side, double (*func)(const Point&),
  //      const NodeList &nlist) const;
    // See element.h for description

//    RVector BndIntFCos (int side, const Surface *surf, const RVector &cntcos,
//        double sigma, double sup, const NodeList &nlist) const;
    // See element.h for description

//    RVector BndIntFCos (int side, const RVector &cntcos, double a,
//        const NodeList &nlist) const;
    // See element.h for description

  
//    RVector BndIntFDelta (int side, const Surface *surf, const RVector &pos,
//        const NodeList &nlist) const;
    // See element.h for description

  //    int GetLocalSubsampleAbsc (const Point *&absc) const;
    // returns abscissae in local coords for numerical integration by uniform
    // subsampling. Return value is the number of samples returned in 'absc'

  //    int GetBndSubsampleAbsc (int side, const Point *&absc) const;
    // abscissae for numerical integration over boundary side in local
    // coordinates of boundary element (dim-1). Use SurfToLocal to convert
    // into local element coordinates

  //    int Intersection (const Point &p1, const Point &p2, Point** pi);
    // creates a list of points where the line defined by p1 and p2 intersects
    // the element (in local coordinates) or starts/ends within the element.
    // Returns the length of the list
  RVector LocalShapeQF (const Point &loc) const;
  RDenseMatrix LocalShapeQD (const Point &loc) const;
  double IntUnitSphereP (const NodeList& nlist, const RVector& P) const;
  double IntUnitSpherePFF(const NodeList& nlist, const int i, const int j, const RVector& P) const;
  double IntUnitSphereFF (const NodeList& nlist, const int i, const int j) const;

  RDenseMatrix LocaltoGlobalMat () const 
  {    return MapLocaltoGlobal;  }

  RDenseMatrix GlobaltoLocalMat () const 
  {    return inverse(MapLocaltoGlobal);  }

  RDenseMatrix FTAMat () const ;


protected:

  //    double ComputeSize (const NodeList &nlist) const;
    // area of triangle in global coordinates

  //    RSymMatrix ComputeIntDD (const NodeList &nlist) const;
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr

  //    void ComputeIntFD (const NodeList &nlist);
    // Calculates the FD integrals over the element:
    // IntFD = Int_el { F_i(r) D_j(r) } dr

  //    RSymMatrix ComputeBndIntFF (const NodeList& nlist) const;
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

private:

  //    double *intfd_0, *intfd_1;
    // IntFD matrix storage

  //    double a0, b0, c0, a1, b1, c1, a2, b2, c2;
    // triangle geometry parameters

#ifdef TRI3_STORE_INTFF
  //    RSymMatrix intff;
    // stores integral over element of product of shape functions
    // Int_el { F_i(r) F_j(r) } dr
    // set by Initialise
#endif
      RDenseMatrix Map2Dto3D;
      RDenseMatrix Map3Dto2D;
      RDenseMatrix MapLocaltoGlobal;
      RVector *V2D;
      RVector uspts[6]; // used by unit sphere integration
};

#endif //!__TRI3D3_H
