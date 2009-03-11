// -*-C++-*-
// ==========================================================================
// Module libfe
// File element.h
// Declaration of classes Element, Element2D, Element3D
//
// Inheritance:
// ------------
// Element ----> Element2D
//         ----> Element3D
// ==========================================================================

#ifndef __ELEMENT_H
#define __ELEMENT_H

#define NCOEFFS 3	// number of coefficients
#define MUA     0	// absorption coeff id
#define P2      1	// `parameter 2' id
#define MUS     1	// alias for P2: scattering coeff
#define KAPPA   1	// alias for P2: diffusion coeff
#define REFIND  2	// refractive index id
#define C_ALL   3	// `all coefficients'

// Element type identifiers
#define ELID_NONE 0
#define ELID_TRI3OLD 1
#define ELID_RCT4 2
#define ELID_TET4 3
#define ELID_WDG6 4
#define ELID_VOX8 5
#define ELID_TRI6 6
#define ELID_TET10 7
#define ELID_TRI6_IP 8
#define ELID_TRI10 9
#define ELID_TRI10_IP 10
#define ELID_TET10_IP 11
#define ELID_WDG18INF 12
#define ELID_QUAD4 13
#define ELID_PIX4 14
#define ELID_TRI3 15
#define ELID_TRI3D3 16
#define ELID_TRI3D6 17
#define ELID_VOX27 18
#define ELID_LINE2D2 19

// Element capability flags
#define ELCAPS_CURVED_BOUNDARY 1 // element can have curved boundaries
#define ELCAPS_SUBSAMPLING 2     // element implements IntFFF and IntFDD by
                                 // subsampling

class Surface;

// ==========================================================================
// class Element
// ==========================================================================
/**
 * \brief Base class for finite element types.
 */
class FELIB Element {

    friend class ElementList;

public:

    /**
     * \brief Creates a new element with no nodes.
     */
    Element ();

    /**
     * \brief Creates a new element as a copy of 'el'
     */
    Element (const Element& el);

    /**
     * \brief Destroys the element.
     */
    virtual ~Element ();

    /**
     * \brief Element initialisation.
     *
     * Calculation of geometric element parameters, including element size and
     * pre-calculation of geometry dependent element matrices.
     * \param nlist Node list, containing node geometry
     * \note Should be called after the node list has been generated, and
     *   whenever the node list changes.
     */
    virtual void Initialise (const NodeList& nlist);

    /**
     * \brief Element setup after mesh initialisation
     *
     * This method is called after all mesh elements have been initialised.
     * It can be used to perform any element initialisation steps that require
     * the rest of the mesh to be initialised.
     */
    virtual void PostInitialisation () {}

    void operator= (const Element& el);
    // basic element assignment

    virtual BYTE Type (void) const = 0;
    // derived classes return the element type id

    virtual unsigned long GetCaps () const = 0;
    // Return element capability flags

    virtual int Dimension (void) const = 0;
    // derived classes return element dimension (2 or 3)

    virtual int nNode (void) const = 0;
    // derived classes return number of nodes per element

    virtual int nSide (void) const = 0;
    // derived classes return number of sides per element

    virtual int nSideNode (int /*side*/) const = 0;
    // derived classes return number of nodes attached to side 'side'

    virtual int SideNode (int /*side*/, int /*node*/) const = 0;
    // derived classes return number of 'node'-th node of side 'side'

    virtual bool IsNode (int node);
    // returns TRUE if node 'node' belongs to the element

    int IsSide (int nn, int *nd);
    // If all nodes in list 'nd' (of size 'nn') belong to one side of
    // the element, then the side number is returned. Else -1 is returned

    virtual bool IsSideNode (int side, int node);
    // returns TRUE if 'node' belongs to 'side' of the element

    bool IsBoundarySide (int side) {
        RANGE_CHECK (side >= 0 && side < nSide());
	return bndside[side];
    }
    // returns TRUE if side 'side' is a boundary side

    bool HasBoundarySide ()
    { return bndel; }
    // returns TRUE if the element contains a boundary side

    bool HasInterfaceSide ()
    { return interfaceel; }
    // returns TRUE if the element contains an interface side

    virtual Point Local (const NodeList& /*nlist*/, const Point& /*glob*/)
	const = 0;
    // abstract; derived classes return the local coordinate corresponding to
    // 'glob'

    virtual Point NodeLocal (int node) const = 0;
    // derived classes return local coordinates of node 'node'

    virtual Point SurfToLocal (int side, const Point &p) const
    { xERROR(Not implemented); return Point(); }

    virtual void MapToSide (int side, Point &loc) const;
    // map Point 'loc' onto side 'side' by moving it along the side's
    // normal

    virtual Point SideCentre (int side) const;
    // returns the centre of side 'side' in local coordinates

    Point Global (const NodeList &nlist, const Point& loc) const;
    // returns the global coordinate corresponding to local coordinate 'loc'

    virtual RDenseMatrix Elgeom (const NodeList& nlist) const;
    // returns a matrix containing the global coordinates of all element nodes

    virtual RVector DirectionCosine (int/*side*/, RDenseMatrix&/*jacin*/) = 0;
    // abstract; derived classes return a vector of direction cosines in global
    // coordinates of the normal to side 'side'.

    virtual const RVector &LNormal (int /*side*/) const = 0;
    // Normal of side 'side' in local coordinates

    virtual double Size() const = 0;
    // return element size (area for 2-D elements, volume for 3-D elements)

    virtual double SideSize (int sd, const NodeList &nlist) const
    { xERROR (Not implemented); return 0; }
    // return size of side 'sd' of the element

    virtual bool LContains (const Point& /*loc*/, bool pad = true) const = 0;
    // abstract; derived classes return TRUE if point 'loc' (in local
    // coordinates of the element) is within the element
    // if pad==true then points very close to the surface are considered
    // to be within the element

    virtual bool GContains (const Point& /*glob*/, const NodeList& /*nlist*/)
	const;
    // abstract; derived classes return TRUE if point 'glob' (in global
    // coordinates) is within the element. nlist is the mesh node coordinate
    // list
  
    virtual int BndSideList (const NodeList& nlist, int *list);
    // creates a list of boundary side numbers in *list
    // and returns the number of boundary sides owned by the element
    // *list has to be allocated before the call to BndSideList

    virtual RVector LocalShapeF (const Point& /*loc*/) const = 0;
    // vector of shape functions u_i(loc) at point 'loc' given in local
    // element coordinates

    virtual RDenseMatrix LocalShapeD (const Point& /*loc*/) const = 0;
    // matrix of shape function derivatives (d u_i)/(d x_j) at point 'loc'
    // given in local element coordinates

    virtual RVector GlobalShapeF (const NodeList &nlist, const Point &glob)
    const { return LocalShapeF (Local (nlist, glob)); }
    // returns the shape functions for each node at point 'glob' given in
    // global coordinates

    virtual RDenseMatrix GlobalShapeD (const NodeList &nlist,
        const Point &glob) const
    { return LocalShapeD (Local (nlist, glob)); }
    // returns the derivatives of shape functions at point 'glob' given in
    // global coordinates

    virtual int QuadRule (int order, const double **wght, const Point **absc)
        const { xERROR(Not implemented); return 0; }
    // returns weights and abscissae of a quadrature rule over the element
    //
    // order  required for
    // -----------------------
    // 1      IntF
    // 2      IntFF, IntD
    // 3      IntFFF, IntFD
    // 4      IntDD
    // etc. (you get the idea)

    virtual double IntF (int i) const = 0;
    // Return integral of shape function over element

    virtual RSymMatrix IntFF () const = 0;
    // Return integral over element of product of shape functions:
    // FF = Int_el { F_i(r) F_j(r) } dr

    virtual double IntFF (int i, int j) const = 0;
    // Return a single element of IntFF

    virtual double IntFFF (int i, int j, int k) const = 0;
    // returns a single element of integral over element of product of three
    // shape functions:
    // IntFFF = Int_el { F_i(r) * F_j(r) * F_k(r) } dr

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

    virtual void IntFFF (double &iii, double &iij, double &ijk) const = 0;
    // returns the values of the FFF tensor for:
    //   all indices equal (iii)
    //   two indices equal (iij)
    //   all indices different (ijk)
#endif

    virtual RSymMatrix IntPFF (const RVector& P) const = 0;
    // Returns integral over element of product of two shape functions and a
    // function P defined in nodal basis:
    // PFF = Int_el { P(r) F_i(r) F_j(r) } dr
    // Vector P contains the nodal solution for the complete mesh

    virtual double IntPFF (int i, int j, const RVector& P) const = 0;
    // returns a single element of IntPFF

    virtual RSymMatrix IntDD () const = 0;
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr

    virtual double IntDD (int i, int j) const = 0;
    // Returns a single element of IntDD

    virtual RVector IntFD (int i, int j) const
    { xERROR(Not implemented); return RVector(); }
    // Returns a single element of IntFD

    virtual double IntFDD (int i, int j, int k) const = 0;
    // returns a single element of integral over element:
    // IntFDD = Int_el { F_i(r) * D_j(r) * D_k(r) } dr

    virtual RSymMatrix IntPDD (const RVector& P) const = 0;
    // Returns integral over element of product of two shape derivatives and a
    // function P defined in nodal basis:
    // PDD = Int_el { P(r) D_i(r) D_j(r) } dr
    // Vector P contains the nodal solution for the complete mesh

    virtual double IntPDD (int i, int j, const RVector &P) const = 0;
    // Returns a single element of IntPDD

    virtual RSymMatrix BndIntFF () const = 0;
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

    virtual double BndIntFFSide (int i, int j,int sd)
    { xERROR(Not implemented); return 0; }
    // Boundary Integral for ONE SIDE OF ELEMENT ONLY

    virtual double BndIntFF (int i, int j) = 0;
    // returns a single element of BndIntFF

    virtual RSymMatrix BndIntPFF (const RVector &P) const = 0;
    // Returns line integral of product of two shape functions and a function
    // P defined in nodal basis along sides of the element which belong to the
    // mesh boundary:
    // lPFF = Int_bnd { P(r) F_i(r) F_j(r) } ds
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

    virtual double BndIntPFF (int i, int j, const RVector &P) const = 0;
    // Returns a single element of BndIntPFF

    // ===============================================================
    // Integrals including mixed derivatives

    virtual double IntFd (int i, int j, int k) const
    { xERROR(Not implemented); return 0.0; }
    // Returns a single element of Int [u_i du_j/dx_k] dr

    virtual double IntPd (const RVector &P, int j, int k) const
    { xERROR(Not implemented); return 0.0; }
    // Returns Int [p(r) du_j/dx_k] dr where p(r) is defined by its nodal
    // basis coefficients P

    virtual RSymMatrix Intdd() const
    { xERROR(Not implemented); return RSymMatrix(); }
    // returns integral over element of mixed shape function derivatives
    // Int [du_i/dx_j du_k/dx_l] dr
    // The dimension of the returned matrix is nd x nd where n is the number
    // of nodes, and d is the dimension of the element (2 or 3). The matrix
    // contains n x n blocks, with each block of dimension d x d

    virtual double IntFdd (int i, int j, int k, int l, int m) const
    { xERROR(Not implemented); return 0; }
    // Int [u_i du_j/dx_l du_k/dx_m] dr

    virtual double IntPdd (const RVector &p, int j, int k, int l, int m) const
    { xERROR(Not implemented); return 0; }
    // Int f(r) du_j/dx_l du_k/dx_m dr
    // where f(r) is given as a nodal vector

  virtual double IntFfd (int i, int j, int k, int l) const
    { xERROR(Not implemented); return 0; }
    // Int [u_i u_j du_k/dx_l] dr

   virtual double IntPfd (const RVector &p, int j, int k, int l) const
    { xERROR(Not implemented here); return 0; }
    // Int f(r) u_j du_k/dx_l dr
    // where f(r) is given as a nodal vector

    //------------- integrate functions on unit sphere --------------
    virtual double IntUnitSphereP (const NodeList& nlist, const RVector& P) const
    { xERROR(Not implemented); return 0.0; }
    // Returns integral of f(r) over unitsphere patch 
    // only implemented for Tri3D3, Tri3D6

     virtual double IntUnitSphereFF (const NodeList& nlist, const int i, const int j) const
    { xERROR(Not implemented); return 0.0; }
    // Returns integral of u_i(r)*u_j(r) over unitsphere patch 
    // only implemented for Tri3D3, Tri3D6

     virtual double IntUnitSpherePFF (const NodeList& nlist, const int i, const int j, const RVector& P) const
    { xERROR(Not implemented); return 0.0; }
    // Returns integral of f(r)* u_i(r)*u_j(r) over unitsphere patch 
    // only implemented for Tri3D3, Tri3D6

    //---------------
    virtual RVector BndIntFX (int side, double (*func)(const Point&),
        const NodeList &nlist) const;
    // Calculates Int [u_i(r) func(r) dr]
    // along side 'side' of the triangle
    // where func is a scalar-valued user supplied function evaluated
    // at a global coordinate r

    virtual RVector BndIntFCos (int side, const Surface *surf,
        const RVector &cntcos, double sigma, double sup,
        const NodeList &nlist) const;
    // Given a surface 'surf' and a cosine profile defined by its
    // parametrised centre 'cntcos', scale 'sigma' and support radius 'sup'
    // return the boundary integral of u_i(xi) * cos(xi) over side sd for all
    // associated nodes i

    virtual RVector BndIntFCos (int side, const RVector &cntcos, double a,
        const NodeList &nlist) const
    { xERROR(Not implemented); return RVector(); }
    // Calculate integral of product of cosine (centered at 'cntcos' and
    // support radius a) with shape functions along 'side'

    virtual RVector BndIntFDelta (int side, const Surface *surf,
        const RVector &pos, const NodeList &nlist) const;
    // Given a surface 'surf' and delta function defined by its parametrised
    // position 'pos', return boundary integral of u_i(xi) * delta(pos-xi)
    // over side 'side' for all associated nodes i

    int GetSubsampleFD (int &n, double *&wght, Point *&absc,
        RVector *&F, RDenseMatrix *&D, const NodeList &nlist) const;
    // returns the abscissae and weights for element subsampling, and
    // the values of the shape functions and derivatives at those points
    // 'n' is the size of the arrays on input (may be changed on output
    // if size was insufficient). return value is the number of subsamples

    int GetBndSubsampleFD (int side, int &n, double *&wght, Point *&absc,
        RVector *&F, RDenseMatrix *&D, const NodeList &nlist) const;
    // returns the abscissae (in global coordinates) and weights for element
    // boundary subsampling along 'side', and the values of the shape functions
    // and derivatives at those points. 'n' is the size of the arrays on input
    // (may be changed on output if size was insufficient). return value is the
    // number of subsamples

    virtual int GetLocalSubsampleAbsc (const Point *&absc) const
    { return 0; }
    // returns abscissae in local coords for numerical integration by uniform
    // subsampling. Return value is the number of samples returned in 'absc'

    virtual int GetBndSubsampleAbsc (int side, const Point *&absc) const
    { return 0; }
    // returns abscissae for numerical integration by uniform
    // subsampling over boundary side in local coordinates of the surface
    // element. Note that abscissae are of dimension dim-1. To convert to
    // local element coordinates, use SurfToLocal.
    // Return value is the number of samples returned in 'absc'

    virtual RDenseMatrix StrainDisplacementMatrix (const Point &glob) const
    { xERROR(Not implemented); return RDenseMatrix(); }
    // Returns the strain displacement matrix of the element at global point
    // 'glob'
    // In 3D, the matrix is of dimension 6 x 3n (n: number of nodes):
    // B = [B1, B2, ... Bn] with
    //
    //      |dNi/dx     0       0  |
    //      |   0    dNi/dy     0  |
    // Bi = |   0       0    dNi/dz|
    //      |dNi/dy  dNi/dx     0  |
    //      |   0    dNi/dz  dNi/dy|
    //      |dNi/dz     0    dNi/dx|

    virtual RDenseMatrix ElasticityStiffnessMatrix (const RDenseMatrix &D)
        const
    { xERROR(Not implemented); return RDenseMatrix(); }
    // Returns the element elasticity stiffness matrix K, given
    // strain matrix D: K = \int_el B^T D B dr
    // where B is strain displacement matrix

    virtual RDenseMatrix ElasticityStiffnessMatrix (double E,
        double nu) const
    { xERROR(Not implemented); return RDenseMatrix(); }
    // Returns the element elasticity stiffness matrix K for isotropic case,
    // given material properties E (Young's modulus) and nu (Poisson's ratio)
    // Strain matrix is calculated as
    //     |d1 d2 d2 0  0  0 |
    //     |d2 d1 d2 0  0  0 |          d1 = E(1-nu)/((1+nu)(1-2nu))
    // D = |d2 d2 d1 0  0  0 |   with   d2 = E nu/((1+nu)(1-2nu))
    //     |0  0  0  d3 0  0 |          d3 = E/(2(1+nu))
    //     |0  0  0  0  d3 0 |
    //     |0  0  0  0  0  d3|

    RDenseMatrix ElasticStrainDisplacement (const RVector &loc,
        const RDenseMatrix &gder) const;
    // Returns strain-displacement matrix for 2D plane elasticity
    // at local position loc, given global shape function derivatives gder
    // See NAG finel routine B2C2
    // only implemented for 2D yet

    virtual RDenseMatrix IsotropicElasticityMatrix (double E, double nu) const
    { xERROR(Not implemented); return RDenseMatrix(); }
    // returns a 6x6 symmetric elasticity matrix for isotropic materials,
    // given Young's modulus E and Poisson's ratio nu

    virtual RVector InitialStrainVector (double E, double nu, const RVector &e0)
    { xERROR(Not implemented); return RVector(); }
    // returns a element strain vector (rhs) for initial nodal strain given
    // in e0

    virtual RVector ThermalExpansionVector (double E, double nu,
        double alpha, double dT)
    { xERROR(Not implemented); return RVector(); }
    // returns element strain vector (rhs) given element thermal expansion
    // coefficient alpha and temperature change dT

    virtual RVector DThermalExpansionVector (double E, double nu)
    { xERROR(Not implemented); return RVector(); }
    // returns derivative of thermal expansion vector with respect to
    // expansion coefficient (assuming temperature change 1)

    virtual int GlobalIntersection (const NodeList &nlist,
	const Point &p1, const Point &p2, Point **list) = 0;
    // same as 'Intersection' but p1 and p2 given in global coords
    // The return list however is still in local coords

    virtual int Intersection (const Point &/*p1*/, const Point &/*p2*/,
        Point** /*pi*/) = 0;
    // abstract; derived classes create a list of points where the line defined
    // by p1 and p2 intersects the element (in local coordinates) or starts/ends
    // within the element. Returns the length of the list

    virtual RDenseMatrix LocaltoGlobalMat () const 
    { xERROR(Not implemented); return RDenseMatrix(); }
    // abstract - local to global matrix

    virtual RDenseMatrix GlobaltoLocalMat () const 
    { xERROR(Not implemented); return RDenseMatrix(); }
    // abstract - global to local matrix

    virtual RDenseMatrix FTAMat () const 
    { xERROR(Not implemented); return RDenseMatrix(); }
    // abstract - transform k-space matrix

    friend std::istream& operator>> (std::istream& i, Element &el);
    friend std::ostream& operator<< (std::ostream& o, const Element &el);
    // standard method for reading/writing an element from/to a stream

    int *Node;
    // list of absolute node numbers for the element's vertices

    //int Region;
    int Region() const {return region;} // region id for segmented meshes
    void SetRegion(int nr) {region = nr;} // set region id 

protected:

    int *sdnbhr;
    // List of side neighbours for the element. A value of -1 indicates no
    // neighbour, i.e. a boundary side. A value of < -1 indicates an invalid
    // (not yet computed) entry.

    bool *bndside;
    // boundary flag for each side of the element. A side is a boundary
    // side if all its nodes are boundary nodes

    bool bndel;
    // boundary flag for the element. The element is a boundary element if
    // it contains one or more boundary sides

    bool interfaceel;
    // interface flag for the element. The element is an interface element if
    // it contains one or more interface sides

    int region;         // region number (-1 = no region)
};


// ==========================================================================
// class Element_Structured
// base class for regular element types (to build structured grids)
//
// Inheritance:
// ------------
// Element
// ---> Element_Structured
// ==========================================================================

class FELIB Element_Structured: public Element {
public:
    Element_Structured (): Element () {}
    Element_Structured (const Element_Structured &el): Element (el) {}
    virtual ~Element_Structured () {}
};


// ==========================================================================
// class Element_Unstructured
// base class for unstructured element types
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
// ==========================================================================

class FELIB Element_Unstructured: public Element {
public:
    Element_Unstructured (): Element () {}
    Element_Unstructured (const Element_Unstructured &el): Element (el) {}
    virtual ~Element_Unstructured () {}
    virtual void Initialise (const NodeList& nlist);
    inline double Size() const { return size; }
    void operator= (const Element_Unstructured& el);

    virtual bool GContains (const Point& glob, const NodeList& nlist)
	const;

    inline RSymMatrix IntDD () const { return intdd; }
    inline double IntDD (int i, int j) const { return intdd.Get(i,j); }
    inline RSymMatrix BndIntFF () const { return intbff; }
    inline double BndIntFF (int i, int j) { return intbff.Get(i,j); }

protected:
    virtual double ComputeSize (const NodeList&) const = 0;
    // computes element size
    // called by `Initialise' to set `size'

    virtual RSymMatrix ComputeIntDD (const NodeList &nlist) const = 0;
    // Returns integral over element of product of shape derivatives:
    // DD = Int_el { D_i(r) D_j(r) } dr
    // called by `Initialise' to set `intdd'


    virtual RSymMatrix ComputeBndIntFF (const NodeList& nlist) const = 0;
    // Returns line integral of product of two shape functions along sides of
    // the element which belong to the mesh boundary
    // If the element does not contain any boundary sides then the returned
    // matrix is empty.

    double size;
    // stores element size (area in 2D, volume in 3D)
    // set by Initialise

    RSymMatrix intdd;
    // stores integral over element of product of shape derivatives
    // Int_el { D_i(r) D_j(r) } dr
    // set by Initialise

    RSymMatrix intbff;
    // stores surface integral over boundary sides of the element of shape
    // functions: Int_bnd { F_i(r) F_j(r) } ds

    Point bbmin, bbmax;
    // bounding box of element in global coordinates
    // set by Initialise
};


// ==========================================================================
// class Element_Structured_2D
//
// Inheritance:
// ------------
// Element
// ---> Element_Structured
//      ---> Element_Structured_2D
// ==========================================================================

class FELIB Element_Structured_2D: public Element_Structured {
public:
    Element_Structured_2D (): Element_Structured () {}
    Element_Structured_2D (const Element_Structured_2D &el)
      : Element_Structured (el) {}
    virtual ~Element_Structured_2D () {}
    virtual int Dimension (void) const { return 2; }
};


// ==========================================================================
// class Element_Structured_3D
//
// Inheritance:
// ------------
// Element
// ---> Element_Structured
//      ---> Element_Structured_3D
// ==========================================================================

class FELIB Element_Structured_3D: public Element_Structured {
public:
    Element_Structured_3D (): Element_Structured () {}
    Element_Structured_3D (const Element_Structured_3D& el)
      : Element_Structured (el) {}
    virtual ~Element_Structured_3D () {}
    virtual int Dimension (void) const { return 3; }
};


// ==========================================================================
// class Element_Unstructured_2D
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_2D
// ==========================================================================

class FELIB Element_Unstructured_2D: public Element_Unstructured {
public:
    Element_Unstructured_2D (): Element_Unstructured () {}
    Element_Unstructured_2D (const Element_Unstructured_2D &el)
      : Element_Unstructured (el) {}
    virtual ~Element_Unstructured_2D () {}
    virtual int Dimension (void) const { return 2; }
};


// ==========================================================================
// class Element_Unstructured_3D
//
// Inheritance:
// ------------
// Element
// ---> Element_Unstructured
//      ---> Element_Unstructured_3D
// ==========================================================================

class FELIB Element_Unstructured_3D: public Element_Unstructured {
public:
    Element_Unstructured_3D (): Element_Unstructured () {}
    Element_Unstructured_3D (const Element_Unstructured_3D& el)
      : Element_Unstructured (el) {}
    virtual ~Element_Unstructured_3D () {}
    virtual int Dimension (void) const { return 3; }

    RDenseMatrix IsotropicElasticityMatrix (double E, double nu) const;
    // returns a 6x6 symmetric elasticity matrix for isotropic materials,
    // given Young's modulus E and Poisson's ratio nu
};

typedef Element *PElement;

// ==========================================================================
// Nonmember functions
// ==========================================================================


#endif // !__ELEMENT_H
