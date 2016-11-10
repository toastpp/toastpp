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

/**
 * \defgroup eltp Element type identifiers
 * Numerical identifiers to distinguish between FEM element types.
 * \sa Element::Type
 */
//@{
#define ELID_NONE 0        ///< undefined element type
#define ELID_TRI3OLD 1     ///< old-style 3-noded triangle
#define ELID_RCT4 2        ///< old-style 3-noded pixel
#define ELID_TET4 3        ///< 4-noded tetrahedron
#define ELID_WDG6 4        ///< 6-noded wedge
#define ELID_VOX8 5        ///< 8-noded voxel
#define ELID_TRI6 6        ///< 6-noded triangle
#define ELID_TET10 7       ///< 10-noded tetrahedron
#define ELID_TRI6_IP 8     ///< 6-noded isoparametric triangle
#define ELID_TRI10 9       ///< 10-noded triangle
#define ELID_TRI10_IP 10   ///< 10-noded isoparametric triangle
#define ELID_TET10_IP 11   ///< 10-noded isoparametric tetrahedron
#define ELID_WDG18INF 12   ///< 18-noded infinite wedge element
#define ELID_QUAD4 13      ///< 4-noded quadrilateral
#define ELID_PIX4 14       ///< 4-noded regular pixel
#define ELID_TRI3 15       ///< 3-noded triangle
#define ELID_TRI3D3 16     ///< 3-noded surface triangle
#define ELID_TRI3D6 17     ///< 6-noded surface triangle
#define ELID_VOX27 18      ///< 27-noded voxel
#define ELID_LINE2D2 19    ///< 2-noded line
//@}

/**
 * \defgroup elcap Element capability flags
 * Bitflags for characterising element properties.
 * \sa Element::GetCaps
 */
//@{
/// element can have curved boundaries
#define ELCAPS_CURVED_BOUNDARY 0x1
/// element implements IntFFF and IntFDD by subsampling
#define ELCAPS_SUBSAMPLING     0x2
//@}

class Surface;
class Element;
class Mesh;

struct ElementSubdivisionData {
    int level;             // subdivision level (>= 0)
    Element *sibling;      // sibling neighbour for dual splits
    bool is_sibling0;      // distinguish the two siblings
};

// ==========================================================================
// class Element
// ==========================================================================
/**
 * \brief Base class for finite element types.
 *
 * This class defines the common interface for all Toast FEM elements, in
 * particular a set of integrals of combinations of shape function,
 * shape function derivatives and nodal functions over the element (Int*).
 */
class FELIB Element {

    friend class ElementList;
    friend class Mesh;

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
     * \brief Create a copy of the element and return a pointer to it
     */
    virtual Element *Copy() = 0;

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
    virtual void PostInitialisation (const NodeList& nlist) {}

    /**
     * \brief %Element assignment.
     * \param el right-hand side element operand
     * \note Allows operations of the form
     *   \code el1 = el2;
     *   \endcode
     * \note This operation is only allowed if both elements are of the same
     *   type (e.g. 4-noded tetrahedra, etc.) It is not possible to change the
     *   type of an element.
     * \note This operation copies the node indices from \e el to \e *this.
     */
    void operator= (const Element& el);
    // basic element assignment

    /**
     * \brief Returns an element type identifier.
     * \return Element type identifier (see \ref eltp).
     */
    virtual BYTE Type () const = 0;

    /**
     * \brief Returns the VTK element type identifier, or 0 if the
     *   element doesn't have a VTK representation.
     */
    virtual BYTE VtkType () const { return 0; }

    /**
     * \brief Returns element capability flags.
     * \return Bitflags for element caps. (see \ref elcap)
     */
    virtual unsigned long GetCaps () const = 0;

    /**
     * \brief Returns the spatial dimension of the element
     * \return 2 for 2-D elements (triangles, pixels etc.), 3 for 3-D elements
     *   (tetrahedra, voxels, etc.)
     */
    virtual int Dimension () const = 0;

    /**
     * \brief Returns the number of nodes associated with the element.
     * \return Number of nodes
     */
    virtual int nNode () const = 0;

    /**
     * \brief Returns the number of element sides
     * \return Number of sides
     * \note For 3-D elements, this returns the number of element faces (4
     *   for tetrahedra, 6, for voxels, etc). For 2-D elements, it returns
     *   the number of boundary segments (3 for triangles, 4 for pixels, etc.)
     */
    virtual int nSide () const = 0;

    /**
     * \brief Returns the number of vertices associated with a side.
     * \param side side index (>= 0)
     * \return Number of vertices connected to the side.
     * \sa nSide, nNode, nSideNode
     */
    virtual int nSideNode (int side) const = 0;

    /**
     * \brief Returns relative node index for a side vertex.
     * \param side side index (>= 0)
     * \param node side node index (>= 0)
     * \return Relative node index
     * \note Returns the node index for one of the vertices of a side.
     * \note \e side must be in the range 0 .. nSide()-1
     * \note \e node must be in the range 0 .. \ref nSideNode (side)
     * \note The return value is in the range 0 .. nNode()-1
     * \sa nNode, nSide, nSideNode
     */
    virtual int SideNode (int side, int node) const = 0;

    /**
     * \brief Checks if a node is part of the element.
     * \param node absolute node index (>= 0)
     * \return \e true if \e node appears in the element's node index list,
     *   \e false otherwise.
     */
    virtual bool IsNode (int node);

    /**
     * \brief Checks if a face defined by vertex nodes is part of the element.
     * \param nn number of node indices in the supplied array
     * \param nd array of absolute node indices of length \e nn
     * \return If all node indices in \e nd are part of the element, and are
     *   associated with a single one of the element's sides, the relative side
     *   index (>= 0) is returned. Otherwise, -1 is returned.
     */
    int IsSide (int nn, int *nd);

    /**
     * \brief Checks if a node is associated with a specific side.
     * \param side side index (>= 0)
     * \param node relative node index (>= 0)
     * \return \e true if \e node belongs to \e side, \e false otherwise.
     * \note side must be in the range 0 .. nSide()-1
     * \note node must be in the range 0 .. nNode()-1
     * \sa nSide, nNode
     */
    virtual bool IsSideNode (int side, int node);

    /**
     * \brief Checks if a side is on the mesh surface.
     * \param side side index (>= 0)
     * \return \e true if the specified side is a boundary side, \e false
     *   otherwise.
     * \note The \ref Initialise method must have been executed prior to
     *   any calls to IsBoundarySide.
     * \sa Initialise, HasBoundarySide
     */
    bool IsBoundarySide (int side) {
        RANGE_CHECK (side >= 0 && side < nSide());
	return bndside[side];
    }

    /**
     * \brief Checks if the element contains any surface sides.
     * \return \e true if the element has at least one boundary side, \e false
     *   otherwise.
     * \note The \ref Initialise method must have been executed prior to
     *   any calls to HasBoundarySide.
     * \sa Initialise, IsBoundarySide
     */
    inline bool HasBoundarySide () const
    { return bndel; }

    /**
     * \brief Checks if the element contains any internal interface sides.
     * \return \e true if the element has at least one internal interface side,
     *   \e false otherwise.
     * \note Interfaces define the boundaries between different mesh regions,
     *   e.g. for the application of boundary conditions.
     * \note The \ref Initialise method must have been executed prior to
     *   any calls to HasInterfaceSide.
     * \sa Initialise, HasBoundarySide
     */
    bool HasInterfaceSide ()
    { return interfaceel; }

    /**
     * \brief Maps a point from global to local element coordinates.
     * \param nlist mesh node list
     * \param glob global point coordinates
     * \return Point coordinates in the element's local reference system.
     * \note Point glob should have dimension Dimension()
     * \note \e glob does not have to be located inside the element.
     * \sa NodeLocal, Global
     */
    virtual Point Local (const NodeList& nlist, const Point& glob)
	const = 0;

    /**
     * \brief Returns the local coordinates of an element node.
     * \param node node index (>= 0)
     * \return Node coordinates in the element's local reference system.
     * \sa Local, Global
     */
    virtual Point NodeLocal (int node) const = 0;

    /**
     * \brief Maps a point from surface coordinates to local element
     *   coordinates.
     * \param side side index (>= 0)
     * \param p surface point coordinates.
     * \return Point coordinates in the element's local reference system.
     * \note Surface point p should have dimension Dimension()-1
     * \sa Local, NodeLocal
     */
    virtual Point SurfToLocal (int side, const Point &p) const
    { ERROR_UNDEF; return Point(); }

    /**
     * \brief Translates a point to an element surface.
     * \param [in] side side index (>= 0)
     * \param [in,out] loc point to be transformed.
     * \note Moves \e loc onto side \e side by translating it along the
     *   side normal.
     */
    virtual void MapToSide (int side, Point &loc) const;

    /**
     * \brief Returns the centre point of a side.
     * \param side side index (>= 0)
     * \return side centre, defined as the barycentre of the associated
     *   nodes.
     */
    virtual Point SideCentre (int side) const;

    /**
     * \brief returns a pointer to the element connected at 'side', or NULL
     *   if this is a boundary side.
     * \param side side index (>= 0)
     * \return Element pointer or NULL
     * \note This function is only supported if Mesh::PopulateNeighbourLists
     *   has been called before.
     */
    Element *SideNeighbour (int side) const;

    /**
     * \brief Returns the list index of the element connected at 'side', or
     *   -1 if this is a boundary side.
     * \param side side index (>= 0)
     * \return Element list index (>- 0) or -1
     * \note This function is only supported if Mesh::PopulateNeighbourLists
     *   has been called before.
     */
    int SideNeighbourIndex (int side) const;

    /**
     * \brief Maps a point from local element coordinates to global
     *   coordinates.
     * \param nlist mesh node list
     * \param loc point in local element coordinates
     * \note Point loc should have dimension Dimension()
     * \note \e loc does not have to be located inside the element.
     * \sa Local, NodeLocal
     */
    Point Global (const NodeList &nlist, const Point& loc) const;

    /**
     * \brief Returns the element's global node coordinates.
     * \param nlist mesh node list
     * \return Global node coordinates in an \ref nNode x \ref Dimension
     *   matrix.
     * \note This method extracts the coordinates of the nodes associated
     *   with the element from \e nlist and returns them as a dense matrix.
     */
    virtual RDenseMatrix Elgeom (const NodeList& nlist) const;

    /**
     * \brief Returns the direction cosines of a side normal.
     * \param side side index (>= 0)
     * \param jacin inverse of Jacobian
     * \return direction cosines of the normal to \e side in global
     *   coordinates.
     * \sa LNormal
     */
    virtual RVector DirectionCosine (int side, RDenseMatrix &jacin) = 0;

    /**
     * \brief Returns a side normal in local coordinates.
     * \param side side index (>= 0)
     * \return Normal to \e side in local element coordinates.
     * \sa DirectionCosine
     */
    virtual const RVector &LNormal (int side) const = 0;

    /**
     * \brief Returns the element size.
     * \return Area (for 2-D elements) or volume (for 3-D elements) of the
     *   element.
     * \sa SideSize
     */
    virtual double Size() const = 0;

    /**
     * \brief Returns the size of an element side.
     * \param side side index (>= 0)
     * \param nlist mesh node list
     * \return Size of the element face (length for 2-D elements, area for
     *   3-D elements).
     * \sa Size
     */
    virtual double SideSize (int side, const NodeList &nlist) const
    { ERROR_UNDEF; return 0; }

    /**
     * Returns the Jacobian matrix of global wrt. local coordinates,
     * evaluated at local point loc
     */
    double Jacobian (const Point &loc, const NodeList *nlist, RDenseMatrix &J)
	const;

    /**
     * Returns the inverse of the Jacobian matrix of global wrt. local
     * coordinates, evaluated at local point loc
     */
    double IJacobian (const Point &loc, const NodeList *nlist,
	RDenseMatrix &IJ) const;

    /**
     * \brief Returns determinant of Jacobian at a given point inside the
     *   element in the local frame.
     * \param loc evaluation point in local coordinates
     * \param nlist node list of the associated mesh
     * \return det(J)
     */
    virtual double DetJ (const Point &loc, const NodeList *nlist = 0) const
    { ERROR_UNDEF; return 0.0; }

    /**
     * \brief Checks if a local point coordinate is inside the element.
     * \param loc point in local coordinates
     * \param pad Apply padding to the element to ensure that boundary points
     *   are considered inside.
     * \return \e true if the point is inside the element, \e false otherwise.
     * \note If \e pad == false, then the status of points exactly on the
     *   element boundary is undefined. They may or may not be considered
     *   inside. Use pad = true to ensure that they are regarded as inside.
     * \sa GContains
     */
    virtual bool LContains (const Point& loc, bool pad = true) const = 0;

    /**
     * \brief Checks if a global point coordinate is inside the element.
     * \param glob point in global coordinates
     * \param nlist mesh node list
     * \return \e true if the point is inside the element, \e false otherwise.
     * \sa LContains
     */
    virtual bool GContains (const Point& glob, const NodeList &nlist) const;
  
    /**
     * \brief Returns a list of boundary sides.
     * \param nlist mesh node list
     * \param address of list to be filled with boundary side indices
     * \return Number of boundary sides in the element.
     * \note This method fills array \e list with side indices (>= 0) of
     *   all boundary sides associated with the element.
     * \note \e list must have been allocated by the caller with sufficient
     *   length to receive all boundary side indices.
     * \bug This method assumes that a side is a boundary side if all its
     *   associated nodes are boundary nodes. This is not always correct, in
     *   particular at corners.
     * \sa IsBoundarySide, HasBoundarySide, IsSideNode
     */
    virtual int BndSideList (const NodeList& nlist, int *list);

    void InitNeighbourSupport ();
    void InitSubdivisionSupport ();

    /**
     * \brief Returns the element's subdivison level
     * \return subdivision level (>= 0)
     * \note By default, all elements have level 0. The subdivision level
     *   can increase during dynamic mesh refinement.
     */
    int SubdivisionLevel () const;

    virtual void Subdivide (Mesh *mesh)
    { ERROR_UNDEF; }

    /**
     * \brief Returns the values of the shape functions at a local point.
     * \param loc point coordinates in the local element system
     * \return Vector (size \ref nNode) of shape function values at the point
     *   for all shape functions associated with the element.
     * \sa GlobalShapeF, LocalShapeD, GlobalShapeD
     */
    virtual RVector LocalShapeF (const Point &loc) const = 0;

    virtual void LocalShapeF (const Point &loc, RVector *fun) const
    { *fun = LocalShapeF (loc); }

    /**
     * \brief Returns the values of the shape function derivatives at a local
     *   point.
     * \param loc point coordinates in the local element system
     * \return Matrix (size \ref nNode x \ref Dimension) of shape function
     *   derivatives \f$ du_i/d x_j \f$ at the point for all shape functions
     *   associated with the element.
     * \sa LocalShapeF, GlobalShapeF, GlobalShapeD
     */
    virtual RDenseMatrix LocalShapeD (const Point &loc) const = 0;

    /**
     * \brief Returns the values of the shape functions at a global point.
     * \param nlist mesh node list
     * \param glob point coordinates in the local element system
     * \return Vector (size \ref nNode) of shape function values at the point
     *   for all shape functions associated with the element.
     * \sa LocalShapeF, LocalShapeD, GlobalShapeD
     */
    virtual RVector GlobalShapeF (const NodeList &nlist, const Point &glob)
    const { return LocalShapeF (Local (nlist, glob)); }

    /**
     * \brief Returns the values of the shape function derivatives at a global
     *   point.
     * \param glob global point coordinates
     * \return Matrix (size \ref nNode x \ref Dimension) of shape function
     *   derivatives \f$ du_i/d x_j \f$ at the point for all shape functions
     *   associated with the element.
     * \sa LocalShapeF, GlobalShapeF, GlobalShapeD
     */
    virtual RDenseMatrix GlobalShapeD (const NodeList &nlist,
        const Point &glob) const;

    /**
     * \brief Returns the weights and abscissae of quadrature rules over the
     *   element.
     * \param order quadrature order (>= 1)
     * \param wght pointer to an array of quadrature weights
     * \param absc pointer to an array of quadrature abscissae
     * \return Number of quadrature points
     * \note The required order of a quadrature rule depends on the type of
     *   integral to be performed. Examples:
     *   - IntF:  order = 1
     *   - IntFF, IntD: order = 2
     *   - intFFF, IntFD: order = 3
     *   - IntDD: order = 4
     *   etc.
     */
    virtual int QuadRule (int order, const double **wght, const Point **absc)
        const { ERROR_UNDEF; return 0; }

    /**
     * \brief Integral of a shape function over the element.
     * \param i node index (range: 0 .. \ref nNode-1)
     * \return value of the integral 
     *   \f[ \int_\Omega u_i(\vec{r}) d\vec{r} \f]
     * \sa IntFF, IntFFF, IntDD, IntFD
     */
    virtual double IntF (int i) const = 0;

    /**
     * \brief Integrals of shape functions over the element.
     * \return vector of integrals
     *   \f[ \int_\Omega u_i(\vec{r}) d\vec{r} \f]
     *   for all element nodes i
     * \sa IntF(int)
     */
    virtual RVector IntF () const;

    /**
     * \brief Integrals of all products of two shape functions over the
     *   element.
     * \return Matrix of size \ref nNode x \ref nNode, containing the integrals
     *   \f[ \int_\Omega u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f]
     * \sa IntFF(int,int)const, IntF, IntFFF, IntDD, IntFD
     */
    virtual RSymMatrix IntFF () const = 0;

    /**
     * \brief Integral of a product of two shape functions over the
     *   element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \return Value of the integral 
     *   \f[ \int_\Omega u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f]
     * \sa IntFF()const, IntF, IntFFF, IntDD, IntFD
     */
    virtual double IntFF (int i, int j) const = 0;

    /**
     * \brief Integral of a product of three shape functions over the
     *   element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \param k third node index (range 0 .. \ref nNode-1)
     * \return Value of the integral 
     *   \f[ \int_\Omega u_i(\vec{r}) u_j(\vec{r}) u_k(\vec{r}) d\vec{r}
     *   \f]
     * \sa IntF, IntFF, IntDD, IntFD
     */
    virtual double IntFFF (int i, int j, int k) const = 0;

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

    virtual void IntFFF (double &iii, double &iij, double &ijk) const = 0;
    // returns the values of the FFF tensor for:
    //   all indices equal (iii)
    //   two indices equal (iij)
    //   all indices different (ijk)
#endif

    /**
     * \brief Integrals of all products of two shape functions and a nodal
     *   function over the element.
     * \param P nodal basis coefficients of a function defined over the mesh
     * \return Matrix of size \ref nNode x \ref nNode, containing the integrals
     *   \f[ \int_\Omega p(\vec{r}) u_i(\vec{r}) u_j(\vec{r}) d\vec{r} =
     *   \sum_k p_k \int_\Omega u_i(\vec{r}) u_j(\vec{r})
     *   u_k(\vec{r}) d\vec{r}  \f]
     * \sa IntPFF(int,int,const RVector&)const, \n
     *   IntF, IntFF, IntFFF, IntFD, IntDD, IntPDD
     */
    virtual RSymMatrix IntPFF (const RVector& P) const = 0;

    /**
     * \brief Integral of a product of two shape functions and a nodal
     *   function over the element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \param P nodal basis coefficients of a function defined over the mesh
     * \return The value of the integral
     *   \f[ \int_\Omega p(\vec{r}) u_i(\vec{r}) u_j(\vec{r}) d\vec{r} =
     *   \sum_k p_k \int_\Omega u_i(\vec{r}) u_j(\vec{r})
     *   u_k(\vec{r}) d\vec{r}  \f]
     * \sa IntPFF(const RVector&)const, \n
     *   IntF, IntFF, IntFFF, IntFD, IntDD, IntPDD
     */
    virtual double IntPFF (int i, int j, const RVector& P) const = 0;

    /**
     * \brief Integrals of all products of two shape function derivatives over
     *   the element.
     * \return Matrix of size \ref nNode x \ref nNode, containing the integrals
     *   \f[ \int_\Omega \nabla u_i(\vec{r}) \nabla u_j(\vec{r}) d\vec{r}
     *   \f]
     * \sa IntDD(int,int)const, \n
     *   IntF, IntFF, IntFFF, IntFD, IntDD, IntPFF, IntPDD
     */
    virtual RSymMatrix IntDD () const = 0;

    /**
     * \brief Integral of a product of two shape function derivatives over
     *   the element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \return Value of the integral
     *   \f[ \int_\Omega \nabla u_i(\vec{r}) \nabla u_j(\vec{r}) d\vec{r}
     *   \f]
     * \sa IntDD()const, \n
     *   IntF, IntFF, IntFFF, IntFD, IntDD, IntPFF, IntPDD
     */
    virtual double IntDD (int i, int j) const = 0;

    /**
     * \brief Integral of a product of a shape function and a shape function
     *   derivative over the element.
     * \param i node index for shape function
     * \param j node index for shape function derivative
     * \return Vector of size \ref Dimension, containing the integral
     *   \f[ \int_\Omega u_i(\vec{r}) \nabla u_j(\vec{r}) d\vec{r} \f]
     * \sa IntDD, IntFDD, IntPDD
     */
    virtual RVector IntFD (int i, int j) const
    { ERROR_UNDEF; return RVector(); }

    /**
     * \brief Integral of a product of a shape function and two shape function
     *   derivatives over the element.
     * \param i node index for shape function
     * \param j first node index for shape function derivative
     * \param k second node index for shape function derivative
     * \return Value of the integral
     *   \f[ \int_\Omega u_i(\vec{r}) \nabla u_j(\vec{r})
     *   \nabla u_k(\vec{r}) d\vec{r} \f]
     * \sa IntDD, IntFD, IntPDD
     */
    virtual double IntFDD (int i, int j, int k) const = 0;

    /**
     * \brief All integrals of products of a nodal function and two shape
     *   function derivatives over the element.
     * \param P nodal function coefficients over the mesh
     * \return Matrix of size \ref nNode x \ref nNode, containing the integrals
     *   \f[ \int_\Omega p(\vec{r}) \nabla u_i(\vec{r})
     *   \nabla u_j(\vec{r}) d\vec{r} =
     *   \sum_k p_k \int_\Omega u_k(\vec{r}) \nabla u_i(\vec{r})
     *   \nabla u_j(\vec{r}) d\vec{r} \f]
     * \sa IntPDD(int,int,const RVector&)const, IntFDD, IntFD, IntPFF
     */
    virtual RSymMatrix IntPDD (const RVector& P) const = 0;

    /**
     * \brief Integrals of a product of a nodal function and two shape
     *   function derivatives over the element.
     * \param i first node index
     * \param j second node index
     * \param P nodal function coefficients over the mesh
     * \return Value of the integral
     *   \f[ \int_\Omega p(\vec{r}) \nabla u_i(\vec{r})
     *   \nabla u_j(\vec{r}) d\vec{r} =
     *   \sum_k p_k \int_\Omega u_k(\vec{r}) \nabla u_i(\vec{r})
     *   \nabla u_j(\vec{r}) d\vec{r} \f]
     * \sa IntPDD(const RVector&)const, IntFDD, IntFD, IntPFF
     */
    virtual double IntPDD (int i, int j, const RVector &P) const = 0;

    /**
     * \brief %Surface integral of a shape function over an element face.
     * \param i node index (range 0 .. \ref nNode-1)
     * \param sd side index (range 0 .. \ref nSide-1)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa SurfIntF(int), BndIntF(int), SurfIntFF(int,int,int)
     */
    virtual double SurfIntF (int i, int sd) const = 0;

    /**
     * \brief %Surface integrals of all shape functions over an element face.
     * \param sd side index (range 0 .. \ref nSide-1)
     * \return Vector of size \ref nNode containing the integrals
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) d\vec{r} \f$
     *   for all element nodes i, where the integration is performed over side
     *   \f$S_{sd}\f$.
     % \sa SurfIntF(int,int), BndIntF()
     */
    virtual RVector SurfIntF (int sd) const;
    
    /**
     * \brief %Surface integral of a shape function over all boundary faces
     *   of the element.
     * \param i node index (range 0 .. \ref nNode-1)
     * \return Value of the integral
     *   \f$ \oint_S u_i(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The returned value is nonzero only if node i is a boundary node.
     */
    virtual double BndIntF (int i) const;
    
    /**
     * \brief %Surface integrals of all shape functions over all boundary
     *   faces of the element.
     * \return Vector of size \ref nNode containing the integrals
     *   \f$ \oint_S u_i(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The returned vector contains nonzero entries at index i only if
     *   node i is a boundary node.
     * \note If the element does not contain boundary sides, the returned
     *   vector is zero.
     * \sa SurfIntF(int), BndIntF(int), BndIntFF
     */
    virtual RVector BndIntF () const;

    /**
     * \brief %Surface integral of a product of two shape functions over an
     *   element side.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \param sd side index (range 0 .. \ref nSide-1)
     * \return Value of the integral
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa BndIntFF()const, BndIntFF(int,int)
     */
    virtual double SurfIntFF (int i, int j, int sd) const = 0;

    /**
     * \brief %Surface integrals of all products of two shape functions over
     *   and element side.
     * \param sd side index (range 0 .. \ref nSide-1)
     * \return Symmetrix matrix of size \ref nNode x \ref nNode, containing
     *   the integrals
     *   \f$ \oint_{S_{sd}} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over side \f$S_{sd}\f$.
     * \sa SurfIntFF(int,int,int)const
     */
    virtual RSymMatrix SurfIntFF (int sd) const;
    
    /**
     * \brief %Surface integrals of all products of two shape functions over
     *   all boundary sides of the element.
     * \return Symmetric matrix of size \ref nNode x \ref nNode, containing
     *   the integrals
     *   \f$ \oint_S u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f$
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The returned matrix contains nonzero entries at (i,j) only if
     *   nodes i and j are both boundary nodes.
     * \note If the element does not contain boundary sides, the returned
     *   matrix is zero.
     * \sa BndIntFF(int,int), SurfIntFF
     */
    virtual RSymMatrix BndIntFF () const;

    /**
     * \brief %Surface integral of a product of two shape functions over
     *   all boundary sides of the element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \return Value of the integral
     *   \f[ \int_{\partial\Omega} u_i(\vec{r}) u_j(\vec{r}) d\vec{r} \f]
     *   where the integration is performed over all sides of the element that
     *   are part of the mesh surface.
     * \note The return value is nonzero only if both nodes \p i and \p j are
     *   boundary nodes.
     * \sa BndIntFF()const, SurfIntFF
     */
    virtual double BndIntFF (int i, int j) const;

    /**
     * \brief Surface integrals of all products of a nodal function and two
     *   shape functions over all boundary sides of the element.
     * \param P nodal function coefficients over the mesh
     * \return Matrix of size \ref nNode x \ref nNode, containing the integrals
     *   \f[ \int_{\partial\Omega} p(\vec{r}) u_i(\vec{r}) u_j(\vec{r})
     *       d\vec{r} = \sum_k p_k \int_{\partial\Omega} u_i(\vec{r})
     *       u_j(\vec{r}) u_k(\vec{r}) d\vec{r} \f]
     * \note If the element does not contain boundary sides, the returned
     *   matrix is zero.
     * \sa BndIntPFF(int,int,const RVector&)const, BndIntFF
     */
    virtual RSymMatrix BndIntPFF (const RVector &P) const = 0;

    /**
     * \brief Surface integrals of a product of a nodal function and two
     *   shape functions over all boundary sides of the element.
     * \param i first node index (range 0 .. \ref nNode-1)
     * \param j second node index (range 0 .. \ref nNode-1)
     * \param P nodal function coefficients over the mesh
     * \return Value of the integral
     *   \f[ \int_{\partial\Omega} p(\vec{r}) u_i(\vec{r}) u_j(\vec{r})
     *       d\vec{r} = \sum_k p_k \int_{\partial\Omega} u_i(\vec{r})
     *       u_j(\vec{r}) u_k(\vec{r}) d\vec{r} \f]
     * \note If the element does not contain boundary sides, the returned
     *   matrix is zero.
     * \sa BndIntPFF(const RVector&)const, BndIntFF
     */
    virtual double BndIntPFF (int i, int j, const RVector &P) const = 0;

    // ===============================================================
    // Integrals including mixed derivatives

    /**
     * \brief Integral of partial shape function derivative over the element
     * \param i node index (range 0 .. \ref nNode-1)
     * \param k coordinate index for the derivative  (range 0 .. \ref
     *   Dimension-1)
     * \return Value of the integral
     *   \f[ \int_\Omega \frac{\partial u_i(\vec{r})}{\partial x_k} d\vec{r}
     *   \f]
     */
    virtual double Intd (int i, int k) const
    { ERROR_UNDEF; return 0.0; }

    /**
     * \brief Integral of the product of a shape function and a partial shape
     *   function derivative over the element.
     * \param i node index for shape function (range 0 .. \ref nNode-1)
     * \param j node index for shape function derivative
     *   (range 0 .. \ref nNode-1)
     * \param k coordinate index for derivative (range 0 .. \ref Dimension-1)
     * \return Value of the integral
     *   \f[ \int_\Omega u_i(\vec{r})
     *   \frac{\partial u_j(\vec{r})}{\partial x_k} d\vec{r} \f]
     */
    virtual double IntFd (int i, int j, int k) const
    { ERROR_UNDEF; return 0.0; }

    /**
     * \brief Integral of the product of a nodal function and a partial shape
     *   function derivative over the element.
     * \param P array of nodal coefficients defining the function
     * \param j node index for shape function derivative
     *   (range 0 .. \ref nNode-1)
     * \param k coordinate index for derivative (range 0 .. \ref Dimension-1)
     * \return Value of the integral
     *   \f[ \sum_i P_i \int_\Omega u_i(\vec{r}) \frac{\partial u_j(\vec{r})}{\partial x_k} d\vec{r} \f]
     */
    virtual double IntPd (const RVector &P, int j, int k) const
    { ERROR_UNDEF; return 0.0; }

    /**
     * \brief Integral of the product of two partial shape function
     *   derivatives over the element.
     * \return Matrix of integrals for all element nodes in each dimension.
     *   \f[ \int_\Omega \frac{\partial u_i(\vec{r})}{\partial x_j} \frac{partial u_k(\vec{r})}{\partial x_l} d\vec{r} \f]
     *   The dimension of the returned matrix is nd x nd, where n is the
     *   number of nodes, and d is the dimension of the element (2 or 3).
     *   The matrix contains n x n blocks, where each block is of
     *   dimension d x d.
     */
    virtual RSymMatrix Intdd() const
    { ERROR_UNDEF; return RSymMatrix(); }

    /**
     * \brief Integral of the product of a shape function and two partial
     *   shape function derivatives over the element.
     * \param i node index for shape function (range 0 .. \ref nNode-1)
     * \param j node index for shape function derivative  (range 0 .. \ref nNode-1)
     * \param k node index for shape function derivative  (range 0 .. \ref nNode-1)
     * \param l coordinate index for derivative (range 0 .. \ref Dimension-1)
     * \param m coordinate index for derivative (range 0 .. \ref Dimension-1)
     * \return Value of the integral
     *   \f[ \int_\Omega u_i(\vec{r}) \frac{\partial u_j(\vec{r})}{\partial x_l} \frac{\partial u_k(\vec{r})}{\partial x_m} \f]
     */
    virtual double IntFdd (int i, int j, int k, int l, int m) const
    { ERROR_UNDEF; return 0; }

    virtual double IntPdd (const RVector &p, int j, int k, int l, int m) const
    { ERROR_UNDEF; return 0; }
    // Int f(r) du_j/dx_l du_k/dx_m dr
    // where f(r) is given as a nodal vector

  virtual double IntFfd (int i, int j, int k, int l) const
    { ERROR_UNDEF; return 0; }
    // Int [u_i u_j du_k/dx_l] dr

   virtual double IntPfd (const RVector &p, int j, int k, int l) const
    { ERROR_UNDEF; return 0; }
    // Int f(r) u_j du_k/dx_l dr
    // where f(r) is given as a nodal vector

    //------------- integrate functions on unit sphere --------------
    virtual double IntUnitSphereP (const NodeList& nlist, const RVector& P) const
    { ERROR_UNDEF; return 0.0; }
    // Returns integral of f(r) over unitsphere patch 
    // only implemented for Tri3D3, Tri3D6

     virtual double IntUnitSphereFF (const NodeList& nlist, const int i, const int j) const
    { ERROR_UNDEF; return 0.0; }
    // Returns integral of u_i(r)*u_j(r) over unitsphere patch 
    // only implemented for Tri3D3, Tri3D6

     virtual double IntUnitSpherePFF (const NodeList& nlist, const int i, const int j, const RVector& P) const
    { ERROR_UNDEF; return 0.0; }
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
    { ERROR_UNDEF; return RVector(); }
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
    { ERROR_UNDEF; return RDenseMatrix(); }
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
    { ERROR_UNDEF; return RDenseMatrix(); }
    // Returns the element elasticity stiffness matrix K, given
    // strain matrix D: K = \int_el B^T D B dr
    // where B is strain displacement matrix

    virtual RDenseMatrix ElasticityStiffnessMatrix (double E,
        double nu) const
    { ERROR_UNDEF; return RDenseMatrix(); }
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
    { ERROR_UNDEF; return RDenseMatrix(); }
    // returns a 6x6 symmetric elasticity matrix for isotropic materials,
    // given Young's modulus E and Poisson's ratio nu

    virtual RVector InitialStrainVector (double E, double nu, const RVector &e0)
    { ERROR_UNDEF; return RVector(); }
    // returns a element strain vector (rhs) for initial nodal strain given
    // in e0

    virtual RVector ThermalExpansionVector (double E, double nu,
        double alpha, double dT)
    { ERROR_UNDEF; return RVector(); }
    // returns element strain vector (rhs) given element thermal expansion
    // coefficient alpha and temperature change dT

    virtual RVector DThermalExpansionVector (double E, double nu)
    { ERROR_UNDEF; return RVector(); }
    // returns derivative of thermal expansion vector with respect to
    // expansion coefficient (assuming temperature change 1)

    /**
     * \brief Return intersection points of a ray with the element surface.
     * \param p1 first ray endpoint (in global frame)
     * \param p2 second ray endpoint (in global frame)
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
     * \note This method assumes global coordinates for the ray endpoints,
     *   but still returns the intersection points in the local element frame.
     * \sa Intersection
     */
    virtual int GlobalIntersection (const NodeList &nlist,
        const Point &p1, const Point &p2, Point *s,
	bool add_endpoints, bool boundary_only);

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
    virtual int Intersection (const Point &p1, const Point &p2,
        Point *s, bool add_endpoints, bool boundary_only) = 0;

    virtual RDenseMatrix LocaltoGlobalMat () const 
    { ERROR_UNDEF; return RDenseMatrix(); }
    // abstract - local to global matrix

    virtual RDenseMatrix GlobaltoLocalMat () const 
    { ERROR_UNDEF; return RDenseMatrix(); }
    // abstract - global to local matrix

    virtual RDenseMatrix FTAMat () const 
    { ERROR_UNDEF; return RDenseMatrix(); }
    // abstract - transform k-space matrix

    friend std::istream& operator>> (std::istream& i, Element &el);
    friend std::ostream& operator<< (std::ostream& o, const Element &el);
    // standard method for reading/writing an element from/to a stream

    int *Node;
    // list of absolute node numbers for the element's vertices

    //int Region;
    inline int Region() const {return region;} // region id for segmented meshes
    inline void SetRegion(int nr) {region = nr;} // set region id 

    virtual void SplitSide (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2)
    { ERROR_UNDEF; }

    virtual void Bisect (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2)
    { ERROR_UNDEF; }

    virtual void MergeAndResplit (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2)

    { ERROR_UNDEF; }

    ElementSubdivisionData *subdivdata;

    Element **sdnbhr;
    // List of side neighbours for the element.
    // A NULL entry indicates no neighbour (boundary side)
    // Note: this list is only created if Mesh::PopulateNeighbourLists is
    // called.

    int *sdnbhridx;
    // list of side neighbour indices for the element
    // A -1 entry indicates no neighbour (boundary side)
    // Note: this list is only created if Mesh::PopulateNeighbourLists is
    // called.

protected:

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

/**
 * \brief Base class for all structured (regular) element types.
 *
 * All structured elements element in a mesh are assumed to have the same
 * shape, size and orientation. Structured element types share global
 * parameters, which makes them more memory-efficient than unstructured types.
 */
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

/**
 * \brief Base class for all unstructured element types.
 *
 * Unstructured elements within a mesh differ in shape and size. Therefore
 * they must store shape parameters on an individual element basis.
 */
class FELIB Element_Unstructured: public Element {
public:
    Element_Unstructured (): Element () {}
    Element_Unstructured (const Element_Unstructured &el): Element (el) {}
    virtual ~Element_Unstructured () {}
    virtual void Initialise (const NodeList& nlist);
    virtual void PostInitialisation (const NodeList& nlist);
    inline double Size() const { return size; }
    void operator= (const Element_Unstructured& el);

    virtual bool GContains (const Point& glob, const NodeList& nlist)
	const;

    inline RSymMatrix IntDD () const { return intdd; }
    inline double IntDD (int i, int j) const { return intdd.Get(i,j); }
    inline RSymMatrix BndIntFF () const { return intbff; }
    inline double BndIntFF (int i, int j) const { return intbff.Get(i,j); }

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

/**
 * \brief Base class for all 2-D structured element types.
 */
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

/**
 * \brief Base class for all 3-D structured element types.
 */
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

/**
 * \brief Base class for all 2-D unstructured element types.
 */
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

/**
 * \brief Base class for all 3-D unstructured element types.
 */
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
