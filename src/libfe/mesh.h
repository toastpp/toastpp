// -*-C++-*-
// ==========================================================================
// Module libfe
// File mesh.h
// Declaration of class Mesh
// ==========================================================================

#ifndef __MESH_H
#define __MESH_H

#define MESH_DIRICHLET BND_DIRICHLET
#define MESH_ROBIN     BND_ROBIN
#define MESH_EXTRAPOL  BND_INTERNAL
#define MESH_ODD       BND_NONE

// system matrix assembly modes
#define ASSEMBLE_FF         0
#define ASSEMBLE_DD         1
#define ASSEMBLE_PFF        2
#define ASSEMBLE_PDD        3
#define ASSEMBLE_BNDPFF     4
#define ASSEMBLE_CFF        5
#define ASSEMBLE_CDD        6
#define ASSEMBLE_iCFF       7
#define ASSEMBLE_iCDD       8
#define ASSEMBLE_PFF_EL     9
#define ASSEMBLE_PDD_EL    10
#define ASSEMBLE_BNDPFF_EL 11
#define ASSEMBLE_BNDFF     12
#define ASSEMBLE_iPFF      13

// rhs assembly modes
#define RHS_P           0
#define RHS_PF          1
#define RHS_BNDPF       2

// ==========================================================================
// Prototypes

FELIB void AddToElMatrix (const Mesh &mesh, int el,
    RGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToElMatrix (const Mesh &mesh, int el,
    FGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToElMatrix (const Mesh &mesh, int el,
    SCGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToElMatrix (const Mesh &mesh, int el,
    CGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    RGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    FGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    SCGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    CGenericSparseMatrix &M, const RVector *coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    CGenericSparseMatrix &M, const double coeff, int mode);

FELIB void AddToSysMatrix (const Mesh &mesh,
    SCGenericSparseMatrix &M, const double coeff, int mode);

FELIB void AddToSysMatrix_elasticity (const Mesh &mesh,
    RGenericSparseMatrix &M, const RVector &modulus,
    const RVector &pratio);

FELIB void AddToRHS_elasticity (const Mesh &mesh, RVector &rhs,
    const RVector *coeff, int mode);

FELIB void AddToRHS_thermal_expansion (const Mesh &mesh,
    RVector &rhs, const RVector &modulus, const RVector &pratio,
    const RVector &thermal_expansion, double deltaT);

FELIB RGenericSparseMatrix *GridMapMatrix (const Mesh &mesh,
    const IVector &gdim, const Point *bbmin, const Point *bbmax,
    const int *elref);

FELIB void AddElasticStrainDisplacementToSysMatrix (const Mesh &mesh,
    RGenericSparseMatrix &M, double E, double nu, int *nf);

FELIB RGenericSparseMatrix *NodeMapMatrix (const Mesh &mesh,
    const IVector &gdim, const Point *bbmin, const Point *bbmax,
    const int *elref);

FELIB RGenericSparseMatrix *NodeMapMatrix2 (const Mesh &mesh,
    const IVector &gdim, const Point *bbmin = 0, const Point *bbmax = 0,
    const int *elref = 0);

FELIB RGenericSparseMatrix *Grid2LinPixMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref = 0);

FELIB RGenericSparseMatrix *LinPix2GridMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref = 0);

FELIB RGenericSparseMatrix *CubicPix2GridMatrix (const IVector &bdim,
    const IVector &gdim, const int *elref = 0);

FELIB void GenerateVoxelPositions (const Mesh &mesh, const IVector &gdim,
     const Point *bbmin, const Point *bbmax, RDenseMatrix &pos);

FELIB void SubsampleLinPixel (const RVector &gf, RVector &bf,
    const IVector &gdim, const IVector &bdim, const int *elref = 0);

FELIB RGenericSparseMatrix *NodeMapMatrix (const Mesh &mesh,
    const IVector &gdim, const Point *bbmin = 0, const Point *bbmax = 0,
    const int *elref = 0);

FELIB RGenericSparseMatrix *GridMapMatrix (const Mesh &mesh,
    const IVector &gdim, const Point *bbmin = 0, const Point *bbmax = 0,
    const int *elref = 0);

FELIB int *GenerateElementPixelRef (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax);

FELIB void RasterLinPixel (RVector &gf, const RVector &bf,
    const IVector &gdim, const IVector &bdim, const int *elref = 0);

FELIB void SubsampleCubPixel (const RVector &gf, RVector &bf,
    const IVector &gdim, const IVector &bdim, const int *elref = 0);

FELIB void RasterCubPixel (RVector &gf, const RVector &bf,
    const IVector &gdim, const IVector &bdim, const int *elref = 0);

FELIB Mesh *Lin2Quad (const Mesh &linmesh);

// ==========================================================================
// class Mesh
// ==========================================================================
/**
 * \brief Finite-element mesh management
 *
 * Defines node coordinates (via NodeList) and element connectivity (via
 * ElementList) and provides general mesh-related methods.
 * This version also defines node-based parameter coefficient sets for
 * optical tomography (via ParameterList), but this should be removed.
 */

class FELIB Mesh {
public:

    NodeList nlist;      // list of mesh nodes
    ElementList elist;   // list of mesh elements

    mutable int **nbhrs;
    // list of edge-adjacent neighbours for each element
    // this is only defined after a call to SetupNeighbourList

    int BndType;
    // either of MESH_DIRICHLET, MESH_ROBIN, MESH_EXTRAPOL, MESH_ODD

    /**
     * \brief Constructs an empty mesh.
     * \note The new mesh has no nodes, elements or parameters.
     * \note The nlist, elist and plist data members must be initialised
     *   to define an actual mesh geometry.
     */
    Mesh ();

    Mesh (const Mesh &mesh);

    /**
     * \brief Destructor
     */
    virtual ~Mesh ();

    /**
     * \brief Initialise the mesh.
     *
     * This method should be called once the mesh geometry (i.e. the
     *   node and element lists) have been defined and after any changes in
     *   these lists.
     * \param mark_boundary If true, this method marks boundary nodes based
     *   on the mesh connectivity defined in the element list. If false, the
     *   boundary flags remain unchanged.
     */
    void Setup (bool mark_boundary=true);

    /**
     * \brief Copy geometry from another mesh.
     * \param mesh Source mesh to be copied.
     * \note This method replaces the current mesh geometry with that of the
     *   argument.
     */
    void Copy (const Mesh &mesh);
  
    /**
     * \brief Mesh dimension
     * \return Returns the mesh dimension (2 or 3)
     */
    int Dimension() const 
    { return (elen() ? elist[0]->Dimension() : 0); }

    /**
     * \brief Number of mesh nodes
     * \return Length of the node list
     */
    int nlen() const { return nlist.Len(); }

    /**
     * \brief Number of mesh elements
     * \return Length of the element list
     */
    int elen() const { return elist.Len(); }

    /**
     * \brief Number of internal (non-boundary) nodes
     * \return Number of internal nodes
     * \note The return value is only valid if boundary nodes have been
     *   identified, e.g. after Setup().
     */
    int ilen() const { return priv_ilen; }

    /**
     * \brief Number of boundary nodes
     * \return Number of boundary nodes
     * \note The return value is only valid if boundary nodes have been
     *   identified, e.g. after Setup().
     */
    int nbnd() const { return priv_nbnd; }

    RDenseMatrix ElGeom (int el) const;
    // returns a matrix containing the global coordinates of all nodes of
    // element el

    double ElSize (int el) const;
    // returns the size of element el (area for 2D, volume for 3D) in global
    // coordinates.

    Point ElCentre (int el) const;
    // returns the centre of element el in global coordinates

    Point ElSideCentre (int el, int sd) const;
    // returns the centre of side sd of element el in global coordinates

    double ElSideSize (int el, int sd) const;
    // returns the size of side 'sd' of element 'el'

    RVector ElDirectionCosine (int el, int sd, Point *p = 0) const;
    // returns the direction cosine of side 'sd' of element 'el' at
    // local point 'p'

    /**
     * \brief Returns the mesh area or volume.
     * \return Mesh area (for 2-D meshes) or volume (for 3-D meshes)
     * \note The return value is the sum of all element areas or volumes.
     */
    double FullSize () const;

    int ElFind (const Point &pt) const;
    // returns the number of the element which contains the global coordinate
    // 'pt', or -1 if no element is found

    int NdFind (const Point &pt, double &dist) const;
    // returns the number of the node closest to pt and its distance

    /**
     * \brief Re-order mesh nodes.
     *
     * This method reorders the mesh nodes such that N[i] <- N[perm[i]].
     * The node list is rearranged, as well as the indices in the element
     * list.
     * \param perm Permutation list (must be of length nlen() and contain
     *   all indices from 0 to nlen()-1)
     * \note After calling this method, any node-based arrays defined for
     *   the original mesh (e.g. nodal parameter lists or fields) must
     *   be updated with the same permutation list to conform to the
     *   modified mesh.
     */
    void Reorder (const IVector &perm);

    /**
     * \brief Remove unused nodes and update the index list
     *
     * If the mesh contains any nodes not referenced by any entry in
     * the index list, those nodes are removed and the index list is
     * updated accordingly.
     * \return true if the mesh was modified, false if no unused nodes
     *   were found.
     * \note Any external object that referenced the original mesh
     *   (nodal arrays, basis mapping matrices, etc.) will be invalid
     *   after this operation.
     */
    bool Shrink ();
    
    double ElDist (int el1, int el2) const;
    // returns the distance between the centres of elements el1 and el2

    virtual void ScaleMesh (double scale);
    // rescale the complete mesh

    virtual void ScaleMesh (const RVector &scale);
    // Anisotropic scaling. 'scale' must have the same dimension as the mesh

    int MaxNodeDiff (int mode=BW_AUTO) const;
    // calculates the maximum node number difference within a single element
    // of the mesh, either using all nodes (mode=BW_TOTAL) or only internal
    // nodes (mode=BW_INTERNAL). By default (mode=BW_AUTO) it uses all nodes
    // for Robin meshes, and internal nodes for all other meshes.
    // To get the semi-bandwith of the system matrix including the diagonal,
    // add 1 to the return value of MaxNodeDiff

    Point NeighbourBarycentre (int node);
    // returns the barycentre position of the neighbours of 'node'.

    void SparseRowStructure (idxtype *&rowptr, idxtype *&colidx, int &nzero) const;
    // generate row and column index lists for a system matrix in
    // compressed row format corresponding to the mesh
    // See SparseLib++ documentation, "Compressed row storage" for format
    // Lists rowptr and colidx are allocated internally and must be freed
    // by the caller when done.

    void NeighbourCount (int *plist, int nnode, bool include_self = false)
        const;
    // returns a list of length intdof or totdof, containing for each node the
    // number of neighbouring nodes, i.e. the number of nodes which share an
    // element with the first node.
    // nnode should be totdof for Robin meshes, or intdof for all others
    // If include_self==TRUE then each entry of plist is incremented by 1

    void SysMatrixStructure (int *nz, int **row_ind, int **col_ind);
    // returns the number of nonzero entries in the system matrix
    // and lists of row and column indices for each nozero entry
    // row_ind and col_ind must be deleted after use

    int FindBoundarySegment (const Point &p, int *n1, int *n2, double *dist1,
	double *dist2) const;
    // Calculate the two neighbouring boundary nodes closest to p, and their
    // distance to p. Return the number of the element owning these nodes
    // Uses physical boundary for extrapolated meshes

    double BoundaryDistance (int node) const;
    // Returns the distance of node `node' from the boundary.
    // (Currently this only returns the distance between `node' and the
    // closest boundary node)

    int BoundaryList (int **bndellist, int **bndsdlist) const;
    // Generate a list of boundary segments. Each segment i is defined by the
    // element number *bndellist[i] and relative side number *bndsdlist[i].
    // return value is the number of boundary segments found.

    /**
     * \brief Returns a list of outer and internal surface faces, based
     *   on element region information.
     * \param [out] idx List of global vertex indices for each surface element
     *   (dimension nsurf x nv, where nsurf is the number of found surface
     *   faces, and nv is the max. number of vertices per surface element)
     */
    int RegionBoundaries (IDenseMatrix &idx);

    bool PullToBoundary (const Point &p, Point &pshift, int &element,
	int &side, double sub = 0) const;
    // calculates the projection pshift of p onto the surface of the nearest
    // boundary segment of the mesh, and returns the element and side number of
    // this segment.
    // if sub > 0 then pshift is not on the surface but at a distance sub
    // below it.
    // if sub < 0 then the distance is calculated internally as 1/mus of the
    // affected element
    // return value is false if no element was found.

    bool PullToBoundary_old (const Point &p, Point &pshift, int &element,
	int &side, double sub = 0) const;
    // This contains the old algorithm which is used whenever the above fails

    void SetupNeighbourList () const;
    // sets up a list of edge-adjacent elements for each element in the mesh
    // and stores it in nbhrs

    void NodeNeighbourList (int **_nnbhrs, int ***_nbhrs) const;
    // generates a list of neighbour nodes for each node in the mesh, and
    // returns it in the 2D array `_nbhrs'. The number of neighbours for each
    // node is returned in the 1D array `_nnbhrs'. Both arrays should be
    // unassigned before the call.

    /**
     * \brief Checks if two mesh elements have a common side.
     * \param [in] el1 element index 1 (>= 0)
     * \param [in] el2 element index 2 (>= 0)
     * \param [out] sd1 pointer to variable receiving the side index of the
     *   connected side in the first element
     * \param [out] sd2 pointer to variable receiving the side index of the
     *   connected side in the second element
     * \return True if the two elements have a common side, false otherwise.
     * \note If sd1 and sd2 are left at their default values (NULL) the
     *   side indices are not returned.
     * \note If the elements are not connected, the *sd1 and *sd2 values are
     *   left unchanged.
     */
    bool ElConnected (int el1, int el2, int *sd1 = NULL, int *sd2 = NULL);

    double Size (Point *centre = 0) const;
    // returns mesh size (= half length of longest bounding box size) and mesh
    // centre if `centre' is supplied.

    void BoundingBox (Point &mmin, Point &mmax, double pad = 0.0) const;
    // returns bounding box of mesh such that for each node N
    // min[i] <= N[i] <= max[i], where i=0..1 or 0..2
    // If pad > 0, the bounding box is enlarged by this margin in all directions

    Point BndIntersect (const Point &pt1, const Point &pt2, int *el = NULL);
    // this calculates the point of intersection of the mesh boundary and the
    // straight line from pt1 to pt2
    // if el != NULL, it gets assigned the index of the intersected element.

    int CheckConsistency () const;
    // performs some tests to check the consistency of the mesh
    // returns 0 if mesh seems OK, or a nonzero errorcode if not
    // aborts on error only if FEM_DEBUG is defined

    void MarkBoundary ();
    // scans through mesh structure and marks nodes on the mesh surface
    // as boundary nodes

    /**
     * \brief Returns the mass matrix for the mesh.
     *
     * Allocates a sparse matrix and fills it by executing the Element::IntFF
     * operation (integral of product of two shape functions) over the
     * elements.
     * \return pointer to dynamically allocated mass matrix
     * \note The user is responsible for deleting the matrix after use.
     */
    RCompRowMatrix *MassMatrix () const;

    friend FELIB void AddToElMatrix (const Mesh &mesh, int el,
        RGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles nodal coefficient vector 'coeff' into element matrix 'M'
    // (only fills entries relevant for the element)

    friend FELIB void AddToElMatrix (const Mesh &mesh, int el,
        FGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles nodal coefficient vector 'coeff' into element matrix 'M'
    // (only fills entries relevant for the element)

    friend FELIB void AddToElMatrix (const Mesh &mesh, int el,
	SCGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles nodal coefficient vector 'coeff' into element matrix 'M'
    // (only fills entries relevant for the element)

    friend FELIB void AddToElMatrix (const Mesh &mesh, int el,
	CGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles nodal coefficient vector 'coeff' into element matrix 'M'
    // (only fills entries relevant for the element)

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        RGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles coefficient vector 'coeff' into system matrix 'M' using
    // the specified mode
    // For modes FF and DD 'coeff' is ignored

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        FGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles coefficient vector 'coeff' into system matrix 'M' using
    // the specified mode
    // For modes FF and DD 'coeff' is ignored

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        SCGenericSparseMatrix &M, const RVector *coeff, int mode);

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        CGenericSparseMatrix &M, const RVector *coeff, int mode);
    // As above but for complex system matrix. Modes can also contain assembly
    // commands for imaginary part

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        CGenericSparseMatrix &M, const double coeff, int mode);

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        SCGenericSparseMatrix &M, const double coeff, int mode);
    // This version adds a component with constant coefficient 'coeff' to
    // the real or imaginary part of system matrix M

    friend FELIB void AddToSysMatrix_elasticity (const Mesh &mesh,
        RGenericSparseMatrix &M, const RVector &modulus,
	const RVector &pratio);

    friend FELIB void AddToRHS_elasticity (const Mesh &mesh, RVector &rhs,
        const RVector *coeff, int mode);
    // Adds a component to an rhs vector, given the specified mode

    friend FELIB void AddToRHS_thermal_expansion (const Mesh &mesh,
        RVector &rhs, const RVector &modulus, const RVector &pratio,
	const RVector &thermal_expansion, double deltaT);

    friend FELIB void AddElasticStrainDisplacementToSysMatrix (const Mesh &mesh,
        RGenericSparseMatrix &M, double E, double nu, int *nf);

    /**
     * \brief Nodal to grid basis mapping
     *
     * Generate a sparse matrix which will map a nodal vector from rectangular
     * area [bbmin,bbmax] to a regular grid of dimensions gdim.
     * \param mesh Mesh defining the unstructured mesh basis
     * \param gdim Grid dimensions (# voxels)
     * \param bbmin min. coordinates of grid bounding box
     * \param bbmax max. coordinates of grid bounding box
     * \param elref element reference list as provided by
     *   GenerateElementPixelRef().
     * \return node->grid mapping matrix
     * \note If bbmin and/or bbmax is not defined then the mesh boundaries are
     *   substituted.
     * \note If elref is not supplied, it is generated on the fly
     * \sa GenerateElementPixelRef, NodeMapMatrix
     */
    friend FELIB RGenericSparseMatrix *GridMapMatrix (const Mesh &mesh,
	const IVector &gdim, const Point *bbmin, const Point *bbmax,
	const int *elref);

    /**
     * \brief Grid to nodal basis mapping
     *
     * Generate a sparse matrix which will map a bitmap to a nodal image
     * gdim contains the bitmap size, [bbmin,bbmax] is the physical area
     * covered by the bitmap.
     * \param mesh Mesh defining the unstructured mesh basis
     * \param gdim Grid dimensions (# voxels)
     * \param bbmin min. coordinates of grid bounding box
     * \param bbmax max. coordinates of grid bounding box
     * \param elref element reference list as provided by
     *   GenerateElementPixelRef().
     * \return grid->node mapping matrix
     * \note If bbmin and/or bbmax is not defined then the
     *   mesh boundaries are substituted. Note if the supplied BB is smaller
     *   than the mesh BB then some nodes will remain undefined after mapping
     * \note If elref is not supplied, it is generated on the fly.
     * \sa GenerateElementPixelRef, GridMapMatrix
     */
    friend FELIB RGenericSparseMatrix *NodeMapMatrix (const Mesh &mesh,
        const IVector &gdim, const Point *bbmin, const Point *bbmax,
	const int *elref);

    friend FELIB RGenericSparseMatrix *NodeMapMatrix2 (const Mesh &mesh,
        const IVector &gdim, const Point *bbmin, const Point *bbmax,
	const int *elref);

    friend FELIB RGenericSparseMatrix *Grid2LinPixMatrix (const IVector &gdim,
        const IVector &bdim, const int *elref);
    // Generate a sparse matrix which maps a bitmap of dimension gdim into
    // a linear pixel basis of dimension bdim

    friend FELIB RGenericSparseMatrix *LinPix2GridMatrix (const IVector &gdim,
        const IVector &bdim, const int *elref);

    friend FELIB RGenericSparseMatrix *CubicPix2GridMatrix (const IVector &bdim,
        const IVector &gdim, const int *elref);

    friend FELIB int *GenerateElementPixelRef (const Mesh &mesh,
        const IVector &gdim, const Point *bbmin, const Point *bbmax);
    // Generate element index list for pixels in grid defined by gdim.
    // This list is required by GridMapMatrix and NodeMapMatrix

    friend FELIB void GenerateVoxelPositions (const Mesh &mesh,
        const IVector &gdim, const Point *bbmin, const Point *bbmax,
	RDenseMatrix &pos);
    // Return the positions of the voxel vertices of a regular grid, defined
    // by (gdim,bbmin,bbmax), in dense matrix pos

    friend FELIB void SubsampleLinPixel (const RVector &gf, RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref);
    // Subsample an image gf from the raster grid of dimension gdim to a
    // (coarser) linear pixel basis of dimension bdim.
    // Both grids are assumed to have identical bounding boxes
    // elref (as returned by GenerateElementPixelRef) contains the mesh support
    // information: elref[i] < 0 indicates that pixel i has no overlap with the
    // mesh. If elref==0 then full support is assumed.

    friend FELIB void RasterLinPixel (RVector &gf, const RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref);
    // Inverse of SubsampleLinPixel: expand linear pixel basis into raster
    // image

    friend FELIB void SubsampleCubPixel (const RVector &gf, RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref);
    // Same as SubsampleLinPixel, but for cubic pixel basis

    friend FELIB void RasterCubPixel (RVector &gf, const RVector &bf,
	const IVector &gdim, const IVector &bdim, const int *elref);
    // Same as RasterLinPixel, but for cubic pixels

    // I/O
    friend FELIB std::istream& operator>> (std::istream& i, Mesh& mesh);
    friend FELIB std::ostream& operator<< (std::ostream& o, Mesh& mesh);
    void put (std::ostream &os, ParameterType p1 = PRM_MUA,
	      ParameterType p2 = PRM_KAPPA, ParameterType p3 = PRM_N);
    void WriteVtk (ostream &os, const RVector &nim);

    /**
     * \brief Write the mesh to file in Gmsh format.
     * \param os output stream
     */
    void WriteGmsh (ostream &os);

    bool trap_load_error;
    // set this to false if mesh input errors should not cause an abort
    // default is true

    Surface *Boundary() const { return boundary; }
    // The parametric (outer) surface of the mesh (or NULL if no surface)

    void SetBoundary (const Surface &_boundary);
    // Set or reset parametric mesh surface

    void PopulateNeighbourLists ();

    void InitSubdivisionSupport ();

    /**
     * \brief Refine an element by subdivision
     * \param el element index (>= 0)
     * \return Refinement level of subdivided elements
     * \note This function may recursively also refine the element's
     *   neighbours, to make sure that refinement levels between neighbours
     *   differ by no more than 1
     */
    int RefineElement (int el);

    RVector *bnd_param;
    // List of boundary parametrisations for all boundary nodes
    // (dimension priv_nbnd). Valid after Setup, and only if boundary != NULL

    int *IndexBnd2Node;
    int *IndexNode2Bnd;
    // index lists to map between absolute node numbers and boundary node
    // numbers. Valid after call to Setup
    // IndexBnd2Node[i] is the absolute node number of the i-th boundary node
    // IndexNode2Bnd[j] is the boundary index of node j, or -1 if j is not
    // a boundary node.

private:
    bool is_set_up;
    // True once Setup has been called

    int priv_ilen;
    //SymMatrix *elk, *elc, *ela;	// list of element matrices

    int priv_nbnd;
    // number of boundary nodes (valid after Setup)

    mutable int lastel_found;
    // contains the element returned by the last call to ElFind
    // The next search starts around this point

    double CalcFullSize () const;
    // returns mesh size by summing up element sizes. This is called by
    // Setup(), and by FullSize() in case Setup() was not invoked

    mutable double fullsize;
    // mesh size (area for 2D meshes, volume for 3D meshes)
    // This is set by Setup() or by the first call to FullSize()

    Surface *boundary;
    // The parametric (outer) surface of the mesh (or NULL if no surface)

    //void SetupElementMatrices (void);
    // generate element matrix lists elc and elk, to be subsequently assembled
    // into system matrices. This function must be called after any
    // manipulations to the mesh topology (eg. Extrapolate) and before system
    // matrix assembly (e.g. before a ForwardSolver is instantiated)
    // The element matrices depend only on mesh topology, not on element
    // mua or kappa (mua and kappa are added only during system matrix
    // assembly). SetupElementMatrices can be called repeatedly to reflect
    // topology changes

    struct BndIntersectParam {
	Point bbmin, bbmax;
    } *intersect_prm;
    // structure to help with ray intersection computation

    BndIntersectParam *ComputeBndIntersectParam();
    
#ifdef TOAST_PARALLEL
    static void Setup_engine (void*,int,int);
#endif
};

    int MakeNodalFreedomArray (int *&nf, int nlen, int dofnod, bool *rest);

    void CreateVoxelMesh (int ex, int ey, int ez, bool *egrid,
			  double dx, double dy, double dz, Mesh &mesh);

#endif // !__MESH_H
