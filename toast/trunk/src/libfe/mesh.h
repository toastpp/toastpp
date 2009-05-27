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

// rhs assembly modes
#define RHS_P           0
#define RHS_PF          1
#define RHS_BNDPF       2

// ==========================================================================
// Prototypes

FELIB RGenericSparseMatrix *Grid2LinPixMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref = 0);

FELIB RGenericSparseMatrix *LinPix2GridMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref = 0);

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
    ParameterList plist; // list of mesh parameters (per node)

    int **nbhrs;
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

    virtual ~Mesh ();
    // destructor

    void Setup();
    // call this once the node list and element list have been assigned to
    // let the mesh perform its setup stuff

    void Copy (const Mesh &mesh);
    // copy `mesh' to *this
  
    int Dimension() const 
    { return (elen() ? elist[0]->Dimension() : 0); }
    // mesh dimensionality

    int nlen() const { return nlist.Len(); }
    // returns size of node list

    int elen() const { return elist.Len(); }
    // returns size of element list

    int ilen() const { return priv_ilen; }
    // number of nodes excluding Dirichlet nodes

    int nbnd() const { return priv_nbnd; }
    // number of boundary nodes (valid after Setup)

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

    double FullSize () const;
    // returns the total area or volume of the mesh (= sum of all elements)

    int ElFind (const Point &pt) const;
    // returns the number of the element which contains the global coordinate
    // 'pt', or -1 if no element is found

    int NdFind (const Point &pt, double &dist) const;
    // returns the number of the node closest to pt and its distance

    void Reorder (int *nindex, int *iindex);
    // reorder node list such that N_new[i] = N_old[nindex[i]]

    double ElDist (int el1, int el2) const;
    // returns the distance between the centres of elements el1 and el2

    virtual void ScaleMesh (double scale);
    // rescale the complete mesh

    virtual void ScaleMesh (const RVector &scale);
    // Anisotropic scaling. 'scale' must have the same dimension as the mesh

    double ParamAverage (const ParameterType prmtp) const;
    // returns the average of parameter `prmtp' over the entire mesh

    int MaxNodeDiff (int mode=BW_AUTO) const;
    // calculates the maximum node number difference within a single element
    // of the mesh, either using all nodes (mode=BW_TOTAL) or only internal
    // nodes (mode=BW_INTERNAL). By default (mode=BW_AUTO) it uses all nodes
    // for Robin meshes, and internal nodes for all other meshes.
    // To get the semi-bandwith of the system matrix including the diagonal,
    // add 1 to the return value of MaxNodeDiff

    Point NeighbourBarycentre (int node);
    // returns the barycentre position of the neighbours of 'node'.

    void SparseRowStructure (int *&rowptr, int *&colidx, int &nzero) const;
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

    void SetupNeighbourList ();
    // sets up a list of edge-adjacent elements for each element in the mesh
    // and stores it in nbhrs

    void NodeNeighbourList (int **_nnbhrs, int ***_nbhrs);
    // generates a list of neighbour nodes for each node in the mesh, and
    // returns it in the 2D array `_nbhrs'. The number of neighbours for each
    // node is returned in the 1D array `_nnbhrs'. Both arrays should be
    // unassigned before the call.

    double Size (Point *centre = 0) const;
    // returns mesh size (= half length of longest bounding box size) and mesh
    // centre if `centre' is supplied.

    void BoundingBox (Point &mmin, Point &mmax) const;
    // returns bounding box of mesh such that for each node N
    // min[i] <= N[i] <= max[i], where i=0..1 or 0..2

    Point BndIntersect (const Point &pt1, const Point &pt2);
    // this calculates the point of intersection of the mesh boundary and the
    // straight line from pt1 to pt2

    void ResetCoeff_homog (ParameterType prmtp, double val)
	{ plist.SetParam (prmtp, val); }
    // resets the coefficient `prmtp' (as defined in param.h) throughout the
    // mesh to homogeneous value `val'

    void ResetCoeff_region (ParameterType prmtp, double val, int region);
    // Resets parameter type `prmtp' to `val' for all nodes belonging to
    // region `region'

    void ResetCoeff_sqrt (ParameterType prmtp, double cnt, double bnd);
    // resets the coefficient 'prmtp' (as defined in param.h) throughout the
    // mesh, as a function of square root of distance from centre, from value
    // 'cnt' at the centre to value 'bnd' at the boundary

    int CheckConsistency () const;
    // performs some tests to check the consistency of the mesh
    // returns 0 if mesh seems OK, or a nonzero errorcode if not
    // aborts on error only if FEM_DEBUG is defined

    void MarkBoundary ();
    // scans through mesh structure and marks nodes on the mesh surface
    // as boundary nodes

    friend FELIB void AddToElMatrix (const Mesh &mesh, int el,
	CGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles nodal coefficient vector 'coeff' into element matrix 'M'
    // (only fills entries relevant for the element)

    friend FELIB void AddToSysMatrix (const Mesh &mesh,
        RGenericSparseMatrix &M, const RVector *coeff, int mode);
    // Assembles coefficient vector 'coeff' into system matrix 'M' using
    // the specified mode
    // For modes FF and DD 'coeff' is ignored

    friend FELIB void AddToSysMatrix (const Mesh &mesh, CGenericSparseMatrix &M,
        const RVector *coeff, int mode);
    // As above but for complex system matrix. Modes can also contain assembly
    // commands for imaginary part

    friend FELIB void AddToSysMatrix (const Mesh &mesh, CGenericSparseMatrix &M,
        const double coeff, int mode);
    // This version adds a component with constant coefficient 'coeff' to
    // the real or imaginary part of system matrix M

    friend FELIB void AddToSysMatrix_elasticity (const Mesh &mesh,
        RGenericSparseMatrix &M, const RVector &modulus,
	const RVector &pratio);

    friend void AddToRHS_elasticity (const Mesh &mesh, RVector &rhs,
        const RVector *coeff, int mode);
    // Adds a component to an rhs vector, given the specified mode

    friend void AddToRHS_thermal_expansion (const Mesh &mesh, RVector &rhs,
        const RVector &modulus, const RVector &pratio,
	const RVector &thermal_expansion, double deltaT);

    friend void AddElasticStrainDisplacementToSysMatrix (const Mesh &mesh,
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

    friend RGenericSparseMatrix *NodeMapMatrix2 (const Mesh &mesh,
        const IVector &gdim, const Point *bbmin = 0, const Point *bbmax = 0,
	const int *elref = 0);

    friend FELIB RGenericSparseMatrix *Grid2LinPixMatrix (const IVector &gdim,
        const IVector &bdim, const int *elref);
    // Generate a sparse matrix which maps a bitmap of dimension gdim into
    // a linear pixel basis of dimension bdim

    friend FELIB RGenericSparseMatrix *LinPix2GridMatrix (const IVector &gdim,
        const IVector &bdim, const int *elref);

    friend RGenericSparseMatrix *CubicPix2GridMatrix (const IVector &bdim,
        const IVector &gdim, const int *elref = 0);

    friend FELIB int *GenerateElementPixelRef (const Mesh &mesh, const IVector &gdim,
        const Point *bbmin, const Point *bbmax);
    // Generate element index list for pixels in grid defined by gdim.
    // This list is required by GridMapMatrix and NodeMapMatrix

    friend FELIB void SubsampleLinPixel (const RVector &gf, RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref);
    // Subsample an image gf from the raster grid of dimension gdim to a
    // (coarser) linear pixel basis of dimension bdim.
    // Both grids are assumed to have identical bounding boxes
    // elref (as returned by GenerateElementPixelRef) contains the mesh support
    // information: elref[i] < 0 indicates that pixel i has no overlap with the
    // mesh. If elref==0 then full support is assumed.

    friend void RasterLinPixel (RVector &gf, const RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref = 0);
    // Inverse of SubsampleLinPixel: expand linear pixel basis into raster
    // image

    friend void SubsampleCubPixel (const RVector &gf, RVector &bf,
        const IVector &gdim, const IVector &bdim, const int *elref = 0);
    // Same as SubsampleLinPixel, but for cubic pixel basis

    friend void RasterCubPixel (RVector &gf, const RVector &bf,
	const IVector &gdim, const IVector &bdim, const int *elref = 0);
    // Same as RasterLinPixel, but for cubic pixels

    // I/O
	friend FELIB std::istream& operator>> (std::istream& i, Mesh& mesh);
	friend FELIB std::ostream& operator<< (std::ostream& o, Mesh& mesh);
	void put (std::ostream &os, ParameterType p1 = PRM_MUA,
	ParameterType p2 = PRM_KAPPA, ParameterType p3 = PRM_N);
  
    bool trap_load_error;
    // set this to false if mesh input errors should not cause an abort
    // default is true

    Surface *Boundary() const { return boundary; }
    // The parametric (outer) surface of the mesh (or NULL if no surface)

    void SetBoundary (const Surface &_boundary);
    // Set or reset parametric mesh surface

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

#ifdef TOAST_PARALLEL
    static void Setup_engine (void*,int,int);
#endif
};

    int MakeNodalFreedomArray (int *&nf, int nlen, int dofnod, bool *rest);

    void CreateVoxelMesh (int ex, int ey, int ez, bool *egrid,
			  double dx, double dy, double dz, Mesh &mesh);

#endif // !__MESH_H
