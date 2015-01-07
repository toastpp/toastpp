// -*-C++-*-

#ifndef __RASTER_H
#define __RASTER_H

/**
 * \file Implements class Raster (base class for image basis mappers).
 */

#include "mathlib.h"
#include "felib.h"

// =========================================================================
/**
 * \brief Base class for mapping between mesh and an independent basis
 *   representation.
 *
 * Derived classes implement specific basis representations.
 */
class STOASTLIB Raster {
public:
    /**
     * \brief Raster constructor.
     * \param _bdim native basis dimensions (the length of this vector is 2 or
     *   3, corresponding to the mesh dimension)
     * \param _gdim high-res grid dimensions (the length of this vector is 2 or
     *   3, corresponding to the mesh dimension)
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0);

    /**
     * \brief Raster destructor.
     */
    virtual ~Raster();

    /**
     * \brief Return problem dimension.
     * \return 2 or 3, depending on problem dimension (corresponds to mesh
     *   dimension).
     */
    const int Dim() const {return dim; }

    /**
     * \brief Return native basis dimensions.
     * \return Vector of length \ref Dim, containing the basis basis
     *   dimensions.
     */
    const IVector &BDim() const { return bdim; }

    /**
     * \brief Return high-res grid dimensions.
     * \return Vector of length \ref Dim, containing the grid dimensions.
     */
    const IVector &GDim() const { return gdim; }

    /**
     * \brief Return size of the grid bounding box.
     * \return Vector of length \ref Dim, containing the physical size of the
     *   grid bounding box.
     * \note The grid bounding box corresponds to the mesh bounding box, unless
     *   an explicit bounding box was specified in the constructor.
     */
    const RVector &GSize() const { return bbsize; }

    /**
     * \brief Return the length of the grid vector representation.
     */
    inline int GLen() const { return glen; }

    /**
     * \brief Return the length of the basis vector representation.
     */
    inline int BLen() const { return blen; }

    /**
     * \brief Return the length of the solution vector representation.
     */
    inline int SLen() const { return slen; }

    /**
     * \brief Return the mesh associated with the basis.
     */
    inline const Mesh &mesh() const { return *meshptr; }

    /**
     * \brief Return the bounding box of the basis support
     * \return 2 x dim dense matrix
     */
    RDenseMatrix BoundingBox () const;

    virtual RDenseMatrix SupportArea (int idx)
    { xERROR("Not implemented"); return RDenseMatrix(); }

    virtual void Refine (int idx)
    { xERROR("Not implemented"); }

    virtual void Refine (int *idx, int nidx)
    { xERROR("Not implemented"); }

    /**
     * \brief Returns a pointer to the element reference array for all grid
     *   points.
     * \return Array of length \ref GLen() containing the element indices
     *   (>= 0) associated with each grid point. Grid points outside the
     *   domain support are indicated by value -1.
     * \sa BElref, GLen
     */
    inline int *Elref() const { return gelref; }

    /**
     * \brief Returns a pointer to the element reference array for all basis
     *   points.
     * \return Array of length \ref BLen() containing the element indices
     *   (>= 0) associated with each basis point. Basis points outside the
     *   domain support are indicated by value -1.
     * \sa Elref, BLen
     */
    inline int *BElref() const { return belref; }

    /**
     * \brief Returns the positions of all basis points in a matrix.
     * \param [out] pos matrix of basis point positions (blen x dim)
     */
    void BasisVoxelPositions (RDenseMatrix &pos) const;

    /**
     * \brief Returns the positions of all solution points in a matrix.
     * \param [out] pos matrix of solution point positions (slen x dim)
     */
    void SolutionVoxelPositions (RDenseMatrix &pos) const;

    virtual void Sample (const RVector &bvec, const IVector &grd, RVector &img)
	const;

    /**
     * \brief Value of basis function b_i at point p
     * Identical to Value_nomask, but returns zero if p is outside the mesh
     * \sa Value_nomask
     */
    virtual double Value (const Point &p, int i, bool is_solidx=true) const;

    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    virtual double Value_nomask (const Point &p, int i, bool is_solidx=true)
	const = 0;

    /**
     * \brief Value of a function, expressed by basis coefficients, at point p.
     * \param p evaluation position
     * \param coeff vector of basis coefficients
     * \param mask mask out domain exterior?
     * \return function value at p
     * \note Computes return value as sum_i b_i(p) coeff(i)
     * \note If coeff.Dim()==blen, the coefficients are taken to be in full
     *   basis expansion. If coeff.Dim()==slen, they are assumed to be in
     *   solution basis, i.e. exclude basis functions without domain overlap.
     */
    double Value (const Point &p, const RVector &coeff, bool mask=true) const;

    /**
     * \brief Map a real-valued field from mesh to grid representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in grid representation
     */
    virtual void Map_MeshToGrid (const RVector &mvec, RVector &gvec) const;

    /**
     * \brief Map a complex-valued field from mesh to grid representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in grid representation
     */
    virtual void Map_MeshToGrid (const CVector &mvec, CVector &gvec) const;

    /**
     * \brief Map a set of fields from mesh to grid representation.
     * \param [in] msol field set in mesh representation
     * \param [out] gsol field set in grid representation
     * \param [in] mapall flag mapping all/active only fields in the set
     */
    virtual void Map_MeshToGrid (const Solution &msol, Solution &gsol,
        bool mapall = false) const;

    /**
     * \brief Map a real-valued field from grid to mesh representation
     * \param [in] gvec field vector in grid representation
     * \param [out] mvec field vector in mesh representation
     */
    virtual void Map_GridToMesh (const RVector &gvec, RVector &mvec) const;

    /**
     * \brief Map a complex-valued field from grid to mesh representation
     * \param [in] gvec field vector in grid representation
     * \param [out] mvec field vector in mesh representation
     */
    virtual void Map_GridToMesh (const CVector &gvec, CVector &mvec) const;

    /**
     * \brief Map a set of fields from grid to mesh representation.
     * \param [in] gsol field set in grid representation
     * \param [out] msol field set in mesh representation
     * \param [in] mapall flag mapping all/active only fields in the set
     */
    virtual void Map_GridToMesh (const Solution &gsol, Solution &msol,
	bool mapall = false) const;

    /**
     * \brief Map a real-valued field from grid to basis representation.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     */
    virtual void Map_GridToBasis (const RVector &gvec, RVector &bvec) const =0;

    /**
     * \brief Map a complex-valued field from grid to basis representation.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     */
    virtual void Map_GridToBasis (const CVector &gvec, CVector &bvec) const =0;

    /**
     * \brief Map a set of fields from grid to basis representation.
     * \param [in] gsol field set in grid representation
     * \param [out] bsol field set in basis representation
     * \param [in] mapall flag mapping all/active only fields in the set
     */
    virtual void Map_GridToBasis (const Solution &gsol, Solution &bsol,
	bool mapall = false) const;

    /**
     * \brief Map a real-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     */
    virtual void Map_BasisToGrid (const RVector &bvec, RVector &gvec) const =0;

    /**
     * \brief Map a complex-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     */
    virtual void Map_BasisToGrid (const CVector &bvec, CVector &gvec) const =0;

    /**
     * \brief Map a set of fields from basis to grid representation.
     * \param [in] bsol field set in basis representation
     * \param [out] gsol field set in grid representation
     * \param [in] mapall flag mapping all/active only fields in the set
     */
    virtual void Map_BasisToGrid (const Solution &bsol, Solution &gsol,
	bool mapall = false) const;

    /**
     * \brief Map a real-valued field from mesh to basis representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in basis representation
     * \default Maps field from mesh to grid representation, and then from grid
     *   to basis representation.
     */
    virtual void Map_MeshToBasis (const RVector &mvec, RVector &bvec) const;

    /**
     * \brief Map a complex-valued field from mesh to basis representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in basis representation
     * \default Maps field from mesh to grid representation, and then from grid
     *   to basis representation.
     */
    virtual void Map_MeshToBasis (const CVector &mvec, CVector &bvec) const;

    /**
     * \brief Map a real-valued field from basis to mesh representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] mvec field vector in mesh representation
     * \default Maps field from basis to grid representation, and then from
     *    grid to mesh representation.
     */
    virtual void Map_BasisToMesh (const RVector &bvec, RVector &mvec) const;

    /**
     * \brief Map a complex-valued field from basis to mesh representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] mvec field vector in mesh representation
     * \default Maps field from basis to grid representation, and then from
     *    grid to mesh representation.
     */
    virtual void Map_BasisToMesh (const CVector &bvec, CVector &mvec) const;

    /**
     * \brief Map a real-valued field from basis to solution representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] svec field vector in solution representation
     * \default Uses the sol2basis mask vector to eliminate basis functions
     *   outside the domain support.
     * \note Derived classes should either initialise the sol2basis vector, or
     *   overload this method. If neither is done, then a default mask based
     *   on the position of each basis function (but disregarding support area)
     *   is used.
     */
    virtual void Map_BasisToSol (const RVector &bvec, RVector &svec) const;

    /**
     * \brief Map a complex-valued field from basis to solution representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] svec field vector in solution representation
     * \default Uses the sol2basis mask vector to eliminate basis functions
     *   outside the domain support.
     * \note Derived classes should either initialise the sol2basis vector, or
     *   overload this method. If neither is done, then a default mask based
     *   on the position of each basis function (but disregarding support area)
     *   is used.
     */
    virtual void Map_BasisToSol (const CVector &bvec, CVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to basis representation.
     * \param [in] svec field vector in solution representation
     * \param [out] bvec field vector in basis representation
     * \default Uses the sol2basis mask vector to copy voxels within the
     *   support of the domain. External voxels are set to zero.
     * \note Derived classes should either initialise the sol2basis vector, or
     *   overload this method. If neither is done, then a default mask based
     *   on the position of each basis function (but disregarding support area)
     *   is used.
     */
    virtual void Map_SolToBasis (const RVector &svec, RVector &bvec) const;

    /**
     * \brief Map a complex-valued field from solution to basis representation.
     * \param [in] svec field vector in solution representation
     * \param [out] bvec field vector in basis representation
     * \default Uses the sol2basis mask vector to copy voxels within the
     *   support of the domain. External voxels are set to zero.
     * \note Derived classes should either initialise the sol2basis vector, or
     *   overload this method. If neither is done, then a default mask based
     *   on the position of each basis function (but disregarding support area)
     *   is used.
     */
    virtual void Map_SolToBasis (const CVector &svec, CVector &bvec) const;

    virtual void Map_SolToBasis (const Solution &ssol, Solution &bsol,
	bool mapall = false) const;

    /**
     * \brief Map a real-valued field from solution to grid representation.
     * \param [in] svec field vector in solution representation
     * \param [out] gvec field vector in grid representation
     * \default Calls Map_SolToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToGrid to map to grid representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_SolToGrid (const RVector &svec, RVector &gvec) const;

    /**
     * \brief Map a complex-valued field from solution to grid representation.
     * \param [in] svec field vector in solution representation
     * \param [out] gvec field vector in grid representation
     * \default Calls Map_SolToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToGrid to map to grid representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_SolToGrid (const CVector &svec, CVector &gvec) const;

    /**
     * \brief Map a real-valued field from grid to solution representation.
     * \param [in] gvec field vector in grid representation
     * \param [out] svec field vector in solution representation
     * \default Calls Map_GridToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToSol to map to solution representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_GridToSol (const RVector &gvec, RVector &svec) const;

    /**
     * \brief Map a complex-valued field from grid to solution representation.
     * \param [in] gvec field vector in grid representation
     * \param [out] svec field vector in solution representation
     * \default Calls Map_GridToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToSol to map to solution representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_GridToSol (const CVector &gvec, CVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to mesh representation.
     * \param [in] svec field vector in solution representation
     * \param [out] mvec field vector in mesh representation
     * \default Calls Map_SolToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToMesh to map to mesh representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_SolToMesh (const RVector &svec, RVector &mvec) const;

    /**
     * \brief Map a complex-valued field from solution to mesh representation.
     * \param [in] svec field vector in solution representation
     * \param [out] mvec field vector in mesh representation
     * \default Calls Map_SolToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToMesh to map to mesh representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_SolToMesh (const CVector &svec, CVector &mvec) const;

    virtual void Map_SolToMesh (const Solution &ssol, Solution &msol,
	bool mapall = false ) const;

    // Mapping from active parameter vector to solution
    virtual void Map_ActiveSolToMesh (const RVector &asol, Solution &msol) const;

    /**
     * \brief Map a real-valued field from mesh to solution representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] svec field vector in solution representation
     * \default Calls Map_MeshToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToSol to map to solution representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_MeshToSol (const RVector &mvec, RVector &svec) const;

    /**
     * \brief Map a complex-valued field from mesh to solution representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] svec field vector in solution representation
     * \default Calls Map_MeshToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToSol to map to solution representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_MeshToSol (const CVector &mvec, CVector &svec) const;

    /**
     * \brief Map a set of fields from mesh to solution representation.
     * \param [in] mvec field set in mesh representation
     * \param [out] svec field set in solution representation
     * \param [in] mapall flag mapping all/active only fields in the set
     * \default Calls Map_MeshToBasis to create an intermediate basis vector,
     *   followed by Map_BasisToSol to map to solution representation.
     * \note Derived classes may overload this method for a direct mapping.
     */
    virtual void Map_MeshToSol (const Solution &msol, Solution &ssol,
        bool mapall = false) const;

    /**
     * \brief Return the mesh->grid transformation matrix.
     */
    inline const RGenericSparseMatrix &Mesh2GridMatrix()  const
    {   return *B; }

    /**
     * \brief Return the mesh->basis transformation matrix.
     * \note Not all basis types may implement the basis->mesh transformation
     *   with an explicit matrix operator. For these classes, calling this
     *   method will throw an error.
     */
    virtual const RGenericSparseMatrix &Mesh2BasisMatrix() const
    {   ERROR_UNDEF;
	static RCompRowMatrix dummy;
	return dummy; }

    /**
     * \brief Return the basis->mesh transformation matrix.
     * \note Not all basis types may implement the basis->mesh transformation
     *   with an explicit matrix operator. For these classes, calling this
     *   method will throw an error.
     */
    virtual const RGenericSparseMatrix &Basis2MeshMatrix() const
    {   ERROR_UNDEF;
	static RCompRowMatrix dummy;
	return dummy; }

    /**
     * \brief Map a linear basis index into a linear solution index.
     * \param basisidx basis index (>= 0)
     * \return Linear solution index (>= 0) or -1 if the specified voxel is not
     *   mapped into the solution.
     * \note For a basis of dimensions nx x ny x nz, the linear basis index for
     *   voxel (i,j,k) is given by idx = i + nx*j + nx*ny*k.
     * \note The returned solution index has all masked voxels removed.
     */
    virtual int GetSolIdx (int basisidx) const;

    /**
     * \brief Map grid indices into linear solution index.
     * \param crd Vector of dimension 2 or 3 with grid indices in each
     *   dimension (>= 0)
     *\ return Linear solution index (>= 0) or -1 if the specified voxel is not
     *   mapped into the solution.
     */
    virtual int GetSolIdx (const IVector &crd) const;

    /**
     * \brief Map a linear solution index into a linear basis index.
     * \param solidx solution index (>= 0)
     * \return Linear basis index (>= 0).
     * \note For a basis of dimensions nx x ny x nz, the linear basis index for
     *   voxel (i,j,k) is given by idx = i + nx*j + nx*ny*k.
     * \note The returned solution index has all masked voxels removed.
     */
    virtual int GetBasisIdx (int solidx) const;

    /**
     * \brief Map a linear basis index into the xy or xyz indices of the
     *   corresponding voxel.
     * \param [in] basisidx linear basis index (>= 0)
     * \param [out] crd Vector of voxel coordinates
     */
    virtual void GetBasisIndices (int basisidx, IVector &crd) const;

    /**
     * \brief Map a linear basis index into the xy or xyz indices of the
     *   corresponding voxel.
     * \param basisidx linear basis index (>= 0)
     * \return Vector of voxel coordinates.
     * \note This version is not threadsafe
     */
    //virtual IVector &GetBasisIndices (int basisidx) const;

    /**
     * \brief Map a linear grid index into the xy or xyz indices of the
     *   corresponding grid voxel.
     * \param grididx linear grid index (>= 0)
     * \return Vector of grid voxel coordinates.
     */
    virtual IVector &GetGridIndices (int grididx) const;

    /**
     * \brief Returns the corresponding solution vector index for a given
     *   basis vector index.
     * \param i basis vector index (0 .. blen-1)
     * \return Solution vector index (0 .. slen-1), or -1 if basis index is
     *   not within the support of the domain.
     */
    inline int Basis2Sol (int i) const { return basis2sol[i]; }

        /**
     * \brief Returns the corresponding basis vector index for a given
     *   solution vector index.
     * \param j solution vector index (0 .. slen-1)
     * \return Basis vector index (0 .. blen-1)
     */
    inline int Sol2Basis (int j) const { return sol2basis[j]; }

    void NeighbourGraph (idxtype *&rowptr, idxtype *&colidx, int &nzero) const;
    void NeighbourGraph (ICompRowMatrix &NG) const;
    IVector NeighbourShift (const ICompRowMatrix &NG, int i, int j) const;
    RVector RNeighbourShift (const ICompRowMatrix &NG, int i, int j) const;

    RVector SmoothImage (const RVector &x, double sd) const;

    RVector * ImageJet(const RVector &x, double sd, bool *iflags) const;

    RVector ImageGradient(const RVector &x, double sd) const;

protected:
    int dim;             ///< problem dimension (2 or 3)
    IVector bdim;        ///< native basis dimensions
    IVector gdim;        ///< high-res grid dimension
    int blen;            ///< length of basis vector (full bb)
    int glen;            ///< length of grid vector (full bb)
    int slen;            ///< length of solution vector (excluding mask)
    RVector bbsize;      ///< physical extents of grid bounding box
    Mesh *meshptr;       ///< mesh pointer
    Point bbmin, bbmax;  ///< grid bounding box
    int *belref;         ///< basis->element reference list (length blen)
    int *gelref;         ///< grid->element reference list (length glen)

    IVector paddim;      ///< dimensions padded to a power of 2
    int padlen;          ///< vector length of padded grid

    /**
     * Vector of length blen with basis support weights.
     * Range: 0 (no mesh support) to 1 (full mesh support)
     * Currently only binary states (0 or 1) are supported.
     */
    RVector bsupport;

    /**
     * \brief basis->solution permutation vector.
     *
     * This vector is of length blen (basis dimension). Each entry contains
     * the index of the corresponding element in the solution vector (range
     * 0 .. slen-1), or -1 if the element is outside the support of the domain.
     * \note By default, this vector is  constructed from the element support
     *   array for the basis voxels. This may need to be redefined by
     *   derived classes to accommodate their basis function support areas.
     */
    IVector basis2sol;

    /**
     * \brief Solution mask vector.
     *
     * This vector is of length slen (solution dimension). Each entry contains
     * the index of the corresponding element in the basis vector (range
     * 0 .. blen-1). Used to map between basis representation (full bounding
     * box) and solution representation (supported voxels only).
     * \note By default, this vector is constructed from the element support
     *   array for the basis voxels. This may need to be redefined by
     *   derived classes to accommodate their basis function support areas.
     */
    IVector sol2basis;

    RGenericSparseMatrix *B;  ///< transformation matrix mesh->grid
    RGenericSparseMatrix *BI; ///< transformation matrix grid->mesh
};

#endif // !__RASTER_H
