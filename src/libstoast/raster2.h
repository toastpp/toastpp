// -*-C++-*-

#ifndef __RASTER2_H
#define __RASTER2_H

/**
 * \file Defines class Raster2. This is the base class for a group of
 *   basis mapper classes which use a least-squares optimisation method
 *   for mapping between mesh basis and their own basis representation.
 *   See papers/basismap_new.tex for details
 */

#include "raster.h"

// =========================================================================
/**
 * \brief Base class for least-squares basis mappers
 */
class STOASTLIB Raster2: public Raster {
public:
    /**
     * \brief Raster factory. Creates a Raster object of the desired type
     *   and returns a pointer to it.
     */
    template<class T>
    static Raster2 *Create (const IVector &_bdim, const IVector &_gdim,
        Mesh *mesh, RDenseMatrix *bb=0, double _map_tol=1e-10)
    {
        Raster2 *raster = new T(_bdim, _gdim, mesh, bb, _map_tol);
	raster->Init();
	return raster;
    }

    /**
     * \brief Basis destructor
     */
    virtual ~Raster2 ();

    virtual RVector Gradient_nomask (const Point &p, int i, bool is_solidx=true)
	const
    { xERROR("Not implemented"); return RVector(); }
    
    /**
     * \brief Map a real-valued coefficient vector from intermediate grid to
     *   intrinsic basis.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     * \note Raster2-based classes don't support intermediate grids, so
     *   this method simply copies gvec to bvec
     */
    void Map_GridToBasis (const RVector &gvec, RVector &bvec) const;

    /**
     * \brief Map a complex-valued coefficient vector from intermediate grid
     *   to intrinsic basis.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     * \note Raster2-based classes don't support intermediate grids, so
     *   this method simply copies gvec to bvec
     */
    void Map_GridToBasis (const CVector &gvec, CVector &bvec) const;

    /**
     * \brief Map a real-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     * \note Raster2-based classes don't support intermediate grids, so
     *   this method simply copies bvec to gvec
     */
    void Map_BasisToGrid (const RVector &bvec, RVector &gvec) const;

    /**
     * \brief Map a complex-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     * \note Raster2-based classes don't support intermediate grids, so
     *   this method simply copies bvec to gvec
     */
    void Map_BasisToGrid (const CVector &bvec, CVector &gvec) const;

    /**
     * \brief Map a real-valued field from mesh to basis representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in basis representation
     */
    void Map_MeshToBasis (const RVector &mvec, RVector &bvec) const;
    void Map_MeshToBasis (const CVector &mvec, CVector &bvec) const;

    /**
     * \brief Map a real-valued field from basis to mesh representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] mvec field vector in mesh representation
     */
    void Map_BasisToMesh (const RVector &bvec, RVector &mvec) const;
    void Map_BasisToMesh (const CVector &bvec, CVector &mvec) const;

    /**
     * \brief Map a real-valued field from basis to solution representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] svec field vector in solution representation
     * \note The 'basis' representation contains the full bounding box,
     *   while the 'solution' representation has the masked voxels removed.
     */
    void Map_BasisToSol (const RVector &bvec, RVector &svec) const;
    void Map_BasisToSol (const CVector &bvec, CVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to basis representation.
     * \param [in] svec field vector in solution representation
     * \param [out] bvec field vector in basis representation
     * \note The 'basis' representation contains the full bounding box,
     *   while the 'solution' representation has the masked voxels removed.
     * \note The masked voxels in bvec are set to zero.
     */
    void Map_SolToBasis (const RVector &svec, RVector &bvec) const;
    void Map_SolToBasis (const CVector &svec, CVector &bvec) const;

    /**
     * \brief Map a real-valued field from mesh to solution representation.
     *   representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] svec field vector in solution representation
     */
    void Map_MeshToSol (const RVector &mvec, RVector &svec) const;
    void Map_MeshToSol (const CVector &mvec, CVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to mesh representation.
     * \param [in] svec field vector in solution representation
     * \param [out] mvec field vector in mesh representation
     */
    void Map_SolToMesh (const RVector &svec, RVector &mvec) const;
    void Map_SolToMesh (const CVector &svec, CVector &mvec) const;

    void Map_MeshToGrid (const RVector &mvec, RVector &gvec) const
    { Map_MeshToBasis (mvec, gvec); }

    void Map_MeshToGrid (const CVector &mvec, CVector &gvec) const
    { Map_MeshToBasis (mvec, gvec); }

    void Map_GridToMesh (const RVector &gvec, RVector &mvec) const
    { Map_BasisToMesh (gvec, mvec); }

    void Map_GridToMesh (const CVector &gvec, CVector &mvec) const
    { Map_BasisToMesh (gvec, mvec); }

    /**
     * \brief Single-element system matrix assembly
     *
     * Assemble single-element contribution for element "el" into global
     * system matrix M, where coefficients (where applicable) are given in
     * pixel basis.
     * \param el element index (>= 0)
     * \param M global system matrix
     * \param pxcoeff pointer to coefficient vector in pixel basis (only
     *   required for integration modes involving a parameter distribution)
     * \param mode integration type index
     */
    virtual void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const = 0;

    const RCompRowMatrix *GetBuu() const { return Buu; }
    const RCompRowMatrix *GetBvv() const { return Bvv; }
    const RCompRowMatrix *GetBuv() const { return Buv; }
    const RCompRowMatrix *GetBvw (const IVector &wdim);

    const RCompRowMatrix *GetDuu () const;
    const RCompRowMatrix *GetDvv () const;
    const RCompRowMatrix *GetDuv () const;
    
    virtual const RCompRowMatrix *GetGuv() const { return 0; }
    virtual const RCompRowMatrix *GetGvv() const { return 0; }

protected:
    /**
     * \brief Basis constructor.
     * \param _bdim basis dimensions (the length of this vector is 2 or 3,
     *   corresponding to the mesh dimension).
     * \param _gdim Not used by this basis type. Must be identical to _bdim.
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \param _map_tol Tolerance limit for least-squares solver in
     *   Map_MeshToBasis and Map_BasisToMesh
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the
     *   basis instance.
     */
    Raster2 (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
       RDenseMatrix *bb=0, double _map_tol=1e-10);

    /**
     * \brief Creates the mass matrix for the class's own basis
     *   representation (Bvv)
     */
    virtual RCompRowMatrix *CreateBvv () const = 0;

    /**
     * \brief Creates the mass matrix for the mixed bases (Buv)
     */
    virtual RCompRowMatrix *CreateBuv () const = 0;

    /**
     * \brief Creates mass matrix for mixed grid bases (Bvw)
     */
    virtual RCompRowMatrix *CreateBvw (const IVector &wdim) const
    { xERROR("Not implemented"); return 0; }

    RCompRowMatrix *CreateDuu () const;

    virtual RCompRowMatrix *CreateDvv () const
    { xERROR("Not implemented"); return 0; }
    
    virtual RCompRowMatrix *CreateDuv () const
    { xERROR("Not implemented"); return 0; }

    /**
     * \brief Creates the matrix for mapping from (rectangular) basis
     *   to (sparse) solution
     */
    RCompRowMatrix *CreateSolMapper () const;

    double map_tol;         ///< tolerance limit for least squares basis map

    RCompRowMatrix *Buu;    ///< mesh mass matrix for mapping
    RCompRowMatrix *Bvv;    ///< pixel mass matrix for mapping
    RCompRowMatrix *Buv;    ///< mixed-basis mapping matrix
    RCompRowMatrix *Bvw;    ///< basis-pixel mixed mapping matrix (on demand)
    IVector Bvw_dim;        ///< dimension for Bvw pixel image

    RPrecon_IC *Buu_precon; ///< preconditioner for Buu
    RPrecon_IC *Bvv_precon; ///< preconditioner for Bvv

    RCompRowMatrix *Buu_Cholesky_L;
    RVector        *Buu_Cholesky_d;
    RCompRowMatrix *Bvv_Cholesky_L;
    RVector        *Bvv_Cholesky_d;

    mutable RCompRowMatrix *Duu; ///< mesh stiffness matrix
    mutable RCompRowMatrix *Dvv; ///< basis stiffness matrix
    mutable RCompRowMatrix *Duv; ///< mixed-basis stiffness matrix
    
    RCompRowMatrix *D;      ///< map basis->solution (limited support)

public:
    /**
     * \brief Polymorphic class initialisation, called by Create().
     */
    virtual void Init();
};

#endif // !__RASTER2_H
