// -*-C++-*-

#ifndef __RASTER_PX2_H
#define __RASTER_PX2_H

/**
 * \file Implements class Raster_Pixel2 (simple pixel/voxel based image
 *   representation).
 *
 * This is a replacement for Raster_Pixel. It uses a more accurate
 * mapping algorithm, but isn't yet fully functional.
 */

#include "raster.h"

void Tetsplit (Mesh *mesh, int el, int cut_orient, double cut_pos);
// exported temporarily

// =========================================================================
/**
 * \brief Pixel (or voxel) basis representation.
 */
class STOASTLIB Raster_Pixel2: public Raster {
public:
    /**
     * \brief Pixel basis constructor.
     * \param _bdim basis dimensions (the length of this vector is 2 or 3,
     *   corresponding to the mesh dimension).
     * \param _gdim Not used by this basis type. Must be identical to _bdim.
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \param _map_tol Tolerance limit for least-squares solver in
     *   Map_MeshToBasis and Map_BasisToMesh
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster_Pixel2 (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0, double _map_tol=1e-10);

    /**
     * \brief Pixel basis destructor.
     */
    ~Raster_Pixel2 ();

    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

    /**
     * \brief Map a real-valued field from grid to basis representation.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     * \note If grid and basis dimensions are identical, this operator just
     *   copies gvec into bvec. Otherwise it uses a linear interpolation
     *   scheme.
     */
    void Map_GridToBasis (const RVector &gvec, RVector &bvec) const;

    /**
     * \brief Map a complex-valued field from grid to basis representation.
     * \param [in] gvec field in grid representation
     * \param [out] bvec field in native basis representation
     */
    void Map_GridToBasis (const CVector &gvec, CVector &bvec) const;

    /**
     * \brief Map a real-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     * \note If grid and basis dimensions are identical, this operator just
     *   copies bvec into gvec. Otherwise it uses a linear interpolation
     *   scheme.
     */
    void Map_BasisToGrid (const RVector &bvec, RVector &gvec) const;

    /**
     * \brief Map a complex-valued field from basis to grid representation.
     * \param [in] bvec field in native basis representation
     * \param [out] gvec field in grid representation
     * \note If grid and basis dimensions are identical, this operator just
     *   copies bvec into gvec. Otherwise it uses a linear interpolation
     *   scheme.
     */
    void Map_BasisToGrid (const CVector &bvec, CVector &gvec) const;

    /**
     * \brief Map a real-valued field from mesh to basis representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] bvec field vector in basis representation
     */
    void Map_MeshToBasis (const RVector &mvec, RVector &bvec) const;

    /**
     * \brief Map a real-valued field from basis to mesh representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] mvec field vector in mesh representation
     */
    void Map_BasisToMesh (const RVector &bvec, RVector &mvec) const;

    /**
     * \brief Map a real-valued field from basis to solution representation.
     * \param [in] bvec field vector in basis representation
     * \param [out] svec field vector in solution representation
     * \note The 'basis' representation contains the full bounding box, while
     *   the 'solution' representation has the masked voxels removed.
     */
    void Map_BasisToSol (const RVector &bvec, RVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to basis representation.
     * \param [in] svec field vector in solution representation
     * \param [out] bvec field vector in basis representation
     * \note The 'basis' representation contains the full bounding box, while
     *   the 'solution' representation has the masked voxels removed.
     * \note The masked voxels in bvec are set to zero.
     */
    void Map_SolToBasis (const RVector &svec, RVector &bvec) const;

    /**
     * \brief Map a real-valued field from mesh to solution representation.
     *   representation.
     * \param [in] mvec field vector in mesh representation
     * \param [out] svec field vector in solution representation
     */
    void Map_MeshToSol (const RVector &mvec, RVector &svec) const;

    /**
     * \brief Map a real-valued field from solution to mesh representation.
     * \param [in] svec field vector in solution representation
     * \param [out] mvec field vector in mesh representation
     */
    void Map_SolToMesh (const RVector &svec, RVector &mvec) const;

    /**
     * \brief Single-element system matrix assembly
     *
     * Assemble single-element contribution for element "el" into global
     * system matrix M, where coefficients (where applicable) are given in pixel
     * basis.
     * \param el element index (>= 0)
     * \param M global system matrix
     * \param pxcoeff pointer to coefficient vector in pixel basis (only
     *   required for integration modes involving a parameter distribution)
     * \param mode integration type index
     */
    void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

private:
    bool grid_is_basis;       ///< grid and basis dimensions identical?
    RCompRowMatrix *D;        ///< map basis->solution (limited support)
    RCompRowMatrix *Buu;      ///< mesh mass matrix for mapping
    RCompRowMatrix *Bvv;      ///< pixel mass matrix for mapping
    RCompRowMatrix *Buv;      ///< mixed mapping matrix

    RPrecon_IC *Buu_precon;   ///< preconditioner for basis->mesh map
    RPrecon_IC *Bvv_precon;   ///< preconditioner for mesh->basis map
    double map_tol;           ///< tolerance limit for least squares basis map

    /**
     * \brief Returns pointer to mixed mapping matrix
     */
    RCompRowMatrix *CreateMixedMassmat () const;
    RCompRowMatrix *CreateMixedMassmat_tri () const;
    RCompRowMatrix *CreateMixedMassmat_tet4 () const;

    RCompRowMatrix *CreatePixelMassmat () const;

    void AddToElMatrix_tri (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

    void AddToElMatrix_tet (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

    int SutherlandHodgman (int el, int xgrid, int ygrid, Point *clip_poly,
        int npoly) const;

#if THREAD_LEVEL==2
    friend void CreateMixedMassmat_tri_pass2_engine (task_data *td);
#endif
};
    
#endif // !__RASTER_PX2_H
