// -*-C++-*-

#ifndef __RASTER_PX2_H
#define __RASTER_PX2_H

/**
 * \file Implements class Raster_Pixel2 (simple bi/tri-linear pixel/voxel
 *   based image representation).
 *
 * This is a replacement for Raster_Pixel. It uses a more accurate
 * mapping algorithm, but isn't yet fully functional.
 */

#include "raster2.h"

void Tetsplit (Mesh *mesh, int el, int cut_orient, double cut_pos);
// exported temporarily

// =========================================================================
/**
 * \brief Pixel (or voxel) basis representation.
 */
class STOASTLIB Raster_Pixel2: public Raster2 {
    friend class Raster2;
    // to make constructor accessible to Raster2 factory

public:
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
    void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

protected:
    /**
     * \brief Pixel basis constructor.
     * \param _bdim basis dimensions (the length of this vector is 2 or 3,
     *   corresponding to the mesh dimension).
     * \param _gdim Not used by this basis type. Must be identical to _bdim.
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \param _map_tol Tolerance limit for least-squares solver in
     *   Map_MeshToBasis and Map_BasisToMesh
     * \note This constructor should not be called directly. Instead, create
     *   a class object via Raster2::Create<Raster_Pixel2>.
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster_Pixel2 (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0, double _map_tol=1e-10);

    RCompRowMatrix *CreateBvv () const;

    RCompRowMatrix *CreateBvw (const IVector &wdim) const;
    RCompRowMatrix *CreateBvw_tri (const IVector &wdim) const;
    RCompRowMatrix *CreateBvw_tet4 (const IVector &wdim) const;

private:
    bool grid_is_basis;       ///< grid and basis dimensions identical?

    /**
     * \brief Returns pointer to mixed mapping matrix
     */
    RCompRowMatrix *CreateBuv () const;
    RCompRowMatrix *CreateBuv_tri () const;
    RCompRowMatrix *CreateBuv_tet4 () const;

    void AddToElMatrix_tri (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

    void AddToElMatrix_tet (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

    int SutherlandHodgman (int el, int xgrid, int ygrid, Point *clip_poly,
        int npoly) const;

#if THREAD_LEVEL==2
    friend void CreateBuv_tri_pass2_engine (task_data *td);
#endif
};
    
#endif // !__RASTER_PX2_H
