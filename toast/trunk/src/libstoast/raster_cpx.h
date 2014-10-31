// -*-C++-*-

#ifndef __RASTER_CPX_H
#define __RASTER_CPX_H

/**
 * \file Defines class Raster_CPixel (piecewise constant pixel/voxel
 *   image basis)
 */

#include "raster2.h"

/**
 * \brief Basis for 2D and 3D image representation with a regular piecewise
 *   constant pixel/voxel grid.
 */

class STOASTLIB Raster_CPixel: public Raster2 {
    friend class Raster2;
    // to make constructor accessible to Raster2 factory

public:
    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;
  
protected:
    /**
     * \brief Pixel basis constructor
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
    Raster_CPixel (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb=0, double _map_tol=1e-10);

    RCompRowMatrix *CreateBasisMassmat () const;
    RCompRowMatrix *CreateMixedMassmat () const;
    RCompRowMatrix *CreateMixedMassmat_tri () const;
    RCompRowMatrix *CreateMixedMassmat_tet4 () const;

private:
    int SutherlandHodgman (int el, int xgrid, int ygrid, Point *clip_poly,
        int npoly) const;

    double elsize;     ///< element volume
};

#endif // !__RASTER_CPX_H
