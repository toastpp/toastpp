// -*-C++-*-

#ifndef __RASTER_CS_H
#define __RASTER_CS_H

/**
 * \file Implements class Raster_CubicSpline (cubic spline blob basis).
 */

#include "raster.h"

// =========================================================================
/**
 * \brief Cubic spline blob basis representation.
 */
class STOASTLIB Raster_CubicSpline: public Raster_Blob {
public:
    /**
     * \brief Cubic spline basis constructor.
     * \param _bdim basis dimensions (the length of this vector is 2 or 3,
     *   corresponding to the mesh dimension).
     * \param _gdim High-res grid dimensions (the length of this vector is 2
     *   or 3, corresponding to the mesh dimension).
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster_CubicSpline (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0);
};

#endif // !__RASTER_CS_H
