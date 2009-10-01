// -*-C++-*-

#ifndef __RASTER_BL_H
#define __RASTER_BL_H

/**
 * \file Implements class Raster_Blob (base class for all radially
 *   symmetric "blob"-type basis expansions).
 */

#include "raster.h"

// =========================================================================
/**
 * \brief Generic blob basis representation.
 */
class STOASTLIB Raster_Blob: public Raster {
public:
    /**
     * \brief Blob basis constructor.
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
    Raster_Blob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0);

    /**
     * \brief Returns the value of a basis function at a given point.
     * \param p point coordinates
     * \param i linear basis index
     */
    virtual double Value (const Point &p, int i) const = 0;
};

#endif // !__RASTER_BL_H
