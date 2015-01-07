// -*-C++-*-

#ifndef __RASTER_CP_H
#define __RASTER_CP_H

/**
 * \file Implements class Raster_CubicPixel (cubic pixel basis)
 */

#include "raster.h"

// =========================================================================
/**
 * \brief Cubic pixel basis representation.
 */
class STOASTLIB Raster_CubicPixel: public Raster {
public:
    /**
     * \brief Cubic pixel basis constructor.
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
    Raster_CubicPixel (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb = 0);

    /**
     * \brief Cubic pixel basis destructor.
     */
    ~Raster_CubicPixel ();

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


protected:
    /**
     * \brief Returns the value of a basis function at a point.
     * \param basisidx basis function index (>= 0)
     * \param p point in global coordinates
     */
    double Value (int basisidx, const Point &p);
    
    /**
     * \brief Returns the value of a basis function at a point
     * \param basisgrd basis function indices for all dimensions
     * \param p point in global coordinates
     */
    double Value (const IVector &basisgrd, const Point &p);

private:
    RVector ibgrid;            ///< inverse basis spacing
    RVector iggrid;            ///< inverse grid spacing
    RGenericSparseMatrix *G;   ///< transformation grid->basis
    RGenericSparseMatrix *GI;  ///< transformation basis->grid
};

#endif // !__RASTER_CP_H
