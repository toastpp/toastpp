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
     * \param _sup blob support radius in units of grid spacing
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster_Blob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, RDenseMatrix *bb = 0);

    virtual ~Raster_Blob();

    /**
     * \brief Returns the value of a basis function at a given point.
     * \param p point coordinates
     * \param i linear basis index
     */
    virtual double Value (const Point &p, int i, bool is_solidx=true) const;

    virtual double Value_nomask (const Point &p, int i, bool is_solidx=true)
	const = 0;

    /**
     * \brief Fill return all basis function values at a given mesh node
     * \param [in] node mesh node index (>= 0)
     * \param [out] nv vector of basis function values
     */
    virtual void NodeValues (int node, RVector &nv) const;
    RVector NodeValues (int node) const
    { RVector tmp(slen); NodeValues (node, tmp); return tmp; }

    void Map_MeshToBasis (const RVector &mvec, RVector &bvec) const;
    void Map_MeshToBasis (const CVector &mvec, CVector &bvec) const;

    void Map_MeshToGrid (const RVector &mvec, RVector &gvec) const;
    void Map_MeshToGrid (const CVector &mvec, CVector &gvec) const;

    void Map_MeshToSol (const RVector &mvec, RVector &svec) const;
    void Map_MeshToSol (const CVector &mvec, CVector &svec) const;

    void Map_BasisToMesh (const RVector &bvec, RVector &mvec) const;
    void Map_BasisToMesh (const CVector &bvec, CVector &mvec) const;

    void Map_GridToMesh (const RVector &gvec, RVector &mvec) const;
    void Map_GridToMesh (const CVector &gvec, CVector &mvec) const;

    void Map_SolToMesh (const RVector &svec, RVector &mvec) const;
    void Map_SolToMesh (const CVector &svec, CVector &mvec) const;

    void Map_GridToBasis (const RVector &gvec, RVector &bvec) const;
    void Map_GridToBasis (const CVector &gvec, CVector &bvec) const;

    void Map_BasisToGrid (const RVector &bvec, RVector &gvec) const;
    void Map_BasisToGrid (const CVector &bvec, CVector &gvec) const;

    void Map_SolToBasis (const RVector &svec, RVector &bvec) const;
    void Map_SolToBasis (const CVector &svec, CVector &bvec) const;

    void Map_BasisToSol (const RVector &bvec, RVector &svec) const;
    void Map_BasisToSol (const CVector &bvec, CVector &svec) const;
    
    void Map_SolToGrid (const RVector &svec, RVector &gvec) const;
    void Map_SolToGrid (const CVector &svec, CVector &gvec) const;

    void Map_GridToSol (const RVector &gvec, RVector &svec) const;
    void Map_GridToSol (const CVector &gvec, CVector &svec) const;

protected:
    double ComputeBasisScale ();

    void ComputeNodeValues ();

    RCompRowMatrix *nodevals;
    // Sparse matrix of size nlen x slen
    // containing nodal basis values for each basis function
    // Initialised by ComputeNodeValues()

    double sup;               // blob support radius in units of grid spacing
    RVector igrid;            // inverse grid spacing in x, y (, z)
    int npad;                 // number of padding layers around bounding box
    Point bbmin_pad, bbmax_pad; // bounding box including padding
};

#endif // !__RASTER_BL_H
