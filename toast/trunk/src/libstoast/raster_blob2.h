// -*-C++-*-

#ifndef __RASTER_BLOB2_H
#define __RASTER_BLOB2_H

#include "raster2.h"

/**
 * \file Implements class Raster_Blob2 (base class for inverse bases that
 *   consist of radially symmetric basis functions with limited support).
 *   The classes derived from this use a least-squares approach for mapping
 *   between FEM and inverse basis (see papers/basismap_new).
 *   The required mixed-basis integrals for mapping matrix Buv are computed
 *   by element subdivision, where the boundary of the radial basis function
 *   is approximated by a polygon.
 */

class STOASTLIB Raster_Blob2: public Raster2 {
public:
    template<class T>
    static Raster2 *Create (const IVector &_bdim, const IVector &_gdim,
        Mesh *mesh, double _sup, double shapeprm, RDenseMatrix *bb=0,
        double _map_tol=1e-10)
    {
        T *raster = new T(_bdim, _gdim, mesh, _sup, shapeprm, bb,
	    _map_tol);
	raster->Init();
	return raster;
    }

    /**
     * \brief Blob basis class constructor.
     * \param _bdim basis grid dimensions
     * \param _gdim intermediate grid dimensions (not used, must be identical
     *    to  _bdim)
     * \param mesh mesh pointer
     * \param _sup blob support radius in units of grid spacing
     * \param shapeprm additional blob shape parameter if required
     * \param bb bounding box
     * \param _map_tol tolerance for linear solver in mapping
     */
    Raster_Blob2 (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, double shapeprm, RDenseMatrix *bb=0,
        double _map_tol=1e-10);

    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

    void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

protected:
    /**
     * \brief Basis function value as a function of radial distance r
     * \param r radial distance in units of grid spacing
     * \return blob basis function value
     */
    virtual double RadValue (double r) const = 0;

    int SutherlandHodgman (int el, double px, double py, Point *clip_poly,
        int npoly) const;

    double sup;    ///< support radius
    double sprm;   ///< shape parameter for basis types that use it
    RVector igrid; ///< inverse grid spacing
private:
    RCompRowMatrix *CreateBasisMassmat () const;
    RCompRowMatrix *CreateBasisMassmat_tri () const;

    /**
     * \brief Returns pointer to mixed mapping matrix
     */
    RCompRowMatrix *CreateMixedMassmat () const;
    RCompRowMatrix *CreateMixedMassmat_tri () const;

    double MC_integral_2D(double basis_dst) const;
};

#endif // !__RASTER_BLOB2_H
