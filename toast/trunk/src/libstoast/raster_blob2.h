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
        Mesh *mesh, double _sup, double shapeprm, double diagscale,
	RDenseMatrix *bb=0, double _map_tol=1e-10, int _npad=0)
    {
        T *raster = new T(_bdim, _gdim, mesh, _sup, shapeprm, diagscale, bb,
	    _map_tol, _npad);
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
	double _sup, double shapeprm, double diagscale, RDenseMatrix *bb=0,
	double _map_tol=1e-10, int _npad=0);

    virtual void Init();

    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

    RVector Gradient_nomask (const Point &p, int i, bool is_solidx=true) const;

    void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

protected:
    /**
     * \brief Basis function value as a function of radial distance r
     * \param r radial distance in units of grid spacing
     * \return blob basis function value
     */
    virtual double RadValue (double r) const = 0;

    virtual double RadGradient (double r) const
    { xERROR("Not implemented"); return 0.0; }
    
    int SutherlandHodgman (int el, double px, double py, Point *clip_poly,
        int npoly) const;

    int SutherlandHodgman (const Point *p0, int np0,
	const Point *p1, int np1,
	Point *clip_poly, const int npoly) const;

    double sup;     ///< support radius
    double sprm;    ///< shape parameter for basis types that use it
    double dgscale; ///< diagonal scaling factor
    double bscale;  ///< global basis scaling factor
    RVector grid;   ///< grid spacing
    RVector igrid;  ///< inverse grid spacing
    int npad;       ///< padding width [number of grid rows/cols] for blob basis
    IVector bdim_pad; ///< basis dimensions including padding
    int blen_pad;   ///< size of padded basis
    
private:
    RCompRowMatrix *CreateBvv () const;
    RCompRowMatrix *CreateBvv_tri () const;
    RCompRowMatrix *CreateBvv_tet4 () const;

    /**
     * \brief Returns pointer to mixed mapping matrix
     */
    RCompRowMatrix *CreateBuv () const;
    RCompRowMatrix *CreateBuv_tri () const;
    RCompRowMatrix *CreateBuv_tet4 () const;

    RCompRowMatrix *CreateDvv () const;
    RCompRowMatrix *CreateDvv_tri () const;

    RCompRowMatrix *CreateDuv () const;
    RCompRowMatrix *CreateDuv_tri () const;
    
    RCompRowMatrix *CreateBvw (const IVector &wdim) const;
    RCompRowMatrix *CreateBvw_tri (const IVector &wdim) const;
    RCompRowMatrix *CreateBvw_tet4 (const IVector &wdim) const;

    RCompRowMatrix *CreateBvw_linw_tri (const IVector &wdim) const;
    
    double MC_integral_2D(double basis_dst) const;
};

#endif // !__RASTER_BLOB2_H
