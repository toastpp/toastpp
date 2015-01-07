// -*-C++-*-

#ifndef __RASTER_GB_H
#define __RASTER_GB_H

/**
 * \file Implements class Raster_GaussBlob (radially symmetric basis
 *   functions with Gaussian profile)
 */

class STOASTLIB Raster_GaussBlob: public Raster_Blob {
public:
    Raster_GaussBlob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sigma, double _sup, RDenseMatrix *bb = 0);

protected:
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

private:
    double sigma, isigma2, a2, scale;
};

#endif // !__RASTER_GB_H
