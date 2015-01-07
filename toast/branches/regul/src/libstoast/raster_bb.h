// -*-C++-*-

#ifndef __RASTER_BB_H
#define __RASTER_BB_H

/**
 * \file Implements class Raster_BesselBlob (radially symmetric basis
 *   functions with Kaiser-Bessel profile)
 */

class STOASTLIB Raster_BesselBlob: public Raster_Blob {
public:
    Raster_BesselBlob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _alpha, double _sup, RDenseMatrix *bb = 0);

protected:
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

private:
    double alpha, a2, scale;
};

#endif // !__RASTER_GB_H
