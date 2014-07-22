// -*-C++-*-

#ifndef __RASTER_RB_H
#define __RASTER_RB_H

/**
 * \file Implements class Raster_RampBlob (radially symmetric basis
 *   functions with linear profile)
 */

class STOASTLIB Raster_RampBlob: public Raster_Blob {
public:
    Raster_RampBlob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, RDenseMatrix *bb = 0);

protected:
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

private:
    double isup, scale;
};

#endif // !__RASTER_RB_H
