// -*-C++-*-

#ifndef __RASTER_BB2_H
#define __RASTER_BB2_H

#include "raster_blob2.h"

// "Bessel-Blob2: radially symmetric blob basis

class STOASTLIB Raster_Blob2_BB: public Raster_Blob2 {
public:
    Raster_Blob2_BB (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, double _alpha, double diagscale,
        RDenseMatrix *bb=0, double _map_tol=1e-10, int _npad=0);

protected:
    double RadValue (double r) const;

private:
    double a2;
};

#endif // !__RASTER_RB2_H
