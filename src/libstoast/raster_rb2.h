// -*-C++-*-

#ifndef __RASTER_RB2_H
#define __RASTER_RB2_H

#include "raster_blob2.h"

class STOASTLIB Raster_Blob2_RB: public Raster_Blob2 {
public:
    Raster_Blob2_RB (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
	double _sup, double shapeprm, double diagscale, RDenseMatrix *bb=0,
	double _map_tol=1e-10, int _npad=0);

protected:
    double RadValue (double r) const;

private:
    double isup; ///< inverse support radius
};

#endif // !__RASTER_RB2_H
