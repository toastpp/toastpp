// -*-C++-*-

#ifndef __RASTER_CS2_H
#define __RASTER_CS2_H

#include "raster_blob2.h"

class STOASTLIB Raster_Blob2_CS: public Raster_Blob2 {
public:
    Raster_Blob2_CS (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
		     RDenseMatrix *bb=0, double _map_tol=1e-10);

protected:
    double RadValue (double r);
};

#endif // !__RASTER_CS2_H
