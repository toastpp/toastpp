// -*-C++-*-

#ifndef __RASTER_GB2_H
#define __RASTER_GB2_H

#include "raster_blob2.h"

/**
 * \file Implements class Raster_Blob2_GB (radially symmetric basis
 *   functions with Gaussian profile)
 */

class STOASTLIB Raster_Blob2_GB: public Raster_Blob2 {
public:
    Raster_Blob2_GB (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
	double _sup, double shapeprm, double diagscale, RDenseMatrix *bb=0,
	double _map_tol=1e-10, int _npad=0);

protected:
    double RadValue (double r) const;

private:
    double a2;
    double isigma2;
};

#endif // !__RASTER_GB2_H
