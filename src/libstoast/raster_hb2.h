// -*-C++-*-

#ifndef __RASTER_HB2_H
#define __RASTER_HB2_H

#include "raster_blob2.h"

/**
 * \file Implements class Raster_Blob2_HB (radially symmetric basis
 *   functions with Hanning profile)
 *
 * b_a(r) = 0.5(1+cos(Pi*r/a))    (a=radius of support)
 * Refs:
 * KM Hanson, GW Wecksung, Appl Opt 24, 4028-4039 (1985)
 */

class STOASTLIB Raster_Blob2_HB: public Raster_Blob2 {
public:
    Raster_Blob2_HB (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
	double _sup, double shapeprm, double diagscale, RDenseMatrix *bb=0,
	double _map_tol=1e-10, int _npad=0);

protected:
    double RadValue (double r) const;

private:
    double fac;
};

#endif // !__RASTER_HB2_H
