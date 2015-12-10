#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_bb2.h"

Raster_Blob2_BB::Raster_Blob2_BB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double _alpha, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, _alpha, diagscale, bb, _map_tol,
		_npad)
{
    a2 = sup*sup;
}

double Raster_Blob2_BB::RadValue (double r) const
{
    if (r >= sup) return 0.0;

    double r2 = r*r;
    double s = sqrt(1.0 - r2/a2);
    return s*s * bessi(2,sprm*s);
}
