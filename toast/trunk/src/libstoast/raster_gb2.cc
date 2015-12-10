#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_gb2.h"

Raster_Blob2_GB::Raster_Blob2_GB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, shapeprm, diagscale, bb, _map_tol,
		_npad)
{
    a2 = sup*sup;
    isigma2 = 1.0/(sprm*sprm);
}

double Raster_Blob2_GB::RadValue (double r) const
{
    if (r >= sup) return 0.0;
    double r2 = r*r;
    return exp(-r2*isigma2);
}
