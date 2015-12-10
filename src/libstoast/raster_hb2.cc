#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_hb2.h"

Raster_Blob2_HB::Raster_Blob2_HB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, shapeprm, diagscale, bb, _map_tol,
		_npad)
{
    fac = Pi/sup;
}

double Raster_Blob2_HB::RadValue (double r) const
{
    if (r >= sup) return 0.0;
    return 1.0 + cos(fac*r);
}
