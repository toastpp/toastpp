#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_rb2.h"

Raster_Blob2_RB::Raster_Blob2_RB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, shapeprm, diagscale, bb, _map_tol,
		_npad)
{
    isup = 1.0/sup;
}

double Raster_Blob2_RB::RadValue (double r) const
{
    return (r >= sup ? 0.0 : 1.0 - r*isup);
}

