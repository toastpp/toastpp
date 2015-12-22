#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_sb2.h"

Raster_Blob2_SB::Raster_Blob2_SB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, shapeprm, diagscale, bb, _map_tol,
		_npad)
{
    int i, j;
    for (i = 0; i <= 4; i++)
	xi[i] = (i-2)*0.5*sup;
    for (i = 0; i <= 4; i++) {
	iomega[i] = 1.0;
	for (j = 0; j <= 4; j++)
	    if (i != j)
		iomega[i] /= xi[i]-xi[j];
    }
}

double Raster_Blob2_SB::RadValue (double r) const
{
    if (r >= sup) return 0.0;

    double x = -r;
    double dx = x-xi[0];
    double d3 = dx*dx*dx;
    double sum = d3*iomega[0];
    if (x > xi[1]) {
	dx = x-xi[1];
	d3 = dx*dx*dx;
	sum += d3*iomega[1];
    }
    sum *= 4.0;
    return sum*bscale;
}

double Raster_Blob2_SB::RadGradient (double r) const
{
    if (r >= sup) return 0.0;

    double x = -r;
    double dx = x-xi[0];
    double d2 = dx*dx;
    double sum = -3.0*iomega[0]*d2;
    if (x > xi[1]) {
	dx = x-xi[1];
	d2 = dx*dx;
	sum += -3.0*iomega[1]*d2;
    }
    sum *= 4.0;
    return sum*bscale;
}
