#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

Raster_BesselBlob::Raster_BesselBlob (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, double _alpha, double _sup,
    RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb), alpha(_alpha)
{
    a2 = sup*sup;
    scale = 1.0/bessi(2,alpha);
    scale *= ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_BesselBlob::Value_nomask (const Point &p, int i) const
{
    RANGE_CHECK (i >= 0 && i < slen);
    int d;
    IVector b;
    GetBasisIndices(GetBasisIdx(i), b);
    RVector r(dim);
    for (d = 0; d < dim; d++) {
	r[d] = (b[d]-npad) + (bbmin[d]-p[d])*igrid[d];
    }
    double r2 = l2normsq(r);
    if (r2 >= a2) return 0.0; // p outside support of basis i
    double s = sqrt (1.0 - r2/a2);
    return scale * s*s * bessi(2,alpha*s);
}
