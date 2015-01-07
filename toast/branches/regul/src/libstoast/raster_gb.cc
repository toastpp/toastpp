#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

Raster_GaussBlob::Raster_GaussBlob (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sigma, double _sup, RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb), sigma(_sigma)
{
    scale = 1.0;
    a2 = sup*sup;
    isigma2 = 1.0/(sigma*sigma);
    scale = ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_GaussBlob::Value_nomask (const Point &p, int i, bool is_solidx)
    const
{
    int ib;
    if (is_solidx) {
	RANGE_CHECK (i >= 0 && i < slen);
	ib = GetBasisIdx(i);
    } else {
	RANGE_CHECK (i >= 0 && i < blen);
	ib = i;
    }
    int d;
    IVector b;
    GetBasisIndices(ib,b);
    RVector r(dim);
    for (d = 0; d < dim; d++) {
	r[d] = (b[d]-npad) + (bbmin[d]-p[d])*igrid[d];
    }
    double r2 = l2normsq(r);
    if (r2 >= a2) return 0.0; // p outside support of basis i
    return scale * exp (-r2*isigma2);
}
