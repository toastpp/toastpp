#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

Raster_HanningBlob::Raster_HanningBlob (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, double _sup, RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb)
{
    a2 = sup*sup;
    fac = M_PI/sup;
    scale = 0.5;
    scale *= ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_HanningBlob::Value_nomask (const Point &p, int i) const
{
    RANGE_CHECK (i >= 0 && i < slen);
    int d;
    IVector b;
    GetBasisIndices(GetBasisIdx(i),b);
    RVector r(dim);
    for (d = 0; d < dim; d++) {
	r[d] = (b[d]-npad) + (bbmin[d]-p[d])*igrid[d];
    }
    double rad = l2norm(r);
    if (rad >= sup) return 0.0; // p outside support of basis i
    return scale * (1.0 + cos(fac*rad));
}
