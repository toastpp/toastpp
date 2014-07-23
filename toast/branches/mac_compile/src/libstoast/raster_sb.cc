#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

Raster_SplineBlob::Raster_SplineBlob (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, double _sup, RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb)
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
    scale = 1.0;
    scale *= ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_SplineBlob::Value_nomask (const Point &p, int i, bool is_solidx)
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
    double rad = l2norm(r);
    if (rad >= sup) return 0.0; // p outside support of basis i
    double x = -rad;
    double dx = x-xi[0];
    double d3 = dx*dx*dx;
    double sum = d3*iomega[0];

    if (x > xi[1]) {
	dx = x-xi[1];
	d3 = dx*dx*dx;
	sum += d3*iomega[1];
    }
    sum *= 4.0;
    return scale * sum;
}
