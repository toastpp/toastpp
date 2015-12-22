#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

// =========================================================================
// class Raster_CPixel
// =========================================================================

Raster_CPixel::Raster_CPixel (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb, double _map_tol)
  : Raster2 (_bdim, _gdim, mesh, bb, _map_tol)
{
    elsize = 1.0;
    for (int i = 0; i < dim; i++)
	elsize *= bbsize[i] / bdim[i];
}

// =========================================================================

double Raster_CPixel::Value_nomask (const Point &p, int i, bool is_solidx)
    const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim[0]*bdim[1]));
    bi -= iz*bdim[0]*bdim[1];
    int iy = bi/bdim[0];
    int ix = bi - iy*bdim[0];

    int px = (int)((p[0]-bbmin[0])/bbsize[0]*bdim[0]);
    int py = (int)((p[1]-bbmin[1])/bbsize[1]*bdim[1]);
    if (px != ix || py != iy) return 0.0;
    if (dim == 3) {
	int pz = (int)((p[2]-bbmin[2])/bbsize[2]*bdim[2]);
	if (pz != iz) return 0.0;
    }
    return 1.0;
}

// =========================================================================

RCompRowMatrix *Raster_CPixel::CreateBvv () const
{
    // construct a sparse diagonal matrix
    idxtype *rowptr = new idxtype[blen+1];
    idxtype *colidx = new idxtype[blen];
    double *val = new double[blen];
    rowptr[0] = 0;
    for (int i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+1;
	colidx[i] = i;
	val[i] = elsize;
    }
    RCompRowMatrix *bvv = new RCompRowMatrix(blen,blen,rowptr,colidx,val);
    delete []rowptr;
    delete []colidx;
    delete []val;

    return bvv;
}

// =========================================================================

RCompRowMatrix *Raster_CPixel::CreateBuv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBuv_tri();
    case ELID_TET4:
	return CreateBuv_tet4();
    default:
	xERROR("Raster_CPixel: Unsupported element type");
	return 0;
    }
}

// ==========================================================================
// Adds contribution from single element to system matrix

void Raster_CPixel::AddToElMatrix (int el,
    RGenericSparseMatrix &M, const RVector *pxcoeff, int mode) const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    //case ELID_TRI6:
    //case ELID_TRI10:
	AddToElMatrix_tri (el, M, pxcoeff, mode);
	break;
    case ELID_TET4:
	//AddToElMatrix_tet (el, M, pxcoeff, mode);
	break;
    default:
	xERROR("Raster_Pixel2: Unsupported element type");
    }
}
