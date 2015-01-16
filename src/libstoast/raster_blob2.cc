#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"

// =========================================================================
// class Raster_Blob2
// =========================================================================

Raster_Blob2::Raster_Blob2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol)
: Raster2 (_bdim, _gdim, mesh, bb, _map_tol), sup(_sup), sprm(shapeprm),
  dgscale(diagscale)
{
    igrid.New(dim);
    for (int i = 0; i < dim; i++)
	igrid[i] = (bdim[i]-1.0)/bbsize[i]; // inverse grid spacing
}

double Raster_Blob2::Value_nomask (const Point &p, int i, bool is_solidx) const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim[0]*bdim[1]));
    bi -= iz*bdim[0]*bdim[1];
    int iy = bi/bdim[0];
    int ix = bi - iy*bdim[0];

    double intvx = bbsize[0]/(double)(bdim[0]-1);
    double intvy = bbsize[1]/(double)(bdim[1]-1);

    double dx = fabs(bbmin[0] + ix*intvx - p[0]);
    double dy = fabs(bbmin[1] + iy*intvy - p[1]);

    double dst2 = dx*dx + dy*dy;;

    if (dim > 2) {
	double intvz = bbsize[2]/(double)(bdim[2]-1);
	double dz = fabs(bbmin[2] + iz*intvz - p[2]);
	dst2 += dz*dz;
    }
    
    return RadValue (sqrt(dst2));
}

// ==========================================================================
// Creates the mixed-basis mass matrix (Buv)

RCompRowMatrix *Raster_Blob2::CreateMixedMassmat () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateMixedMassmat_tri();
    case ELID_TET4:
	return CreateMixedMassmat_tet4();
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

RCompRowMatrix *Raster_Blob2::CreateBasisMassmat () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBasisMassmat_tri();
    case ELID_TET4:
	return CreateBasisMassmat_tet4();
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

void Raster_Blob2::AddToElMatrix (int el, RGenericSparseMatrix &M,
    const RVector *pxcoeff, int mode) const
{
    // TODO
}

