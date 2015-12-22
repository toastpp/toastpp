#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"

// =========================================================================
// class Raster_Blob2
// =========================================================================

Raster_Blob2::Raster_Blob2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double shapeprm, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
    : Raster2 (_bdim, _gdim, mesh, bb, _map_tol), sup(_sup), sprm(shapeprm),
	       dgscale(diagscale), npad(_npad)
{
    bdim_pad = bdim + 2*npad;
    blen_pad = prod(bdim_pad);
    
    grid.New(dim);
    igrid.New(dim);
    for (int i = 0; i < dim; i++) {
        grid[i] = bbsize[i]/(bdim[i]-1.0); // grid spacing
	igrid[i] = 1.0/grid[i];            // inverse grid spacing
    }
}

void Raster_Blob2::Init ()
{
    int i, j, k;
    double brescale;
    
    // determine scaling factor so that sum of basis coefficients approx 1
    bscale = 1.0;
    brescale = 0.0;
    IVector nnbr(dim);
    for (i = 0; i < dim; i++)
        nnbr[i] = (int)ceil(sup*igrid[i]);
    int nnbrz = (dim < 3 ? 0 : nnbr[2]);
    for (k = -nnbrz; k <= nnbrz; k++) {
        double dz = (k ? k/igrid[2] : 0.0);
        for (j = -nnbr[1]; j <= nnbr[1]; j++) {
	    double dy = j/igrid[1];
	    for (i = -nnbr[0]; i <= nnbr[0]; i++) {
	        double dx = i/igrid[0];
		double dst = sqrt(dx*dx + dy*dy + dz*dz);
		brescale += RadValue(dst);
	    }
	}
    }
    bscale = 1.0/brescale;
    // Note: we need to compute bscale before calling Raster2::Init, so that
    // the correct bscale is used for calculating the mapping matrices

    Raster2::Init();
}

double Raster_Blob2::Value_nomask (const Point &p, int i, bool is_solidx) const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim_pad[0]*bdim_pad[1]));
    bi -= iz*bdim_pad[0]*bdim_pad[1];
    int iy = bi/bdim_pad[0];
    int ix = bi - iy*bdim_pad[0];

    ix -= npad;
    iy -= npad;

    double dx = bbmin[0] + ix*grid[0] - p[0];
    double dy = bbmin[1] + iy*grid[1] - p[1];

    double dst2 = dx*dx + dy*dy;;

    if (dim > 2) {
	iz -= npad;
	double dz = bbmin[2] + iz*grid[2] - p[2];
	dst2 += dz*dz;
    }
    
    return RadValue (sqrt(dst2));
}

// ==========================================================================

RVector Raster_Blob2::Gradient_nomask (const Point &p, int i, bool is_solidx)
    const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim_pad[0]*bdim_pad[1]));
    bi -= iz*bdim_pad[0]*bdim_pad[1];
    int iy = bi/bdim_pad[0];
    int ix = bi - iy*bdim_pad[0];

    ix -= npad;
    iy -= npad;

    double dx = p[0] - (bbmin[0] + ix*grid[0]);
    double dy = p[1] - (bbmin[1] + iy*grid[1]);
    double dz = 0.0;
    double dst2 = dx*dx + dy*dy;
    double phi = atan2(dy,dx);
    double theta = 0.0;

    if (dim > 2) {
	iz -= npad;
	dz = p[2] - (bbmin[2] + iz*grid[2]);
	dst2 += dz*dz;
    }
    double r = sqrt(dst2);
    if (dim > 2 && r)
	theta = acos(dz/r);
    
    double rg = bscale * RadGradient (r);
    RVector g(dim);
    g[0] = cos(phi)*rg;
    g[1] = sin(phi)*rg;
    if (dim > 2) {
	double sint = sin(theta);
	g[0] *= sint;
	g[1] *= sint;
	g[2] = cos(theta)*rg;
    }
    return g;
}

// ==========================================================================
// Creates the mixed-basis mass matrix (Buv)

RCompRowMatrix *Raster_Blob2::CreateBuv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBuv_tri();
    case ELID_TET4:
	return CreateBuv_tet4();
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

// ==========================================================================
// Creates the mixed-basis stiffness matrix (Duv)

RCompRowMatrix *Raster_Blob2::CreateBvv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBvv_tri();
    case ELID_TET4:
	return CreateBvv_tet4();
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

RCompRowMatrix *Raster_Blob2::CreateDvv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateDvv_tri();
    case ELID_TET4:
	//return CreateDvv_tet4();
	return 0;
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

RCompRowMatrix *Raster_Blob2::CreateDuv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateDuv_tri();
    case ELID_TET4:
	//return CreateBuv_tet4();
	return 0;
    default:
	xERROR("Raster_Blob2: Unsupported element type");
	return 0;
    }
}

RCompRowMatrix *Raster_Blob2::CreateBvw (const IVector &pxdim) const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBvw_tri (pxdim);
    case ELID_TET4:
	return CreateBvw_tet4 (pxdim);
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

