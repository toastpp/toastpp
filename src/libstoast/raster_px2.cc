#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "timing.h"

using namespace std;

// =========================================================================
// class Raster_Pixel2 (v.2)
// =========================================================================

Raster_Pixel2::Raster_Pixel2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb, double _map_tol)
: Raster2 (_bdim, _gdim, mesh, bb, _map_tol)
{
    if (toastVerbosity > 0) {
        cout << "--> Type............"
	     << (mesh->Dimension() == 2 ? "Bi":"Tri") << "-linear" << endl;
    }
}

// =========================================================================

Raster_Pixel2::~Raster_Pixel2 ()
{
}

// =========================================================================

double Raster_Pixel2::Value_nomask (const Point &p, int i, bool is_solidx) const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim[0]*bdim[1]));
    bi -= iz*bdim[0]*bdim[1];
    int iy = bi/bdim[0];
    int ix = bi - iy*bdim[0];

    double intvx = bbsize[0]/(double)(bdim[0]-1);
    double intvy = bbsize[1]/(double)(bdim[1]-1);

    double dx = fabs(bbmin[0] + ix*intvx - p[0]);
    double vx = max(0.0, 1.0-dx/intvx);
    double dy = fabs(bbmin[1] + iy*intvy - p[1]);
    double vy = max(0.0, 1.0-dy/intvy);
    double v = vx*vy;

    if (dim > 2) {
	double intvz = bbsize[2]/(double)(bdim[2]-1);
	double dz = fabs(bbmin[2] + iz*intvz - p[2]);
	double vz = max(0.0, 1.0-dz/intvz);
	v *= vz;
    }
    return v;
}

// =========================================================================
// Creates the pixel-pixel mass matrix (Bvv)

RCompRowMatrix *Raster_Pixel2::CreateBvv () const
{
    // neighbour stencils
    static const int nstencil2 = 9;
    static const int nstencil3 = 27;

    // integral weights on unit square for 2D problems
    static const double w2[3] = {
	1.0/9.0,  // diagonal
	1.0/18.0, // edge neighbours
	1.0/36.0  // diagonal neighbours
    }; 
    // integral weights on unit cube for 3D problems
    static const double w3[4] = {
	1.0/27.0,  // diagonal
	1.0/54.0,  // edge neighbours
	1.0/108.0, // side diagonal neighbours
	1.0/216.0  // volume diagonal neighbour
    };

    int i, j, k, k0, k1, r, c, nstencil;
    int c0i, c0j, c0k, c1i, c1j, c1k;
    double scale = 1.0;
    RVector d(dim);
    for (i = 0; i < dim; i++) {
	d[i] = bbsize[i]/(bdim[i]-1.0);
	scale = scale*d[i];
    }

    nstencil = (dim==2 ? nstencil2 : nstencil3);
    int *stencil = new int[nstencil];
    if (dim == 2) {
	for (j = 0; j < 3; j++)
	    for (i = 0; i < 3; i++) {
		stencil[i+j*3] = (i-1) + (j-1)*bdim[0];
	    }
    } else {
	for (k = 0; k < 3; k++)
	    for (j = 0; j < 3; j++)
		for (i = 0; i < 3; i++) {
		    stencil[i+(j+k*3)*3] =
			(i-1) + ((j-1) + (k-1)*bdim[1])*bdim[0];
		}
    }

    // evaluate sparse matrix structure
    int *nrowentry = new int[blen];
    int nentry = 0;
    for (r = 0; r < blen; r++) {
	nrowentry[r] = 0;
	for (i = 0; i < nstencil; i++) {
	    c = r + stencil[i];
	    if (c >= 0 && c < blen) nrowentry[r]++;
	}
	nentry += nrowentry[r];
    }

    int *rowptr = new int[blen+1];
    int *colidx = new int[nentry];
	
    rowptr[0] = 0;
    for (r = i = 0; r < blen; r++) {
	rowptr[r+1] = rowptr[r] + nrowentry[r];
	for (j = 0; j < nstencil; j++) {
	    c = r + stencil[j];
	    if (c >= 0 && c < blen)
		colidx[i++] = c;
	}
    }
 
    RCompRowMatrix *M = new RCompRowMatrix (blen, blen, rowptr, colidx);
    delete []rowptr;
    delete []colidx;
    delete []stencil;
    delete []nrowentry;

    // fill matrix
    if (dim == 2) { // 2D case
	for (j = 0; j < bdim[1]-1; j++) {
	    for (i = 0; i < bdim[0]-1; i++) {
		int base = i + j*bdim[0];
		for (k1 = 0; k1 < 4; k1++) {
		    int idx1 = base + k1%2 + (k1/2)*bdim[0];
		    for (k0 = 0; k0 < 4; k0++) {
			int idx0 = base + k0%2 + (k0/2)*bdim[0];
			double v = (k0 == k1 ? w2[0] :
				    k0+k1==3 ? w2[2] : w2[1])*scale;
			(*M)(idx0,idx1) += v;
		    }
		}
	    }
	}
    } else { // 3D case
	int stride_j = bdim[0];
	int stride_k = bdim[0]*bdim[1];
	for (k = 0; k < bdim[2]-1; k++) {
	    for (j = 0; j < bdim[1]-1; j++) {
		for (i = 0; i < bdim[0]-1; i++) {
		    int base = i + j*stride_j + k*stride_k;
		    for (k1 = 0; k1 < 8; k1++) {
			c1i = k1%2;
			c1j = (k1/2)%2;
			c1k = k1/4;
			int idx1 = base + c1i + c1j*stride_j + c1k*stride_k;
			for (k0 = 0; k0 < 8; k0++) {
			    c0i = k0%2;
			    c0j = (k0/2)%2;
			    c0k = k0/4;
			    int idx0 = base + c0i + c0j*stride_j +
				c0k*stride_k;
			    int rank = abs(c0i-c1i) + abs(c0j-c1j) +
				abs(c0k-c1k);
			    (*M)(idx0,idx1) += w3[rank]*scale;
			}
		    }
		}
	    }
	}
    }
    
    return M;
}

// ==========================================================================
// Creates the basis-pixel mass matrix (Bvw) for mapping from a raster image
// with piecewise constant pixel values to the basis

RCompRowMatrix *Raster_Pixel2::CreateBvw (const IVector &wdim)
    const
{
    xASSERT(wdim.Dim() == dim, "Invalid length of pixel dimension vector");

    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
	return CreateBvw_tri (wdim);
    case ELID_TET4:
	return CreateBvw_tet4 (wdim);
    default:
	xERROR("Raster_Pixel2: unsupported element type");
	return 0;
    }
}

// ==========================================================================
// Creates the mixed-basis mass matrix (Buv)

RCompRowMatrix *Raster_Pixel2::CreateBuv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBuv_tri();
    case ELID_TET4:
	return CreateBuv_tet4();
    default:
	xERROR("Raster_Pixel2: Unsupported element type");
	return 0;
    }
}

// ==========================================================================
// Adds contribution from single element to system matrix

void Raster_Pixel2::AddToElMatrix (int el,
    RGenericSparseMatrix &M, const RVector *pxcoeff, int mode) const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    //case ELID_TRI6:
    //case ELID_TRI10:
	AddToElMatrix_tri (el, M, pxcoeff, mode);
	break;
    case ELID_TET4:
	AddToElMatrix_tet (el, M, pxcoeff, mode);
	break;
    default:
	xERROR("Raster_Pixel2: Unsupported element type");
    }
}
