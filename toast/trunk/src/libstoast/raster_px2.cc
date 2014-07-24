#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tri_qr.h"

using namespace std;

// ==========================================================================
// class Raster_Pixel2 (v.2)

Raster_Pixel2::Raster_Pixel2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb)
    : Raster (_bdim, _gdim, mesh, bb)
{
    int i, j;

    xASSERT(bdim==gdim,
	    "This basis type doesn't support intemediate grid basis");

    // Compute the matrices for the least squares mapping between
    // mesh and pixel basis
    Buu = meshptr->MassMatrix();
    Bvv = CreatePixelMassmat();
    Buv = CreateMixedMassmat();

    // formulate basis->solution mapping in sparse matrix
    idxtype *rowptr = new idxtype[slen+1];
    idxtype *colidx = new idxtype[slen];
    double *val = new double[slen];
    for (i = 0; i <= slen; i++) rowptr[i] = i; // each row has one entry
    for (i = 0; i < slen; i++) val[i] = 1.0;
    for (i = j = 0; i < blen; i++)
        if (bsupport[i] > 0.0) colidx[j++] = i;
    D = new RCompRowMatrix (slen, blen, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;

    if (toastVerbosity > 0) {
        cout << "--> Type............" << (mesh->Dimension() == 2 ? "Bi":"Tri")
	     << "-linear" << endl;
    }
}

// ==========================================================================

Raster_Pixel2::~Raster_Pixel2 ()
{
    delete D;
    delete Buu;
    delete Bvv;
    delete Buv;
}

// ==========================================================================

double Raster_Pixel2::Value_nomask (const Point &p, int i, bool is_solidx) const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim[0]*bdim[1]));
    bi -= iz*bdim[0]*bdim[1];
    int iy = bi/bdim[0];
    int ix = bi - iy*bdim[0];

    double intvx = (bbmax[0]-bbmin[0])/(double)(bdim[0]-1);
    double intvy = (bbmax[1]-bbmin[1])/(double)(bdim[1]-1);

    double dx = fabs(bbmin[0] + ix*intvx - p[0]);
    double vx = max(0.0, 1.0-dx/intvx);
    double dy = fabs(bbmin[1] + iy*intvy - p[1]);
    double vy = max(0.0, 1.0-dy/intvy);
    double v = vx*vy;

    if (dim > 2) {
	double intvz = (bbmax[2]-bbmin[2])/(double)(bdim[2]-1);
	double dz = fabs(bbmin[2] + iz*intvz - p[2]);
	double vz = max(0.0, 1.0-dz/intvz);
	v *= vz;
    }
    return v;
}

// ==========================================================================

void Raster_Pixel2::Map_GridToBasis (const RVector &gvec, RVector &bvec) const
{
    bvec = gvec; // NOP
}

// ==========================================================================

void Raster_Pixel2::Map_GridToBasis (const CVector &gvec, CVector &bvec) const
{
    bvec = gvec; // NOP
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToGrid (const RVector &bvec, RVector &gvec) const
{
    gvec = bvec; // NOP
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToGrid (const CVector &bvec, CVector &gvec) const
{
    gvec = bvec; // NOP
}

// ==========================================================================

void Raster_Pixel2::Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
{
    double tol = 1e-10;
    PCG (*Bvv, ATx(*Buv,mvec), bvec, tol);
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    double tol = 1e-10;
    PCG (*Buu, Ax(*Buv,bvec), mvec, tol);
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    D->Ax (bvec, svec);
}

// ==========================================================================

void Raster_Pixel2::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    D->ATx (svec, bvec);
}

// ==========================================================================

void Raster_Pixel2::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    RVector bvec(blen);
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster_Pixel2::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    RVector bvec(blen);
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}

// ==========================================================================
// Sutherland-Hodgman polygon clipping algorithm
// el: element index
// xgrid, ygrid: x and y basis indices for clipping cell
// clip_poly:    pointer to vertex array defining the clipped polygon on output
// npoly:        length of clip_poly list
// return value: number of vertices returned (or -1 if npoly is not sufficiently large)

typedef struct { double x, y; } vec_t, *vec;

static inline double dot(vec a, vec b)
{
    return a->x * b->x + a->y * b->y;
}

static inline double cross(vec a, vec b)
{
    return a->x * b->y - a->y * b->x;
}
 
static inline vec vsub(vec a, vec b, vec res)
{
    res->x = a->x - b->x;
    res->y = a->y - b->y;
    return res;
}
 
/* tells if vec c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
static int left_of(vec a, vec b, vec c)
{
    vec_t tmp1, tmp2;
    double x;
    vsub(b, a, &tmp1);
    vsub(c, b, &tmp2);
    x = cross(&tmp1, &tmp2);
    return x < 0 ? -1 : x > 0;
}
 
static int line_sect(vec x0, vec x1, vec y0, vec y1, vec res)
{
    vec_t dx, dy, d;
    vsub(x1, x0, &dx);
    vsub(y1, y0, &dy);
    vsub(x0, y0, &d);
    /* x0 + a dx = y0 + b dy ->
       x0 X dx = y0 X dx + b dy X dx ->
       b = (x0 - y0) X dx / (dy X dx) */
    double dyx = cross(&dy, &dx);
    if (!dyx) return 0;
    dyx = cross(&d, &dx) / dyx;
    if (dyx <= 0 || dyx >= 1) return 0;
    
    res->x = y0->x + dyx * dy.x;
    res->y = y0->y + dyx * dy.y;
    return 1;
}
 
/* === polygon stuff === */
typedef struct { int len, alloc; vec v; } poly_t, *poly;
 
static poly poly_new()
{
    return (poly)calloc(1, sizeof(poly_t));
}
 
static void poly_free(poly p)
{
    free(p->v);
    free(p);
}

static void poly_append(poly p, vec v)
{
    if (p->len >= p->alloc) {
	p->alloc *= 2;
	if (!p->alloc) p->alloc = 4;
	p->v = (vec)realloc(p->v, sizeof(vec_t) * p->alloc);
    }
    p->v[p->len++] = *v;
}

static void poly_reset(poly p)
{
    p->len = 0;
}
 
/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
*/
static int poly_winding(poly p)
{
	return left_of(p->v, p->v + 1, p->v + 2);
}
 
static void poly_edge_clip(poly sub, vec x0, vec x1, int left, poly res)
{
	int i, side0, side1;
	vec_t tmp;
	vec v0 = sub->v + sub->len - 1, v1;
	res->len = 0;
 
	side0 = left_of(x0, x1, v0);
	if (side0 != -left) poly_append(res, v0);
 
	for (i = 0; i < sub->len; i++) {
		v1 = sub->v + i;
		side1 = left_of(x0, x1, v1);
		if (side0 + side1 == 0 && side0)
			/* last point and current straddle the edge */
			if (line_sect(x0, x1, v0, v1, &tmp))
				poly_append(res, &tmp);
		if (i == sub->len - 1) break;
		if (side1 != -left) poly_append(res, v1);
		v0 = v1;
		side0 = side1;
	}
}
 
static poly poly_clip(poly sub, poly clip)
{
	int i;
	static poly p1 = poly_new(), p2 = poly_new();
	poly tmp;
	poly_reset(p1); poly_reset(p2);

	int dir = poly_winding(clip);
	poly_edge_clip(sub, clip->v + clip->len - 1, clip->v, dir, p2);
	for (i = 0; i < clip->len - 1; i++) {
		tmp = p2; p2 = p1; p1 = tmp;
		if(p1->len == 0) {
			p2->len = 0;
			break;
		}
		poly_edge_clip(p1, clip->v + i, clip->v + i + 1, dir, p2);
	}
 
	//poly_free(p1);
	return p2;
}
 
int Raster_Pixel2::SutherlandHodgman (int el, int xgrid, int ygrid,
    Point *clip_poly, int npoly) const
{
    int i,j;

    static int nlocal = 10;
    static Point *list = new Point[nlocal];
    if (nlocal < npoly) {
	delete []list;
	list = new Point[nlocal=npoly];
    }

    double xmin = bbmin[0] +
	(bbmax[0]-bbmin[0])*(double)xgrid/(double)(bdim[0]-1);
    double xmax = bbmin[0] +
	(bbmax[0]-bbmin[0])*(double)(xgrid+1)/(double)(bdim[0]-1);
    double ymin = bbmin[1] +
	(bbmax[1]-bbmin[1])*(double)ygrid/(double)(bdim[1]-1);
    double ymax = bbmin[1] +
	(bbmax[1]-bbmin[1])*(double)(ygrid+1)/(double)(bdim[1]-1);

    Element *pel = meshptr->elist[el];
    int ni, nv = pel->nNode();
    if (nv > npoly) return -1;

    vec_t tri[3];
    for (i = 0; i < 3; i++) {
	tri[i].x = meshptr->nlist[pel->Node[i]][0];
	tri[i].y = meshptr->nlist[pel->Node[i]][1];
    }
    poly_t subject = {3, 0, tri};

    vec_t clip[4];
    clip[0].x = clip[3].x = xmin;
    clip[1].x = clip[2].x = xmax;
    clip[0].y = clip[1].y = ymin;
    clip[2].y = clip[3].y = ymax;
    poly_t clipper = {4, 0, clip};

    poly res = poly_clip (&subject, &clipper);
    if (res->len > npoly) return -1;

    for (i = 0; i < res->len; i++) {
	clip_poly[i].New(2);
	clip_poly[i][0] = res->v[i].x;
	clip_poly[i][1] = res->v[i].y;
    }
    return res->len;
}

// ==========================================================================
// Creates a mixed mass matrix 

// Calculates overlaps between a triangle and a pixel, subdivides the triangle
// and performs the integral over the sub-triangles using a numerical
// quadrature rule

RCompRowMatrix *Raster_Pixel2::CreateMixedMassmat () const
{
    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, jj, idx_i, idx_j;
    int imin, imax, jmin, jmax;
    double b;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);

    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    // pass 1: determine matrix fill structure
    int *nimin = new int[n];
    int *nimax = new int[n];
    int *njmin = new int[n];
    int *njmax = new int[n];
    for (i = 0; i < n; i++) {
	nimin[i] = bdim[0];
	njmin[i] = bdim[1];
	nimax[i] = njmax[i] = -1;
    }
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	}
	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
	for (i = 0; i < pel->nNode(); i++) {
	    int nidx = pel->Node[i];
	    if (imin < nimin[nidx]) nimin[nidx] = imin;
	    if (imax > nimax[nidx]) nimax[nidx] = imax;
	    if (jmin < njmin[nidx]) njmin[nidx] = jmin;
	    if (jmax > njmax[nidx]) njmax[nidx] = jmax;
	}
    }

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (r = 0; r < n; r++) {
	int nentry = (nimax[r]-nimin[r]+2)*(njmax[r]-njmin[r]+2);
	rowptr[r+1] = rowptr[r]+nentry;
    }
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (r = k = 0; r < n; r++) {
	for (j = njmin[r]; j <= njmax[r]+1; j++) {
	    for (i = nimin[r]; i <= nimax[r]+1; i++) {
		colidx[k++] = i + j*bdim[0];
	    }
	}
    }

    RCompRowMatrix *Buv = new RCompRowMatrix (n, blen, rowptr, colidx);
    delete []rowptr;
    delete []colidx;
    delete []nimin;
    delete []nimax;
    delete []njmin;
    delete []njmax;

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	xASSERT(pel->Type() == ELID_TRI3, "Currently only implemented for 3-noded triangles");

	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	}

	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));

	// perform subdivision for every intersecting pixel
	const int npoly = 10;
	Point poly[npoly];
	Triangle3 t3;
	NodeList n3(3);
	for (i = 0; i < 3; i++) {
	    t3.Node[i] = i;
	    n3[i].New(2);
	}
	int nv;
	RVector fun(3);
	double v;
	for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++) {
		nv = SutherlandHodgman (el, i, j, poly, npoly);

		// split into nv-2 triangles
		for (k = 0; k < nv-2; k++) {
		    n3[0][0] = poly[0][0];
		    n3[0][1] = poly[0][1];
		    n3[1][0] = poly[k+1][0];
		    n3[1][1] = poly[k+1][1];
		    n3[2][0] = poly[k+2][0];
		    n3[2][1] = poly[k+2][1];
		    t3.Initialise(n3);
		    double djac = t3.Size()*2.0;
		    
		    // map quadrature points into global frame
		    for (m = 0; m < np; m++) {
			Point glob = t3.Global(n3, absc[m]);
			Point loc = pel->Local(meshptr->nlist, glob);
			fun = pel->LocalShapeF (loc);
			v = wght[m] * djac;
			for (jj = 0; jj < 4; jj++) {
			    idx_j = i + jj%2 + (j+jj/2)*bdim[0];
			    b = Value_nomask (glob, idx_j, false);
			    for (ii = 0; ii < 3; ii++) {
				idx_i = pel->Node[ii];
				(*Buv)(idx_i, idx_j) += v*b*fun[ii];
			    }
			}
		    }
		}
	    }
    }

    Buv->Shrink();
    return Buv;
}

// ==========================================================================
// Creates the pixel-pixel mass matrix (Bvv)

RCompRowMatrix *Raster_Pixel2::CreatePixelMassmat () const
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
    // WARNING: Currently only supports 2D
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
    
    return M;
}
