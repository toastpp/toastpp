#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tri_qr.h"
#include "tet_qr.h"

using namespace std;

// ==========================================================================
// class Raster_Pixel2 (v.2)

Raster_Pixel2::Raster_Pixel2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb, double _map_tol)
: Raster (_bdim, _gdim, mesh, bb), map_tol(_map_tol)
{
    int i, j;

    xASSERT(bdim==gdim,
	    "This basis type doesn't support intemediate grid basis");

    // Compute the matrices for the least squares mapping between
    // mesh and pixel basis
    Buu = meshptr->MassMatrix();
    Bvv = CreatePixelMassmat();
    Buv = CreateMixedMassmat();

    ofstream ofs1("Buu.dat");
    Buu->ExportRCV(ofs1);
    ofstream ofs2("Buv.dat");
    Buv->ExportRCV(ofs2);
    ofstream ofs3("Bvv.dat");
    Bvv->ExportRCV(ofs3);

    Buu_precon = new RPrecon_IC; Buu_precon->Reset (Buu);
    Bvv_precon = new RPrecon_IC; Bvv_precon->Reset (Bvv);

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
    delete Buu_precon;
    delete Bvv_precon;
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
    double tol = map_tol;
    int nit = PCG (*Bvv, ATx(*Buv,mvec), bvec, tol, Bvv_precon);
    //std::cerr << "Map_MeshToBasis: niter=" << nit << std:: endl;
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    double tol = map_tol;
    int nit = PCG (*Buu, Ax(*Buv,bvec), mvec, tol, Buu_precon);
    //std::cerr << "Map_BasisToMesh: niter=" << nit << std:: endl;
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
    // for now
    if (meshptr->elist[0]->Type() == ELID_TET4)
	return CreateMixedMassmat_tet4();

    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, jj, idx_i, idx_j;
    int imin, imax, jmin, jmax;
    double b, djac;

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
	//xASSERT(pel->Type() == ELID_TRI3, "Currently only implemented for 3-noded triangles");

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
	RVector fun;
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
		    //double djac = t3.Size()*2.0;
		    
		    // map quadrature points into global frame
		    for (m = 0; m < np; m++) {
			djac = t3.DetJ(absc[m], &n3);
			Point glob = t3.Global(n3, absc[m]);
			Point loc = pel->Local(meshptr->nlist, glob);
			fun = pel->LocalShapeF (loc);
			v = wght[m] * djac;
			for (jj = 0; jj < 4; jj++) {
			    idx_j = i + jj%2 + (j+jj/2)*bdim[0];
			    b = Value_nomask (glob, idx_j, false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
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




// Mesh class with pre-allocated element and node lists, so we don't
// need to dynamically grow the lists too often
// nlen_used and elen_used are the length of the lists actually in use

class BufMesh: public Mesh
{
public:
    BufMesh(): Mesh()
    { nlen_used = elen_used = 0; }

    void SubSetup()
    {
	for (int el=0; el < elen_used; el++)
	    elist[el]->Initialise(nlist);
	for (int el=0; el < elen_used; el++)
	    elist[el]->PostInitialisation(nlist);
    }

    int nlen_used, elen_used;
};

RCompRowMatrix *Raster_Pixel2::CreateMixedMassmat_tet4 () const
{
    void Tetsplit (BufMesh *mesh, int el, int cut_orient, double cut_pos);

    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, jj, idx_i, idx_j;
    int imin, imax, jmin, jmax, kmin, kmax;
    RVector fun;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    double xmin, xmax, ymin, ymax, zmin, zmax, djac, b, v;

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_4_14 (&wght, &absc);

    // pass 1: determine matrix fill structure
    int *nimin = new int[n];
    int *nimax = new int[n];
    int *njmin = new int[n];
    int *njmax = new int[n];
    int *nkmin = new int[n];
    int *nkmax = new int[n];
    for (i = 0; i < n; i++) {
	nimin[i] = bdim[0];
	njmin[i] = bdim[1];
	nkmin[i] = bdim[2];
	nimax[i] = njmax[i] = nkmax[i] = -1;
    }
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	double ezmin = meshptr->nlist[pel->Node[0]][2];
	double ezmax = meshptr->nlist[pel->Node[0]][2];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	    ezmin = min (ezmin, meshptr->nlist[pel->Node[j]][2]);
	    ezmax = max (ezmax, meshptr->nlist[pel->Node[j]][2]);
	}
	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor((bdim[2]-1) * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-2, (int)floor((bdim[2]-1) * (ezmax-bbmin[2])/zrange));
	for (i = 0; i < pel->nNode(); i++) {
	    int nidx = pel->Node[i];
	    if (imin < nimin[nidx]) nimin[nidx] = imin;
	    if (imax > nimax[nidx]) nimax[nidx] = imax;
	    if (jmin < njmin[nidx]) njmin[nidx] = jmin;
	    if (jmax > njmax[nidx]) njmax[nidx] = jmax;
	    if (kmin < nkmin[nidx]) nkmin[nidx] = kmin;
	    if (kmax > nkmax[nidx]) nkmax[nidx] = kmax;
	}
    }

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (r = 0; r < n; r++) {
	int nentry = (nimax[r]-nimin[r]+2)*(njmax[r]-njmin[r]+2)*(nkmax[r]-nkmin[r]+2);
	rowptr[r+1] = rowptr[r]+nentry;
    }
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (r = m = 0; r < n; r++) {
	for (k = nkmin[r]; k <= nkmax[r]+1; k++) {
	    for (j = njmin[r]; j <= njmax[r]+1; j++) {
		for (i = nimin[r]; i <= nimax[r]+1; i++) {
		    colidx[m++] = i + j*bdim[0] + k*bdim[0]*bdim[1];
		}
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
    delete []nkmin;
    delete []nkmax;

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
	std::cout << "Buv: processing el " << el << " of " << nel << std::endl;
	Element *pel = meshptr->elist[el];
	xASSERT(pel->Type() == ELID_TET4, "Currently only implemented for 4-noded tetrahedra");

	double orig_size = pel->Size();

	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	double ezmin = meshptr->nlist[pel->Node[0]][2];
	double ezmax = meshptr->nlist[pel->Node[0]][2];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	    ezmin = min (ezmin, meshptr->nlist[pel->Node[j]][2]);
	    ezmax = max (ezmax, meshptr->nlist[pel->Node[j]][2]);
	}

	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor((bdim[2]-1) * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-2, (int)floor((bdim[2]-1) * (ezmax-bbmin[2])/zrange));

	// Create a mesh containing this single element
	for (i = 0; i < 4; i++) {
	    submesh.elist[0]->Node[i] = i;
	    submesh.nlist[i] = meshptr->nlist[pel->Node[i]];
	}
	submesh.elist[0]->SetRegion(0);
	submesh.elen_used = 1;
	submesh.nlen_used = 4;

#ifdef UNDEF
	// DEBUG
	submesh.SubSetup();
	double orig_size = submesh.CalcFullSize();
	//ofstream ofs1("dbg1.msh");
	//ofs1 << submesh << std::endl;
	std::cout << "Initial element size: " << orig_size << std::endl;
#endif

	// perform subdivisions along x-cutting planes
	for (i = imin; i < imax; i++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[0] + (i+1)*dx;
	    for (m = 0; m < elen; m++) {
		int subreg = submesh.elist[m]->Region() & 0x3FF;
		if (subreg >= i-imin)
		    Tetsplit (&submesh, m, 0, cut_pos);
	    }
	}
	    
#ifdef UNDEF
	// DEBUG
	//ofstream ofs2("dbg2.msh");
	//ofs2 << submesh << std::endl;
	submesh.SubSetup();
	double sz = 0.0;
	std::cout << "Subdiv element sizes: " << std::endl;
	for (i = 0; i < submesh.elen_used; i++) {
	    double s = submesh.elist[i]->Size();
	    std::cout << s << std::endl;
	    sz += s;
	    if (s <= 0)
		std::cerr << "**** x-split: El size(" << el << ',' << i
			  << ")=" << s << std::endl;
	}
	std::cout << "Total: " << sz << std::endl;
#endif

	// perform subdivisions along y-cutting planes
	for (j = jmin; j < jmax; j++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[1] + (j+1)*dy;
	    for (m = 0; m < elen; m++) {
		int subreg = (submesh.elist[m]->Region() >> 10) & 0x3FF;
		if (subreg >= j-jmin)
		    Tetsplit (&submesh, m, 1, cut_pos);
	    }
	}
	    
#ifdef UNDEF
	// DEBUG
	//ofstream ofs3("dbg3.msh");
	//ofs3 << submesh << std::endl;
	submesh.SubSetup();
	sz = 0.0;
	std::cout << "Subdiv element sizes: " << std::endl;
	for (i = 0; i < submesh.elen_used; i++) {
	    double s = submesh.elist[i]->Size();
	    std::cout << s << std::endl;
	    sz += s;
	    if (s <= 0)
		std::cerr << "**** y-split: El size(" << el << ',' << i
			  << ")=" << s << std::endl;
	}
	std::cout << "Total: " << sz << std::endl;
#endif

	// perform subdivisions along z-cutting planes
	for (k = kmin; k < kmax; k++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[2] + (k+1)*dz;
	    for (m = 0; m < elen; m++) {
		int subreg = (submesh.elist[m]->Region() >> 20) & 0x3FF;
		if (subreg >= k-kmin)
		    Tetsplit (&submesh, m, 2, cut_pos);
	    }
	}

#ifdef UNDEF
	// DEBUG
	//ofstream ofs4("dbg4.msh");
	//ofs4 << submesh << std::endl;
	submesh.SubSetup();
	sz = 0.0;
	std::cout << "Subdiv element sizes: " << std::endl;
	for (i = 0; i < submesh.elen_used; i++) {
	    double s = submesh.elist[i]->Size();
	    std::cout << s << std::endl;
	    sz += s;
	    if (s <= 0)
		std::cerr << "**** z-split: El size(" << el << ',' << i
			  << ")=" << s << std::endl;
	}
	std::cout << "Total: " << sz << std::endl;
#endif

	// now perform quadrature on all sub-elements
	submesh.SubSetup();
	int nskipped = 0;
	for (i = 0; i < submesh.elen_used; i++) {
	    Element *psel = submesh.elist[i];
	    if (psel->Size() <= orig_size*1e-10) {
		nskipped++;
		continue;
	    }
	    // TEMPORARY: skip tiny subelements

	    // find the voxel cell containing the sub-tet
	    int reg = psel->Region();
	    int ci = imin + reg & 0x3FF;
	    int cj = jmin + (reg >> 10) & 0x3FF;
	    int ck = kmin + (reg >> 20) & 0x3FF;
	    int cellidx = ci + (cj + ck*bdim[1])*bdim[0];
	    // map quadrature points into global frame
	    for (m = 0; m < np; m++) {
		djac = psel->DetJ(absc[m], &submesh.nlist);
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		v = wght[m] * djac;
		// loop over the vertices of the voxel cell
		for (jj = 0; jj < 8; jj++) {
		    idx_j = cellidx + jj%2 + ((jj/2)%2)*bdim[0] +
			(jj/4)*bdim[0]*bdim[1];
		    b = Value_nomask (glob, idx_j, false);
		    for (ii = 0; ii < fun.Dim(); ii++) {
			idx_i = pel->Node[ii];
			(*Buv)(idx_i, idx_j) += v*b*fun[ii];
		    }
		}
	    }
	    
	}
	std::cout << "skipped: " << nskipped << " of " << submesh.elen_used
		  << " sub-elements" << std::endl;
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



// =========================================================================
// =========================================================================
// Tetrahedra splitting code. This is preliminary and should eventually
// be shifted elsewhere (e.g. to the individual element classes)
// =========================================================================
// =========================================================================

// Calculate the intersection point pi of a line defined by two points p0
// and p1, and a cutting plane parallel to two coordinate axes
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// cut_pos is the position of the intersection of the cutting plane with
// the corresponding axis
// It is assumed that there is an intersection, i.e. the line is not
// parallel to the plane
inline void Intersect_line_plane (int cut_orient, double cut_pos,
			   const Point &p0, const Point &p1, Point &pi)
{
    double s = (cut_pos-p0[cut_orient])/(p1[cut_orient]-p0[cut_orient]);
    for (int i = 0; i < 3; i++)
	pi[i] = (i == cut_orient ? cut_pos : p0[i] + (p1[i]-p0[i])*s);
}


// Check if a point is left of a cutting plane
inline bool Point_is_left (const Point &p, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    return (p[cut_orient] < cut_pos-eps);
}


// Check if a point is right of a cutting plane
inline bool Point_is_right (const Point &p, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    return (p[cut_orient] > cut_pos+eps);
}


// Check if a tetrahedron has a portion left of a cutting plane
bool Tet_extends_left (Mesh *mesh, int el, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Check if a tetrahedron has a portion right of a cutting plane
bool Tet_extends_right (Mesh *mesh, int el, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Calculate the intersection points of the edges of tetrahedron el
// with a cutting plane parallel to two axes.
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// cut_pos is the position of the intersection of the cutting plane with
// the corresponding axis
// Return value is the number of intersection points found (0, 3, or 4)
// The points are written to the array pointed to by isect (must be allocated to
// at least 4 points of dimension 3)
int CalcIntersections_tet(Mesh *mesh, int el, int cut_orient, double cut_pos, Point *isect)
{
    const double eps = 1e-8;
    int i, j, k, nleft, nright, nsingle;
    Element *pel = mesh->elist[el];

    for (i = nleft = nright = 0; i < 4; i++) {
	Point &pt = mesh->nlist[pel->Node[i]];
	if      (Point_is_left  (pt, cut_orient, cut_pos)) nleft++;
	else if (Point_is_right (pt, cut_orient, cut_pos)) nright++;
   }
    if (!nleft || !nright) return 0; // no intersection

    if (nleft == 1 || nright == 1) { // triangular intersection

	// find the single node
	if (nleft == 1) {
	    for (i = 0; i < 4; i++)
		if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		    nsingle = i; break;
		}
	} else {
	    for (i = 0; i < 4; i++)
		if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		    nsingle = i; break;
		}
	}
	Point &psingle = mesh->nlist[pel->Node[nsingle]];
	double vsingle = psingle[cut_orient];

	// find the intersection points on the 3 edges connected to psingle
	for (i = j = 0; i < 4; i++) {
	    if (i == nsingle) continue;
	    Point &p = mesh->nlist[pel->Node[i]];
	    Intersect_line_plane (cut_orient, cut_pos, psingle, p, isect[j++]);
	}
	return 3;

    } else { // quadrilateral intersection

	int nleft[2], nright[2];
	for (i = j = k = 0; i < 4; i++) {
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
		nleft[j++] = i;
	    else
		nright[k++] = i;
	}

	// re-arrange a bit for simpler computation
	if (nleft[0] != 0 && nleft[1] != 0) {
	    nleft[0] = 0;
	    nleft[1] = (nright[0] == 0 ? nright[1]:nright[0]);
	    nright[0] = (nleft[1] == 1 ? 2:1);
	    nright[1] = (nleft[1] == 3 ? 2:3);
	} else {
	    if (nleft[0] != 0) {
		nleft[1] = nleft[0];
		nleft[0] = 0;
	    }
	    if (nright[0] > nright[1]) {
		int tmp = nright[0]; nright[0] = nright[1]; nright[1] = tmp;
	    }
	}
	

	static const int tet_edge[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	// the 6 edges of the tetrahedron, defined by the node indices

	for (i = j = 0; i < 6; i++) {
	    if ((tet_edge[i][0] == nleft[0] && tet_edge[i][1] == nleft[1]) ||
		(tet_edge[i][0] == nright[0] && tet_edge[i][1] == nright[1]))
		continue; // we skip the 2 edges with nodes on the same side of the cutting plane
	    Intersect_line_plane (cut_orient, cut_pos,
				  mesh->nlist[pel->Node[tet_edge[i][0]]],
				  mesh->nlist[pel->Node[tet_edge[i][1]]],
				  isect[j++]);
	}
	return 4;

    }
}


// Split a tetrahedron that has one node on one side of the cutting plane,
// the other 3 on the other side, into 8 subtets.
// 'inode' is the index of the single node (0..3)
// isect is a list of 3 global points defining the intersections of the
// tet edges with the cutting plane
// On exit, the returned mesh has modified element el, added the 7 additional
// tetrahedra, and appended the additional nodes to its node list
// If mod_el_idx is set, it should point to an integer array of length >= 8,
// and will receive the indices of the 8 modified/added tetrahedra
// Return value is the number of modified/added elements (8)

int Tetsplit_1_3 (BufMesh *mesh, int el, int inode, const Point *isect,
		   int *mod_el_idx=0)
{
    int i, j;

    Element *pel[8];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	for (i = 0; i < 7; i++)
	    mod_el_idx[i+1] = mesh->elen_used+i;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+7 > ebuf)
	for (i = 0; i < mesh->elen_used+7-ebuf; i++)
	    mesh->elist.Append(new Tetrahedron4);
    
    for (i = 1; i < 8; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+6 > nbuf) {
	nlist.Append(mesh->nlen_used+6-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 6;

    int nidx[10];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    for (i = 0; i < 6; i++)
	nidx[i+4] = nlen+i;

    switch (inode) {
    case 0:
	nlist[nlen+0] = isect[0];
	nlist[nlen+1] = isect[1];
	nlist[nlen+2] = isect[2];
	nlist[nlen+3] = (nlist[nidx[1]]+nlist[nidx[2]]) * 0.5;
	nlist[nlen+4] = (nlist[nidx[1]]+nlist[nidx[3]]) * 0.5;
	nlist[nlen+5] = (nlist[nidx[2]]+nlist[nidx[3]]) * 0.5;
	break;
    case 1:
	nlist[nlen+0] = isect[0];
	nlist[nlen+1] = (nlist[nidx[0]]+nlist[nidx[2]]) * 0.5;
	nlist[nlen+2] = (nlist[nidx[0]]+nlist[nidx[3]]) * 0.5;
	nlist[nlen+3] = isect[1];
	nlist[nlen+4] = isect[2];
	nlist[nlen+5] = (nlist[nidx[2]]+nlist[nidx[3]]) * 0.5;
	break;
    case 2:
	nlist[nlen+0] = (nlist[nidx[0]]+nlist[nidx[1]]) * 0.5;
	nlist[nlen+1] = isect[0];
	nlist[nlen+2] = (nlist[nidx[0]]+nlist[nidx[3]]) * 0.5;
	nlist[nlen+3] = isect[1];
	nlist[nlen+4] = (nlist[nidx[1]]+nlist[nidx[3]]) * 0.5;
	nlist[nlen+5] = isect[2];
	break;
    case 3:
	nlist[nlen+0] = (nlist[nidx[0]]+nlist[nidx[1]]) * 0.5;
	nlist[nlen+1] = (nlist[nidx[0]]+nlist[nidx[2]]) * 0.5;
	nlist[nlen+2] = isect[0];
	nlist[nlen+3] = (nlist[nidx[1]]+nlist[nidx[2]]) * 0.5;
	nlist[nlen+4] = isect[1];
	nlist[nlen+5] = isect[2];
	break;
    }

    int local_nidx[8][4] = {
	{0,4,5,6},{4,1,7,8},{5,7,2,9},{6,8,9,3},
	{6,4,5,8},{4,7,5,8},{5,7,9,8},{5,9,6,8}
    };
    for (i = 0; i < 8; i++)
	for (j = 0; j < 4; j++)
	    pel[i]->Node[j] = nidx[local_nidx[i][j]];
	    
    return 8;
}


// Split a tetrahedron that has two nodes on either side of the cutting plane,
// into 6 subtets.
// 'inode' is the list of two node indices on the same side of the cutting
// plane.
// isect is a list of 4 global points defining the intersections of the
// tet edges with the cutting plane
// On exit, the returned mesh has modified element el, added the 5 additional
// tetrahedra, and appended the additional nodes to its node list
// If mod_el_idx is set, it should point to an integer array of length >= 6,
// and will receive the indices of the 6 modified/added tetrahedra
// Return value is the number of modified/added elements (6)

int Tetsplit_2_2 (BufMesh *mesh, int el, int *inode, const Point *isect,
		   int *mod_el_idx=0)
{
    int i, j, k;
    int onode[2];
    for (i = j = 0; i < 4; i++) {
	for (k = 0; k < 2; k++)
	    if (inode[k] == i) break;
	if (k == 2)
	    onode[j++] = i;
    }
	
    Element *pel[6];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");
    dASSERT(inode[0] == 0, "Inconsistent node order");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	for (i = 0; i < 5; i++)
	    mod_el_idx[i+1] = mesh->elen_used+i;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+5 > ebuf)
	for (i = 0; i < mesh->elen_used+7-ebuf; i++)
	    mesh->elist.Append (new Tetrahedron4);

    for (i = 1; i < 6; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+4 > nbuf) {
	nlist.Append(mesh->nlen_used+4-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 4;

    int nidx[8];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    for (i = 0; i < 4; i++)
	nidx[i+4] = nlen+i;

    for (i = 0; i < 4; i++)
	nlist[nlen+i] = isect[i];

    int local_nidx_a[6][4] = { // split 01-23
	{5,7,4,3},{4,7,6,2},{3,7,4,2},
	{5,7,0,4},{4,7,1,6},{7,1,0,4}};
    int local_nidx_b[6][4] = { // split 02-13
	{4,6,5,1},{5,6,7,3},{1,6,5,3},
	{4,6,0,5},{5,6,2,7},{6,2,0,5}};
    int local_nidx_c[6][4] = { // split 03-12
	{5,7,4,2},{4,7,6,1},{2,7,4,1},
	{5,7,0,4},{4,7,3,6},{7,3,0,4}};
    int (&local_nidx)[6][4] = (inode[1] == 1 ? local_nidx_a :
			       inode[1] == 2 ? local_nidx_b :
			       local_nidx_c);
#ifdef UNDEF
    int local_nidx_a[6][4] = {
	{inode[0],inode[1],6,7},{5,inode[0],4,7},{4,inode[0],6,7},
	{onode[1],5,4,7},{onode[1],4,onode[0],7},{4,6,onode[0],7}
    };
    int local_nidx_b[6][4] = {
	{inode[0],inode[1],7,6},{4,inode[0],5,6},{5,inode[0],7,6},
	{onode[1],4,5,6},{onode[1],6,onode[0],4},{5,7,onode[1],6}
    };

    bool case_a = (inode[0]!=0 || inode[1]!=3) && (inode[0]!=1 || inode[1]!=3);
    int (&local_nidx)[6][4] = (case_a ? local_nidx_a : local_nidx_b);
#endif

    for (i = 0; i < 6; i++)
	for (j = 0; j < 4; j++)
	    pel[i]->Node[j] = nidx[local_nidx[i][j]];

    return 6;
}


// split a tetrahedron in the mesh along a cutting plane.
// cut_orient: orientation of the cutting plane (0: x=const, 1: y=const, 2: z=const)
// cut_pos: position of the cutting plane along the corresponding orientation
// The voxel cell for each of the resulting sub-tetrahedra is encoded in their
// region label. Tetrahedra left of the cutting plane retain their region value.
// Tetrahedra right of the cutting plane get their region index along the
// corresponding slice orientation incremented by one.
// The i,j,k cell grid indices are encoded in the region index as follows:
// reg = i + j<<10 + k<<20.
// i,j,k can thus have values between 0 and 1023.
// i,j,k are indices for the sub-cell grid bounding the given element, rather
// than absolute values

void Tetsplit (BufMesh *mesh, int el, int cut_orient, double cut_pos)
{
    Element *pel = mesh->elist[el];
    int regidx = pel->Region();

    // compute the intersection points
    int i, nisect;
    int elidx[8];
    static Point isect[4];
    if (!isect[0].Dim())
	for (i = 0; i < 4; i++) isect[i].New(3);
    nisect = CalcIntersections_tet (mesh, el, cut_orient, cut_pos, isect);

    if (!nisect) { // no intersections: tet is either fully left or fully right

	if (Tet_extends_right (mesh, el, cut_orient, cut_pos)) // increment slice index
	    pel->SetRegion(regidx + (1 << (cut_orient*10)));

    } else if (nisect == 3) {

	// find the single node to be cut
	int singlend, leftnd, rightnd, nleft, nright;
	for (i = nleft = nright = 0; i < 4; i++)
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		nleft++; leftnd = i;
	    } else if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		nright++; rightnd = i;
	    }
	singlend = (nleft == 1 ? leftnd : rightnd);

	int nelidx = Tetsplit_1_3 (mesh, el, singlend, isect, elidx);

	// modify slice indices for all sub-tetrahedra that ended up right of the cutting plane
	for (i = 0; i < nelidx; i++)
	    if (Tet_extends_right (mesh, elidx[i], cut_orient, cut_pos))
		mesh->elist[elidx[i]]->SetRegion (regidx + (1 << cut_orient*10));
	    else
		mesh->elist[elidx[i]]->SetRegion (regidx);

	
    } else if (nisect == 4) {

	// find the two nodes each left and right of the cutting plane
	int leftnd[2], rightnd[2], nleft, nright;
	for (i = nleft = nright = 0; i < 4; i++) {
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
		leftnd[nleft++] = i;
	    else
		rightnd[nright++] = i;
	}

	int *grp0 = (leftnd[0] == 0 || leftnd[1] == 0 ? leftnd : rightnd);
	int nelidx = Tetsplit_2_2 (mesh, el, grp0, isect, elidx);

	// modify slice indices for all sub-tetrahedra that ended up right of the cutting plane
	for (i = 0; i < nelidx; i++)
	    if (Tet_extends_right (mesh, elidx[i], cut_orient, cut_pos))
		mesh->elist[elidx[i]]->SetRegion (regidx + (1 << cut_orient*10));
	    else
		mesh->elist[elidx[i]]->SetRegion (regidx);
	
    }
}
