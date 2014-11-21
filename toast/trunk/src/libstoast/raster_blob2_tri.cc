// Implementation of those routines in Raster_Blob2 specific to triangular
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"
#include "tri_qr.h"
#include "sh.h"

// MC integral of two blob basis functions, with distance basis_dst between
// their centres (basis_dst in units of grid spacing)
double Raster_Blob2::MC_integral_2D(double basis_dst) const
{
    if (basis_dst >= 2.0*sup) return 0.0;

    double px, py, r1, r2, v1, v2, sum = 0.0;
    const int nsample = 100000;
    for (int i = 0; i < nsample; i++) {
	px = (drand48()-0.5)*2.0*sup;
	py = (drand48()-0.5)*2.0*sup;
	r1 = sqrt(px*px+py*py);
	v1 = RadValue(r1);
	if (basis_dst) {
	    px = basis_dst-px;
	    r2 = sqrt(px*px+py*py);
	    v2 = RadValue(r2);
	} else v2 = v1;
	sum += v1*v2;
    }
    return sum/(double)nsample*4.0*sup*sup;
}

RCompRowMatrix *Raster_Blob2::CreateBasisMassmat_tri () const
{
    int i, j, k, r, c;
    double dx, dy, dst;
    RVector d(dim);
    for (i = 0; i < dim; i++) {
	d[i] = bbsize[i]/(bdim[i]-1.0);
    }

    // neighbour stencils
    int range = (int)ceil(2.0*sup/d[0]); // assumes isotropic grid spacing
    int nstencil_x = range*2+1;
    int nstencil_y = nstencil_x;
    int nstencil = nstencil_x * nstencil_y;

    std::cerr << "range=" << range << std::endl;

    double *w = new double[nstencil];    // store distance-dependent weights
    double *wdst = new double[nstencil]; // stored distances
    int nw = 0;                          // number of stored weights

    int *stencil = new int[nstencil];        // relative stencil indices
    double *stencilw = new double[nstencil]; // stencil weights
    for (j = 0; j < nstencil_y; j++) {
	dy = (double)(j-range)*d[1];
	for (i = 0; i < nstencil_x; i++) {
	    dx = (double)(i-range)*d[0];
	    dst = sqrt(dx*dx+dy*dy);
	    for (k = 0; k < nw; k++) { // check if weight is available
		if (fabs(wdst[k]-dst) < 1e-10)
		    break;
	    }
	    if (k == nw) { // compute weight
		wdst[nw] = dst;
		w[nw] = MC_integral_2D (dst);
		nw++;
	    }
	    stencilw[i+j*nstencil_x] = w[k];
	    stencil[i+j*nstencil_x] = (i-range) + (j-range)*bdim[0];
	}
    }
    delete []w;
    delete []wdst;

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
    double *val = new double[nentry];

    rowptr[0] = 0;
    for (r = j = 0; r < blen; r++) {
	rowptr[r+1] = rowptr[r] + nrowentry[r];
	for (i = 0; i < nstencil; i++) {
	    c = r + stencil[i];
	    if (c >= 0 && c < blen) {
		colidx[j] = c;
		val[j] = stencilw[i];
		j++;
	    }
	}
    }
    xASSERT(j == nentry, "Inconsistency!");
    delete []stencil;
    delete []stencilw;

    RCompRowMatrix *Bvv = new RCompRowMatrix (blen, blen, rowptr, colidx, val);
    delete []rowptr;
    delete []colidx;
    delete []val;
    delete []nrowentry;

    Bvv->Shrink();
    return Bvv;
}

RCompRowMatrix *Raster_Blob2::CreateMixedMassmat_tri () const
{
    int i, j, k, ii, jj, nv, m, el, nd, idx_i, idx_j;
    int nel = meshptr->elen(), n = meshptr->nlen();
    bool intersect;
    double rx, ry, djac, v;
    RVector fun;

    // grid range and spacing (assumes regular arrangement of blob
    // basis functions)
    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double radlimit2 = sup*sup;
    
    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    int *npx = new int[n];
    for (i = 0; i < n; i++) npx[i] = 0;

    // pass 1: determine matrix fill structure
    double px, py;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = 0; j < bdim[1]; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = 0; i < bdim[0]; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++)
			npx[pel->Node[k]]++;
		}
	    }
	}
    }

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (i = 0; i < n; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (i = 0; i < n; i++)
	npx[i] = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = 0; j < bdim[1]; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = 0; i < bdim[0]; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++) {
			nd = pel->Node[k];
			colidx[rowptr[nd]+npx[nd]++] = i + j*bdim[0];
		    }
		}
	    }
	}
    }

    RCompRowMatrix *buv = new RCompRowMatrix (n, blen, rowptr, colidx);
    delete []rowptr;
    delete []colidx;
    delete []npx;

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	const int npoly = 10;
	Point poly[npoly];
	Triangle3 t3;
	NodeList n3(3);
	for (i = 0; i < 3; i++) {
	    t3.Node[i] = i;
	    n3[i].New(2);
	}

	for (j = 0; j < bdim[1]; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = 0; i < bdim[0]; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    nv = SutherlandHodgman (el, px, py, poly, npoly);
		    idx_j = i + j*bdim[0];

		    // split into nv-2 triangles
		    for (k = 0; k < nv-2; k++) {
			n3[0][0] = poly[0][0];
			n3[0][1] = poly[0][1];
			n3[1][0] = poly[k+1][0];
			n3[1][1] = poly[k+1][1];
			n3[2][0] = poly[k+2][0];
			n3[2][1] = poly[k+2][1];
			t3.Initialise(n3);

			for (m = 0; m < np; m++) {
			    djac = t3.DetJ(absc[m], &n3);
			    Point glob = t3.Global(n3, absc[m]);
			    Point loc = pel->Local(meshptr->nlist, glob);
			    fun = pel->LocalShapeF (loc);
			    v = wght[m] * djac * Value_nomask(glob,idx_j,false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				(*buv)(idx_i, idx_j) += v*fun[ii];
			    }
			}
		    }
		}    
	    }
	}
    }

    buv->Shrink();
    return buv;
}

int Raster_Blob2::SutherlandHodgman (int el, double px, double py,
    Point *clip_poly, int npoly) const
{
    int i;

    Point *list = new Point[npoly];

    Element *pel = meshptr->elist[el];
    int ni, nv = pel->nNode();
    if (nv > npoly) return -1;

    vec_t tri[3];
    for (i = 0; i < 3; i++) {
	tri[i].x = meshptr->nlist[pel->Node[i]][0];
	tri[i].y = meshptr->nlist[pel->Node[i]][1];
    }
    poly_t subject = {3, 0, tri};

    const int nseg = 16; // segments for approximating support circle
    vec_t clip[nseg];
    for (i = 0; i < nseg; i++) {
	double phi = (double)i/(double)nseg * Pi2;
	clip[i].x = cos(phi)*sup + px;
	clip[i].y = sin(phi)*sup + py;
    }
    poly_t clipper = {nseg, 0, clip};

    poly res = poly_clip (&subject, &clipper);
    int len = res->len;
    if (len > npoly) {
	len = -1;
    } else {
	for (i = 0; i < len; i++) {
	    clip_poly[i].New(2);
	    clip_poly[i][0] = res->v[i].x;
	    clip_poly[i][1] = res->v[i].y;
	}
    }
    poly_free(res);
    delete []list;
    return len;
}
