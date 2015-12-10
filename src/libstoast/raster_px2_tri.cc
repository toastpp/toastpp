// Implementation of those routines in Raster_Pixel2 specific to triangular
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tri_qr.h"
#include "sh.h"

// =========================================================================
// Sutherland-Hodgman polygon clipping algorithm
// el:           element index
// xgrid, ygrid: x and y basis indices for clipping cell
// clip_poly:    pointer to vertex array defining the clipped polygon on
//               output
// npoly:        length of clip_poly list
// return value: number of vertices returned (or -1 if npoly is not
//               sufficiently large)

int Raster_Pixel2::SutherlandHodgman (int el, int xgrid, int ygrid,
    Point *clip_poly, int npoly) const
{
    int i,j;

    Point *list = new Point[npoly];

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

// ==========================================================================
// Create the mixed-basis mass matrix Buv by subdividing triangular elements
// at pixel boundaries
// Calculates overlaps between a triangle and a pixel, subdivides the triangle
// and performs the integral over the sub-triangles using a numerical
// quadrature rule

// --------------------------------------------------------------------------
// Engine for parallel computation

#if THREAD_LEVEL==2
struct CREATEMIXEDMASSMAT_TRI_PASS2_THREADDATA {
    const Raster_Pixel2 *raster;
    RCompRowMatrix *Buv;
    Mesh *mesh;
    const IVector *bdim;
    const Point *bbmin;
    const Point *bbmax;
};

void CreateBuv_tri_pass2_engine (task_data *td)
{
    int itask = td->proc;
    int ntask = td->np;
    CREATEMIXEDMASSMAT_TRI_PASS2_THREADDATA *thdata =
	(CREATEMIXEDMASSMAT_TRI_PASS2_THREADDATA*)td->data;
    Mesh *mesh = thdata->mesh;
    int nel = mesh->elen();
    int el0 = (itask*nel)/ntask;
    int el1 = ((itask+1)*nel)/ntask;
    const Raster_Pixel2 *raster = thdata->raster;
    RCompRowMatrix *Buv = thdata->Buv;
    const IVector &bdim = *thdata->bdim;
    const Point &bbmin = *thdata->bbmin;
    const Point &bbmax = *thdata->bbmax;
    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];

    const idxtype *rowptr, *colidx;
    Buv->GetSparseStructure (&rowptr, &colidx);
    RCompRowMatrix Buv_local(Buv->nRows(), Buv->nCols(), rowptr, colidx);

    int i, j, k, m, el, ii, jj, idx_i, idx_j, imin, imax, jmin, jmax;
    double b, djac;

    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    for (el = el0; el < el1; el++) {
	Element *pel = mesh->elist[el];

	// element bounding box
	double exmin = mesh->nlist[pel->Node[0]][0];
	double exmax = mesh->nlist[pel->Node[0]][0];
	double eymin = mesh->nlist[pel->Node[0]][1];
	double eymax = mesh->nlist[pel->Node[0]][1];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, mesh->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, mesh->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, mesh->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, mesh->nlist[pel->Node[j]][1]);
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
		nv = raster->SutherlandHodgman (el, i, j, poly, npoly);

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
			Point loc = pel->Local(mesh->nlist, glob);
			fun = pel->LocalShapeF (loc);
			v = wght[m] * djac;
			for (jj = 0; jj < 4; jj++) {
			    idx_j = i + jj%2 + (j+jj/2)*bdim[0];
			    b = raster->Value_nomask (glob, idx_j, false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				Buv_local(idx_i, idx_j) += v*b*fun[ii];
			    }
			}
		    }
		}
	    }
    }

    // now add thread contribution into global matrix
    double *v = Buv->ValPtr();
    const double *vloc = Buv_local.ValPtr();
    int nz = Buv->nVal();
    Task::UserMutex_lock();
    for (i = 0; i < nz; i++)
	v[i] += vloc[i];
    Task::UserMutex_unlock();
}
#endif // THREAD_LEVEL

RCompRowMatrix *Raster_Pixel2::CreateBuv_tri () const
{
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
#if THREAD_LEVEL==2 // call parallel engine
    static CREATEMIXEDMASSMAT_TRI_PASS2_THREADDATA thdata;
    thdata.raster = this;
    thdata.Buv = Buv;
    thdata.mesh = meshptr;
    thdata.bdim = &bdim;
    thdata.bbmin = &bbmin;
    thdata.bbmax = &bbmax;
    Task::Multiprocess (CreateBuv_tri_pass2_engine, &thdata);
#else // !THREAD_LEVEL: serial computation
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
#endif

    Buv->Shrink();
    return Buv;
}

// ==========================================================================
// Creates the basis-pixel mass matrix (Bvw) for mapping from a raster image
// with piecewise constant pixel values to the basis

RCompRowMatrix *Raster_Pixel2::CreateBvw_tri (const IVector &wdim)
    const
{
    const int sublen = 1000;
    int i, j, k, ii, jj, si, sj, pi, pj, idx_i, idx_j, idx;
    double x0, x1, y0, y1, xcnt, ycnt;
    double px0, px1, py0, py1, v;
    double xmin, xmax, ymin, ymax;
    bool intersect;
    int *npx = new int[blen];
    for (i = 0; i < blen; i++) npx[i] = 0;
    int plen = wdim[0]*wdim[1];
    RVector psize(2), bsize(2);
    Point p(2);
    for (i = 0; i < 2; i++) {
	psize[i] = bbsize[i]/wdim[i];
	bsize[i] = bbsize[i]/(bdim[i]-1.0);
    }
    double dx = bsize[0]/sublen;
    double dy = bsize[1]/sublen;

    // pass 1: evaluate sparse matrix structure
    for (j = 0; j < bdim[1]; j++) {
	ycnt = bsize[1]*j;
	y0 = max(0.0, ycnt-bsize[1]);
	y1 = min(bbsize[1], ycnt+bsize[1]);
	for (i = 0; i < bdim[0]; i++) {
	    idx_i = i + j*bdim[0];
	    xcnt = bsize[0]*i;
	    x0 = max(0.0, xcnt-bsize[0]);
	    x1 = min(bbsize[0], xcnt+bsize[0]);
	    for (pj = 0; pj < wdim[1]; pj++) {
		py0 = psize[1]*pj;
		py1 = py0 + psize[1];
		for (pi = 0; pi < wdim[0]; pi++) {
		    px0 = psize[0]*pi;
		    px1 = px0 + psize[0];
		    if (px0 <= x1 && px1 >= x0 && py0 <= y1 && py1 >= y0)
			npx[idx_i]++;
		}
	    }
	}
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen];
    int *colidx = new int[nz];
    double *scale = new double[nz];
    
    // pass 2: sparsity pattern
    for (j = 0; j < bdim[1]; j++) {
	ycnt = bsize[1]*j;
	y0 = max(0.0, ycnt-bsize[1]);
	y1 = min(bbsize[1], ycnt+bsize[1]);
	for (i = 0; i < bdim[0]; i++) {
	    idx_i = i + j*bdim[0];
	    xcnt = bsize[0]*i;
	    x0 = max(0.0, xcnt-bsize[0]);
	    x1 = min(bbsize[0], xcnt+bsize[0]);	    
	    for (pj = 0; pj < wdim[1]; pj++) {
		py0 = psize[1]*pj;
		py1 = py0 + psize[1];
		for (pi = 0; pi < wdim[0]; pi++) {
		    px0 = psize[0]*pi;
		    px1 = px0 + psize[0];
		    idx_j = pi + pj*wdim[0];
		    if (px0 <= x1 && px1 >= x0 && py0 <= y1 && py1 >= y0) {
			xmin = max(x0,px0);
			xmax = min(x1,px1);
			ymin = max(y0,py0);
			ymax = min(y1,py1);
			idx = rowptr[idx_i]+npx[idx_i]++;
			colidx[idx] = idx_j;
			scale[idx] = (xmax-xmin)*(ymax-ymin);
		    }
		}
	    }
	}
    }

    RCompRowMatrix *bvw = new RCompRowMatrix (blen, plen, rowptr, colidx);
    ICompRowMatrix *count = new ICompRowMatrix (blen, plen, rowptr, colidx);
    
    delete []npx;
    delete []rowptr;
    delete []colidx;

    // pass 3: fill the matrix
    for (j = 0; j < bdim[1]-1; j++) {
	y0 = bsize[1]*j;
	y1 = y0 + bsize[1];
	for (i = 0; i < bdim[0]-1; i++) {
	    x0 = bsize[0]*i;
	    x1 = x0 + bsize[0];
		
	    for (sj = 0; sj < sublen; sj++) {
		p[1] = y0 + (sj+0.5)*dy;
		pj = (int)(p[1]/psize[1]);
		p[1] += bbmin[1];
		for (si = 0; si < sublen; si++) {
		    p[0] = x0 + (si+0.5)*dx;
		    pi = (int)(p[0]/psize[0]);
		    p[0] += bbmin[0];
		    idx_j = pi + pj*wdim[0];
		    for (jj = 0; jj < 2; jj++) {
			for (ii = 0; ii < 2; ii++) {
			    idx_i = (i+ii) + (j+jj)*bdim[0];
			    v = Value_nomask(p, idx_i, false);
			    (*bvw)(idx_i, idx_j) += v;
			    (*count)(idx_i, idx_j) += 1;
			}
		    }
		}
	    }
	}
    }

    double *pbvw = bvw->ValPtr();
    int *pcount = count->ValPtr();

    for (i = 0; i < nz; i++) {
	if (pcount[i]) {
	    pbvw[i] *= scale[i]/pcount[i];
	}
    }
    delete []scale;
    
    return bvw;
}

// ==========================================================================
// Assemble single-element contribution for element "el" into global
// system matrix M, where coefficients (where applicable) are given in pixel
// basis. "mode" is the integration type.

void Raster_Pixel2::AddToElMatrix_tri (int el, RGenericSparseMatrix &M,
    const RVector *pxcoeff, int mode) const
{
    // All "parameter-free" integrals can be done in the standard way
    if (mode != ASSEMBLE_PFF && mode != ASSEMBLE_PDD) {
	::AddToElMatrix (*meshptr, el, M, pxcoeff, mode);
	return;
    }

    dASSERT(pxcoeff, "AddToElMatrix: requires a parameter pointer");

    Element *pel = meshptr->elist[el];
    int nnode = pel->nNode();
    int i, j, k, m, ii, jj, kk, idx_i, idx_j, idx_k, nv;
    int imin, imax, jmin, jmax;
    double b, v;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];

    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    // SPLIT ELEMENT

    // compute bounding box
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
    RVector fun;

    for (i = imin; i <= imax; i++) {
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
		    
		switch (mode) {
		case ASSEMBLE_PFF:
		    for (m = 0; m < np; m++) {
			v = t3.DetJ(absc[m], &n3) * wght[m];
			Point glob = t3.Global(n3, absc[m]);
			Point loc = pel->Local(meshptr->nlist, glob);
			fun = pel->LocalShapeF (loc);
			for (kk = 0; kk < 4; kk++) {
			    idx_k = i + kk%2 + (j+kk/2)*bdim[0];
			    b = v * (*pxcoeff)[idx_k] *
				Value_nomask (glob, idx_k, false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				for (jj = 0; jj < fun.Dim(); jj++) {
				    idx_j = pel->Node[jj];
				    M(idx_i,idx_j) += b * fun[ii] * fun[jj];
				}
			    }
			}
		    }
		    break;
		case ASSEMBLE_PDD: {
		    RDenseMatrix IJ(2,2), der(2,3);
		    double dd;

		    // assumes that the derivative is constant over the element
		    RVector vtxfun[3];
		    for (ii = 0; ii < 3; ii++)
			vtxfun[ii] = pel->GlobalShapeF(meshptr->nlist,
						       n3[ii]);
		    for (ii = 0; ii < 3; ii++) {
			der(0,ii) = vtxfun[1][ii]-vtxfun[0][ii];
			der(1,ii) = vtxfun[2][ii]-vtxfun[0][ii];
		    }

		    for (m = 0; m < np; m++) {
			Point glob = t3.Global(n3, absc[m]);
			Point loc = pel->Local(meshptr->nlist, glob);
			fun = pel->LocalShapeF (loc);
			v = t3.IJacobian (absc[m], &n3, IJ)
			    * wght[m];
			for (kk = 0; kk < 4; kk++) {
			    idx_k = i + kk%2 + (j+kk/2)*bdim[0];
			    b = v * (*pxcoeff)[idx_k] *
				Value_nomask (glob, idx_k, false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				for (jj = 0; jj < fun.Dim(); jj++) {
				    idx_j = pel->Node[jj];
				    dd =
					((IJ(0,0)*der(0,ii)+IJ(0,1)*der(1,ii)) *
					 (IJ(0,0)*der(0,jj)+IJ(0,1)*der(1,jj)) +
					 (IJ(1,0)*der(0,ii)+IJ(1,1)*der(1,ii)) *
					 (IJ(1,0)*der(0,jj)+IJ(1,1)*der(1,jj)));
				    M(idx_i,idx_j) += b*dd;
				}
			    }
			}
		    }
		    } break;
		default:
		    xERROR("Integration type not supported");
		}
	    }
	}
    }
}
