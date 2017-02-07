// Implementation of those routines in Raster_CPixel specific to triangular
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tri_qr.h"
#include "sh.h"

double TriIntF (double *x, double *y, double *z);

RCompRowMatrix *Raster_CPixel::CreateBuv_tri () const
{
    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, idx_i, idx_j;
    int imin, imax, jmin, jmax;
    double djac;

    double xrange = bbsize[0];
    double yrange = bbsize[1];
    
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
	imin = max (0, (int)floor(bdim[0] * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-1, (int)floor(bdim[0]*(exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor(bdim[1] * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-1, (int)floor(bdim[1]*(eymax-bbmin[1])/yrange));
	for (i = 0; i < pel->nNode(); i++) {
	    int nidx = pel->Node[i];
	    if (imin < nimin[nidx]) nimin[nidx] = imin;
	    if (imax > nimax[nidx]) nimax[nidx] = imax;
	    if (jmin < njmin[nidx]) njmin[nidx] = jmin;
	    if (jmax > njmax[nidx]) njmax[nidx] = jmax;
	}
    }

    int *rowptr = new int [n+1];
    rowptr[0] = 0;
    for (r = 0; r < n; r++) {
	int nentry = (nimax[r]-nimin[r]+1)*(njmax[r]-njmin[r]+1);
	rowptr[r+1] = rowptr[r]+nentry;
    }
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (r = k = 0; r < n; r++) {
	for (j = njmin[r]; j <= njmax[r]; j++) {
	    for (i = nimin[r]; i <= nimax[r]; i++) {
		colidx[k++] = i + j*bdim[0];
	    }
	}
    }
    RCompRowMatrix *buv = new RCompRowMatrix (n, blen, rowptr, colidx);

    delete []rowptr;
    delete []colidx;
    delete []nimin;
    delete []nimax;
    delete []njmin;
    delete []njmax;

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	bool numerical_quadrature = (pel->Type() != ELID_TRI3);
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
	imin = max (0, (int)floor(bdim[0] * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-1, (int)floor(bdim[0]*(exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor(bdim[1] * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-1, (int)floor(bdim[1]*(eymax-bbmin[1])/yrange));

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
		    
		    if (t3.Size() <= 0.0) continue;

		    if (numerical_quadrature) {
			// map quadrature points into global frame
			for (m = 0; m < np; m++) {
			    djac = t3.DetJ(absc[m], &n3);
			    Point glob = t3.Global(n3, absc[m]);
			    Point loc = pel->Local(meshptr->nlist, glob);
			    fun = pel->LocalShapeF (loc);
			    v = wght[m] * djac;
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				(*buv)(idx_i, idx_j) += v*fun[ii];
			    }
			}

		    } else { // exact integral for linear shape funct

			RVector nfun[3];
			double x[3], y[3], z[3];
			for (m = 0; m < 3; m++) {
			    // loop over nodes of subdivided triangle
			    Point loc = pel->Local(meshptr->nlist, n3[m]);
			    nfun[m] = pel->LocalShapeF(loc);
			    x[m] = n3[m][0];
			    y[m] = n3[m][1];
			}
			for (ii = 0; ii < 3; ii++) {
			    // loop over nodes of base triangle
			    idx_i = pel->Node[ii];
			    for (m = 0; m < 3; m++) {
				z[m] = nfun[m][ii];
			    }
			    (*buv)(idx_i, idx_j) += TriIntF(x,y,z);
			}
		    }
		}
	    }
    }

    buv->Shrink();
    return buv;
}


int Raster_CPixel::SutherlandHodgman (int el, int xgrid, int ygrid,
    Point *clip_poly, int npoly) const
{
    int i;

    Point *list = new Point[npoly];

    double xmin = bbmin[0] +
	bbsize[0]*(double)xgrid/(double)bdim[0];
    double xmax = bbmin[0] +
	bbsize[0]*(double)(xgrid+1)/(double)bdim[0];
    double ymin = bbmin[1] +
	bbsize[1]*(double)ygrid/(double)bdim[1];
    double ymax = bbmin[1] +
	bbsize[1]*(double)(ygrid+1)/(double)bdim[1];

    Element *pel = meshptr->elist[el];
    int nv = pel->nNode();
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
// Assemble single-element contribution for element "el" into global
// system matrix M, where coefficients (where applicable) are given in pixel
// basis. "mode" is the integration type.

void Raster_CPixel::AddToElMatrix_tri (int el, RGenericSparseMatrix &M,
    const RVector *pxcoeff, int mode) const
{
    // All "parameter-free" integrals can be done in the standard way
    if (mode != ASSEMBLE_PFF && mode != ASSEMBLE_PDD) {
	::AddToElMatrix (*meshptr, el, M, pxcoeff, mode);
	return;
    }

    dASSERT(pxcoeff, "AddToElMatrix: requires a parameter pointer");

    Element *pel = meshptr->elist[el];
    int i, j, k, m, ii, jj, idx_i, idx_j, idx_k, nv;
    int imin, imax, jmin, jmax;
    double b, v;

    double xrange = bbsize[0];
    double yrange = bbsize[1];

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
    imin = max (0, (int)floor(bdim[0] * (exmin-bbmin[0])/xrange));
    imax = min (bdim[0]-1, (int)floor(bdim[0] * (exmax-bbmin[0])/xrange));
    jmin = max (0, (int)floor(bdim[1] * (eymin-bbmin[1])/yrange));
    jmax = min (bdim[1]-1, (int)floor(bdim[1] * (eymax-bbmin[1])/yrange));

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
	    idx_k = i + j*bdim[0];

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
			b = v * (*pxcoeff)[idx_k];
			Point glob = t3.Global(n3, absc[m]);
			Point loc = pel->Local(meshptr->nlist, glob);
			fun = pel->LocalShapeF (loc);
			for (ii = 0; ii < fun.Dim(); ii++) {
			    idx_i = pel->Node[ii];
			    for (jj = 0; jj < fun.Dim(); jj++) {
				idx_j = pel->Node[jj];
				M(idx_i,idx_j) += b * fun[ii] * fun[jj];
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
			b = v * (*pxcoeff)[idx_k];
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
		    } break;
		default:
		    xERROR("Integration type not supported");
		}
	    }
	}
    }
}

double TriIntF (double *x, double *y, double *z)
{
    return (x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1])) *
	(z[0]+z[1]+z[2]) / 6.0;
}
