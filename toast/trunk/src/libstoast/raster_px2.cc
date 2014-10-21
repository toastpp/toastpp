#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "timing.h"

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
    tic();
    Buu = meshptr->MassMatrix();
    double t_uu = toc();
    tic();
    Bvv = CreatePixelMassmat();
    double t_vv = toc();
    tic();
    Buv = CreateMixedMassmat();
    double t_uv = toc();

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
}

// ==========================================================================

void Raster_Pixel2::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    double tol = map_tol;
    int nit = PCG (*Buu, Ax(*Buv,bvec), mvec, tol, Buu_precon);
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
// Creates a mixed mass matrix 

#if THREAD_LEVEL==2
struct CREATEMIXEDMASSMAT_PASS2_THREADDATA {
    const Raster_Pixel2 *raster;
    RCompRowMatrix *Buv;
    const Point *bbmin;
    const Point *bbmax;
};

void CreateMixedMassmat_tri_pass2_engine (task_data *td)
{
    int itask = td->proc;
    int ntask = td->np;
    CREATEMIXEDMASSMAT_PASS2_THREADDATA *thdata =
	(CREATEMIXEDMASSMAT_PASS2_THREADDATA*)td->data;
    const Raster_Pixel2 *raster = thdata->raster;
    const Mesh &mesh = raster->mesh();
    const IVector &bdim = raster->BDim();
    int nel = mesh.elen();
    int el0 = (itask*nel)/ntask;
    int el1 = ((itask+1)*nel)/ntask;
    RCompRowMatrix *Buv = thdata->Buv;
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
	Element *pel = mesh.elist[el];

	// element bounding box
	double exmin = mesh.nlist[pel->Node[0]][0];
	double exmax = mesh.nlist[pel->Node[0]][0];
	double eymin = mesh.nlist[pel->Node[0]][1];
	double eymax = mesh.nlist[pel->Node[0]][1];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, mesh.nlist[pel->Node[j]][0]);
	    exmax = max (exmax, mesh.nlist[pel->Node[j]][0]);
	    eymin = min (eymin, mesh.nlist[pel->Node[j]][1]);
	    eymax = max (eymax, mesh.nlist[pel->Node[j]][1]);
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
			Point loc = pel->Local(mesh.nlist, glob);
			pel->LocalShapeF (loc, &fun);
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
    double *bval = Buv->ValPtr();
    const double *blocval = Buv_local.ValPtr();
    int nz = Buv->nVal();
    Task::UserMutex_lock();
    for (i = 0; i < nz; i++)
	bval[i] += blocval[i];
    Task::UserMutex_unlock();
}
#endif

RCompRowMatrix *Raster_Pixel2::CreateMixedMassmat () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateMixedMassmat_tri();
    case ELID_TET4:
	return CreateMixedMassmat_tet4();
    default:
	xERROR("Raster_Pixel2: Unsupported element type");
	return 0;
    }
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

// ==========================================================================
// Creates a mixed mass matrix 

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
