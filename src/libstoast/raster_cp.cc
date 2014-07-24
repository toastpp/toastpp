#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "timing.h"

using namespace std;

// ==========================================================================
// class Raster_CubicPixel

Raster_CubicPixel::Raster_CubicPixel (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, RDenseMatrix *bb)
    : Raster (_bdim, _gdim, mesh, bb)
{
    int i, j, k;
    Point p(dim);

    // inverse basis spacing
    ibgrid.New (dim);
    for (int i = 0; i < dim; i++) 
	ibgrid[i] = (bdim[i]-1)/(bbmax[i]-bbmin[i]);

    // inverse grid spacing
    iggrid.New (dim);
    for (int i = 0; i < dim; i++) 
	iggrid[i] = (gdim[i]-1)/(bbmax[i]-bbmin[i]);

    // set up grid->basis mapping
    G = new RCompRowMatrix (blen, glen);
    RVector grow(glen);

    int *g_rownz = new int[blen];
    int **g_colidx = new int*[blen];
    double **g_val = new double*[blen];
    int nz, nz_tot = 0;
    double v, gsum;

    tic();

    double **swght = new double*[dim];
    for (i = 0; i < dim; i++) {
	swght[i] = new double[gdim[i]];
    }
    IVector imin(3), imax(3), basisgrd;
    int grd[3];
    int jx, jy, jz, zofs, yofs;

    for (i = 0; i < blen; i++) {
	nz = 0;
	gsum = 0.0;
	grow.Clear();
	GetBasisIndices (i, basisgrd);
	imin = 100000; imax = -1;
	if (dim < 3) imin[2] = imax[2] = 0;
	// find grid support range for this basis function
	for (k = 0; k < dim; k++) {
	    for (j = 0; j < gdim[k]; j++) {
		double dr = fabs (basisgrd[k] - (j/iggrid[k])*ibgrid[k]);
		if (fabs(dr) < 2.0) {
		    double dm2 = dr-2.0;
		    double dm1 = dr-1.0;
		    double dp1 = dr+1.0;
		    double dp2 = dr+2.0;
		    swght[k][j] = (fabs(dm2*dm2*dm2)
				   - 4.0*fabs(dm1*dm1*dm1)
				   + 6.0*fabs(dr*dr*dr)
				   - 4.0*fabs(dp1*dp1*dp1)
				   + fabs(dp2*dp2*dp2)) / 12.0;
		    if (j < imin[k]) imin[k] = j;
		    if (j > imax[k]) imax[k] = j;
		}
	    }
	}
	for (jz = imin[2]; jz <= imax[2]; jz++) {
	    zofs = gdim[0]*gdim[1]*jz;
	    grd[2] = jz;
	    for (jy = imin[1]; jy <= imax[1]; jy++) {
		yofs = zofs + jy*gdim[0];
		grd[1] = jy;
		for (jx = imin[0]; jx <= imax[0]; jx++) {
		    j = yofs + jx;
		    grd[0] = jx;

		    if (gelref[j] >= 0) {
			v = 1.0;
			for (k = 0; k < dim; k++)
			    v *= swght[k][grd[k]];
			if (v > 0) {
			    grow[j] = v;
			    gsum += v;
			    nz++;
			}
		    }
		}
	    }
	}

	if ((g_rownz[i] = nz)) {
	    g_colidx[i] = new int[nz];
	    g_val[i] = new double[nz];
	    nz_tot += nz;
	    for (j = k = 0; j < glen; j++) {
		if (grow[j]) {
		    g_colidx[i][k] = j;
		    g_val[i][k] = grow[j]/gsum;
		    k++;
		}
	    }
	}
    }
    idxtype *rowptr = new idxtype[blen+1];
    idxtype *colidx = new idxtype[nz_tot];
    double *val = new double[nz_tot];
    rowptr[0] = 0;
    for (i = j = 0; i < blen; i++) {
	if (g_rownz[i]) {
	    memcpy (colidx+j, g_colidx[i], g_rownz[i]*sizeof(idxtype));
	    memcpy (val+j, g_val[i], g_rownz[i]*sizeof(double));
	    j += g_rownz[i];
	    delete []g_val[i];
	    delete []g_colidx[i];
	}
	rowptr[i+1] = rowptr[i] + g_rownz[i];
    }
    delete []g_val;
    delete []g_colidx;
    delete []g_rownz;
    ((RCompRowMatrix*)G)->Initialise(rowptr, colidx, val);
    delete []rowptr;
    delete []colidx;
    delete []val;    

    // set up basis->grid mapping
    GI = new RCompRowMatrix (glen, blen);
    RVector brow(blen);
    tic();

    int *gi_rownz = new int[glen];
    int **gi_colidx = new int*[glen];
    double **gi_val = new double*[glen];
    nz_tot = 0;

    for (j = 0; j < glen; j++) {
	if (gelref[j] >= 0) {
	    nz = 0;
	    IVector &grd = GetGridIndices (j);
	    for (k = 0; k < dim; k++)
		p[k] = bbmin[k] + grd[k]/iggrid[k];
	    for (i = 0; i < blen; i++) {
		brow[i] = Value (i, p);
		if (brow[i]) nz++;   
	    }
	    gi_rownz[j] = nz;
	    gi_colidx[j] = new int[nz];
	    gi_val[j] = new double[nz];
	    nz_tot += nz;
	    for (i = k = 0; i < blen; i++) {
		if (brow[i]) {
		    gi_colidx[j][k] = i;
		    gi_val[j][k] = brow[i];
		    k++;
		}
	    }
	} else {
	    gi_rownz[j] = 0;
	}
    }
    rowptr = new idxtype[glen+1];
    colidx = new idxtype[nz_tot];
    val = new double[nz_tot];
    rowptr[0] = 0;
    for (j = i = 0; j < glen; j++) {
	if (gi_rownz[j]) {
	    memcpy (colidx+i, gi_colidx[j], gi_rownz[j]*sizeof(idxtype));
	    memcpy (val+i, gi_val[j], gi_rownz[j]*sizeof(double));
	    i += gi_rownz[j];
	    delete []gi_val[j];
	    delete []gi_colidx[j];
	}
	rowptr[j+1] = rowptr[j] + gi_rownz[j];
    }
    delete []gi_val;
    delete []gi_colidx;
    delete []gi_rownz;
    ((RCompRowMatrix*)GI)->Initialise(rowptr, colidx, val);
    delete []rowptr;
    delete []colidx;
    delete []val;

    for (i = 0; i < dim; i++)
	delete []swght[i];
    delete []swght;

    // set up basis->solution mapping
    basis2sol.New (blen);
    rowptr = ((RCompRowMatrix*)G)->rowptr;
    for (i = slen = 0; i < blen; i++) {
	basis2sol[i] = (rowptr[i+1] > rowptr[i] ? slen++ : -1);
    }
    sol2basis.New (slen);
    for (i = slen = 0; i < blen; i++) {
	if (rowptr[i+1] > rowptr[i]) sol2basis[slen++] = i;
    }

    if (toastVerbosity > 0) {
        cout << "--> Type............" << (mesh->Dimension() == 2 ? "Bi":"Tri")
	     << "-cubic" << endl;
    }
}

Raster_CubicPixel::~Raster_CubicPixel ()
{
    delete G;
    delete GI;
}

double Raster_CubicPixel::Value (int basisidx, const Point &p)
{
    int i;
    IVector grd;
    GetBasisIndices (basisidx, grd);
    
    double dr, s, val = 1.0;
    for (i = 0; i < dim; i++) {
	dr = fabs (grd[i]/*-1*/ + (bbmin[i]-p[i])*ibgrid[i]);
	if (dr >= 2.0) return 0.0; // point outside basis function support
	double dm2 = dr-2.0;
	double dm1 = dr-1.0;
	double dp1 = dr+1.0;
	double dp2 = dr+2.0;
	s = (fabs(dm2*dm2*dm2) - 4.0*fabs(dm1*dm1*dm1) + 6.0*fabs(dr*dr*dr)
	  - 4.0*fabs(dp1*dp1*dp1) + fabs(dp2*dp2*dp2)) / 12.0;
	val *= s;
    }
    return val;
}

// ==========================================================================

double Raster_CubicPixel::Value (const IVector &basisgrd, const Point &p)
{
    int i;
    static double ddr[3];

    for (i = 0; i < dim; i++) {
	ddr[i] = fabs (basisgrd[i]/*-1*/ + (bbmin[i]-p[i])*ibgrid[i]);
	if (ddr[i] >= 2.0) return 0.0; // point outside basis function support
    }

    double dr, s, val = 1.0;
    for (i = 0; i < dim; i++) {
	dr = ddr[i];
	double dm2 = dr-2.0;
	double dm1 = dr-1.0;
	double dp1 = dr+1.0;
	double dp2 = dr+2.0;
	s = (fabs(dm2*dm2*dm2) - 4.0*fabs(dm1*dm1*dm1) + 6.0*fabs(dr*dr*dr)
	  - 4.0*fabs(dp1*dp1*dp1) + fabs(dp2*dp2*dp2)) / 12.0;
	val *= s;
    }
    return val;    
}

// ==========================================================================

double Raster_CubicPixel::Value_nomask (const Point &p, int i, bool is_solidx)
    const
{
    // to be done
    return 0.0;
}

// ==========================================================================

void Raster_CubicPixel::Map_GridToBasis (const RVector &gvec, RVector &bvec)
    const
{
    G->Ax (gvec, bvec);
}

// ==========================================================================

void Raster_CubicPixel::Map_GridToBasis (const CVector &gvec, CVector &bvec)
    const
{
    ((RCompRowMatrix*)G)->Ax_cplx (gvec, bvec);
}

// ==========================================================================

void Raster_CubicPixel::Map_BasisToGrid (const RVector &bvec, RVector &gvec)
    const
{
    GI->Ax (bvec, gvec);
}

// ==========================================================================

void Raster_CubicPixel::Map_BasisToGrid (const CVector &bvec, CVector &gvec)
    const
{
    ((RCompRowMatrix*)GI)->Ax_cplx (bvec, gvec);
}
