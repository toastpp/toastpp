#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

using namespace std;

// ==========================================================================
// class Raster_Blob

Raster_Blob::Raster_Blob (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, RDenseMatrix *bb)
    : Raster (_bdim, _bdim, mesh, bb), sup(_sup)
      // for now, we ignore the intermediate gdim grid
{
    const double eps = 1e-6;
    int i, j, k, d, idx, nd;

    nodevals = 0;

    // inverse grid spacing
    igrid.New(dim);
    for (d = 0; d < dim; d++)
	igrid[d] = (bdim[d]-1.0)/(bbmax[d]-bbmin[d]);

    // update grid size to accommodate padding around bounding box
    npad = (int)ceil (sup-eps); // number of padding layers
    for (d = 0, blen = 1; d < dim; d++)
	blen *= (bdim[d] += 2*npad);
    int nx = bdim[0], ny = bdim[1], nz = (dim > 2 ? bdim[2]:1);

    for (d = 0; d < dim; d++) gdim[d] = bdim[d]; // for now
    glen = blen; // for now

    // bounding box including padding
    bbmin_pad = bbmin;
    bbmax_pad = bbmax;
    for (d = 0; d < dim; d++) {
    	bbmin_pad[d] -= npad/igrid[d];
    	bbmax_pad[d] += npad/igrid[d];
    }
    delete []belref;
    delete []gelref;
    belref = GenerateElementPixelRef (*meshptr, bdim, &bbmin_pad, &bbmax_pad);
    gelref = GenerateElementPixelRef (*meshptr, gdim, &bbmin_pad, &bbmax_pad);

    // re-calculate the basis2sol and sol2basis permutation arrays
    basis2sol.New(blen);
    RVector r(3);
    double dx, rad2, sup2 = sup*sup;

    for (k = slen = idx = 0; k < nz; k++) {
	r[2] = (double)k;
	for (j = 0; j < ny; j++) {
	    r[1] = (double)j;
	    for (i = 0; i < nx; i++) {
		r[0] = (double)i;
		basis2sol[idx] = -1;
		for (nd = 0; nd < mesh->nlen(); nd++) {
		    for (d = 0, rad2 = 0.0; d < dim; d++) {
			dx = (mesh->nlist[nd][d]-bbmin[d])*igrid[d] +
			    npad - r[d];
			rad2 += dx*dx;
		    }
		    if (rad2 < sup2) {
			basis2sol[idx] = slen++;
			break;
		    }
		}
		idx++;
	    }
	}
    }
    sol2basis.New(slen);
    for (i = idx = 0; i < blen; i++)
	if (basis2sol[i] >= 0)
	    sol2basis[idx++] = i;
}

Raster_Blob::~Raster_Blob ()
{
    if (nodevals)
	delete nodevals;
}

void Raster_Blob::ComputeNodeValues ()
{
    int i, j, idx, nz;
    int nlen = meshptr->nlen();
    RVector nv(slen);

    if (nodevals) delete nodevals;
    nodevals = 0;

    // pass 1: generate row ptr list
    int *rowptr = new int[nlen+1];
    rowptr[0] = 0;
    for (i = 0; i < nlen; i++) {
	NodeValues (i, nv);
	for (j = nz = 0; j < slen; j++)
	    if (nv[j]) nz++;
	rowptr[i+1] = rowptr[i] + nz;
    }

    // pass 2: generate col idx and data lists
    int *colidx = new int[rowptr[nlen]];
    double *data = new double[rowptr[nlen]];
    for (i = idx = 0; i < nlen; i++) {
	NodeValues (i, nv);
	for (j = 0; j < slen; j++) 
	    if (nv[j]) {
		colidx[idx] = j;
		data[idx] = nv[j];
		idx++;
	    }
    }

    nodevals = new RCompRowMatrix (nlen, slen, rowptr, colidx, data);

    delete []rowptr;
    delete []colidx;
    delete []data;
}

double Raster_Blob::Value (const Point &p, int i, bool is_solidx) const
{
    double v = Value_nomask (p, i, is_solidx);
    return (v && meshptr->ElFind (p) >= 0 ? v : 0.0);
}

void Raster_Blob::NodeValues (int node, RVector &nv) const
{
    if (nodevals) {
	nv = nodevals->Row (node);
    } else {
	for (int i = 0; i < slen; i++)
	    nv[i] = Value_nomask (meshptr->nlist[node], i);
    }
}

void Raster_Blob::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    int i;
    if (svec.Dim() != slen) svec.New(slen);
    RVector nv(slen), fsum(slen), sum(slen);
    for (i = 0; i < meshptr->nlen(); i++) {
	NodeValues (i, nv);
	fsum += nv * mvec[i];
	sum += nv;
    }
    for (i = 0; i < slen; i++)
	svec[i] = (sum[i] ? fsum[i]/sum[i] : 0.0);
}

void Raster_Blob::Map_MeshToSol (const CVector &mvec, CVector &svec) const
{
    int i, j;
    if (svec.Dim() != slen) svec.New(slen);
    RVector nv(slen), sum(slen);
    CVector fsum(slen);
    for (i = 0; i < meshptr->nlen(); i++) {
	NodeValues (i, nv);
	for (j = 0; j < slen; j++)
	    fsum[j] += mvec[i] * nv[j];
	sum += nv;
    }
    for (i = 0; i < slen; i++)
	svec[i] = (sum[i] ? fsum[i]/sum[i] : 0.0);
}

void Raster_Blob::Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
{
    int i, j;
    if (bvec.Dim() != blen) bvec.New(blen);
    RVector nv(slen), fsum(slen), sum(slen);
    for (i = 0; i < meshptr->nlen(); i++) {
	NodeValues (i, nv);
	fsum += nv * mvec[i];
	sum += nv;
    }
    for (j = 0; j < blen; j++) {
	i = basis2sol[j];
	bvec[j] = (i >= 0 && sum[i] ? fsum[i]/sum[i] : 0.0);
    }
}

void Raster_Blob::Map_MeshToBasis (const CVector &mvec, CVector &bvec) const
{
    int i, j;
    if (bvec.Dim() != blen) bvec.New(blen);
    RVector nv(slen), sum(slen);
    CVector fsum(slen);
    for (i = 0; i < meshptr->nlen(); i++) {
	NodeValues (i, nv);
	for (j = 0; j < slen; j++)
	    fsum[j] += mvec[i] * nv[j];
	sum += nv;
    }
    for (j = 0; j < blen; j++) {
	i = basis2sol[j];
	bvec[j] = (i >= 0 && sum[i] ? fsum[i]/sum[i] : 0.0);
    }
}

void Raster_Blob::Map_MeshToGrid (const RVector &mvec, RVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_MeshToBasis (mvec, gvec);
}

void Raster_Blob::Map_MeshToGrid (const CVector &mvec, CVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_MeshToBasis (mvec, gvec);
}

void Raster_Blob::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    int i, j;
    if (mvec.Dim() != meshptr->nlen()) mvec.New(meshptr->nlen());
    for (i = 0; i < meshptr->nlen(); i++) {
	RVector nv = NodeValues(i);
	double sum = 0.0;
	for (j = 0; j < slen; j++)
	    sum += bvec[sol2basis[j]] * nv[j];
	mvec[i] = sum;
    }
}

void Raster_Blob::Map_BasisToMesh (const CVector &bvec, CVector &mvec) const
{
    int i, j;
    if (mvec.Dim() != meshptr->nlen()) mvec.New(meshptr->nlen());
    for (i = 0; i < meshptr->nlen(); i++) {
	RVector nv = NodeValues(i);
	std::complex<double> sum = std::complex<double>(0,0);
	for (j = 0; j < slen; j++)
	    sum += bvec[sol2basis[j]] * nv[j];
	mvec[i] = sum;
    }
}

void Raster_Blob::Map_GridToMesh (const RVector &gvec, RVector &mvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_BasisToMesh (gvec, mvec);
}

void Raster_Blob::Map_GridToMesh (const CVector &gvec, CVector &mvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_BasisToMesh (gvec, mvec);
}

void Raster_Blob::Map_SolToMesh (const RVector &bvec, RVector &mvec) const
{
    int i, j;
    if (mvec.Dim() != meshptr->nlen()) mvec.New(meshptr->nlen());
    for (i = 0; i < meshptr->nlen(); i++) {
	RVector nv = NodeValues(i);
	double sum = 0.0;
	for (j = 0; j < slen; j++)
	    sum += bvec[j] * nv[j];
	mvec[i] = sum;
    }
}

void Raster_Blob::Map_SolToMesh (const CVector &bvec, CVector &mvec) const
{
    int i, j;
    if (mvec.Dim() != meshptr->nlen()) mvec.New(meshptr->nlen());
    for (i = 0; i < meshptr->nlen(); i++) {
	RVector nv = NodeValues(i);
	std::complex<double> sum = std::complex<double>(0,0);
	for (j = 0; j < slen; j++)
	    sum += bvec[j] * nv[j];
	mvec[i] = sum;
    }
}

void Raster_Blob::Map_GridToBasis (const RVector &gvec, RVector &bvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    bvec = gvec;
}

void Raster_Blob::Map_GridToBasis (const CVector &gvec, CVector &bvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    bvec = gvec;
}

void Raster_Blob::Map_BasisToGrid (const RVector &bvec, RVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    gvec = bvec;
}

void Raster_Blob::Map_BasisToGrid (const CVector &bvec, CVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    gvec = bvec;
}

void Raster_Blob::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    if (bvec.Dim() != blen) bvec.New(blen);
    bvec.Clear();
    for (int i = 0; i < slen; i++)
	bvec[sol2basis[i]] = svec[i];
}

void Raster_Blob::Map_SolToBasis (const CVector &svec, CVector &bvec) const
{
    if (bvec.Dim() != blen) bvec.New(blen);
    bvec.Clear();
    for (int i = 0; i < slen; i++)
	bvec[sol2basis[i]] = svec[i];
}

void Raster_Blob::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    if (svec.Dim() != slen) svec.New(slen);
    for (int i = 0; i < slen; i++)
	svec[i] = bvec[sol2basis[i]];
}

void Raster_Blob::Map_BasisToSol (const CVector &bvec, CVector &svec) const
{
    if (svec.Dim() != slen) svec.New(slen);
    for (int i = 0; i < slen; i++)
	svec[i] = bvec[sol2basis[i]];
}

void Raster_Blob::Map_SolToGrid (const RVector &svec, RVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_SolToBasis (svec, gvec);
}

void Raster_Blob::Map_SolToGrid (const CVector &svec, CVector &gvec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_SolToBasis (svec, gvec);
}

void Raster_Blob::Map_GridToSol (const RVector &gvec, RVector &svec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_BasisToSol (gvec, svec);
}

void Raster_Blob::Map_GridToSol (const CVector &gvec, CVector &svec) const
{
    // for now, we ignore the intermediate basis and assume basis=grid
    Map_BasisToSol (gvec, svec);
}

double Raster_Blob::ComputeBasisScale ()
{
    RVector bunit (slen);
    RVector nf (meshptr->nlen());
    bunit = 1.0;
    Map_SolToMesh (bunit, nf);
    return 1.0/mean(nf);
}
