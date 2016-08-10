#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

// ==========================================================================
// STENCIL2D and STENCIL3D: (used by PCG_PRECON_SPARSEJTJ only): number of
// neighbours to include in the JTJ sparsity graph

#define STENCIL2D 5  // valid entries: 5, 9, 13
#define STENCIL3D 7  // valid entries: 7, 19, 27, 33

using namespace std;

static int POW2[] = {1,2,4,8,16,32,64,128,256,512,1024,2096,4192};

// ==========================================================================
// class Raster

Raster::Raster (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
    RDenseMatrix *bb)
{
    int i, j;

    // set up dimensions
    bdim = _bdim;      // native basis dimension
    gdim = _gdim;      // high-res grid dimension
    dim = bdim.Dim();
    meshptr = mesh;
    xASSERT(dim == gdim.Dim(),
	    "Basis and grid have incompatible dimensions");
    xASSERT(dim == meshptr->Dimension(),
	    "Basis and mesh have incompatible dimensions");

    // set up bounding box
    if (bb) {
	xASSERT(bb->nCols() == dim && bb->nRows() == 2,
		"Invalid bounding box dimensions");
	bbmin = bb->Row(0);
	bbmax = bb->Row(1);
    } else {
      meshptr->BoundingBox (bbmin, bbmax, 0 /*-1e-8*/);
    }
    bbsize = bbmax-bbmin;
    Point bbmin_pad = bbmin /*- 2e8*/;
    Point bbmax_pad = bbmax /*+ 2e8*/;

    for (i = 0, blen = glen = 1; i < dim; i++) {
	blen *= bdim[i];
	glen *= gdim[i];
    }

    belref = GenerateElementPixelRef (*meshptr, bdim, &bbmin, &bbmax);
    gelref = GenerateElementPixelRef (*meshptr, gdim, &bbmin, &bbmax);

    RVector s(glen);
    bsupport.New (blen);
    for (i = 0; i < glen; i++)
	if (gelref[i] >= 0) s[i] = 1.0;

    // map mask into user basis
    if (glen == blen) {
	bsupport = s;
    } else {
	SubsampleLinPixel (s, bsupport, gdim, bdim, 0);
    }

    // calculate default basis->solution mapping index list
    // (default is to use the basis voxel support array)

    basis2sol.New (blen);    
    for (i = slen = 0; i < blen; i++) {
        if (bsupport[i] > 0.0) basis2sol[i] = slen++;
	else                   basis2sol[i] = -1;
    }
    sol2basis.New (slen);
    for (i = slen = 0; i < blen; i++) {
	if (bsupport[i] > 0.0) sol2basis[slen++] = i;
    }

    // set up transformations between mesh and high-res grid
    B     = GridMapMatrix (*meshptr, gdim, &bbmin, &bbmax, gelref);
    BI    = NodeMapMatrix (*meshptr, gdim, &bbmin_pad, &bbmax_pad, gelref);
    ((RCompRowMatrix*)BI)->Shrink();

    // set up power-2 padding
    paddim.New(dim);
    for (i = 0, padlen = 1; i < dim; i++) {
	j = 0;
	while (POW2[j] < gdim[i])  // find power of 2 greate than dimension
	    j++;
	paddim[i] = POW2[j];
	padlen *= paddim[i];
    }

    if (toastVerbosity > 0) {
        cout << "Basis:" << endl;
	cout << "--> Grid size......." << blen << " [";
	for (i = 0; i < dim; i++)
	    cout << bdim[i] << (i==dim-1 ? ']':'x');
	cout << endl;
	cout << "--> Sol. size......." << slen << endl;
	cout << "--> Mesh size......." << meshptr->nlen() << endl;
	if (bb) {
	    cout << "--> Grid bounding box:" << endl;
	    for (i = 0; i < 2; i++) {
	        for (j = 0; j < 2; j++)
		    cout << "  " << bb->Get(i,j);
		cout << endl;
	    }
	}
    }
}

// ==========================================================================

Raster::~Raster()
{
    delete []belref;
    delete []gelref;
    delete B;
    delete BI;
}

// ==========================================================================

RDenseMatrix Raster::BoundingBox () const
{
    RDenseMatrix bb(2,dim);
    for (int i = 0; i < dim; i++) {
	bb(0,i) = bbmin[i];
	bb(1,i) = bbmax[i];
    }
    return bb;
}

// ==========================================================================

void Raster::BasisVoxelPositions (RDenseMatrix &pos) const
{
    GenerateVoxelPositions (*meshptr, bdim, &bbmin, &bbmax, pos);
}

// ==========================================================================

void Raster::SolutionVoxelPositions (RDenseMatrix &pos) const
{
    RDenseMatrix bpos;
    GenerateVoxelPositions (*meshptr, bdim, &bbmin, &bbmax, bpos);
    pos.New (slen,dim);
    int i, j, k;
    for (j = k = 0; j < blen; j++) {
        if (basis2sol[j] >= 0) {
	    for (i = 0; i < dim; i++)
	        pos(k,i) = bpos(j,i);
	    k++;
	}
    }
}

// ==========================================================================

double Raster::Value (const Point &p, int i, bool is_solidx) const
{
    double v = Value_nomask (p, i, is_solidx);
    return (v && meshptr->ElFind(p) >= 0 ? v : 0.0);
}

// ==========================================================================

double Raster::Value (const Point &p, const RVector &coeff, bool mask) const
{
    bool is_solidx = true;
    int len = coeff.Dim();
    if (len == blen)
	is_solidx = false;
    else if (len != slen)
	xERROR("Invalid length of coefficient vector.");

    double sum = 0.0;
    if (!mask || meshptr->ElFind(p) >= 0) {
	for (int i = 0; i < len; i++)
	    sum += Value_nomask(p, i, is_solidx) * coeff[i];
    }
    return sum;
}

// ==========================================================================

void Raster::Map_MeshToGrid (const RVector &mvec, RVector &gvec) const
{
    B->Ax (mvec, gvec);
}

// ==========================================================================

void Raster::Map_MeshToGrid (const CVector &mvec, CVector &gvec) const
{
    ((RCompRowMatrix*)B)->Ax_cplx (mvec, gvec);
}

// ==========================================================================

void Raster::Map_MeshToGrid (const Solution &msol, Solution &gsol, bool mapall)
    const
{
    for (int i = 0; i < msol.nParam(); i++)
        if (mapall || msol.active[i])
  	    Map_MeshToGrid (msol.param[i], gsol.param[i]);
}

// ==========================================================================

void Raster::Map_GridToMesh (const RVector &gvec, RVector &mvec) const
{
    BI->Ax (gvec, mvec);
}

// ==========================================================================

void Raster::Map_GridToMesh (const CVector &gvec, CVector &mvec) const
{
    ((RCompRowMatrix*)BI)->Ax_cplx (gvec, mvec);
}

// ==========================================================================

void Raster::Map_GridToMesh (const Solution &gsol, Solution &msol, bool mapall)
    const
{
    for (int i = 0; i < gsol.nParam(); i++)
        if (mapall || gsol.active[i])
	    Map_GridToMesh (gsol.param[i], msol.param[i]);
}

// ==========================================================================

void Raster::Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
{
    RVector gvec(glen);
    Map_MeshToGrid (mvec, gvec);
    Map_GridToBasis (gvec, bvec);
}

// ==========================================================================

void Raster::Map_MeshToBasis (const CVector &mvec, CVector &bvec) const
{
    CVector gvec(glen);
    Map_MeshToGrid (mvec, gvec);
    Map_GridToBasis (gvec, bvec);
}

// ==========================================================================

void Raster::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    RVector gvec(glen);
    Map_BasisToGrid (bvec, gvec);
    Map_GridToMesh (gvec, mvec);
}

// ==========================================================================

void Raster::Map_BasisToMesh (const CVector &bvec, CVector &mvec) const
{
    CVector gvec(glen);
    Map_BasisToGrid (bvec, gvec);
    Map_GridToMesh (gvec, mvec);
}

// ==========================================================================

void Raster::Map_GridToBasis (const Solution &gsol, Solution &bsol,
    bool mapall) const
{
    for (int i = 0; i < gsol.nParam(); i++)
        if (mapall || gsol.active[i])
	    Map_GridToBasis (gsol.param[i], bsol.param[i]);
}

// ==========================================================================

void Raster::Map_BasisToGrid (const Solution &bsol, Solution &gsol,
    bool mapall) const
{
    for (int i = 0; i < bsol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_BasisToGrid (bsol.param[i], gsol.param[i]);
}

// ==========================================================================

void Raster::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    if (slen == blen) {
	svec = bvec;  // identity operation
    } else {
	if (svec.Dim() != slen) svec.New (slen);
	for (int i = 0; i < slen; i++)
	    svec[i] = bvec[sol2basis[i]];
    }
}

// ==========================================================================

void Raster::Map_BasisToSol (const CVector &bvec, CVector &svec) const
{
    if (slen == blen) {
	svec = bvec;  // identity operation
    } else {
	if (svec.Dim() != slen) svec.New (slen);
	for (int i = 0; i < slen; i++)
	    svec[i] = bvec[sol2basis[i]];
    }
}

// ==========================================================================

void Raster::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    if (slen == blen) {
	bvec = svec;  // identity operation
    } else {
	if (bvec.Dim() != blen) bvec.New (blen);
	else                    bvec.Clear();
	for (int i = 0; i < slen; i++)
	    bvec[sol2basis[i]] = svec[i];
    }
}

// ==========================================================================

void Raster::Map_SolToBasis (const CVector &svec, CVector &bvec) const
{
    if (slen == blen) {
	bvec = svec;  // identity operation
    } else {
	if (bvec.Dim() != blen) bvec.New (blen);
	else                    bvec.Clear();
	for (int i = 0; i < slen; i++)
	    bvec[sol2basis[i]] = svec[i];
    }
}

// ==========================================================================

void Raster::Map_SolToBasis (const Solution &ssol, Solution &bsol,
    bool mapall) const
{
    for (int i = 0; i < ssol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_SolToBasis (ssol.param[i], bsol.param[i]);
}

// ==========================================================================

void Raster::Map_SolToGrid (const RVector &svec, RVector &gvec) const
{
    RVector bvec;
    Map_SolToBasis (svec, bvec);
    Map_BasisToGrid (bvec, gvec);
}

// ==========================================================================

void Raster::Map_SolToGrid (const CVector &svec, CVector &gvec) const
{
    CVector bvec;
    Map_SolToBasis (svec, bvec);
    Map_BasisToGrid (bvec, gvec);
}

// ==========================================================================

void Raster::Map_GridToSol (const RVector &gvec, RVector &svec) const
{
    RVector bvec;
    Map_GridToBasis (gvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster::Map_GridToSol (const CVector &gvec, CVector &svec) const
{
    CVector bvec;
    Map_GridToBasis (gvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    RVector bvec;
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}

// ==========================================================================

void Raster::Map_SolToMesh (const CVector &svec, CVector &mvec) const
{
    CVector bvec;
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}

// ==========================================================================

void Raster::Map_SolToMesh (const Solution &ssol, Solution &msol,
    bool mapall) const
{
    for (int i = 0; i < ssol.nParam(); i++)
        if (mapall || ssol.active[i])
	    Map_SolToMesh (ssol.param[i], msol.param[i]);
}

// ==========================================================================

void Raster::Map_ActiveSolToMesh (const RVector &asol, Solution &msol) const
{
    int i, j = 0;
    for (i = 0; i < msol.nParam(); i++)
        if (msol.active[i]) {
	    RVector prm (asol, slen*j++, slen);
	    Map_SolToMesh (prm, msol.param[i]);
	}
}

// ==========================================================================

void Raster::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    RVector bvec;
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster::Map_MeshToSol (const CVector &mvec, CVector &svec) const
{
    CVector bvec;
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster::Map_MeshToSol (const Solution &msol, Solution &ssol, bool mapall)
    const
{
    for (int i = 0; i < msol.nParam(); i++)
        if (mapall || msol.active[i])
	    Map_MeshToSol (msol.param[i], ssol.param[i]);    
}

// ==========================================================================

void Raster::Sample (const RVector &svec, const IVector &grd, RVector &img)
    const
{
    dASSERT(svec.Dim() == slen, "Invalid basis vector length");
    dASSERT(grd.Dim() == dim, "Invalid grid dimension");
    int i, j, k, s, idx;
    int nx = grd[0];
    int ny = grd[1];
    int nz = (dim >= 3 ? grd[2] : 1);
    int ilen = nx*ny*nz;
    if (img.Dim() != ilen)
	img.New(ilen);
    double v;
    double dx = (bbmax[0]-bbmin[0])/(grd[0]-1.0);
    double dy = (bbmax[1]-bbmin[1])/(grd[1]-1.0);
    double dz = (dim >= 3 ? (bbmax[2]-bbmin[2])/(grd[2]-1.0) : 0.0);
    Point p(dim);

    for (k = idx = 0; k < nz; k++) {
	if (dim >= 3) p[2] = bbmin[2] + k*dz;
	for (j = 0; j < ny; j++) {
	    p[1] = bbmin[1] + j*dy;
	    for (i = 0; i < nx; i++) {
		p[0] = bbmin[0] + i*dx;
		for (s = 0, v = 0.0; s < slen; s++)
		    v += Value(p, s) * svec[s];
		img[idx++] = v;
	    }
	}
    }
}

// ==========================================================================

int Raster::GetSolIdx (int basisidx) const
{
    dASSERT (basisidx >= 0 && basisidx < blen, "Argument 1 index out of range");
    return basis2sol[basisidx];
}

// ==========================================================================

int Raster::GetSolIdx (const IVector &crd) const
{
    dASSERT (crd.Dim() == dim, "Argument 1 unexpected size");
    if (dim == 2)
	return basis2sol[crd[0] + bdim[0]*crd[1]];
    else
	return basis2sol[crd[0] + bdim[0]*(crd[1] + bdim[1]*crd[2])];
}

// ==========================================================================

int Raster::GetBasisIdx (int solidx) const
{
    dASSERT (solidx >= 0 && solidx < slen, "Argument 1 index out of range");
    return sol2basis[solidx];
}

// ==========================================================================

void Raster::GetBasisIndices (int basisidx, IVector &crd) const
{
    if (crd.Dim() != dim) crd.New(dim);
    if (dim == 3) {
	crd[2] = basisidx / (bdim[0]*bdim[1]);
	basisidx -= crd[2]*bdim[0]*bdim[1];
    }
    crd[0] = basisidx % bdim[0];
    crd[1] = basisidx / bdim[0];
}

// ==========================================================================
// Warning: this version is not threadsafe
#ifdef UNDEF
IVector &Raster::GetBasisIndices (int basisidx) const
{
    static IVector crd;
    if (crd.Dim() != dim) crd.New(dim);
    if (dim == 3) {
	crd[2] = basisidx / (bdim[0]*bdim[1]);
	basisidx -= crd[2]*bdim[0]*bdim[1];
    }
    crd[0] = basisidx % bdim[0];
    crd[1] = basisidx / bdim[0];
    return crd;
}
#endif
// ==========================================================================

IVector &Raster::GetGridIndices (int grididx) const
{
    static IVector crd;
    if (crd.Dim() != dim) crd.New(dim);
    if (dim == 3) {
	crd[2] = grididx / (gdim[0]*gdim[1]);
	grididx -= crd[2]*gdim[0]*gdim[1];
    }
    crd[0] = grididx % gdim[0];
    crd[1] = grididx / gdim[0];
    return crd;
}

// ==========================================================================

void Raster::NeighbourGraph (idxtype *&rowptr, idxtype *&colidx, int &nzero) const
{
    // Produces a graph in the solution basis (i.e. supported voxels only)

    //cout << "Setting up neighbour graph" << endl;

#if STENCIL3D == 7
    const int stencil3D_len = 7;
    int stencil3D[7][3] = {
        {0,0,-1},{0,-1,0},{-1,0,0},{0,0,0},{1,0,0},{0,1,0},{0,0,1}
    };
#elif STENCIL3D == 19
    const int stencil3D_len = 19;
    int stencil3D[19][3] = {
        {0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{0,1,-1},
	{-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0},
	{-1,1,0},{0,1,0},{1,1,0},
	{0,-1,1},{-1,0,1},{0,0,1},{1,0,1},{0,1,1}
    };
#elif STENCIL3D == 27
    const int stencil3D_len = 27;
    int stencil3D[27][3] = {
        {-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},
	{-1,1,-1},{0,1,-1},{1,1,-1},
        {-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0},
	{-1,1,0},{0,1,0},{1,1,0},
        {-1,-1,1},{0,-1,1},{1,-1,1},{-1,0,1},{0,0,1},{1,0,1},
	{-1,1,1},{0,1,1},{1,1,1}
    };
#elif STENCIL3D == 33
    const int stencil3D_len = 33;
    int stencil3D[33][3] = {
        {0,0,-2},
        {-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},
	{-1,1,-1},{0,1,-1},{1,1,-1},{0,-2,0},
        {-1,-1,0},{0,-1,0},{1,-1,0},{-2,0,0},{-1,0,0},{0,0,0},{1,0,0},
	{2,0,0},{-1,1,0},{0,1,0},{1,1,0},{0,2,0},
        {-1,-1,1},{0,-1,1},{1,-1,1},{-1,0,1},{0,0,1},{1,0,1},
	{-1,1,1},{0,1,1},{1,1,1},{0,0,2}
    };
#else
    cerr << "Error: STENCIL3D not defined" << endl;
    exit (1);
#endif

#if STENCIL2D == 5
    const int stencil2D_len = 5;
    int stencil2D[5][2] = {
        {0,-1},{-1,0},{0,0},{1,0},{0,1}
    };
#elif STENCIL2D == 9
    const int stencil2D_len = 9;
    int stencil2D[9][2] = {
        {-1,-1},{0,-1},{1,-1},{-1,0},{0,0},{1,0},{-1,1},{0,1},{1,1}
    };
#elif STENCIL2D == 13
    const int stencil2D_len = 13;
    int stencil2D[13][2] = {
        {0,-2},{-1,-1},{0,-1},{1,-1},{-2,0},{-1,0},{0,0},{1,0},{2,0},
	{-1,1},{0,1},{1,1},{0,2}
    };
#else
    cerr << "Error: STENCIL2D not defined" << endl;
    exit (1);
#endif

    idxtype *browptr = new idxtype[blen+1];
    idxtype *bcolidx;
    browptr[0] = 0;
    int bnzero = 0;

    int nx = bdim[0], ny = bdim[1], nz = (dim>2 ? bdim[2]:1);
    int i, j, k, r, c, s, nr;
    int xofs, yofs, zofs;

    // neighbour position offsets
    int dy = bdim[0], dz = bdim[0]*bdim[1];

    // step 1: build row pointer list
    if (dim == 2) {

        for (j = r = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        nr = 0;
		for (s = 0; s < stencil2D_len; s++) {
		    yofs = j + stencil2D[s][1];
		    xofs = i + stencil2D[s][0];
		    if (yofs >= 0 && yofs < ny &&
			xofs >= 0 && xofs < nx) nr++;
		}
		bnzero += nr;
		browptr[++r] = bnzero;
	    }
	}

	bcolidx = new idxtype[bnzero];
	for (j = c = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        for (s = 0; s < stencil2D_len; s++) {
		    yofs = j + stencil2D[s][1];
		    xofs = i + stencil2D[s][0];
		    if (yofs >= 0 && yofs < ny &&
			xofs >= 0 && xofs < nx)
		      bcolidx[c++] = xofs + yofs*dy;
		}
	    }
	}

    } else { // dim == 3

        for (k = r = 0; k < nz; k++) {
	    for (j = 0; j < ny; j++) {
	        for (i = 0; i < nx; i++) {
		    nr = 0;
		    for (s = 0; s < stencil3D_len; s++) {
		        zofs = k + stencil3D[s][2];
			yofs = j + stencil3D[s][1];
			xofs = i + stencil3D[s][0];
			if (zofs >= 0 && zofs < nz &&
			    yofs >= 0 && yofs < ny &&
			    xofs >= 0 && xofs < nx) nr++;
		    }
		    bnzero += nr;
		    browptr[++r] = bnzero;
		}
	    }
	}

	bcolidx = new idxtype[bnzero];
	for (k = c = 0; k < nz; k++) {
	    for (j = 0; j < ny; j++) {
	        for (i = 0; i < nx; i++) {
		    for (s = 0; s < stencil3D_len; s++) {
		        zofs = k + stencil3D[s][2];
			yofs = j + stencil3D[s][1];
			xofs = i + stencil3D[s][0];
			if (zofs >= 0 && zofs < nz &&
			    yofs >= 0 && yofs < ny &&
			    xofs >= 0 && xofs < nx) 
			  bcolidx[c++] = xofs + yofs*dy + zofs*dz;
		    }
		}
	    }
	}

    }

    // now we need to remove rows and columns with unsupported
    // raster points, i.e. transform basis -> solution
    nzero = 0;
    rowptr = new idxtype[slen+1];
    rowptr[0] = 0;
    for (i = r = 0; i < blen; i++) {
        if (basis2sol[i] >= 0) {
  	    nr = 0;
	    for (j = browptr[i]; j < browptr[i+1]; j++) {
	        c = bcolidx[j];
		if (basis2sol[c] >= 0) nr++;
	    }
	    nzero += nr;
	    rowptr[++r] = nzero;
	}
    }
    colidx = new idxtype[nzero];
    for (i = k = r = 0; i < blen; i++) {
        if (basis2sol[i] >= 0) {
	    for (j = browptr[i]; j < browptr[i+1]; j++) {
	        c = bcolidx[j];
		if (basis2sol[c] >= 0) colidx[k++] = basis2sol[c];
	    }
	}
    }
    delete []browptr;
    delete []bcolidx;

    //cout << nzero << " entries in neighbour graph found" << endl;
}

// ==========================================================================

void Raster::NeighbourGraph (ICompRowMatrix &NG) const
{
    // Produces a graph in the solution basis (i.e. supported voxels only)
    // Returns graph structure in a sparse integer matrix

    // This version additionally builds a connectivity list and returns it
    // in array 'connect'. This can be used, together
    // with rowptr and colidx, to create a sparse integer matrix with entries
    // that encode the neighbour direction
    // byte 0: -1/0/+1 for x-offset
    // byte 1: -1/0/+1 for y-offset
    // byte 2: -1/0/+1 for z-offset
    // byte 3 (and higher): not used

    //cout << "Setting up neighbour graph" << endl;

#if STENCIL3D == 7
    const int stencil3D_len = 7;
    int stencil3D[7][3] = {
        {0,0,-1},{0,-1,0},{-1,0,0},{0,0,0},{1,0,0},{0,1,0},{0,0,1}
    };
#elif STENCIL3D == 19
    const int stencil3D_len = 19;
    int stencil3D[19][3] = {
        {0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{0,1,-1},
	{-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0},
	{-1,1,0},{0,1,0},{1,1,0},
	{0,-1,1},{-1,0,1},{0,0,1},{1,0,1},{0,1,1}
    };
#elif STENCIL3D == 27
    const int stencil3D_len = 27;
    int stencil3D[27][3] = {
        {-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},
	{-1,1,-1},{0,1,-1},{1,1,-1},
        {-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0},
	{-1,1,0},{0,1,0},{1,1,0},
        {-1,-1,1},{0,-1,1},{1,-1,1},{-1,0,1},{0,0,1},{1,0,1},
	{-1,1,1},{0,1,1},{1,1,1}
    };
#elif STENCIL3D == 33
    const int stencil3D_len = 33;
    int stencil3D[33][3] = {
        {0,0,-2},
        {-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},
	{-1,1,-1},{0,1,-1},{1,1,-1},{0,-2,0},
        {-1,-1,0},{0,-1,0},{1,-1,0},{-2,0,0},{-1,0,0},{0,0,0},{1,0,0},
	{2,0,0},{-1,1,0},{0,1,0},{1,1,0},{0,2,0},
        {-1,-1,1},{0,-1,1},{1,-1,1},{-1,0,1},{0,0,1},{1,0,1},
	{-1,1,1},{0,1,1},{1,1,1},{0,0,2}
    };
#else
    cerr << "Error: STENCIL3D not defined" << endl;
    exit (1);
#endif

#if STENCIL2D == 5
    const int stencil2D_len = 5;
    int stencil2D[5][2] = {
        {0,-1},{-1,0},{0,0},{1,0},{0,1}
    };
#elif STENCIL2D == 9
    const int stencil2D_len = 9;
    int stencil2D[9][2] = {
        {-1,-1},{0,-1},{1,-1},{-1,0},{0,0},{1,0},{-1,1},{0,1},{1,1}
    };
#elif STENCIL2D == 13
    const int stencil2D_len = 13;
    int stencil2D[13][2] = {
        {0,-2},{-1,-1},{0,-1},{1,-1},{-2,0},{-1,0},{0,0},{1,0},{2,0},
	{-1,1},{0,1},{1,1},{0,2}
    };
#else
    cerr << "Error: STENCIL2D not defined" << endl;
    exit (1);
#endif

    idxtype *browptr = new idxtype[glen+1];
    idxtype *bcolidx;
    idxtype *rowptr, *colidx;
	int nzero, *bconnect, *connect;
    browptr[0] = 0;
    int bnzero = 0;

    int nx = gdim[0], ny = gdim[1], nz = (dim>2 ? gdim[2]:1);
    int i, j, k, r, c, s, nr;
    int xofs, yofs, zofs;

    // neighbour position offsets
    int dx = 1, dy = gdim[0], dz = gdim[0]*gdim[1];

    // step 1: build row pointer list
    if (dim == 2) {

        for (j = r = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        nr = 0;
		for (s = 0; s < stencil2D_len; s++) {
		    yofs = j + stencil2D[s][1];
		    xofs = i + stencil2D[s][0];
		    if (yofs >= 0 && yofs < ny &&
			xofs >= 0 && xofs < nx) nr++;
		}
		bnzero += nr;
		browptr[++r] = bnzero;
	    }
	}

	bcolidx = new idxtype[bnzero];
	bconnect = new int[bnzero];
	for (j = c = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        for (s = 0; s < stencil2D_len; s++) {
		    yofs = j + stencil2D[s][1];
		    xofs = i + stencil2D[s][0];
		    if (yofs >= 0 && yofs < ny &&
		      xofs >= 0 && xofs < nx) {
			bcolidx[c] = xofs + yofs*dy;
			bconnect[c] = 
			    (stencil2D[s][0] & 0xFF) |
			    ((stencil2D[s][1] & 0xFF) << 8);
			c++;
		    }
		}
	    }
	}

    } else { // dim == 3

        for (k = r = 0; k < nz; k++) {
	    for (j = 0; j < ny; j++) {
	        for (i = 0; i < nx; i++) {
		    nr = 0;
		    for (s = 0; s < stencil3D_len; s++) {
		        zofs = k + stencil3D[s][2];
			yofs = j + stencil3D[s][1];
			xofs = i + stencil3D[s][0];
			if (zofs >= 0 && zofs < nz &&
			    yofs >= 0 && yofs < ny &&
			    xofs >= 0 && xofs < nx) nr++;
		    }
		    bnzero += nr;
		    browptr[++r] = bnzero;
		}
	    }
	}

	bcolidx = new idxtype[bnzero];
	bconnect = new int[bnzero];
	for (k = c = 0; k < nz; k++) {
	    for (j = 0; j < ny; j++) {
	        for (i = 0; i < nx; i++) {
		    for (s = 0; s < stencil3D_len; s++) {
		        zofs = k + stencil3D[s][2];
			yofs = j + stencil3D[s][1];
			xofs = i + stencil3D[s][0];
			if (zofs >= 0 && zofs < nz &&
			    yofs >= 0 && yofs < ny &&
			    xofs >= 0 && xofs < nx) {
			    bcolidx[c] = xofs + yofs*dy + zofs*dz;
			    bconnect[c] =
				(stencil3D[s][0] & 0xFF) |
				((stencil3D[s][1] & 0xFF) << 8) |
				((stencil3D[s][2] & 0xFF) << 16);
			    c++;
			}
		    }
		}
	    }
	}

    }

    // now we need to remove rows and columns with unsupported
    // raster points, i.e. transform basis -> solution
    nzero = 0;
    rowptr = new idxtype[slen+1];
    rowptr[0] = 0;
    for (i = r = 0; i < glen; i++) {
        if (basis2sol[i] >= 0) {
  	    nr = 0;
	    for (j = browptr[i]; j < browptr[i+1]; j++) {
	        c = bcolidx[j];
		if (basis2sol[c] >= 0) nr++;
	    }
	    nzero += nr;
	    rowptr[++r] = nzero;
	}
    }
    colidx = new idxtype[nzero];
    connect = new int[nzero];
    for (i = k = r = 0; i < glen; i++) {
        if (basis2sol[i] >= 0) {
	    for (j = browptr[i]; j < browptr[i+1]; j++) {
	        c = bcolidx[j];
		if (basis2sol[c] >= 0) {
		    colidx[k] = basis2sol[c];
		    connect[k] = bconnect[j];
		    k++;
		}
	    }
	}
    }

    NG.New(slen, slen);
    NG.Initialise (rowptr, colidx, connect);

    delete []browptr;
    delete []bcolidx;
    delete []bconnect;
    delete []rowptr;
    delete []colidx;
    delete []connect;
    //cout << nzero << " entries in neighbour graph found" << endl;
}

// ==========================================================================

IVector Raster::NeighbourShift (const ICompRowMatrix &NG, int i, int j) const
{
    int d;
    if (!NG.Exists (i,j))
	return IVector();  // invalid index -> return empty vector

    int val = NG.Get(i,j);
    //    cout << val;
    IVector shift(dim);
    for (d = 0; d < dim; d++) {
	char c = (char)((val >> (d*8)) & 0xFF);
	shift[d] = (int)c;
    }
    return shift;
}

// ==========================================================================

RVector Raster::RNeighbourShift (const ICompRowMatrix &NG, int i, int j) const
{
    IVector ishift;
    ishift = NeighbourShift(NG, i, j);
    //    cout << ishift;
    RVector shift(ishift.Dim());
    for(int d =0; d < ishift.Dim(); d ++)
      shift[d] = ishift[d];
    //    cout << shift;

    return shift;
}

// ==========================================================================

RVector Raster::ImageGradient(const RVector &x, double sd) const
{
    RVector gn;
    int ngim = (dim ==3 ? 10 : 6);
    //	bool *iflag = new bool [dimflag ==3 ? 10 : 6];
    bool *iflag =new bool[ngim]; 
    for (int ng = 0; ng < ngim; ng++)
	iflag[ng] = false;
    iflag[1] = iflag[2] = true;
    if(dim == 3)
	iflag[3] = true;
    RVector *gx = ImageJet(x,sd,iflag);
    gn = gx[1]*gx[1] + gx[2]*gx[2];
    if(dim == 3)
	gn = gn + gx[3]*gx[3];
    gn =  sqrt(gn);
    delete []gx;
    delete []iflag;
    return gn;
}

// ==========================================================================


#define RESCALE_W0
RVector Raster::SmoothImage(const RVector &x, double sd) const
{
    // smooth image by sd in each dimension
    double max = -1;
    RVector sx(x);

    for(int k = 0; k < sx.Dim(); k++)
        max = (sx[k] > max ? sx[k] : max);
    LOGOUT("Max solbasis (on input) = %f",max);
    RVector gsx(glen);            // image in grid coordinates
    //    cout << "smoothing image \n";
    Map_SolToGrid(sx,gsx);
    //    cout << "Map_SolToGrid done\n";
    for(int k = 0; k < gsx.Dim(); k++)
          max = (gsx[k] > max ? gsx[k] : max);
    //    cout << "Max gridbasis (before FFT)= " << max << endl;
    float *data = new float[2*padlen]; 
    for(int id = 0; id < 2*padlen; id++)
      data[id] = 0;
    if(padlen == glen) {
      for(int k = 0; k < glen; k++) {
	  data[2*k] = (float)gsx[k];  // load up real part of data
	  //	  data[2*k+1] = 0;    // load up imag part of data
      }
    }
    else {
      int xstart  =(paddim[0]-gdim[0])/2;
      int ystart  =(paddim[1]-gdim[1])/2;
      int zstart; // cannot scope it inside switch
    switch(dim) {
    case 2 :
      cout << "x,y start : " << xstart << " " << ystart << endl;
      for(int iy = ystart, k=0; iy < ystart + gdim[1]; iy ++) {
	int offy = iy*paddim[0];
	for(int ix = xstart; ix < xstart + gdim[0]; ix ++) {
	  data[2*(offy + ix)]  = (float)gsx[k++]; // real
	}
	//	cout << k << " "; cout.flush();
      }
      //      cout << endl;
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      zstart  =(paddim[2]-gdim[2])/2;

      for(int iz = zstart, k = 0; iz < zstart + gdim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       for(int iy = ystart; iy < ystart + gdim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	for(int ix = xstart; ix <  xstart + gdim[0]; ix ++) {
	  data[2*(offtot + ix)]  = (float)gsx[k++]; // real 
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }

    }
    int *nn = new int[dim];
    float *sdd = new float[dim];
    float *w0 = new float[dim];
    RVector *filter  = new RVector[dim];

    //    cout << "calling FFT dimensions ";  
    int ncount = 1;
    for(int idim = 0; idim < dim; idim++) {

	  nn[idim] = paddim[dim-idim-1];  // image dimensions
#ifdef RESCALE_W0
	  w0[idim] = (float)(2*Pi/paddim[0]);
	  sdd[idim] = (float)sd; // CAREFUL! RESCALING SD!
#else
	  w0[idim] = 2*M_PI/paddim[idim];
	  sdd[idim] = sd; 
#endif
	  /* OLD VERSION
	  nn[idim] = paddim[dim-idim-1];  // image dimensions
	  w0[idim] = 2*M_PI/paddim[idim];
	  sdd[idim] = sd; // *gdim[idim]/gdim[0]; // CAREFUL! RESCALIN SD!
	  */

	  filter[idim].New(paddim[idim]);
	  //       	  cout << "nn " << nn[idim] << " w0 "<< w0[idim] << " " << " sd "<< sdd[idim] << " ";
	  ncount *= paddim[idim];
    }
    fourn(data-1,nn-1,dim,1); // forward FFT
    //    cout << "\nFFT done\n";


    /* now construct filters in Fourier domain */
    for(int idim = 0; idim < dim; idim++) {
      for(int j = 0; j <= paddim[idim]/2 ; j++)
	filter[idim][j] = exp(-SQR(j*w0[idim] *sdd[idim])/2);
      for(int j =  1 + paddim[idim]/2 ; j < paddim[idim] ; j++)
	filter[idim][j] = filter[idim][paddim[idim]-j];
    }
    // filter it
    switch(dim) {
    case 2 :
      for(int iy = 0; iy < paddim[1]; iy ++) {
	int offy = iy*paddim[0];
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  data[2*(offy + ix)]  *= (float)(filter[0][ix]*filter[1][iy]); // real
	  data[2*(offy + ix)+1] *= (float)(filter[0][ix]*filter[1][iy]); // imag
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      for(int iz = 0; iz < paddim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       for(int iy = 0; iy < paddim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  data[2*(offtot + ix)]   *= (float)(filter[0][ix]*filter[1][iy]*filter[2][iz]); // real
	  data[2*(offtot + ix)+1] *= (float)(filter[0][ix]*filter[1][iy]*filter[2][iz]); // imag
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }
    fourn(data-1,nn-1,dim,-1); // inverse FFT
      
    if(padlen == glen) {
      for(int k = 0; k < glen; k++)
//        gsx[k] = sqrt(SQR(data[2*k]) + SQR(data[2*k+1]) );     // abs FFT
        gsx[k] = data[2*k]/ncount;     // real FFT
    }
    else {
      int xstart  =(paddim[0]-gdim[0])/2;
      int ystart  =(paddim[1]-gdim[1])/2;
      int zstart;
    switch(dim) {
    case 2 :

      for(int iy = ystart, k=0; iy < ystart + gdim[1]; iy ++) {
	int offy = iy*paddim[0];
	for(int ix = xstart; ix < xstart + gdim[0]; ix ++) {
	  gsx[k++] =data[2*(offy + ix)]/ncount ; // real
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      zstart  =(paddim[2]-gdim[2])/2;

      for(int iz = zstart, k = 0; iz < zstart + gdim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       for(int iy = ystart; iy < ystart + gdim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	for(int ix = xstart; ix <  xstart + gdim[0]; ix ++) {
	  gsx[k++] =data[2*(offtot + ix)]/ncount  ; // real 
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }
    }
    for(int k = 0; k < gsx.Dim(); k++)
        max = (gsx[k] > max ? gsx[k] : max);
    //    cout << "Max gridbasis (after FFT)= " << max << endl;
    //    cout << "\nCalling Map_GridToSol\n";
    Map_GridToSol(gsx,sx);
    //	cout << "Map_GridToSol done\n";
    max = -1;
    for(int k = 0 ; k < sx.Dim(); k++)
        max = (sx[k] > max ? sx[k] : max);
    LOGOUT("Max solbasis (after smoothing) = %f", max);
    delete [] data;
    delete [] nn;
    delete [] w0;
    delete [] filter;
    return sx;
}

RVector *Raster::ImageJet(const RVector &x, double sd, bool *iflags) const
{
  /*********************************************************************/
  /* returned values are :
     2D : 0              - smoothed image 
          1,2            - smoothed x,y derivatives
          3,4            - smoothed xx,yy derivatives
	  5              - smoothed xy derivative
     3D : 0              - smoothed image 
          1,2,3          - smoothed x,y,z derivatives
          4,5,6          - smoothed xx,yy,zz derivatives
	  7,8,9          - smoothed xy,xz,yz derivatives

     iflags is a set of Booleans. The corresponding image is only
     constructed if the flag is set.
  */
  static int D2D[2][6] = { {0, 1, 0, 2, 0, 1 }, 
		           {0, 0, 1, 0, 2, 1}};
  static int D3D[3][10] = { {0, 1, 0, 0, 2, 0, 0, 1, 1, 0 },
		            {0, 0, 1, 0, 0, 2, 0, 1, 0, 1},
		            {0, 0, 0, 1, 0, 0, 2, 0, 1, 1} };
    int ngim = (dim == 3 ? 10 : 6); // the number of images to be output
    // smooth image by sd in each dimension
    double max = -1;
    RVector *gx = new RVector [ngim];

    for(int k = 0; k < x.Dim(); k++)
        max = (x[k] > max ? x[k] : max);
    LOGOUT("Max solbasis (on input) = %f",max);
    RVector gsx(glen);            // image in grid coordinates
    //    cout << "smoothing image \n";
    Map_SolToGrid(x,gsx);
    //    cout << "Map_SolToGrid done\n";
    for(int k = 0; k < gsx.Dim(); k++)
          max = (gsx[k] > max ? gsx[k] : max);
    //    cout << "Max gridbasis (before FFT)= " << max << endl;

    float *data = new float[2*padlen]; 
    for(int id = 0; id < 2*padlen; id++)
      data[id] = 0;
    if(padlen == glen) {
      for(int k = 0; k < glen; k++) {
	  data[2*k] = (float)gsx[k]; // load up real part of data
	  data[2*k+1] = 0;    // load up imag part of data
      }
    }
    else {
      int xstart  =(paddim[0]-gdim[0])/2;
      int ystart  =(paddim[1]-gdim[1])/2;
      int zstart; // cannot scope it inside switch
    switch(dim) {
    case 2 :

      cout << "x,y start : " << xstart << " " << ystart << endl;
      for(int iy = ystart, k=0; iy < ystart + gdim[1]; iy ++) {
	int offy = iy*paddim[0];
	for(int ix = xstart; ix < xstart + gdim[0]; ix ++) {
	  data[2*(offy + ix)]  = (float)gsx[k++]; // real
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      zstart  =(paddim[2]-gdim[2])/2;

      for(int iz = zstart, k = 0; iz < zstart + gdim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       for(int iy = ystart; iy < ystart + gdim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	for(int ix = xstart; ix <  xstart + gdim[0]; ix ++) {
	  data[2*(offtot + ix)]  = (float)gsx[k++]; // real 
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }

    }

    int *nn = new int[dim];
    float *sdd = new float[dim];
    float *w0 = new float[dim];
    RVector *filter  = new RVector[dim];
    RVector *d1  = new RVector[dim]; // 1st derivative in 1D

    //    cout << "calling FFT dimensions ";  
    int ncount = 1;
    for(int idim = 0; idim < dim; idim++) {
	  nn[idim] = paddim[dim-idim-1];  // image dimensions
#ifdef RESCALE_W0
	  w0[idim] = (float)(2*Pi/paddim[0]);
	  sdd[idim] = (float)sd; // CAREFUL! RESCALING SD!
#else
	  w0[idim] = 2*M_PI/paddim[idim];
	  sdd[idim] = sd; 
#endif
	  filter[idim].New(paddim[idim]);
	  d1[idim].New(paddim[idim]);
	  //       	  cout << "nn " << nn[idim] << " w0 "<< w0[idim] << " " << " sd "<< sdd[idim] << " ";
	  ncount *= paddim[idim];
    }
    fourn(data-1,nn-1,dim,1); // forward FFT

    float *data0 = new float[2*padlen]; // need a copy
    for (int n =0; n < 2*padlen; n++)
      data0[n] = data[n];

    //    cout << "\nFFT done\n";


    /* now construct filters in Fourier domain */
    for(int idim = 0; idim < dim; idim++) {
      for(int j = 0; j <= paddim[idim]/2 ; j++) {
	filter[idim][j] = (sdd[idim]==0 ? 1 : /* sdd[idim]* */
		    exp(-SQR(j*w0[idim] *sdd[idim])/2));// /sqrt(2*M_PI);
	d1[idim][j] = j*w0[idim]; // recall it is imaginary...
      }
      for(int j =  1+paddim[idim]/2 ; j < paddim[idim] ; j++) {
	filter[idim][j] = filter[idim][paddim[idim]-j];
	d1[idim][j] = (j-paddim[idim])*w0[idim];
      }
      //      cout << filter[idim] << endl;
    }
 
    // filter it
    std::complex<double> *dz = new std::complex<double>[dim];

    for(int ng = 0; ng < ngim; ng++) {
      if(!iflags[ng]) // skip if not set
	continue;
      //      cout << "Image Jet " << ng << " size " << x.Dim() << " scale " << sd << endl;
      gx[ng].New(x.Dim()); // output images are in solbasis

    switch(dim) {
    case 2 :
      for(int iy = 0; iy < paddim[1]; iy ++) {
	int offy = iy*paddim[0];
	dz[1] = std::complex<double>(0.0,d1[1][iy]);
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  dz[0] = std::complex<double>(0.0,d1[0][ix]);

	  std::complex<double> z(data0[2*(offy + ix)],data0[2*(offy + ix)+1]);
	  //	  cout << z << " " << (filter[0][ix]*filter[1][iy]) << " ";
	  z *= std::complex<double>(filter[0][ix]*filter[1][iy],0.0); // no toast::complex*real??
	  //       	  cout << z << " ";
	  for(int idim = 0; idim < dim; idim++) {
	    for(int p = 0; p < D2D[idim][ng] ; p++)
	      z *= dz[idim];
	  }
	  //	  cout << dz[0] << " " << dz[1] << " ";
	  //       	  cout << z << " "; 
	  data[2*(offy + ix)]  = (float)z.real();
	  data[2*(offy + ix)+1] = (float)z.imag();
	  //	  cout << data[2*(offy + ix)] << " " << data[2*(offy + ix)+1] << endl;
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      for(int iz = 0; iz < paddim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       dz[2] = std::complex<double>(0.0,d1[2][iz]);
       for(int iy = 0; iy < paddim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	dz[1] = std::complex<double>(0.0,d1[1][iy]);
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  dz[0] = std::complex<double>(0.0,d1[0][ix]);

	  std::complex<double> z(data0[2*(offtot + ix)],data0[2*(offtot + ix)+1]);
	  //	  cout << z << " " << (filter[0][ix]*filter[1][iy]) << " ";
	  z *= std::complex<double>(filter[0][ix]*filter[1][iy]*filter[2][iz],0.0); // no toast::complex*real??
	  //       	  cout << z << " ";
	  for(int idim = 0; idim < dim; idim++) {
	    for(int p = 0; p < D3D[idim][ng] ; p++)
	      z *= dz[idim];
	  }
	
	  data[2*(offtot + ix)]  = (float)z.real();
	  data[2*(offtot + ix)+1] = (float)z.imag();
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }
    fourn(data-1,nn-1,dim,-1); // inverse FFT
      

    if(padlen == glen) {
      for(int k = 0; k < glen; k++)
//        gsx[k] = sqrt(SQR(data[2*k]) + SQR(data[2*k+1]) );     // abs FFT
        gsx[k] = data[2*k]/ncount;     // real FFT
    }
    else {
      int xstart  =(paddim[0]-gdim[0])/2;
      int ystart  =(paddim[1]-gdim[1])/2;
      int zstart;
    switch(dim) {
    case 2 :

      for(int iy = ystart, k=0; iy < ystart + gdim[1]; iy ++) {
	int offy = iy*paddim[0];
	for(int ix = xstart; ix < xstart + gdim[0]; ix ++) {
	  gsx[k++] =data[2*(offy + ix)]/ncount ; // real
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      zstart  =(paddim[2]-gdim[2])/2;

      for(int iz = zstart, k = 0; iz < zstart + gdim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       for(int iy = ystart; iy < ystart + gdim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	for(int ix = xstart; ix <  xstart + gdim[0]; ix ++) {
	  gsx[k++] =data[2*(offtot + ix)]/ncount  ; // real 
	}
       }
      }
      break;
    default :
      cerr << "Cannot do smoothing with dimension " << dim << endl;
    }
    }

    for(int k = 0; k < gsx.Dim(); k++)
        max = (gsx[k] > max ? gsx[k] : max);
    //    cout << "Max gridbasis (after FFT)= " << max << endl;
    //    cout << "\nCalling Map_GridToSol\n";
    Map_GridToSol(gsx,gx[ng]);
    //	cout << "Map_GridToSol done\n";
    max = -1;
    for(int k = 0 ; k < gx[ng].Dim(); k++)
        max = (gx[ng][k] > max ? gx[ng][k] : max);
    //    cout << "Max solbasis (after smoothing) = " << max << endl;
    } // end of images loop
    LOGOUT("Max solbasis (after smoothing) = %f", max);
    delete []data;
    delete []data0;
    delete []nn;
    delete []w0;
    delete []filter;
    delete []d1;
    delete []dz;
    return gx;
#ifdef NOT_DONE_YET

#endif
}
