#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "timing.h"

using namespace std;
using namespace toast;

// ==========================================================================
// STENCIL2D and STENCIL3D: (used by PCG_PRECON_SPARSEJTJ only): number of
// neighbours to include in the JTJ sparsity graph

#define STENCIL2D 5  // valid entries: 5, 9, 13
#define STENCIL3D 7  // valid entries: 7, 19, 27, 33

// ==========================================================================
// class Raster

Raster::Raster (IVector &_gdim, Mesh *mesh)
{
    int i, j;

    bCoarseBasis = false; // no lo-res solution basis supported
    dim = _gdim.Dim();
    meshptr = mesh;
    xASSERT(dim == meshptr->Dimension(),
	    Raster and mesh have incompatible dimensions);
    meshptr->BoundingBox (bbmin, bbmax);
    gdim.New(dim); gdim = _gdim;
    bdim.New(dim); bdim = _gdim;
    paddim.New(dim);
    gsize.New(dim);

    for (i = 0, glen = 1, padlen=1; i < dim; i++) {
        j = 0; 
	while(POW2[j] < gdim[i]) // find power of 2 greater than dimension
	  j++;  // do nothing
	paddim[i] =POW2[j];
        glen *= gdim[i];
	padlen  *= paddim[i];
	gsize[i] = bbmax[i]-bbmin[i];
    }
    blen = glen;

    //tic();
    elref = GenerateElementPixelRef (*meshptr, gdim, &bbmin, &bbmax);
    //cout << "GenerateElementPixelRef: t=" << toc() << endl;
    //tic();
    B     = GridMapMatrix (*meshptr, gdim, &bbmin, &bbmax, elref);
    ((RCompRowMatrix*)B)->Shrink();
    //cout << "GridMapMatrix:           t=" << toc() << endl;
    //tic();

#if MAP_GRIDTOMESH == MGM_PSEUDOINVERSE
    if (B->nRows() > B->nCols()) {
	cout << "Over-determined formulation" << endl;
	BB = transp(*(RCompRowMatrix*)B) * *(RCompRowMatrix*)B;
	BBover = true;
    } else {
	cout << "Under-determined formulation" << endl;
	BB = *(RCompRowMatrix*)B * transp(*(RCompRowMatrix*)B);
	BBover = false;
    }
    RVector diag = BB.Diag();
    double lambda = l2norm(diag)*1e-2;
    int n = BB.nRows();
    int *drowptr = new int[n+1];
    int *dcolidx = new int[n];
    double *dval = new double[n];
    for (i = 0; i <= n; i++) drowptr[i] = i;
    for (i = 0; i <  n; i++) dcolidx[i] = i;
    for (i = 0; i <  n; i++) dval[i] = lambda;
    RCompRowMatrix dg (n, n, drowptr, dcolidx, dval);
    delete []drowptr;
    delete []dcolidx;
    delete []dval;
    BB += dg;
    pBTB = new RPrecon_Diag;
    pBTB->Reset (&BB);
#else
    pBTB = 0;
#endif

#if MAP_GRIDTOMESH == MGM_TRANSPOSE
    BI = GridMapMatrix (*meshptr, gdim, &bbmin, &bbmax, elref);
    BI->Transpone();
#else
    BI = NodeMapMatrix (*meshptr, gdim, &bbmin, &bbmax, elref);
    //BI = NodeMapMatrix2 (*meshptr, gdim, &bbmin, &bbmax, elref);
    ((RCompRowMatrix*)BI)->Shrink();
#endif
    //cout << "NodeMapMatrix:           t=" << toc() << endl;

    //tic();
    // calculate mesh support weighting factors for all pixels in
    // user basis.
    RVector s(glen);
    bsupport.New (glen);
    // s is the support mask for the hires pixel grid
    for (i = 0; i < glen; i++)
        if (elref[i] >= 0) s[i] = 1.0;

    // map mask into user basis
    //SubsampleLinPixel (s, bsupport, gdim, bdim, 0);
    bsupport = s; // for now only binary mask is supported

    // calculate basis->solution mapping index list
    basis2sol.New (glen);    
    for (i = slen = 0; i < glen; i++) {
        if (bsupport[i] > 0.0) basis2sol[i] = slen++;
	else                   basis2sol[i] = -1;
    }
    sol2basis.New (slen);
    for (i = slen = 0; i < glen; i++) {
	if (bsupport[i] > 0.0) sol2basis[slen++] = i;
    }

    // formulate basis->solution mapping in sparse matrix
    int *rowptr = new int[slen+1];
    int *colidx = new int[slen];
    double *val = new double[slen];
    for (i = 0; i <= slen; i++) rowptr[i] = i; // each row has one entry
    for (i = 0; i < slen; i++) val[i] = 1.0;
    for (i = j = 0; i < glen; i++)
        if (bsupport[i] > 0.0) colidx[j++] = i;
    D = new RCompRowMatrix (slen, glen, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;
    //cout << "Raster constructor: t=" << toc() << endl;
}

Raster::Raster (IVector &_gdim, IVector &_bdim, Mesh *mesh)
{
    int i, j;

    bCoarseBasis = true; // this version supports a lo-res solution basis
    dim     = _bdim.Dim();
    meshptr = mesh;
    xASSERT(dim == _gdim.Dim(),
	    Raster and user basis have incompatible dimensions);
    xASSERT(dim == meshptr->Dimension(),
	    Raster and mesh have incompatible dimensions);

    meshptr->BoundingBox (bbmin, bbmax);
    gdim.New(dim); gdim = _gdim;
    bdim.New(dim); bdim = _bdim;
    paddim.New(dim);
    gsize.New(dim);

    for (i = 0, glen = blen = 1, padlen=1; i < dim; i++) {
        j = 0; 
	while(POW2[j] < gdim[i]) // find power of 2 greater than dimension
	  j++;  // do nothing
	paddim[i] =POW2[j];
        glen *= gdim[i];
	padlen  *= paddim[i];
	blen *= bdim[i];
	gsize[i] = bbmax[i]-bbmin[i];
    }

    tic();
    elref = GenerateElementPixelRef (*meshptr, gdim, &bbmin, &bbmax);
    belref = GenerateElementPixelRef (*meshptr, bdim, &bbmin, &bbmax);
    cout << "GenerateElementPixelRef: t=" << toc() << endl;
    tic();
    B     = GridMapMatrix (*meshptr, gdim, &bbmin, &bbmax, elref);
    ((RCompRowMatrix*)B)->Shrink();
    cout << "GridMapMatrix:           t=" << toc() << endl;
    tic();
#if MAP_GRIDTOMESH == MGM_PSEUDOINVERSE
    xERROR(Raster: Not implemented: mapping via pseudo-inverse);
#endif
#if MAP_GRIDTOMESH == MGM_TRANSPOSE
    xERROR(Raster: Not implemented: mapping via transpose);
#endif
    pBTB = 0;
    BI    = NodeMapMatrix (*meshptr, gdim, &bbmin, &bbmax, elref);
    ((RCompRowMatrix*)BI)->Shrink();
    cout << "NodeMapMatrix:           t=" << toc() << endl;
    tic();
    B2    = Grid2LinPixMatrix (gdim, bdim, elref);
    ((RCompRowMatrix*)B2)->Shrink();
    cout << "Grid2LinPixMatrix:       t=" << toc() << endl;
    tic();
    B2I   = LinPix2GridMatrix (gdim, bdim, elref);
    ((RCompRowMatrix*)B2I)->Shrink();

#ifdef UNDEF
    // to be incorporated
    delete B2I;
    B2I = CubicPix2GridMatrix (bdim, gdim, elref);
    delete B2;
    B2 = new RCompRowMatrix (*(RCompRowMatrix*)B2I);
    B2->Transpone();

    RCompRowMatrix tmp1 = *(RCompRowMatrix*)B2 * *(RCompRowMatrix*)B2I;
#endif

    cout << "LinPix2GridMatrix:       t=" << toc() << endl;
    tic();
    C     = new RCompRowMatrix;
    cout << "Multiplying B2(" << B2->nRows() << ','
	 << B2->nCols() << ',' << B2->nVal() << ')' << endl;
    cout << "with         B(" << B->nRows() << ','
	 << B->nCols() << ',' << B->nVal() << ')' << endl;
    ((RCompRowMatrix*)B2)->AB (*(RCompRowMatrix*)B, *C);
    ((RCompRowMatrix*)C)->Shrink();
    cout << "Result       C(" << C->nRows() << ','
	 << C->nCols() << ',' << C->nVal() << ')' << endl;
    cout << "SparseMatrixMult:        t=" << toc() << endl;
    tic();
    CI    = new RCompRowMatrix;
    cout << "Multiplying BI(" << BI->nRows() << ','
	 << BI->nCols() << ',' << BI->nVal() << ')' << endl;
    cout << "with       B2I(" << B2I->nRows() << ','
	 << B2I->nCols() << ',' << B2I->nVal() << ')' << endl;
    ((RCompRowMatrix*)BI)->AB (*(RCompRowMatrix*)B2I, *CI);
    ((RCompRowMatrix*)CI)->Shrink();
    cout << "Result      CI(" << CI->nRows() << ','
	 << CI->nCols() << ',' << CI->nVal() << ')' << endl;
    cout << "SparseMatrixMult:        t=" << toc() << endl;

    tic();
    // calculate mesh support weighting factors for all pixels in
    // user basis.
    RVector s(glen);
    bsupport.New (blen);
    // s is the support mask for the hires pixel grid
    for (i = 0; i < glen; i++)
        if (elref[i] >= 0) s[i] = 1.0;

    // map mask into user basis
    SubsampleLinPixel (s, bsupport, gdim, bdim, 0);

    // calculate basis->solution mapping index list
    basis2sol.New (blen);    
    for (i = slen = 0; i < blen; i++) {
        if (bsupport[i] > 0.0) basis2sol[i] = slen++;
	else                   basis2sol[i] = -1;
    }
    sol2basis.New (slen);
    for (i = slen = 0; i < blen; i++) {
	if (bsupport[i] > 0.0) sol2basis[slen++] = i;
    }

    // formulate basis->solution mapping in sparse matrix
    int *rowptr = new int[slen+1];
    int *colidx = new int[slen];
    double *val = new double[slen];
    for (i = 0; i <= slen; i++) rowptr[i] = i; // each row has one entry
    for (i = 0; i < slen; i++) val[i] = 1.0;
    for (i = j = 0; i < blen; i++)
        if (bsupport[i] > 0.0) colidx[j++] = i;
    D = new RCompRowMatrix (slen, blen, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;
    cout << "Raster constructor: t=" << toc() << endl;
}

Raster::~Raster ()
{
    delete []elref;
    delete B;
    delete BI;
    delete D;
    if (pBTB) delete pBTB;

    if (bCoarseBasis) {
	delete []belref;
	delete B2;
	delete B2I;
	delete C;
	delete CI;
    }
}

void Raster::GetPixelCoords (int i, IVector &pix) const
{
    int idx = sol2basis[i];
    if (dim > 2) {
	pix[2] = idx / (bdim[0]*bdim[1]);
	idx -= pix[2] * bdim[0]*bdim[1];
    }
    pix[1] = idx / bdim[0];
    pix[0] = idx % bdim[0];
}

const RGenericSparseMatrix &Raster::Mesh2BasisMatrix() const
{
    if (bCoarseBasis) return *C;
    else return Mesh2GridMatrix();
}

const RGenericSparseMatrix &Raster::Basis2MeshMatrix() const
{
    if (bCoarseBasis) return *CI;
    else return Grid2MeshMatrix();
}

void Raster::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
  
    RVector bvec(blen);
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

void Raster::Map_GridToSol (const RVector &gvec, RVector &svec) const
{
    if (bCoarseBasis) {
	RVector bvec(blen);
	Map_GridToBasis (gvec, bvec);
	Map_BasisToSol (bvec, svec);
    } else {
	Map_BasisToSol (gvec, svec);
    }
}

void Raster::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    RVector bvec(blen);
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}

void Raster::Map_GridToMeshPI (const RVector &gvec, RVector &mvec) const
{
    // map grid->mesh by pseudo-inverse of the mesh->grid matrix
    double tol = 1e-10;
    int maxit = 100000;
    if (BBover) {
	RVector b = B->ATx (gvec);
	int n = CG (BB, b, mvec, tol, pBTB, maxit);
	if (n < maxit) { // converged
	    cerr << "G->M: niter=" << n << ", res=" << tol << endl;
	} else {
	    cerr << "No convergence, use direct mapping\n";
	    BI->Ax (gvec, mvec);
	}
    } else {
	RVector m(gvec.Dim());
	int n = CG (BB, gvec, m, tol, pBTB, maxit);
	if (n < maxit) { // converged
	    mvec = B->ATx (m);
	    cerr << "G->M: niter=" << n << ", res=" << tol << endl;
	} else {
	    cerr << "No convergence\n";
	    //cerr << "No convergence, use direct mapping\n";
	    //BI->Ax (gvec, mvec);
	}
    }
    RVector tmp;
    B->Ax(mvec,tmp);
    cerr << "residual: " << l2norm(tmp-gvec)/l2norm(gvec) << endl;
}

void Raster::Map_MeshToSol (const CVector &mvec, CVector &svec) const
{
    RVector tmp(svec.Dim());
    Map_MeshToSol (Re(mvec), tmp); SetReal (svec, tmp);
    Map_MeshToSol (Im(mvec), tmp); SetImag (svec, tmp);
}

void Raster::Map_SolToMesh (const CVector &svec, CVector &mvec) const
{
    RVector tmp(mvec.Dim());
    Map_SolToMesh (Re(svec), tmp); SetReal (mvec, tmp);
    Map_SolToMesh (Im(svec), tmp); SetImag (mvec, tmp);
}

void Raster::Map_GridToSol (const CVector &gvec, CVector &svec) const
{
    CVector bvec(blen);
    Map_GridToBasis (gvec, bvec);
    Map_BasisToSol (bvec, svec);
}

void Raster::Map_MeshToGrid (const Solution &msol, Solution &gsol, bool mapall)
    const
{
    for (int i = 0; i < msol.nParam(); i++)
        if (mapall || msol.active[i])
  	    Map_MeshToGrid (msol.param[i], gsol.param[i]);
}

void Raster::Map_GridToMesh (const Solution &gsol, Solution &msol, bool mapall)
    const
{
    for (int i = 0; i < gsol.nParam(); i++)
        if (mapall || gsol.active[i])
	    Map_GridToMesh (gsol.param[i], msol.param[i]);
}

void Raster::Map_GridToBasis (const Solution &gsol, Solution &bsol,
    bool mapall) const
{
    for (int i = 0; i < gsol.nParam(); i++)
        if (mapall || gsol.active[i])
	    Map_GridToBasis (gsol.param[i], bsol.param[i]);
}

void Raster::Map_BasisToGrid (const Solution &bsol, Solution &gsol,
    bool mapall) const
{
    for (int i = 0; i < bsol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_BasisToGrid (bsol.param[i], gsol.param[i]);
}

void Raster::Map_MeshToBasis (const Solution &msol, Solution &bsol,
    bool mapall) const
{
    for (int i = 0; i < msol.nParam(); i++)
        if (mapall || msol.active[i])
	    Map_MeshToBasis (msol.param[i], bsol.param[i]);
}

void Raster::Map_BasisToMesh (const Solution &bsol, Solution &msol,
    bool mapall) const
{
    for (int i = 0; i < bsol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_BasisToMesh (bsol.param[i], msol.param[i]);
}

void Raster::Map_BasisToSol (const Solution &bsol, Solution &ssol,
    bool mapall) const
{
    for (int i = 0; i < bsol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_BasisToSol (bsol.param[i], ssol.param[i]);
}

void Raster::Map_SolToBasis (const Solution &ssol, Solution &bsol,
    bool mapall) const
{
    for (int i = 0; i < ssol.nParam(); i++)
        if (mapall || bsol.active[i])
	    Map_SolToBasis (ssol.param[i], bsol.param[i]);
}

void Raster::Map_MeshToSol (const Solution &msol, Solution &ssol,
    bool mapall) const
{
    for (int i = 0; i < msol.nParam(); i++)
        if (mapall || msol.active[i])
	    Map_MeshToSol (msol.param[i], ssol.param[i]);
}

void Raster::Map_SolToMesh (const Solution &ssol, Solution &msol,
    bool mapall) const
{
    for (int i = 0; i < ssol.nParam(); i++)
        if (mapall || ssol.active[i])
	    Map_SolToMesh (ssol.param[i], msol.param[i]);
}

void Raster::Map_ActiveSolToMesh (const RVector &asol, Solution &msol) const
{
    int i, j = 0;
    for (i = 0; i < msol.nParam(); i++)
        if (msol.active[i]) {
	    RVector prm (asol, slen*j++, slen);
	    Map_SolToMesh (prm, msol.param[i]);
	}
}

void Raster::NeighbourGraph (int *&rowptr, int *&colidx, int &nzero) const
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

    int *browptr = new int[blen+1];
    int *bcolidx;
    browptr[0] = 0;
    int bnzero = 0;

    int nx = bdim[0], ny = bdim[1], nz = (dim>2 ? bdim[2]:1);
    int i, j, k, r, c, s, nr;
    int xofs, yofs, zofs;

    // neighbour position offsets
    int dx = 1, dy = bdim[0], dz = bdim[0]*bdim[1];
    int dx1 = -dx;
    int dx2 =  dx;
    int dy1 = -dy;
    int dy2 =  dy;
    int dz1 = -dz;
    int dz2 =  dz;

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

	bcolidx = new int[bnzero];
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

	bcolidx = new int[bnzero];
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
    rowptr = new int[slen+1];
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
    colidx = new int[nzero];
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

    int *browptr = new int[glen+1];
    int *bcolidx, *bconnect;
    int *rowptr, *colidx, *connect, nzero;
    browptr[0] = 0;
    int bnzero = 0;

    int nx = gdim[0], ny = gdim[1], nz = (dim>2 ? gdim[2]:1);
    int i, j, k, r, c, s, nr;
    int xofs, yofs, zofs;

    // neighbour position offsets
    int dx = 1, dy = gdim[0], dz = gdim[0]*gdim[1];
    int dx1 = -dx;
    int dx2 =  dx;
    int dy1 = -dy;
    int dy2 =  dy;
    int dz1 = -dz;
    int dz2 =  dz;

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

	bcolidx = new int[bnzero];
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

	bcolidx = new int[bnzero];
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
    rowptr = new int[slen+1];
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
    colidx = new int[nzero];
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

void Raster::MapGraph2Sol (const int *browptr, const int *bcolidx, const int bnzero, int *&srowptr, int *&scolidx, int &snzero) const
{
  /* convert the graph on basis map to graph on solution map
  */
    cout << "Setting up solution neighbour graph" << endl;

    srowptr = new int[slen+1];
    srowptr[0] = snzero = 0;

    int i, r, c, row, col, i0, i1;

    // step 1: build row pointer list
    for (row = r = 0; row < blen; row++) {
      i0 = browptr[row]; i1 = browptr[row+1];
      if(basis2sol[row] < 0)               // pixel is not in support
	continue;
      for(i = i0 ; i < i1; i++)   {        // look at neighbours.
	    col = bcolidx[i];
	    if(basis2sol[col] < 0)         // neighbour is not in support
	       continue;
	    snzero++;
      }
      srowptr[++r] = snzero;
    }
    // step 2: build column index list
    scolidx = new int[snzero];
    for (row = r = c=0; row < blen; row++) {
      i0 = browptr[row]; i1 = browptr[row+1];
      if(basis2sol[row] < 0)               // pixel is not in support
	continue;
//      cout << "columns " << srowptr[r] << "-" << srowptr[r+1] -1<< ": ";
      for(i = i0 ; i < i1; i++)   {        // look at neighbours.
	    col = bcolidx[i];
	    if(basis2sol[col] < 0)         // neighbour is not in support
	       continue;
//	    cout << "(" << row << "," << col << ")->("<< r << "," << basis2sol[col] << ");";
	    scolidx[c++] = basis2sol[col];
      }
      r++;
//      cout << endl;
    }
    cout << snzero << " entries in solution neighbour graph found" << endl;
}

void Raster::BuildRegularRestriction (const IVector &gdim, RCompRowMatrix &R)
{
    // Builds restriction matrix from grid gdim to gdim_1,
    // where gdim[i] = 2 gdim_1[i] - 1
    // thus gdim[i] = 2k-1, i=0,1[,2] for some k is required

    int i, j, k, m, idx;
    int x1 = gdim[0], y1 = gdim[1], z1 = (gdim.Dim() > 2 ? gdim[2] : 1);
    int x2 = (x1-1)/2, y2 = (y1-1)/2, z2 = (z1 > 1 ? (z1-1)/2 : 1);
    int len1 = x1*y1*z1, len2 = x2*y2*z2;
    int stridex = 1;
    int stridey = x1;
    int stridez = x1*y1;

    R.New (len2, len1);
    RVector row(len1);
    
    for (k = m = 0; k < z2; k++) {
        for (j = 0; j < y2; j++) {
	    for (i = 0; i < x2; i++) {
	        // build stencil
	        row.Clear();
	        // 1. in-plane
	        idx = k*stridez*2 + j*stridey*2 + i*2;
		row[idx] = 1.0; // node itself
		if (i > 0)
		    row[idx-stridex] = 0.5; // left neighbour
		if (i < x2-1)
		    row[idx+stridex] = 0.5; // right neighbour
		if (j > 0) {
		    row[idx-stridey] = 0.5; // bottom neighbour
		    if (i > 0)
		        row[idx-stridey-stridex] = 0.25; // bottom left
		    if (i < x2-1)
		        row[idx-stridey+stridex] = 0.25; // bottom right
		}
		if (j < y2-1) {
		    row[idx+stridey] = 0.5; // top neighbour
		    if (i > 0)
		        row[idx+stridey-stridex] = 0.25; // top left
		    if (i < x2-1)
		        row[idx+stridey+stridex] = 0.25; // top right
		}
		// 2. plane below
		if (k > 0) {
		    row[idx-stridez] = 0.5;
		    if (i > 0)    row[idx-stridez-stridex] = 0.25;
		    if (i < x2-1) row[idx-stridez+stridex] = 0.25;
		    if (j > 0) {
		        row[idx-stridez-stridey] = 0.25;
			if (i > 0)    row[idx-stridez-stridey-stridex] = 0.125;
			if (i < x2-1) row[idx-stridez-stridey+stridex] = 0.125;
		    }
		    if (j < y2-1) {
		        row[idx-stridez+stridey] = 0.25;
			if (i > 0)    row[idx-stridez+stridey-stridex] = 0.125;
			if (i < x2-1) row[idx-stridex+stridey+stridez] = 0.125;
		    }
		}
		// 3. plane above
		if (k < z2-1) {
		    row[idx+stridez] = 0.5;
		    if (i > 0)    row[idx+stridez-stridex] = 0.25;
		    if (i < x2-1) row[idx+stridez+stridex] = 0.25;
		    if (j > 0) {
		        row[idx+stridez-stridey] = 0.25;
			if (i > 0)    row[idx+stridez-stridey-stridex] = 0.125;
			if (i < x2-1) row[idx+stridez-stridey+stridex] = 0.125;
		    }
		    if (j < y2-1) {
		        row[idx+stridez+stridey] = 0.25;
			if (i > 0)    row[idx+stridez+stridey-stridex] = 0.125;
			if (i < x2-1) row[idx+stridez+stridey+stridex] = 0.125;
		    }
		}
		R.SetRow (m++, row); // assemble into restriction matrix
	    }
	}
    }
}

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


#define RESCALE_W0
RVector Raster::SmoothImage(const RVector &x, double sd) const
{
    // smooth image by sd in each dimension
    double max = -1;
    RVector sx(x);

    for(int k = 0; k < sx.Dim(); k++)
        max = (sx[k] > max ? sx[k] : max);
    LOGOUT_1PRM("Max solbasis (on input)= ",max);
    //    cout << "Max solbasis (on input) = " << max << endl;
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
    //    cout << "Max solbasis (after smoothing) = " << max << endl;
    LOGOUT_1PRM("Max solbasis (after smoothing) = ", max);
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
    LOGOUT_1PRM("Max solbasis (on input)= ",max);
    //    cout << "Max solbasis (on input) = " << max << endl;
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
    complex *dz = new complex[dim];

    for(int ng = 0; ng < ngim; ng++) {
      if(!iflags[ng]) // skip if not set
	continue;
      //      cout << "Image Jet " << ng << " size " << x.Dim() << " scale " << sd << endl;
      gx[ng].New(x.Dim()); // output images are in solbasis

    switch(dim) {
    case 2 :
      for(int iy = 0; iy < paddim[1]; iy ++) {
	int offy = iy*paddim[0];
	dz[1] = complex(0.0,d1[1][iy]);
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  dz[0] = complex(0.0,d1[0][ix]);

	  complex z(data0[2*(offy + ix)],data0[2*(offy + ix)+1]);
	  //	  cout << z << " " << (filter[0][ix]*filter[1][iy]) << " ";
	  z *= complex(filter[0][ix]*filter[1][iy],0.0); // no complex*real??
	  //       	  cout << z << " ";
	  for(int idim = 0; idim < dim; idim++) {
	    for(int p = 0; p < D2D[idim][ng] ; p++)
	      z *= dz[idim];
	  }
	  //	  cout << dz[0] << " " << dz[1] << " ";
	  //       	  cout << z << " "; 
	  data[2*(offy + ix)]  = (float)re(z);
	  data[2*(offy + ix)+1] = (float)im(z);
	  //	  cout << data[2*(offy + ix)] << " " << data[2*(offy + ix)+1] << endl;
	}
      }
      break;
    case 3 :
      cout <<"Dimension 3 !\n";
      for(int iz = 0; iz < paddim[2]; iz ++) {
       int offz = iz*paddim[0]*paddim[1];
       dz[2] = complex(0.0,d1[2][iz]);
       for(int iy = 0; iy < paddim[1]; iy ++) {
	int offtot = offz + iy*paddim[0];
	dz[1] = complex(0.0,d1[1][iy]);
	for(int ix = 0; ix < paddim[0]; ix ++) {
	  dz[0] = complex(0.0,d1[0][ix]);

	  complex z(data0[2*(offtot + ix)],data0[2*(offtot + ix)+1]);
	  //	  cout << z << " " << (filter[0][ix]*filter[1][iy]) << " ";
	  z *= complex(filter[0][ix]*filter[1][iy]*filter[2][iz],0.0); // no complex*real??
	  //       	  cout << z << " ";
	  for(int idim = 0; idim < dim; idim++) {
	    for(int p = 0; p < D3D[idim][ng] ; p++)
	      z *= dz[idim];
	  }
	
	  data[2*(offtot + ix)]  = (float)re(z);
	  data[2*(offtot + ix)+1] = (float)im(z);
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
    LOGOUT_1PRM("Max solbasis (after smoothing) = ", max);
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
