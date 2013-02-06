#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

using namespace std;
using namespace toast;

// ==========================================================================
// class Raster_Pixel

Raster_Pixel::Raster_Pixel (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb)
    : Raster (_bdim, _gdim, mesh, bb)
{
    int i, j;

    if (bdim == gdim) {
	// basis and grid are identical
	G  = NULL;
	GI = NULL;
	C  = B;
	CI = BI;
    } else {
	// set up mapping between high-res grid and native pixel basis
	G  = Grid2LinPixMatrix (gdim, bdim, gelref);
	((RCompRowMatrix*)G)->Shrink();
	GI = LinPix2GridMatrix (gdim, bdim, gelref);
	((RCompRowMatrix*)GI)->Shrink();
	// set up mapping between mesh and native pixel basis
	C = new RCompRowMatrix;
	((RCompRowMatrix*)G)->AB (*(RCompRowMatrix*)B, *(RCompRowMatrix*)C);
	((RCompRowMatrix*)C)->Shrink();
	CI = new RCompRowMatrix;
	((RCompRowMatrix*)BI)->AB (*(RCompRowMatrix*)GI, *(RCompRowMatrix*)CI);
	((RCompRowMatrix*)CI)->Shrink();
    }

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

Raster_Pixel::~Raster_Pixel ()
{
    if (G)  {
	delete G;
	delete C;
    }
    if (GI) {
	delete GI;
	delete CI;
    }
    delete D;
}

// ==========================================================================

void Raster_Pixel::Map_GridToBasis (const RVector &gvec, RVector &bvec) const
{
    if (G) G->Ax (gvec, bvec);
    else   bvec = gvec;
}

// ==========================================================================

void Raster_Pixel::Map_GridToBasis (const CVector &gvec, CVector &bvec) const
{
    if (G) ((RCompRowMatrix*)G)->Ax_cplx (gvec, bvec);
    else   bvec = gvec;
}

// ==========================================================================

void Raster_Pixel::Map_BasisToGrid (const RVector &bvec, RVector &gvec) const
{
    if (GI) GI->Ax (bvec, gvec);
    else    gvec = bvec;
}

// ==========================================================================

void Raster_Pixel::Map_BasisToGrid (const CVector &bvec, CVector &gvec) const
{
    if (GI) ((RCompRowMatrix*)GI)->Ax_cplx (bvec, gvec);
    else    gvec = bvec;
}

// ==========================================================================

void Raster_Pixel::Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
{
    C->Ax (mvec, bvec);
}

// ==========================================================================

void Raster_Pixel::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    CI->Ax (bvec, mvec);
}

// ==========================================================================

void Raster_Pixel::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    D->Ax (bvec, svec);
}

// ==========================================================================

void Raster_Pixel::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    D->ATx (svec, bvec);
}

// ==========================================================================

void Raster_Pixel::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    RVector bvec(blen);
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// ==========================================================================

void Raster_Pixel::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    RVector bvec(blen);
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}
