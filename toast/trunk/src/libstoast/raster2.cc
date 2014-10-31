#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

// =========================================================================
// class Raster2
// =========================================================================

Raster2::Raster2 (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb, double _map_tol)
: Raster (_bdim, _gdim, mesh, bb), map_tol(_map_tol)
{
    xASSERT(bdim==gdim,
	    "This basis type doesn't support intemediate grid basis");

    Buu = 0;
    Bvv = 0;
    Buv = 0;
    Buu_precon = 0;
    Bvv_precon = 0;
    D = 0;
}

// =========================================================================

Raster2::~Raster2 ()
{
    if (Buu) delete Buu;
    if (Bvv) delete Bvv;
    if (Buv) delete Buv;
    if (Buu_precon) delete Buu_precon;
    if (Bvv_precon) delete Bvv_precon;
    if (D) delete D;
}

// =========================================================================

void Raster2::Init ()
{
    // mass matrix for mesh basis
    Buu = meshptr->MassMatrix();
    // mass matrix for intrinsic basis
    Bvv = CreateBasisMassmat();
    // mass matrix for the mixed bases
    Buv = CreateMixedMassmat();

    // basis->solution mapper
    D = CreateSolMapper();

    // preconditioners
    Buu_precon = new RPrecon_IC; Buu_precon->Reset (Buu);
    Bvv_precon = new RPrecon_IC; Bvv_precon->Reset (Bvv);
}

// =========================================================================

void Raster2::Map_GridToBasis (const RVector &gvec, RVector &bvec) const
{
    bvec = gvec; // NOP
}

// =========================================================================

void Raster2::Map_GridToBasis (const CVector &gvec, CVector &bvec) const
{
    bvec = gvec; // NOP
}

// =========================================================================

void Raster2::Map_BasisToGrid (const RVector &bvec, RVector &gvec) const
{
    gvec = bvec; // NOP
}

// =========================================================================

void Raster2::Map_BasisToGrid (const CVector &bvec, CVector &gvec) const
{
    gvec = bvec; // NOP
}

// =========================================================================

void Raster2::Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
{
    double tol = map_tol;
    int nit = PCG (*Bvv, ATx(*Buv,mvec), bvec, tol, Bvv_precon);
}

// =========================================================================

void Raster2::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    double tol = map_tol;
    int nit = PCG (*Buu, Ax(*Buv,bvec), mvec, tol, Buu_precon);
}

// =========================================================================

void Raster2::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    D->Ax (bvec, svec);
}

// =========================================================================

void Raster2::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    D->ATx (svec, bvec);
}

// =========================================================================

void Raster2::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    RVector bvec(blen);
    Map_MeshToBasis (mvec, bvec);
    Map_BasisToSol (bvec, svec);
}

// =========================================================================

void Raster2::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    RVector bvec(blen);
    Map_SolToBasis (svec, bvec);
    Map_BasisToMesh (bvec, mvec);
}

// =========================================================================

RCompRowMatrix *Raster2::CreateSolMapper () const
{
    // maybe this should go into the base class!

    int i, j;

    // formulate basis->solution mapping in sparse matrix
    idxtype *rowptr = new idxtype[slen+1];
    idxtype *colidx = new idxtype[slen];
    double *val = new double[slen];
    for (i = 0; i <= slen; i++) rowptr[i] = i; // each row has one entry
    for (i = 0; i < slen; i++) val[i] = 1.0;
    for (i = j = 0; i < blen; i++)
        if (bsupport[i] > 0.0) colidx[j++] = i;
    RCompRowMatrix *d = new RCompRowMatrix (slen, blen, rowptr, colidx,
	val);

    delete []rowptr;
    delete []colidx;
    delete []val;

    return d;
}
