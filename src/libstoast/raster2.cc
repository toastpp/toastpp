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
    Bvw = 0;
    Buu_precon = 0;
    Bvv_precon = 0;
    Buu_Cholesky_L = 0;
    Buu_Cholesky_d = 0;
    Bvv_Cholesky_L = 0;
    Bvv_Cholesky_d = 0;
    D = 0;
}

// =========================================================================

Raster2::~Raster2 ()
{
    if (Buu) delete Buu;
    if (Bvv) delete Bvv;
    if (Buv) delete Buv;
    if (Bvw) {
        delete Bvw;
	Bvw = 0;
    }
    if (Buu_precon) {
	delete Buu_precon;
    }
    if (Bvv_precon) {
	delete Bvv_precon;
    }
    if (Buu_Cholesky_L) {
	delete Buu_Cholesky_L;
    }
    if (Buu_Cholesky_d) {
	delete Buu_Cholesky_d;
    }
    if (Bvv_Cholesky_L) {
	delete Bvv_Cholesky_L;
    }
    if (Bvv_Cholesky_d) {
	delete Bvv_Cholesky_d;
    }
    if (D) {
	delete D;
    }
}

// =========================================================================

void Raster2::Init ()
{
    if (Buu) {
	delete Buu;
    }
    if (Bvv) {
	delete Bvv;
    }
    if (Buv) {
	delete Buv;
    }
    if (Buu_precon) {
	delete Buu_precon;
    }
    if (Bvv_precon) {
	delete Bvv_precon;
    }
    if (Buu_Cholesky_L) {
	delete Buu_Cholesky_L;
    }
    if (Buu_Cholesky_d) {
	delete Buu_Cholesky_d;
    }
    if (Bvv_Cholesky_L) {
	delete Bvv_Cholesky_L;
    }
    if (Bvv_Cholesky_d) {
	delete Bvv_Cholesky_d;
    }
    if (D) {
	delete D;
    }

    // mass matrix for mesh basis
    std::cerr << "enter MassMatrix" << std::endl;
    Buu = meshptr->MassMatrix();
    std::cerr << "back from MassMatrix" << std::endl;

    // mass matrix for intrinsic basis
    Bvv = CreateBasisMassmat();
    std::cerr << "back from BasisMatrix" << std::endl;

    // mass matrix for the mixed bases
    Buv = CreateMixedMassmat();

    // basis->solution mapper
    D = CreateSolMapper();

    // preconditioners
    if (map_tol) {
	Buu_precon = new RPrecon_IC; Buu_precon->Reset (Buu);
	Bvv_precon = new RPrecon_IC; Bvv_precon->Reset (Bvv);
	Buu_Cholesky_L = 0;
	Buu_Cholesky_d = 0;
	Bvv_Cholesky_L = 0;
	Bvv_Cholesky_d = 0;
    } else {
	Buu_precon = 0;
	Bvv_precon = 0;
	int nlen = meshptr->nlen();
	int *rowptr, *colidx;
	Buu->SymbolicCholeskyFactorize (rowptr, colidx);
	Buu_Cholesky_L = new RCompRowMatrix(nlen, nlen, rowptr, colidx);
	delete []rowptr;
	delete []colidx;
	Buu_Cholesky_d = new RVector(nlen);
	Bvv->SymbolicCholeskyFactorize (rowptr, colidx);
	Bvv_Cholesky_L = new RCompRowMatrix(blen, blen, rowptr, colidx);
	delete []rowptr;
	delete []colidx;
	Bvv_Cholesky_d = new RVector(blen);
	CholeskyFactorize (*Buu, *Buu_Cholesky_L, *Buu_Cholesky_d, true);
	CholeskyFactorize (*Bvv, *Bvv_Cholesky_L, *Bvv_Cholesky_d, true);
    }

    glen = slen = blen;
}

// =========================================================================

const RCompRowMatrix *Raster2::GetBvw (const IVector &wdim)
{
    if (Bvw) {
        if (Bvw_dim == wdim) return Bvw;
	else delete Bvw;
    }
    Bvw = CreateBasisPixelMassmat (wdim);
    return Bvw;
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
    if (map_tol) {
	double tol = map_tol;
	int nit = PCG (*Bvv, ATx(*Buv,mvec), bvec, tol, Bvv_precon);
    } else {
	CholeskySolve (*Bvv_Cholesky_L, *Bvv_Cholesky_d, ATx(*Buv,mvec), bvec);
    }
}

// =========================================================================

void Raster2::Map_MeshToBasis (const CVector &mvec, CVector &bvec) const
{
    int nit = 0;
    RVector bvec_tmp(bvec.Dim());
    RVector mvec_tmp(mvec.Dim());

    mvec_tmp = Re(mvec);
    if (map_tol) {
	double tol = map_tol;
	nit += PCG (*Bvv, ATx(*Buv,mvec_tmp), bvec_tmp, tol, Bvv_precon);
    } else {
	CholeskySolve (*Bvv_Cholesky_L, *Bvv_Cholesky_d, ATx(*Buv,mvec_tmp), bvec_tmp);
    }
    SetReal(bvec, bvec_tmp);

    mvec_tmp = Im(mvec);
    if (map_tol) {
	double tol = map_tol;
	nit += PCG (*Bvv, ATx(*Buv,mvec_tmp), bvec_tmp, tol, Bvv_precon);
    } else {
	CholeskySolve (*Bvv_Cholesky_L, *Bvv_Cholesky_d, ATx(*Buv,mvec_tmp), bvec_tmp);
    }
    SetImag(bvec, bvec_tmp);
}

// =========================================================================

void Raster2::Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
{
    if (map_tol) {
	double tol = map_tol;
	int nit = PCG (*Buu, Ax(*Buv,bvec), mvec, tol, Buu_precon);
    } else {
	CholeskySolve (*Buu_Cholesky_L, *Buu_Cholesky_d, Ax(*Buv,bvec), mvec);
    }
}

// =========================================================================

void Raster2::Map_BasisToMesh (const CVector &bvec, CVector &mvec) const
{
    int nit = 0;
    RVector bvec_tmp(bvec.Dim());
    RVector mvec_tmp(mvec.Dim());

    bvec_tmp = Re(bvec);
    if (map_tol) {
	double tol = map_tol;
	nit += PCG (*Buu, Ax(*Buv,bvec_tmp), mvec_tmp, tol, Buu_precon);
    } else {
	CholeskySolve (*Buu_Cholesky_L, *Buu_Cholesky_d, Ax(*Buv,bvec_tmp), mvec_tmp);
    }
    SetReal(mvec, mvec_tmp);

    bvec_tmp = Im(bvec);
    if (map_tol) {
	double tol = map_tol;
	nit += PCG (*Buu, Ax(*Buv,bvec_tmp), mvec_tmp, tol, Buu_precon);
    } else {
	CholeskySolve (*Buu_Cholesky_L, *Buu_Cholesky_d, Ax(*Buv,bvec_tmp), mvec_tmp);
    }
    SetImag(mvec, mvec_tmp);
}

// =========================================================================

void Raster2::Map_BasisToSol (const RVector &bvec, RVector &svec) const
{
    svec = bvec;
}

// =========================================================================

void Raster2::Map_BasisToSol (const CVector &bvec, CVector &svec) const
{
    svec = bvec;
}

// =========================================================================

void Raster2::Map_SolToBasis (const RVector &svec, RVector &bvec) const
{
    bvec = svec;
}

// =========================================================================

void Raster2::Map_SolToBasis (const CVector &svec, CVector &bvec) const
{
    bvec = svec;
}

// =========================================================================

void Raster2::Map_MeshToSol (const RVector &mvec, RVector &svec) const
{
    Map_MeshToBasis (mvec, svec);
}

// =========================================================================

void Raster2::Map_MeshToSol (const CVector &mvec, CVector &svec) const
{
    Map_MeshToBasis (mvec, svec);
}

// =========================================================================

void Raster2::Map_SolToMesh (const RVector &svec, RVector &mvec) const
{
    Map_BasisToMesh (svec, mvec);
}

// =========================================================================

void Raster2::Map_SolToMesh (const CVector &svec, CVector &mvec) const
{
    Map_BasisToMesh (svec, mvec);
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
