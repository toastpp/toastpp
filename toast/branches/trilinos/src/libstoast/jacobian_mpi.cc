// =========================================================================
// jacobian_mpi.cc (libstoast)
// Implementation of Jacobian calculation (distributed version)
// =========================================================================

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "dnsmatrix_mpi.h"
#include "jacobian_mpi.h"

using namespace toast;

// ============================================================================
// Local prototypes

void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg);

void GenerateJacobian_cw_mua_grid (const Raster *raster, const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrixMPI &J);

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are calculated from mvec if required.

void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg)
{
    if (raster) {
	GenerateJacobian_grid (raster, mesh, mvec, dphi, aphi, dscale,
			       Jmod, Jarg);
    } else {
	//GenerateJacobian_mesh (mesh, mvec, dphi, aphi, dscale, Jmod, Jarg);
    }
}

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are passed in as arguments.

void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg)
{
    if (raster) {
	//GenerateJacobian_grid (raster, mesh, dphi, aphi, proj, dscale,
	//		       Jmod, Jarg);
    } else {
	//GenerateJacobian_mesh (mesh, dphi, aphi, proj, dscale, Jmod, Jarg);
    }
}

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are calculated from mvec if required.

void GenerateJacobian_cw_mua (const Raster *raster, const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrixMPI &J)
{
    if (raster) {
	GenerateJacobian_cw_mua_grid (raster, mesh, mvec, dphi, aphi, dscale,
            J);
    } else {
	//GenerateJacobian_cw_mua_mesh (mesh, mvec, dphi, aphi, dscale, J);
    }
}

// ============================================================================
// Generate Jacobian matrix in grid basis defined by 'raster'
// Projections are calculated from mvec if required

void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg)
{
    int i, j, k, jj, idx, dim, glen, slen, nQ, nM;
    bool q_ok;
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    dim  = raster->Dim();
    glen = raster->GLen();
    slen = raster->SLen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;

    int r0, r1;
    Jmod.MPIRange (&r0, &r1);
    // Row range for current process: This is assumed to be the same for Jarg

    CVector pmdf(glen*2);                    // a single row of J
    CVector pmdf_mua(pmdf, 0,    glen);      // absorption part
    CVector pmdf_kap(pmdf, glen, glen);      // diffusion part
    CVector pmdf_basis(slen);                // PMDF mapped into solution basis

    CVector proj;
    CVector cdfield(glen);
    CVector *cafield = new CVector[nM];      // array of all adjoint fields
    CVector *cdfield_grad = new CVector[dim];// forward field gradient
    CVector *cafield_grad = new CVector[dim];// adjoint field gradient

    for (i = idx = 0; i < nQ; i++) {
	q_ok = false;
	for (j = jj = 0; j < nM; j++) {
	    if (!mesh->Connected(i,j)) continue;
	    if (idx >= r0 && idx < r1) { // row in range of process
		if (!q_ok) {
		    raster->Map_MeshToGrid (dphi[i], cdfield);
		    ImageGradient (gdim, gsize, cdfield, cdfield_grad,
				   raster->Elref());
		    proj = ProjectSingle (mesh, i, mvec, dphi[i]);
		    q_ok = true;
		}
		if (!cafield[j].Dim()) {
		    cafield[j].New(glen);
		    raster->Map_MeshToGrid (aphi[j], cafield[j]);
		}

		ImageGradient (gdim, gsize, cafield[j], cafield_grad,
			       raster->Elref());

		pmdf_mua = PMDF_mua (cdfield, cafield[j]);
		pmdf_kap = PMDF_kap (cdfield_grad, cafield_grad, dim);

		if (dscale == DATA_LOG)   // map to log data
		    pmdf = PMDF_log (pmdf, proj[jj]);

		// map into solution basis: 1. absorption
		raster->Map_GridToSol (pmdf_mua, pmdf_basis);
		for (k = 0; k < slen; k++) {
		    Jmod(idx,k)      = pmdf_basis[k].re;
		    Jarg(idx,k)      = pmdf_basis[k].im;
		}
		// map into solution basis: 2. diffusion
		raster->Map_GridToSol (pmdf_kap, pmdf_basis);
		for (k = 0; k < slen; k++) {
		    Jmod(idx,k+slen) = pmdf_basis[k].re;
		    Jarg(idx,k+slen) = pmdf_basis[k].im;
		}

	    }
	    idx++;
	    jj++;
	}
    }

    delete []cdfield_grad;
}


// ============================================================================
// Generate Jacobian matrix in grid basis defined by 'raster'
// Projections are calculated from mvec if required

void GenerateJacobian_cw_mua_grid (const Raster *raster, const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrixMPI &J)
{
    int i, j, jj, k, idx, dim, nQ, nM, nQM, slen, glen;
    bool q_ok;
    const IVector &gdim = raster->GDim();
    const IVector &bdim = raster->BDim();
    const RVector &gsize = raster->GSize();
    dim  = raster->Dim();
    glen = raster->GLen();
    slen = raster->SLen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;

    int r0, r1;
    J.MPIRange (&r0, &r1); // Row range for current process

    RVector pmdf(glen);                      // a single row of J
    RVector pmdf_basis(slen);                // PMDF mapped into solution basis

    RVector proj;
    RVector cdfield(glen);                   // a single forward field
    RVector *cafield = new RVector[nM];      // array of all adjoint fields

    for (i = idx = 0; i < nQ; i++) {
	q_ok = false;

	for (j = jj = 0; j < nM; j++) {
	    if (!mesh->Connected (i,j)) continue;
	    if (idx >= r0 && idx < r1) { // row in range of process
		if (!q_ok) {
		    // resample field and field gradient to fine grid
		    raster->Map_MeshToGrid (dphi[i], cdfield);
		    proj = ProjectSingle (mesh, i, mvec, dphi[i]);
		    q_ok = true;
		}
		if (!cafield[j].Dim()) {
		    cafield[j].New(glen);
		    raster->Map_MeshToGrid (aphi[j], cafield[j]);
		}

		pmdf = PMDF_mua (cdfield, cafield[j]);

		if (dscale == DATA_LOG)   // map to log data
		    pmdf = PMDF_log (pmdf, proj[jj]);

		// map into solution basis
		raster->Map_GridToSol (pmdf, pmdf_basis);
		for (k = 0; k < slen; k++)
		    J(idx,k) = pmdf_basis[k];
	    }
	    idx++;
	    jj++;
	}
    }
    delete []cafield;
}


