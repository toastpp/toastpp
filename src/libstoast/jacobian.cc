// =========================================================================
// jacobian.cc (libstoast)
// Implementation of Jacobian calculation
// =========================================================================

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "timing.h"

// ============================================================================
// Local prototypes

// CW case Jacobians
void GenerateJacobian_cw_grid (const Raster *raster, const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua, RDenseMatrix *Jkap);

void GenerateJacobian_cw_mesh (const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua, RDenseMatrix *Jkap, bool elbasis);

// Frequency-domain case Jacobians
void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J);

void GenerateJacobian_mesh (const QMMesh *mesh, const CCompRowMatrix &mvec,
    const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J, bool elbasis);

void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J);

void GenerateJacobian_mesh (const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J, bool elbasis);

template<class T>
TVector<T> IntFG (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g);

template<class T>
TVector<T> IntGradFGradG (const Mesh &mesh,
    const TVector<T> &f, const TVector<T> &g);

template<class T>
TVector<T> IntFG_el (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g);

template<class T>
TVector<T> IntGradFGradG_el (const Mesh &mesh,
    const TVector<T> &f, const TVector<T> &g);

// ============================================================================
// PMDF calculations

// absorption PMDF
template<class T>
inline TVector<T> PMDF_mua (const TVector<T> &dphi, const TVector<T> &aphi)
{
    return -dphi*aphi;
}

// diffusion PMDF
template<class T>
inline TVector<T> PMDF_kap (const TVector<T> *Gdphi, const TVector<T> *Gaphi,
    int dim)
{
    TVector<T> pmdf(Gdphi[0].Dim());
    for (int j = 0; j < dim; j++)
	pmdf -= Gdphi[j]*Gaphi[j];
    return pmdf;
}

// convert to log data PMDF for measurement m
template<class T>
inline TVector<T> PMDF_log (const TVector<T> &pmdf, T &m)
{
    return pmdf/m;
}


// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are calculated from mvec if required.

void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J)
{
    if (toastVerbosity > 0)
	std::cout << "Jacobian:" << std::endl;

    if (raster == RASTER_NDBASIS)
	GenerateJacobian_mesh (mesh, mvec, dphi, aphi, dscale, J, false);
    else if (raster == RASTER_ELBASIS)
	GenerateJacobian_mesh (mesh, mvec, dphi, aphi, dscale, J, true);
    else
	GenerateJacobian_grid (raster, mesh, mvec, dphi, aphi, dscale, J);

    if (toastVerbosity > 0)
	std::cout << "--> Dimensions......" << J.nRows() << 'x' << J.nCols()
		  << std::endl;
}

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are passed in as arguments.

void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J)
{
    if (toastVerbosity > 0)
	std::cout << "Jacobian:" << std::endl;

    if (raster == RASTER_NDBASIS)
	GenerateJacobian_mesh (mesh, dphi, aphi, proj, dscale, J, false);
    else if (raster == RASTER_ELBASIS)
	GenerateJacobian_mesh (mesh, dphi, aphi, proj, dscale, J, true);
    else
	GenerateJacobian_grid (raster, mesh, dphi, aphi, proj, dscale, J);

    if (toastVerbosity > 0)
	std::cout << "--> Dimensions......" << J.nRows() << 'x' << J.nCols()
		  << std::endl;
}

// ============================================================================
// Generate raw Jacobian matrix in grid basis (if raster is provided) or mesh
// basis

void GenerateJacobian_cw (const Raster *raster, const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua, RDenseMatrix *Jkap)
{
    if (toastVerbosity > 0)
	std::cout << "Jacobian:" << std::endl;

    if (raster == RASTER_NDBASIS)
	GenerateJacobian_cw_mesh (mesh, dphi, aphi, Jmua, Jkap, false);
    else if (raster == RASTER_ELBASIS)
	GenerateJacobian_cw_mesh (mesh, dphi, aphi, Jmua, Jkap, true);
    else
	GenerateJacobian_cw_grid (raster, mesh, dphi, aphi, Jmua, Jkap);

    if (toastVerbosity > 0) {
	if (Jmua)
	    std::cout << "--> Dim(Jmua)......." << Jmua->nRows() << 'x'
		      << Jmua->nCols() << std::endl;
	if (Jkap)
	    std::cout << "--> Dim(Jkap)......." << Jkap->nRows() << 'x'
		      << Jkap->nCols() << std::endl;
    }
}

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are calculated from mvec if required.

void GenerateJacobian_cw (const Raster *raster, const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrix *Jmua, RDenseMatrix *Jkap)
{
    // Raw Jacobian
    GenerateJacobian_cw (raster, mesh, dphi, aphi, Jmua, Jkap);

    // Map Jacobian to log data: dy/dx -> dlog y/dx
    if (dscale == DATA_LOG) {
	RVector iproj = inv (ProjectAll (mesh, mvec, dphi, DATA_LIN));
	if (Jmua) Jmua->RowScale (iproj);
	if (Jkap) Jkap->RowScale (iproj);
    }

    if (toastVerbosity > 0)
	std::cout << "--> Data scale......"
		  << (dscale == DATA_LOG ? "LOG":"LIN") << std::endl;
}

// ============================================================================
// Generate Jacobian matrix in grid basis (if raster is provided) or mesh basis
// Projections are passed in as arguments.

void GenerateJacobian_cw (const Raster *raster, const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi, const RVector *proj,
    DataScale dscale, RDenseMatrix *Jmua, RDenseMatrix *Jkap)
{
    // Raw Jacobian
    GenerateJacobian_cw (raster, mesh, dphi, aphi, Jmua, Jkap);

    // Map Jacobian to log data: dy/dx -> dlog y/dx
    if (dscale == DATA_LOG) {
	RVector iproj = inv (*proj);
	if (Jmua) Jmua->RowScale (iproj);
	if (Jkap) Jkap->RowScale (iproj);
    }
}

// ============================================================================
// Generate Jacobian matrix in grid basis defined by 'raster'
// Projections are calculated from mvec if required

void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J)
{
    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Grid" << std::endl;

    int i, j, jj, k, idx, dim, nQ, nM, nQM, slen, glen;
    const IVector &gdim = raster->GDim();
    const IVector &bdim = raster->BDim();
    const RVector &gsize = raster->GSize();
    dim  = raster->Dim();
    glen = raster->GLen();
    slen = raster->SLen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;

    CVector pmdf(glen*2);                    // a single row of J
    CVector pmdf_mua(pmdf, 0,    glen);      // absorption part
    CVector pmdf_kap(pmdf, glen, glen);      // diffusion part
    CVector pmdf_basis(slen);                // PMDF mapped into solution basis

    CVector proj;
    CVector cdfield(glen);                   // a single forward field
    CVector *cafield = new CVector[nM];      // array of all adjoint fields
    CVector *cdfield_grad = new CVector[dim];// forward field gradient
    CVector *cafield_grad = new CVector[dim];// adjoint field gradient

    // resample all measurement fields to fine grid
    //cout << "Allocating " << glen*nM*8*2 << " bytes for cafield" << endl;
    for (i = 0; i < nM; i++) {
        cafield[i].New (glen);
	raster->Map_MeshToGrid (aphi[i], cafield[i]);
    }

    for (i = idx = 0; i < nQ; i++) {

        // resample field and field gradient for source i to fine grid
	raster->Map_MeshToGrid (dphi[i], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, raster->Elref());

	proj = ProjectSingle (mesh, i, mvec, dphi[i], DATA_LIN);

	for (j = jj = 0; j < nM; j++) {
	    if (!mesh->Connected (i,j)) continue;

	    // measurement field gradient
	    ImageGradient (gdim, gsize, cafield[j], cafield_grad,
	        raster->Elref());

	    pmdf_mua = PMDF_mua (cdfield, cafield[j]);
	    pmdf_kap = PMDF_kap (cdfield_grad, cafield_grad, dim);

	    if (dscale == DATA_LOG)   // map to log data
		pmdf = PMDF_log (pmdf, proj[jj]);

	    // map into solution basis: 1. absorption
	    raster->Map_GridToSol (pmdf_mua, pmdf_basis);
	    for (k = 0; k < slen; k++) {
	        J(idx,k)          = pmdf_basis[k].real();
	        J(idx+nQM,k)      = pmdf_basis[k].imag();
	    }
	    // map into solution basis: 2. diffusion
	    raster->Map_GridToSol (pmdf_kap, pmdf_basis);
	    for (k = 0; k < slen; k++) {
	        J(idx,k+slen)     = pmdf_basis[k].real();
	        J(idx+nQM,k+slen) = pmdf_basis[k].imag();
	    }

	    idx++;
	    jj++;
	}
    }

    // scale with voxel size
    double sz = 1.0;
    for (i = 0; i < dim; i++)
	sz *= gsize[i] / (bdim[i]-1.0);
    J *= sz;

    delete []cafield;
    delete []cdfield_grad;
    delete []cafield_grad;    
}


// ============================================================================
// Generate Jacobian matrix in grid basis defined by 'raster'
// Projections are passed as argument

void GenerateJacobian_grid (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J)
{
    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Grid" << std::endl;

    int i, j, jj, k, idx, dim, nQ, nM, nQM, slen, glen;
    const IVector &gdim = raster->GDim();
    const IVector &bdim = raster->BDim();
    const RVector &gsize = raster->GSize();
    dim  = raster->Dim();
    glen = raster->GLen();
    slen = raster->SLen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;

    CVector pmdf(glen*2);                    // a single row of J
    CVector pmdf_mua(pmdf, 0,    glen);      // absorption part
    CVector pmdf_kap(pmdf, glen, glen);      // diffusion part
    CVector pmdf_basis(slen);                // PMDF mapped into solution basis

    CVector cdfield(glen);                   // a single forward field
    CVector *cafield = new CVector[nM];      // array of all adjoint fields
    CVector *cdfield_grad = new CVector[dim];// forward field gradient
    CVector *cafield_grad = new CVector[dim];// adjoint field gradient

    // resample all measurement fields to fine grid
    //cout << "Allocating " << glen*nM*8*2 << " bytes for cafield" << endl;
    for (i = 0; i < nM; i++) {
        cafield[i].New (glen);
	raster->Map_MeshToGrid (aphi[i], cafield[i]);
    }

    for (i = idx = 0; i < nQ; i++) {

        // resample field and field gradient for source i to fine grid
	raster->Map_MeshToGrid (dphi[i], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, raster->Elref());

	for (j = jj = 0; j < nM; j++) {
	    if (!mesh->Connected (i,j)) continue;

	    // measurement field gradient
	    ImageGradient (gdim, gsize, cafield[j], cafield_grad,
	        raster->Elref());

	    pmdf_mua = PMDF_mua (cdfield, cafield[j]);
	    pmdf_kap = PMDF_kap (cdfield_grad, cafield_grad, dim);

	    if (dscale == DATA_LOG)   // map to log data
		pmdf = PMDF_log (pmdf, (*proj)[idx]);

	    // map into solution basis: 1. absorption
	    raster->Map_GridToSol (pmdf_mua, pmdf_basis);
	    for (k = 0; k < slen; k++) {
	        J(idx,k)          = pmdf_basis[k].real();
	        J(idx+nQM,k)      = pmdf_basis[k].imag();
	    }
	    // map into solution basis: 2. diffusion
	    raster->Map_GridToSol (pmdf_kap, pmdf_basis);
	    for (k = 0; k < slen; k++) {
	        J(idx,k+slen)     = pmdf_basis[k].real();
	        J(idx+nQM,k+slen) = pmdf_basis[k].imag();
	    }

	    idx++;
	    jj++;
	}
    }
    delete []cafield;
    delete []cdfield_grad;
    delete []cafield_grad;    

    // scale with voxel size
    double sz = 1.0;
    for (i = 0; i < dim; i++)
	sz *= gsize[i] / (bdim[i]-1.0);
    J *= sz;
}



// ============================================================================
// This version doesn't use base mapping and works directly on the mesh basis

void GenerateJacobian_mesh (const QMMesh *mesh, const CCompRowMatrix &mvec,
    const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J, bool elbasis)
{
    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Mesh ("
		  << (elbasis ? "elements":"nodal") << std::endl;

    int i, j, jj, k, n, nsol, idx, nQ, nM, nQM;
    n    = mesh->nlen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    nsol = (elbasis ? mesh->elen() : n);

    CVector pmdf (nsol*2);
    CVector pmdf_mua (pmdf, 0, nsol);
    CVector pmdf_kap (pmdf, n, nsol);

    CVector proj;
    CVector cdfield(n);

	if (elbasis && !mesh->nbhrs)
		mesh->SetupNeighbourList();

    for (i = idx = 0; i < nQ; i++) {

	proj = ProjectSingle (mesh, i, mvec, dphi[i], DATA_LIN);

	for (j = jj = 0; j < nM; j++) {

	    if (!mesh->Connected (i,j)) continue;
	    if (elbasis) {
		pmdf_mua = IntFG_el (*mesh, dphi[i], aphi[j]);
		pmdf_kap = IntGradFGradG_el (*mesh, dphi[i], aphi[j]);
	    } else {
		pmdf_mua = IntFG (*mesh, dphi[i], aphi[j]);
		pmdf_kap = IntGradFGradG (*mesh, dphi[i], aphi[j]);
	    }

	    if (dscale == DATA_LOG)   // map to log data
		pmdf = PMDF_log (pmdf, proj[jj]);

	    // map into solution basis
	    for (k = 0; k < nsol; k++) {
	        J(idx,k)          = -pmdf_mua[k].real();
	        J(idx+nQM,k)      = -pmdf_mua[k].imag();
	        J(idx,k+nsol)     = -pmdf_kap[k].real();
	        J(idx+nQM,k+nsol) = -pmdf_kap[k].imag();
	    }

	    idx++;
	    jj++;
	}
    }
}

// ============================================================================
// This version doesn't use base mapping and works directly on the mesh basis
// The projection vector is passed in rather than calculated from mvec

void GenerateJacobian_mesh (const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J, bool elbasis)
{
    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Mesh ("
		  << (elbasis ? "elements":"nodal") << std::endl;

    int i, j, jj, k, n, nsol, idx, nQ, nM, nQM;
    n    = mesh->nlen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    nsol = (elbasis ? mesh->elen() : n);

    CVector pmdf (nsol*2);
    CVector pmdf_mua (pmdf, 0, nsol);
    CVector pmdf_kap (pmdf, nsol, nsol);

    //CVector proj;
    //CVector cdfield(n);

	if (elbasis && !mesh->nbhrs)
		mesh->SetupNeighbourList();

    for (i = idx = 0; i < nQ; i++) {

	//proj.New (mesh->nQMref[i]);
	//Project (*mesh, i, mvec, dphi[i], proj);

	for (j = jj = 0; j < nM; j++) {

	    if (!mesh->Connected (i,j)) continue;
	    if (elbasis) {
            pmdf_mua = IntFG_el (*mesh, dphi[i], aphi[j]);
		    pmdf_kap = IntGradFGradG_el (*mesh, dphi[i], aphi[j]);
	    } else {
            pmdf_mua = IntFG (*mesh, dphi[i], aphi[j]);
            pmdf_kap = IntGradFGradG (*mesh, dphi[i], aphi[j]);
	    }

	    if (dscale == DATA_LOG)   // map to log data
		pmdf = PMDF_log (pmdf, (*proj)[idx]);

		// map into solution basis
	    for (k = 0; k < nsol; k++) {
	        J(idx,k)          = -pmdf_mua[k].real();
	        J(idx+nQM,k)      = -pmdf_mua[k].imag();
	        J(idx,k+nsol)     = -pmdf_kap[k].real();
	        J(idx+nQM,k+nsol) = -pmdf_kap[k].imag();
	    }

	    idx++;
	    jj++;
	}
    }
}


// ============================================================================
// Generate Jacobian matrix in grid basis defined by 'raster'
// This computes a 'raw' Jacobian dy/dx without any data or parameter scaling
// (y: cw intensity, x: mua and/or kappa)

#if THREAD_LEVEL==2

struct GENJAC_CW_GRID_THREADDATA {
    const Raster *raster;
    const QMMesh *mesh;
    const RVector *dphi;
    const RVector *aphi;
    RDenseMatrix *Jmua;
    RDenseMatrix *Jkap;
};

void GenerateJacobian_cw_grid_engine (task_data *td)
{
    int itask = td->proc;
    int ntask = td->np;
    GENJAC_CW_GRID_THREADDATA *thdata = (GENJAC_CW_GRID_THREADDATA*)td->data;
    int nq  = thdata->mesh->nQ;
    int nm  = thdata->mesh->nM;
    int q0  = (itask*nq)/ntask;
    int q1  = ((itask+1)*nq)/ntask;
    int q, m;
    const IVector &gdim = thdata->raster->GDim();
    const RVector &gsize = thdata->raster->GSize();
    int dim  = thdata->raster->Dim();
    int glen = thdata->raster->GLen();
    int slen = thdata->raster->SLen();
    RVector cdfield(glen);                   // a single forward field
    RVector *cafield = new RVector[nm];      // array of all adjoint fields
    RVector *cdfield_grad = 0;
    RVector **cafield_grad = 0;
    RVector pmdf(glen);                      // a single row of J
    RVector pmdf_basis(slen);                // PMDF mapped into solution basis

    double *Jmua_ptr = 0;
    double *Jkap_ptr = 0;
    if (thdata->Jmua) {
	Jmua_ptr = thdata->Jmua->ValPtr()+thdata->mesh->Qofs[q0]*slen;
    }
    if (thdata->Jkap) {
	Jkap_ptr = thdata->Jkap->ValPtr()+thdata->mesh->Qofs[q0]*slen;
	cdfield_grad = new RVector[dim];
	cafield_grad = new RVector*[nm];
    }

    // resample all measurement fields to fine grid
    // note: this should be done in a separate parallel kernel
    for (m = 0; m < nm; m++) {
        cafield[m].New (glen);
	thdata->raster->Map_MeshToGrid (thdata->aphi[m], cafield[m]);

	if (thdata->Jkap) {
	    cafield_grad[m] = new RVector[dim];
	    ImageGradient (gdim, gsize, cafield[m], cafield_grad[m]);
	}
    }

    for (q = q0; q < q1; q++) {
        // resample field for source i to fine grid
	thdata->raster->Map_MeshToGrid (thdata->dphi[q], cdfield);

	// generate gradient
	if (thdata->Jkap)
	    ImageGradient (gdim, gsize, cdfield, cdfield_grad,
                thdata->raster->Elref());

	// loop over detectors
	for (m = 0; m < nm; m++) {
	    if (!thdata->mesh->Connected (q,m)) continue;

	    // absorption component
	    if (thdata->Jmua) {
		pmdf = PMDF_mua (cdfield, cafield[m]);

		// map into solution basis
		thdata->raster->Map_GridToSol (pmdf, pmdf_basis);
		memcpy (Jmua_ptr, pmdf_basis.data_buffer(),
			slen*sizeof(double));
		Jmua_ptr += slen;
	    }

	    // diffusion component
	    if (thdata->Jkap) {
		pmdf = PMDF_kap (cdfield_grad, cafield_grad[m], dim);

		// map into solution basis
		thdata->raster->Map_GridToSol (pmdf, pmdf_basis);
		memcpy (Jkap_ptr, pmdf_basis.data_buffer(),
			slen*sizeof(double));
		Jkap_ptr += slen;
	    }
	}
    }
    
    delete []cafield;
    if (cdfield_grad)
	delete []cdfield_grad;
    if (cafield_grad) {
	for (m = 0; m < nm; m++)
	    delete []cafield_grad[m];
	delete []cafield_grad;
    }
}

#endif // THREAD_LEVEL

void GenerateJacobian_cw_grid (const Raster *raster, const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua, RDenseMatrix *Jkap)
{
    int i, dim, glen, slen;
    const IVector &bdim = raster->BDim();
    const RVector &gsize = raster->GSize();
    dim  = raster->Dim();
    glen = raster->GLen();
    slen = raster->SLen();

    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Grid" << std::endl;

    double t0 = walltic();

#if THREAD_LEVEL==2
    int nqm  = mesh->nQM;
    if (Jmua) Jmua->New (nqm, slen);
    if (Jkap) Jkap->New (nqm, slen);
    static GENJAC_CW_GRID_THREADDATA thdata;
    thdata.raster = raster;
    thdata.mesh   = mesh;
    thdata.dphi   = dphi;
    thdata.aphi   = aphi;
    thdata.Jmua   = Jmua;
    thdata.Jkap   = Jkap;
    Task::Multiprocess (GenerateJacobian_cw_grid_engine, &thdata);
#else
    const IVector &gdim = raster->GDim();
    int j, nQ, nM, nQM;
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;

    RVector pmdf(glen);                      // a single row of J
    RVector pmdf_basis(slen);                // PMDF mapped into solution basis

    RVector cdfield(glen);                   // a single forward field
    RVector *cafield = new RVector[nM];      // array of all adjoint fields
    RVector *cdfield_grad = 0;
    RVector **cafield_grad = 0;

    double *Jmua_ptr = 0;
    double *Jkap_ptr = 0;
    if (Jmua) {
	Jmua->New (nQM, slen);
	Jmua_ptr = Jmua->ValPtr();
    }
    if (Jkap) {
	Jkap->New (nQM, slen);
	Jkap_ptr = Jkap->ValPtr();
	cdfield_grad = new RVector[dim];
	cafield_grad = new RVector*[nM];
    }

    // resample all measurement fields to fine grid
    for (i = 0; i < nM; i++) {
        cafield[i].New (glen);
	raster->Map_MeshToGrid (aphi[i], cafield[i]);
	if (Jkap) {
	    cafield_grad[i] = new RVector[dim];
	    ImageGradient (gdim, gsize, cafield[i], cafield_grad[i]);
	}
    }

    for (i = 0; i < nQ; i++) {

        // resample field for source i to fine grid
	raster->Map_MeshToGrid (dphi[i], cdfield);
	// generate gradient
	if (Jkap)
	    ImageGradient (gdim, gsize, cdfield, cdfield_grad,
                raster->Elref());

	// loop over detectors
	for (j = 0; j < nM; j++) {
	    if (!mesh->Connected (i,j)) continue;

	    // absorption component
	    if (Jmua) {
		pmdf = PMDF_mua (cdfield, cafield[j]);
		// map into solution basis
		raster->Map_GridToSol (pmdf, pmdf_basis);
		memcpy (Jmua_ptr, pmdf_basis.data_buffer(),
			slen*sizeof(double));
		Jmua_ptr += slen;
	    }

	    // diffusion component
	    if (Jkap) {
		pmdf = PMDF_kap (cdfield_grad, cafield_grad[j], dim);

		// map into solution basis
		raster->Map_GridToSol (pmdf, pmdf_basis);
		memcpy (Jkap_ptr, pmdf_basis.data_buffer(),
			slen*sizeof(double));
		Jkap_ptr += slen;
	    }
	}
    }
    
    delete []cafield;
    if (cdfield_grad)
	delete []cdfield_grad;
    if (cafield_grad) {
	for (i = 0; i < nM; i++)
	    delete []cafield_grad[i];
	delete []cafield_grad;
    }
#endif // !THREAD_LEVEL

    // scale with voxel size
    double sz = 1.0;
    for (i = 0; i < dim; i++)
	sz *= gsize[i] / (bdim[i]-1.0);
    if (Jmua) *Jmua *= sz;
    if (Jkap) *Jkap *= sz;

    std::cerr << "T(Jacobian) = " << walltoc(t0) << std::endl;
}


// ============================================================================
// This version doesn't use base mapping and works directly on the mesh basis
// This computes a 'raw' Jacobian dy/dx without any data or parameter scaling
// (y: cw intensity, x: mua and/or kappa)

void GenerateJacobian_cw_mesh (const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua, RDenseMatrix *Jkap, bool elbasis)
{
    if (toastVerbosity > 0)
	std::cout << "--> Basis...........Mesh ("
		  << (elbasis ? "elements":"nodal") << ')' << std::endl;

    int i, j, n, nsol, nQ, nM, nQM;
    n    = mesh->nlen();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    nsol = (elbasis ? mesh->elen() : n);

    RVector pmdf (nsol);

    double *Jmua_ptr = 0;
    double *Jkap_ptr = 0;

    if (Jmua) {
	Jmua->New(nQM, nsol);
	Jmua_ptr = Jmua->ValPtr();
    }
    if (Jkap) {
	Jkap->New(nQM, nsol);
	Jkap_ptr = Jkap->ValPtr();
    }

	if (elbasis && !mesh->nbhrs)
		mesh->SetupNeighbourList();

    for (i = 0; i < nQ; i++) {

	// loop over detectors
	for (j = 0; j < nM; j++) {

	    if (!mesh->Connected (i,j)) continue;

	    // absorption component
	    if (Jmua) {
		if (elbasis)
		    pmdf = -IntFG_el (*mesh, dphi[i], aphi[j]);
		else
		    pmdf = -IntFG (*mesh, dphi[i], aphi[j]);
		
		// copy into Jacobian
		memcpy (Jmua_ptr, pmdf.data_buffer(), nsol*sizeof(double));
		Jmua_ptr += nsol;
	    }

	    // diffusion component
	    if (Jkap) {
		if (elbasis)
		    pmdf = -IntGradFGradG_el (*mesh, dphi[i], aphi[j]);
		else
		    pmdf = -IntGradFGradG (*mesh, dphi[i], aphi[j]);

		// copy into Jacobian
		memcpy (Jkap_ptr, pmdf.data_buffer(), nsol*sizeof(double));
		Jkap_ptr += nsol;
	    }
	}
    }
}

// ============================================================================

template<class T>
TVector<T> IntFG (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, k, nj, nk, bs;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
        pel   = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
	        nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFFF(i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFFF(i,j,k);
		}
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}

// ============================================================================

template<class T>
TVector<T> IntGradFGradG (const Mesh &mesh,
    const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, k, nj, nk, bs;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
		nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFDD (i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFDD (i,j,k);
		}
		// we exploit the fact that IntFDD(i,j,k) is symmetric in
		// j and k: IntFDD(i,j,k) = IntFDD(i,k,j), so that two terms
		// can be combined in each pass of the inner (k) loop
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}


// ============================================================================
// Compute the product of two nodal field on (piecewise constant)
// element basis

template<class T>
TVector<T> IntFG_el (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, ni, nj;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.elen());

	for (el = 0; el < mesh.elen(); el++) {
		pel   = mesh.elist[el];
		nnode = pel->nNode();
		node  = pel->Node;
		sum = (T)0;

		for (i = 0; i < nnode; i++) {
			ni = node[i];
			for (j = 0; j < nnode; j++) {
				nj = node[j];
				sum += (f[ni] * g[nj]) * pel->IntFF(i,j);
			}
		}
	    tmp[el] += sum;
    }
    return tmp;
}


// ============================================================================

template<class T>
TVector<T> IntGradFGradG_el (const Mesh &mesh,
    const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, ni, nj;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.elen());

    for (el = 0; el < mesh.elen(); el++) {
		pel = mesh.elist[el];
		nnode = pel->nNode();
		node  = pel->Node;
		sum = (T)0;

		for (i = 0; i < nnode; i++) {
			ni = node[i];
			for (j = 0; j < nnode; j++) {
				nj = node[j];
				sum += (f[ni] * g[nj]) * pel->IntDD (i,j);
			}
		}
	    tmp[el] += sum;
    }
    return tmp;
}


// ============================================================================

template<class T>
STOASTLIB void ImageGradient (const IVector &gdim, const RVector &gsize,
    const TVector<T> &im, TVector<T> *grad, const int *mask)
{
    // this should be done by Fourier transform
    int x, y, z, idx, dim = gdim.Dim();
    int nx = gdim[0], ny = gdim[1], nz = (dim >= 3 ? gdim[2]:1);
    int n = nx*ny*nz;

    int *lmask;
    if (mask) {
	lmask = (int*)mask;
    } else {
	lmask = new int[n];
	memset (lmask, 0, n*sizeof(int));
    }

    // x gradient
    double dx = gsize[0]/nx, ix = 1.0/dx, i2x = 0.5*ix;
    TVector<T> &gradx = grad[0];
    gradx.New (n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt   = (lmask[idx] >= 0);
		bool bleft  = (x > 0 && (lmask[idx-1] >= 0));
		bool bright = (x < nx-1 && (lmask[idx+1] >= 0));
		if (bleft && bright) {
		    gradx[idx] = (im[idx+1]-im[idx-1]) * i2x;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradx[idx] = 0.0;
		} else if (bleft) {
		    gradx[idx] = (im[idx]-im[idx-1]) * ix;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (bright) {
		    gradx[idx] = (im[idx+1]-im[idx]) * ix;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else {
		    gradx[idx] = (T)0;
		}
		idx++;
	    }
	}
    }

    // y gradient
    double dy = gsize[1]/ny, iy = 1.0/dy, i2y = 0.5*iy;
    TVector<T> &grady = grad[1];
    grady.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt  = (lmask[idx] >= 0);
		bool bup   = (y > 0 && (lmask[idx-nx] >= 0));
		bool bdown = (y < ny-1 && (lmask[idx+nx] >= 0));
		if (bup && bdown) {
		    grady[idx] = (im[idx+nx]-im[idx-nx]) * i2y;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    grady[idx] = 0.0;
		} else if (bup) {
		    grady[idx] = (im[idx]-im[idx-nx]) * iy;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (bdown) {
		    grady[idx] = (im[idx+nx]-im[idx]) * iy;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else {
		    grady[idx] = (T)0;
		}
		idx++;
	    }
	}
    }
    if (dim < 3) return;

    // z gradient
    double dz = gsize[2]/nz, iz = 1.0/dz, i2z = 0.5*iz;
    int stridez = nx*ny;
    TVector<T> &gradz = grad[2];
    gradz.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt = (lmask[idx] >= 0);
		bool bfront = (z > 0 && (lmask[idx-stridez] >= 0));
		bool bback  = (z < nz-1 && (lmask[idx+stridez] >= 0));
	        if (bfront && bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx-stridez]) * i2z;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradz[idx] = 0.0;
		} else if (bfront) {
		    gradz[idx] = (im[idx]-im[idx-stridez]) * iz;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else if (bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx]) * iz;
		    //INC_FLOPS_ADD(1);
		    //INC_FLOPS_MUL(1);
		} else {
		    gradz[idx] = (T)0;
		}
		idx++;
	    }
	}
    }

    if (!mask) delete []lmask;
}

// instantiation
template STOASTLIB void ImageGradient (const IVector &gdim,
    const RVector &gsize, const RVector &im, RVector *grad, const int *mask);
template STOASTLIB void ImageGradient (const IVector &gdim,
    const RVector &gsize, const CVector &im, CVector *grad, const int *mask);

template RVector IntFG (const Mesh &mesh,
    const RVector &f, const RVector &g);
template CVector IntFG (const Mesh &mesh,
    const CVector &f, const CVector &g);

template RVector IntGradFGradG (const Mesh &mesh,
    const RVector &f, const RVector &g);
template CVector IntGradFGradG (const Mesh &mesh,
    const CVector &f, const CVector &g);

template RVector IntFG_el (const Mesh &mesh,
    const RVector &f, const RVector &g);
template CVector IntFG_el (const Mesh &mesh,
    const CVector &f, const CVector &g);

template RVector IntGradFGradG_el (const Mesh &mesh,
    const RVector &f, const RVector &g);
template CVector IntGradFGradG_el (const Mesh &mesh,
    const CVector &f, const CVector &g);
