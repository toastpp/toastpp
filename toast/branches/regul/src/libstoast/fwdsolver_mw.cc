#define __FWDSOLVER_MW_CC
#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "fwdsolver_mw.h"
#include "timing.h"

using namespace std;

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, LSOLVER linsolver,
    double tol)
    : TFwdSolver<T> (mesh, linsolver, tol)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, char *solver, double tol)
    : TFwdSolver<T> (mesh, solver, tol)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, ParamParser &pp)
    : TFwdSolver<T> (mesh, pp)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::~TFwdSolverMW ()
{
    Cleanup();
}

// =========================================================================

template<class T>
TVector<T> TFwdSolverMW<T>::ProjectAll_wavel (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    const QMMesh *mesh = this->meshptr;
    int i, n = mesh->nlen(), nq = mesh->nQ;
    int nofwavel = sol.nofwavel;
    int w;

    TVector<T> proj(nofwavel*mesh->nQM);

#ifndef TOAST_MPI
    TVector<T> *phi = new TVector<T>[nq];
    for (i = 0; i < nq; i++) phi[i].New(n);

    for (w = 0; w < nofwavel; w++) {
	this->Reset (*sol.swsol[w], omega);
	this->CalcFields (qvec, phi);
	TVector<T> proj_w (proj, w*mesh->nQM, mesh->nQM);
	proj_w = this->ProjectAll (mvec, phi, scl);
    }
    delete []phi;
#else
    int nfwd = nofwavel*nq;    // number of sources over all wavelengths

    int nfwd_part = (nfwd+sze-1)/sze;       // number of sources per processor
    int wq0 = min(nfwd,rnk*nfwd_part);    // low source index
    int wq1 = min(nfwd,(rnk+1)*nfwd_part);// high source index+1
    int r, qi, ofs, len, w_old = -1;
    TVector<T> phi(n);

    if (!qidx) { // first call: initialise counters and offsets
	qidx = new int[nq];           // data offset for each source
	projall_count = new int[sze]; // projection values per processor
	projall_ofs   = new int[sze]; // projection offset for each processor

	for (i = 0; i < nq; i++)
	    qidx[i] = (!i ? 0 : qidx[i-1] + mesh->nQMref[i-1]);

	for (r = 0; r < sze; r++) {
	    int wwq0 = min(nfwd,r*nfwd_part);
	    int wwq1 = min(nfwd,(r+1)*nfwd_part);
	    projall_count[r] = 0;
	    if (wwq0 < wwq1) {
		for (i = wwq0; i < wwq1; i++) {
		    qi = i%nq;
		    projall_count[r] += mesh->nQMref[qi];
		}
	    }
	    projall_ofs[r] = (!r ? 0 : projall_ofs[r-1] + projall_count[r-1]);
	}
    }

    if (wq0 < wq1) {
	for (i = wq0; i < wq1; i++) {
	    w   = i/nq;                   // wavelength
	    qi  = i%nq;                   // source
	    ofs = w*mesh->nQM + qidx[qi]; // data offset
	    len = mesh->nQMref[qi];       // data length
	    if (w != w_old) {             // initialise wavelength
		this->Reset (*sol.swsol[w], omega);
		w_old = w;
	    }
	    CalcField(qvec.Row(qi), phi);
	    TVector<T> proj_i(proj, ofs, len);
	    proj_i = ProjectSingle (mesh, qi, mvec, phi);
	}
    }

    MPI_Allgatherv (MPI_IN_PLACE, 0, mpitp, proj.data_buffer(),
		    projall_count, projall_ofs, mpitp, MPI_COMM_WORLD);

    if (scl == DATA_DEFAULT) scl = this->dscale;
    if (scl != DATA_LIN) proj = log(proj);
#endif

    return proj;
}

// =========================================================================

template<>
RVector RFwdSolverMW::ProjectAll_wavel_real (const RCompRowMatrix &qvec,
    const RCompRowMatrix &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    // real case: nothing to do
    return ProjectAll_wavel (qvec, mvec, sol, omega, scl);
}

template<>
RVector CFwdSolverMW::ProjectAll_wavel_real (const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    // complex case: split real and imaginary parts
    CVector cproj = ProjectAll_wavel (qvec, mvec, sol, omega, scl);

    const QMMesh *mesh = this->meshptr;
    int nofwavel = sol.nofwavel;
    int nqm = mesh->nQM;
    int n = cproj.Dim();
    int i;
    RVector proj(n*2);
    for (i = 0; i < nofwavel; i++) {
        CVector cproj_i(cproj, i*nqm, nqm);
	RVector proj_re (proj, i*2*nqm, nqm);
	RVector proj_im (proj, (i*2+1)*nqm, nqm);
	proj_re = Re(cproj_i);
	proj_im = Im(cproj_i);
    }
    return proj;
}

// =========================================================================

template<class T>
void TFwdSolverMW<T>::Setup ()
{
#ifdef TOAST_MPI
    sze   = TMPI::Size();
    rnk   = TMPI::Rank();
    mpitp = TMPI::MPIType<T>();
    projall_count = 0;
    projall_ofs   = 0;
    qidx          = 0;
#endif // TOAST_MPI
}

// =========================================================================

template<class T>
void TFwdSolverMW<T>::Cleanup ()
{
#ifdef TOAST_MPI
    if (projall_count) {
	delete []projall_count;
	projall_count = 0;
    }
    if (projall_ofs) {
	delete []projall_ofs;
	projall_ofs = 0;
    }
    if (qidx) {
	delete []qidx;
	qidx = 0;
    }
#endif // !TOAST_MPI
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class STOASTLIB TFwdSolverMW<double>;
template class STOASTLIB TFwdSolverMW<std::complex<double> >;

#endif // NEED_EXPLICIT_INSTANTIATION
