#define __FWDSOLVER_CC
#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "fwdsolver_zslu.h"
#include "fwdsolver_cslu.h"
#ifdef MPI_FWDSOLVER
#include "toast_mpi.h"
#endif
#ifdef USE_CUDA_FLOAT
#include "toastcuda.h"
#endif

using namespace std;

// =========================================================================

template<class T>
TFwdSolver<T>::TFwdSolver (const QMMesh *mesh, LSOLVER linsolver, double tol)
{
    Setup ();
    meshptr = mesh;
    solvertp = linsolver;
    iterative_tol = tol;
}

// =========================================================================

template<class T>
TFwdSolver<T>::TFwdSolver (const QMMesh *mesh, const char *solver, double tol, int nth)
{
    Setup (nth);
    meshptr = mesh;
    SetLinSolver (solver, tol);
}

// =========================================================================

template<class T>
TFwdSolver<T>::TFwdSolver (const QMMesh *mesh, ParamParser &pp)
{
    Setup ();
    meshptr = mesh;
    ReadParams (pp);
}

// =========================================================================

template<class T>
TFwdSolver<T>::~TFwdSolver ()
{
    if (F)      delete F;
    if (FL)     delete FL;
    if (Fd)     delete Fd;
    if (B)      delete B;
    if (precon) delete precon;
    if (pphi)   delete []pphi;

    DeleteType ();
    CLEANUP_MPI();
}

// =========================================================================

template<class T>
void TFwdSolver<T>::Setup (int nth)
{
    // set initial and default parameters
    SuperLU = 0;
    dscale = DATA_LIN;
    solvertp = LSOLVER_ITERATIVE;
    precontp = PRECON_NULL;
    iterative_tol = 1e-6;
    iterative_maxit = 0;
    unwrap_phase = false;
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;
    pphi = 0;

    SetupType (nth);
    SETUP_MPI();
}

// =========================================================================

template<>
void TFwdSolver<std::complex<double> >::SetupType (int nth)
{
    SuperLU = new ZSuperLU (nth);
}

template<>
void TFwdSolver<std::complex<float> >::SetupType (int nth)
{
    SuperLU = new CSuperLU ();
}

template<class T>
void TFwdSolver<T>::SetupType (int nth)
{
}

// =========================================================================

template<>
void TFwdSolver<std::complex<double> >::DeleteType ()
{
    delete (ZSuperLU*)SuperLU;
}

template<>
void TFwdSolver<std::complex<float> >::DeleteType ()
{
    delete (CSuperLU*)SuperLU;
}

template<class T>
void TFwdSolver<T>::DeleteType ()
{
}

// =========================================================================

template<class T>
void TFwdSolver<T>::SetDataScaling (DataScale scl)
{
    xASSERT(scl != DATA_DEFAULT, "Invalid input argument");
    dscale = scl;
}

template<class T>
DataScale TFwdSolver<T>::GetDataScaling () const
{
    return dscale;
}

// =========================================================================

template<class T>
void TFwdSolver<T>::SetPrecon (PreconType type)
{
    if (solvertp == LSOLVER_DIRECT) return; // no preconditioner here
    if (type == precontp) return;           // nothing to do
    if (precon) delete precon;
    precontp = type;
    switch (precontp) {
    case PRECON_DIAG: precon = new TPrecon_Diag<T>; break;
    case PRECON_ICH:  precon = new TPrecon_IC<T>;   break;
    default:          precon = new TPrecon_Null<T>; break;
    }
}

template<>
void TFwdSolver<float>::Allocate ()
{
    int *rowptr, *colidx, nzero;
    int n = meshptr->nlen();

    // allocate system matrix
    meshptr->SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
#ifdef MPI_FWDSOLVER
    F = new FCompRowMatrixMPI (n, n, rowptr, colidx);
#else
    F = new FCompRowMatrix (n, n, rowptr, colidx);
#endif
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
	F->SymbolicCholeskyFactorize (rowptr, colidx);
	if (FL) delete FL;
	FL = new FCompRowMatrix (n, n, rowptr, colidx);
	if (Fd) delete Fd;
	Fd = new FVector (n);
	delete []rowptr;
	delete []colidx;
    } else {
	if (precon) delete precon;
	switch (precontp) {
	case PRECON_DIAG: precon = new FPrecon_Diag; break;
	case PRECON_ICH:  precon = new FPrecon_IC;   break;
	default:          precon = new FPrecon_Null; break;
	}
    }
}

template<>
void TFwdSolver<double>::Allocate ()
{
    idxtype *rowptr, *colidx;
    int nzero;
    int n = meshptr->nlen();

    // allocate system matrix
    meshptr->SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
#ifdef MPI_FWDSOLVER
    F = new RCompRowMatrixMPI (n, n, rowptr, colidx);
#else
    F = new RCompRowMatrix (n, n, rowptr, colidx);
#endif
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
	F->SymbolicCholeskyFactorize (rowptr, colidx);
	if (FL) delete FL;
	FL = new RCompRowMatrix (n, n, rowptr, colidx);
	if (Fd) delete Fd;
	Fd = new RVector (n);
	delete []rowptr;
	delete []colidx;
    } else {
	if (precon) delete precon;
	switch (precontp) {
	case PRECON_DIAG: precon = new RPrecon_Diag; break;
	case PRECON_ICH:  precon = new RPrecon_IC;   break;
	default:          precon = new RPrecon_Null; break;
	}
    }
}

template<>
void TFwdSolver<std::complex<double> >::Allocate ()
{
    idxtype *rowptr, *colidx;
	int nzero;
    int n = meshptr->nlen();

    // allocate system matrix
    meshptr->SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
#ifdef MPI_FWDSOLVER
    F = new CCompRowMatrixMPI (n, n, rowptr, colidx);
#else
    F = new CCompRowMatrix (n, n, rowptr, colidx);
#endif
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
	//lu_data.Setup(F);
	((ZSuperLU*)SuperLU)->Reset (F);
    } else {
	if (precon) delete precon;
	switch (precontp) {
	case PRECON_DIAG: precon = new CPrecon_Diag; break;
	case PRECON_ICH:  precon = new CPrecon_IC;   break;
	// NOTE: IC preconditioner shouldn't really work for symmetric
	// complex problem, but it seems to anyway
	default:          precon = new CPrecon_Null; break;
	}
    }
}

template<>
void TFwdSolver<std::complex<float> >::Allocate ()
{
    int *rowptr, *colidx, nzero;
    int n = meshptr->nlen();

    // allocate system matrix
    meshptr->SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
#ifdef MPI_FWDSOLVER
    F = new SCCompRowMatrixMPI (n, n, rowptr, colidx);
#else
    F = new SCCompRowMatrix (n, n, rowptr, colidx);
#endif
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
        ((CSuperLU*)SuperLU)->Reset (F);
	// lu_data.Setup(meshptr->nlen(), F);
    } else {
	if (precon) delete precon;
	switch (precontp) {
	case PRECON_DIAG: precon = new SCPrecon_Diag; break;
	case PRECON_ICH:  precon = new SCPrecon_IC;   break;
	default:          precon = new SCPrecon_Null; break;
	}
    }
}

template<>
void TFwdSolver<float>::AssembleSystemMatrix (const Solution &sol,
    double omega, bool elbasis)
{
	xASSERT (meshptr, "Mesh reference not defined."); 
    xASSERT(omega==0, "Nonzero omega parameter not allowed here");
	// real version

    RVector prm;

    // To improve accuracy, we assemble the system matrix in double
    // precision, and map it to single precision after assembly

    RCompRowMatrix FF (F->nRows(), F->nCols(), F->rowptr, F->colidx);
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, FF, &prm,
                    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, FF, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, FF, &prm,
		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);

    int i, nz = F->nVal();
    float *fval = F->ValPtr();
    double *dval = FF.ValPtr();
    for (i = 0; i < nz; i++) *fval++ = (float)*dval++;
 
#ifdef UNDEF
    F->Zero();
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
#endif
}

template<>
void TFwdSolver<double>::AssembleSystemMatrix (const Solution &sol,
    double omega, bool elbasis)
{
    xASSERT (meshptr, "Mesh reference not defined."); 
    xASSERT(omega==0, "Nonzero omega parameter not allowed here");
    // real version

    RVector prm;

    F->Zero();
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
}

template<>
void TFwdSolver<std::complex<float> >::AssembleSystemMatrix (
    const Solution &sol, double omega, bool elbasis)
{
    xASSERT (meshptr, "Mesh reference not defined."); 

    // complex version
    RVector prm;

    // To improve accuracy, we assemble the system matrix in double
    // precision, and map it to single precision after assembly

    CCompRowMatrix FF (F->nRows(), F->nCols(), F->rowptr, F->colidx);
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, FF, &prm,
		    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, FF, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, FF, &prm,
		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
    AddToSysMatrix (*meshptr, FF, omega, ASSEMBLE_iCFF);

    int i, nz = F->nVal();
    std::complex<float> *sval = F->ValPtr();
    std::complex<double> *cval = FF.ValPtr();
    for (i = 0; i < nz; i++)
        sval[i] = (std::complex<float>)cval[i];

#ifdef UNDEF
    F->Zero();
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, *F, &prm,
    		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
    AddToSysMatrix (*meshptr, *F, omega, ASSEMBLE_iCFF);
#endif
}

template<>
void TFwdSolver<std::complex<double> >::AssembleSystemMatrix (
    const Solution &sol, double omega, bool elbasis)
{
	xASSERT (meshptr, "Mesh reference not defined."); 

	// complex version
    RVector prm;

    F->Zero();
    prm = sol.GetParam (OT_CMUA);
    AddToSysMatrix (*meshptr, *F, &prm,
        elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    AddToSysMatrix (*meshptr, *F, &prm,
        elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    AddToSysMatrix (*meshptr, *F, &prm,
        elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
    AddToSysMatrix (*meshptr, *F, omega, ASSEMBLE_iCFF);
}

template<class T>
void TFwdSolver<T>::AssembleSystemMatrixComponent (const RVector &prm, int type)
{
	xASSERT (meshptr, "Mesh reference not defined."); 

	F->Zero();
    AddToSysMatrix (*meshptr, *F, &prm, type);
}

template<class T>
void TFwdSolver<T>::AssembleMassMatrix (const Mesh *mesh)
{
    if (!mesh) mesh = meshptr;
	xASSERT (meshptr, "Mesh reference not defined."); 

    if (!B) { // allocate on the fly
	idxtype *rowptr, *colidx;
	int nzero, n = mesh->nlen();
	mesh->SparseRowStructure (rowptr, colidx, nzero);
	B = new TCompRowMatrix<T> (n, n, rowptr, colidx);
	delete []rowptr;
	delete []colidx;
    } else {
	B->Zero();
    }

    AddToSysMatrix (*mesh, *B, (RVector*)0, ASSEMBLE_FF);
}

template<>
void TFwdSolver<float>::Reset (const Solution &sol, double omega, bool elbasis)
{
    // real version
    AssembleSystemMatrix (sol, omega, elbasis);
    if (solvertp == LSOLVER_DIRECT)
	CholeskyFactorize (*F, *FL, *Fd, true);
    else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<double>::Reset (const Solution &sol, double omega, bool elbasis)
{
    // real version
    AssembleSystemMatrix (sol, omega, elbasis);
    if (solvertp == LSOLVER_DIRECT)
	CholeskyFactorize (*F, *FL, *Fd, true);
    else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<std::complex<float> >::Reset (const Solution &sol,
    double omega, bool elbasis)
{
    // single complex version
    AssembleSystemMatrix (sol, omega, elbasis);
    if (solvertp == LSOLVER_DIRECT) {
	//lu_data.Setup (F);
	((CSuperLU*)SuperLU)->Reset (F);
    } else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<std::complex<double> >::Reset (const Solution &sol,
    double omega, bool elbasis)
{
    // complex version
    AssembleSystemMatrix (sol, omega, elbasis);
    if (solvertp == LSOLVER_DIRECT)
	//lu_data.Setup (F);
	((ZSuperLU*)SuperLU)->Reset (F);
    else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<float>::CalcField (const TVector<float> &qvec,
    TVector<float> &phi, IterativeSolverResult *res, int th) const
{
    // calculate the (real) photon density field for a given (real)
    // source distribution. Use only if data type is INTENSITY

    if (solvertp == LSOLVER_DIRECT) {
	CholeskySolve (*FL, *Fd, qvec, phi);
    } else {
	double tol = iterative_tol;
	int it = IterativeSolve (*F, qvec, phi, tol, precon, iterative_maxit);
	if (res) {
	    res->it_count = it;
	    res->rel_err = tol;
	}
    }
}

template<>
void TFwdSolver<double>::CalcField (const TVector<double> &qvec,
    TVector<double> &phi, IterativeSolverResult *res, int th) const
{
    // calculate the (real) photon density field for a given (real)
    // source distribution. Use only if data type is INTENSITY

    if (solvertp == LSOLVER_DIRECT) {
	CholeskySolve (*FL, *Fd, qvec, phi);
    } else {
	double tol = iterative_tol;
	int it = IterativeSolve (*F, qvec, phi, tol, precon, iterative_maxit);
	if (res) {
	    res->it_count = it;
	    res->rel_err = tol;
	}
    }
}

template<>
void TFwdSolver<std::complex<float> >::CalcField (
    const TVector<std::complex<float> > &qvec,
    TVector<std::complex<float> > &cphi, IterativeSolverResult *res, int th) const
{
    if (solvertp == LSOLVER_DIRECT) {
        ((CSuperLU*)SuperLU)->CalcField (qvec, cphi, res);
    } else {
        double tol = iterative_tol;
	int it = IterativeSolve (*F, qvec, cphi, tol, precon, iterative_maxit);
	if (res) {
	    res->it_count = it;
	    res->rel_err = tol;
	}
    }
}

template<>
void TFwdSolver<std::complex<double> >::CalcField (
    const TVector<std::complex<double> > &qvec,
    TVector<std::complex<double> > &cphi, IterativeSolverResult *res, int th) const
{
    // calculate the complex field for a given source distribution

#ifdef DO_PROFILE
    times (&tm);
    clock_t time0 = tm.tms_utime;
#endif
    if (solvertp == LSOLVER_DIRECT) {
	((ZSuperLU*)SuperLU)->CalcField (qvec, cphi, res, th);
        //SuperMatrix B, X;
	//int n = meshptr->nlen();
	//
	//doublecomplex *rhsbuf = (doublecomplex*)qvec.data_buffer();
	//doublecomplex *xbuf   = (doublecomplex*)cphi.data_buffer();
	//zCreate_Dense_Matrix (&B, n, 1, rhsbuf, n, SLU_DN, SLU_Z,SLU_GE);
	//zCreate_Dense_Matrix (&X, n, 1, xbuf, n, SLU_DN, SLU_Z, SLU_GE);
	//
	//lu_data.Solve (&B, &X);
	//
	//Destroy_SuperMatrix_Store (&B);
	//Destroy_SuperMatrix_Store (&X);
    } else {
        double tol = iterative_tol;
	int it = IterativeSolve (*F, qvec, cphi, tol, precon, iterative_maxit);
	if (res) {
	    res->it_count = it;
	    res->rel_err = tol;
	}
    }
#ifdef DO_PROFILE
    times (&tm);
    solver_time += (double)(tm.tms_utime-time0)/(double)HZ;
#endif
    //cerr << "Solve time = " << solver_time << endl;
}

// ==========================================================================
// CalcFields

#if THREAD_LEVEL==2

template<class T>
struct CALCFIELDS_THREADDATA {
    TCompRowMatrix<T> *F;
    const TCompRowMatrix<T> *qvec;
    TVector<T> *phi;
    TPreconditioner<T> *precon;
    double tol;
    int maxit;
    IterativeSolverResult *res;
};

template<class T>
void CalcFields_engine (task_data *td)
{
    int itask = td->proc;
    int ntask = td->np;
    CALCFIELDS_THREADDATA<T> *thdata = (CALCFIELDS_THREADDATA<T>*)td->data;
    int nq = thdata->qvec->nRows();
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int q, it, itsum = 0;
    double tol, err_max = 0.0;

    for (q = q0; q < q1; q++) {
        tol = thdata->tol;
	it = IterativeSolve (*thdata->F, thdata->qvec->Row(q), thdata->phi[q],
			tol, thdata->precon, thdata->maxit);
	itsum += it;
	if (tol > err_max) err_max = tol;
    }
    if (thdata->res) {
        Task::UserMutex_lock();
	thdata->res->it_count += itsum;
	if (err_max > thdata->res->rel_err) thdata->res->rel_err = err_max;
	Task::UserMutex_unlock();
    }
}
#ifdef NEED_EXPLICIT_INSTANTIATION
template void CalcFields_engine<double> (task_data *td);
template void CalcFields_engine<float> (task_data *td);
template void CalcFields_engine<std::complex<double> > (task_data *td);
template void CalcFields_engine<std::complex<float> > (task_data *td);
#endif

#endif // THREAD_LEVEL==2

template<>
void TFwdSolver<std::complex<double> >::CalcFields (const CCompRowMatrix &qvec,
    CVector *phi, IterativeSolverResult *res) const
{
    // calculate the fields for all sources
    static IterativeSolverResult s_res_single;
    IterativeSolverResult *res_single = 0;
    if (res) {
        res_single = &s_res_single;
	res->it_count = 0;
	res->rel_err = 0.0;
    }

    if (solvertp == LSOLVER_DIRECT) {
	((ZSuperLU*)SuperLU)->CalcFields (qvec, phi, res);
    //for (int i = 0; i < nq; i++) {
	//    CalcField (qvec.Row(i), phi[i], res);
	//}
    } else {
#if THREAD_LEVEL==2

        //dASSERT(g_tpool, ThreadPool not initialised);
        static CALCFIELDS_THREADDATA<std::complex<double> > thdata;
	thdata.F      = F;
	thdata.qvec   = &qvec;
	thdata.phi    = phi;
	thdata.precon = precon;
	thdata.tol    = iterative_tol;
	thdata.maxit  = iterative_maxit;
	thdata.res    = res;
#ifdef UNDEF
    int nq = qvec.nRows();
	g_tpool->ProcessSequence (CalcFields_engine<std::complex<double> >, &thdata,
				  0, nq, 1);
#endif
	Task::Multiprocess (CalcFields_engine<std::complex<double> >, &thdata);
#else
        int nq = qvec.nRows();
        CVector *qv = new CVector[nq];
	for (int i = 0; i < nq; i++) qv[i] = qvec.Row(i);
        IterativeSolve (*F, qv, phi, nq, iterative_tol, iterative_maxit,
            precon, res);
	delete []qv;
#endif // THREAD_LEVEL==2
    }
}

template<class T>
void TFwdSolver<T>::CalcFields (const TCompRowMatrix<T> &qvec,
    TVector<T> *phi, IterativeSolverResult *res) const
{
    // calculate the fields for all sources
    static IterativeSolverResult s_res_single;
    IterativeSolverResult *res_single = 0;
    if (res) {
        res_single = &s_res_single;
	res->it_count = 0;
	res->rel_err = 0.0;
    }

#ifndef MPI_FWDSOLVER
    int i, nq = qvec.nRows();

    if (solvertp == LSOLVER_DIRECT) {
        LOGOUT1_INIT_PROGRESSBAR ("CalcFields", 50, nq);
	for (i = 0; i < nq; i++) {
	    CalcField (qvec.Row(i), phi[i], res_single);
	    if (res) {
	        if (res_single->it_count > res->it_count)
		    res->it_count = res_single->it_count;
		if (res_single->rel_err > res->rel_err)
		    res->rel_err = res_single->rel_err;
	    }
	    LOGOUT1_PROGRESS(i);
	}
    } else {
#if THREAD_LEVEL==2
        //dASSERT(g_tpool, ThreadPool not initialised);
        static CALCFIELDS_THREADDATA<T> thdata;
	thdata.F      = F;
	thdata.qvec   = &qvec;
	thdata.phi    = phi;
	thdata.precon = precon;
	thdata.tol    = iterative_tol;
	thdata.maxit  = iterative_maxit;
	thdata.res    = res;
#ifdef UNDEF
	g_tpool->ProcessSequence (CalcFields_engine<T>, &thdata, 0, nq, 1);
#endif
	Task::Multiprocess (CalcFields_engine<T>, &thdata);
#else
        TVector<T> *qv = new TVector<T>[nq];
	for (i = 0; i < nq; i++) qv[i] = qvec.Row(i);
        IterativeSolve (*F, qv, phi, nq, iterative_tol, iterative_maxit,
            precon, res);
	delete []qv;
#endif
    }

#else

    int i, q0, q1, n = meshptr->nlen();
    int nq = qvec.nRows();

    // store all projections in a single vector, so we can gather the results
    TVector<T> phi_tot(nq*n);
    TVector<T> *phi_single = new TVector<T>[nq];
    for (i = 0; i < nq; i++)
	phi_single[i].Relink (phi_tot, i*n, n);

    // distributed field calculation
    CalcFields_proc (qvec, phi_single);

    // Synchronise: collect fields from all processes
    int *ofs = new int[sze];
    int *len = new int[sze];
    for (i = 0; i < sze; i++) {
	ofs[i] = Q0[i]*n;
	len[i] = (Q1[i]-Q0[i])*n;
    }
    MPI_Allgatherv (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, phi_tot.data_buffer(),
		    len, ofs, mpitp, MPI_COMM_WORLD);

    for (i = 0; i < nq; i++)
	phi[i] = phi_single[i];

    delete []ofs;
    delete []len;
    delete []phi_single;

#endif
}

// =========================================================================

template<class T>
void TFwdSolver<T>::CalcFields (int q0, int q1, const TCompRowMatrix<T> &qvec,
    TVector<T> *phi, IterativeSolverResult *res) const
{
    // calculate fields for sources q0 <= q < q1
    static IterativeSolverResult s_res_single;
    IterativeSolverResult *res_single = 0;
    if (res) {
        res_single = &s_res_single;
	res->it_count = 0;
	res->rel_err = 0.0;
    }

    for (int i = q0; i < q1; i++) {
        CalcField (qvec.Row(i), phi[i], res_single);
	if (res) {
	    if (res_single->it_count > res->it_count)
	        res->it_count = res_single->it_count;
	    if (res_single->rel_err > res->rel_err)
	        res->rel_err = res_single->rel_err;
	}
    }
}

// =========================================================================

template<class T>
TVector<T> TFwdSolver<T>::ProjectSingle (int q, const TCompRowMatrix<T> &mvec,
    const TVector<T> &phi, DataScale scl) const
{
	xASSERT (meshptr, "Mesh reference not defined."); 

	if (scl == DATA_DEFAULT) scl = dscale;
    return ::ProjectSingle (meshptr, q, mvec, phi, scl);
}

// =========================================================================

template<class T>
TVector<T> TFwdSolver<T>::ProjectAll (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl)
{
	xASSERT (meshptr, "Mesh reference not defined."); 

	if (scl == DATA_DEFAULT) scl = dscale;

#ifndef MPI_FWDSOLVER

    return ::ProjectAll (meshptr, mvec, phi, scl);

#else

    // distributed calculation of projections
    TVector<T> proj (meshptr->nQM);
    ProjectAll_proc (mvec, phi, scl, proj);

    // synchronise the projection vector
    int i;
    int *ofs = new int[sze];
    int *len = new int[sze];
    for (i = 0; i < sze; i++) {
	ofs[i] = meshptr->Qofs[Q0[i]];
	len[i] = (Q1[i] == nQ ? meshptr->nQM : meshptr->Qofs[Q1[i]])-ofs[i];
    }
    MPI_Allgatherv (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, proj.data_buffer(),
		    len, ofs, mpitp, MPI_COMM_WORLD);

    delete []ofs;
    delete []len;
    return proj;

#endif // MPI_FWDSOLVER
}

// =========================================================================

template<class T>
TVector<T> TFwdSolver<T>::ProjectAll (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    const QMMesh *mesh = MeshPtr();
    int i, n = mesh->nlen(), nq = mesh->nQ;

    if (!pphi) { // allocate persistent workspace
	pphi = new TVector<T>[nq];
	for (i = 0; i < nq; i++) pphi[i].New(n);
    } else
	for (i = 0; i < nq; i++) pphi[i].Clear();

    Reset (sol, omega);
    // for MPI, this is a bottleneck - how to do in parallel?

#ifndef MPI_FWDSOLVER
    CalcFields (qvec, pphi);
#else
    CalcFields_proc (qvec, pphi);
#endif // MPI_FWDSOLVER

    return ProjectAll (mvec, pphi, scl);
}

// =========================================================================

template<>
STOASTLIB RVector TFwdSolver<double>::UnfoldComplex (const RVector &vec)
   const
{
    // nothing to do for real case
    return vec;
}

template<>
STOASTLIB RVector TFwdSolver<std::complex<double> >::UnfoldComplex (
   const CVector &vec) const
{
    int n = vec.Dim();
    RVector rvec(n*2);
    RVector rvec_r(rvec,0,n); rvec_r = Re(vec);
    RVector rvec_i(rvec,n,n); rvec_i = Im(vec);
    return rvec;
}

template<class T>
RVector TFwdSolver<T>::UnfoldComplex (const TVector<T> &vec)
   const
{
	xERROR("Not implemented");
	return RVector();
}

template<>
STOASTLIB FVector TFwdSolver<float>::UnfoldSComplex (const FVector &vec)
   const
{
    // nothing to do for real case
    return vec;
}

template<>
STOASTLIB FVector TFwdSolver<std::complex<float> >::UnfoldSComplex (
   const SCVector &vec) const
{
    int n = vec.Dim();
    FVector rvec(n*2);
    FVector rvec_r(rvec,0,n); rvec_r = Re(vec);
    FVector rvec_i(rvec,n,n); rvec_i = Im(vec);
    return rvec;
}

// =========================================================================

template<>
STOASTLIB RVector TFwdSolver<double>::ProjectAll_real (const RCompRowMatrix &mvec,
    const RVector *phi, DataScale scl)
{
    return ProjectAll (mvec, phi, scl);
}

template<>
STOASTLIB RVector TFwdSolver<std::complex<double> >::ProjectAll_real (
    const CCompRowMatrix &mvec, const CVector *phi, DataScale scl)
{
    if (!unwrap_phase || scl == DATA_LIN) { // nothing to unwrap
	return UnfoldComplex (ProjectAll (mvec, phi, scl));
    } else {
	RVector proj(meshptr->nQM*2);
#ifdef UNDEF
	//RVector mproj (proj, 0, meshptr->nQM);
	//RVector pproj (proj, meshptr->nQM, meshptr->nQM);
	int i, q, m, len, mofs = 0, pofs = meshptr->nQM, nq = meshptr->nQ;
	double scale;
	for (q = 0; q < nq; q++) {
	    len = meshptr->nQMref[q];
	    CVector lphi = log(phi[q]);
	    RVector lnamp = Re(lphi);
	    RVector phase = Im(lphi);

	    NimPhaseUnwrap (meshptr, phase, meshptr->Q[q]);

	    for (i = 0; i < len; i++) {
		m = meshptr->QMref[q][i];
		// log amplitude projection
		toast::complex c = log (dot (phi[q], mvec.Row(m)));
		proj[mofs++] = c.re;

		// phase projection
		RVector mv_r = Re(mvec.Row(m));
		// not sure what to do if mvec has imaginary component
		scale = sum(mv_r);
		proj[pofs++] = dot (phase, mv_r) / scale;

	    }
	}
#else
	int i;
	proj = UnfoldComplex (ProjectAll (mvec, phi, scl));
	RVector pproj (proj, meshptr->nQM, meshptr->nQM);
	for (i = 0; i < meshptr->nQM; i++)
	    if (pproj[i] > 0.0)
		pproj[i] -= Pi2;
#endif

	return proj;
    }
}

template<class T>
RVector TFwdSolver<T>::ProjectAll_real (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl)
{
    xERROR("Not implemented");
	return RVector();
}

// =========================================================================

template<>
STOASTLIB FVector TFwdSolver<float>::ProjectAll_singlereal (const FCompRowMatrix &mvec,
    const FVector *phi, DataScale scl)
{
    return ProjectAll (mvec, phi, scl);
}

template<>
STOASTLIB FVector TFwdSolver<std::complex<float> >::ProjectAll_singlereal (
    const SCCompRowMatrix &mvec, const SCVector *phi, DataScale scl)
{
    return UnfoldSComplex (ProjectAll (mvec, phi, scl));
}

template<class T>
FVector TFwdSolver<T>::ProjectAll_singlereal (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl)
{
    xERROR("Not implemented");
	return FVector();
}

// =========================================================================

template<>
STOASTLIB RVector TFwdSolver<double>::ProjectAll_real (const RCompRowMatrix &qvec,
    const RCompRowMatrix &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    return ProjectAll (qvec, mvec, sol, omega, scl);
}

template<>
STOASTLIB RVector TFwdSolver<std::complex<double> >::ProjectAll_real (
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const Solution &sol, double omega, DataScale scl)
{
    return UnfoldComplex (ProjectAll (qvec, mvec, sol, omega, scl));
}

template<class T>
RVector TFwdSolver<T>::ProjectAll_real (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    xERROR("Not implemented");
	return RVector();
}

// =========================================================================

template<>
STOASTLIB FVector TFwdSolver<float>::ProjectAll_singlereal (const FCompRowMatrix &qvec,
    const FCompRowMatrix &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    return ProjectAll (qvec, mvec, sol, omega, scl);
}

template<>
STOASTLIB FVector TFwdSolver<std::complex<float> >::ProjectAll_singlereal (
    const SCCompRowMatrix &qvec, const SCCompRowMatrix &mvec,
    const Solution &sol, double omega, DataScale scl)
{
    return UnfoldSComplex (ProjectAll (qvec, mvec, sol, omega, scl));
}

template<class T>
FVector TFwdSolver<T>::ProjectAll_singlereal (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    xERROR("Not implemented");
	return FVector();
}

// =========================================================================

template<class T>
void TFwdSolver<T>::ReadParams (ParamParser &pp)
{
    char cbuf[256];
    solvertp = LSOLVER_UNDEFINED;
    dscale = DATA_DEFAULT;

    if (pp.GetString ("DATASCALE", cbuf)) {
	if (!strcasecmp (cbuf, "LIN"))
	    dscale = DATA_LIN;
	else if (!strcasecmp (cbuf, "LOG"))
	    dscale = DATA_LOG;
    }
    while (dscale == DATA_DEFAULT) {
	int cmd;
	cout << "\nSelect data scaling:\n";
	cout << "(1) Linear\n";
	cout << "(2) Logarithmic\n";
	cout << "[1|2] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: dscale = DATA_LIN;
	        break;
	case 2: dscale = DATA_LOG;
	        break;
	}
    }

    if (pp.GetString ("LINSOLVER", cbuf)) {
        if (!strcasecmp (cbuf, "DIRECT"))
	    solvertp = LSOLVER_DIRECT;
	else if (!strcasecmp (cbuf, "CG"))
	    solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_CG;
	else if (!strcasecmp (cbuf, "BICGSTAB"))
  	    solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_BICGSTAB;
	else if (!strcasecmp (cbuf, "BICG"))
	    solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_BICG;
	else if (!strcasecmp (cbuf, "GMRES"))
	    solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_GMRES;
	else if (!strcasecmp (cbuf, "GAUSSSEIDEL"))
	    solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_GAUSSSEIDEL;
    }
    while (solvertp == LSOLVER_UNDEFINED) {
        int cmd;
        cout << "\nSelect linear solver:\n";
	cout << "(1) Direct\n";
	cout << "(2) CG (Conjugate gradient)\n";
	cout << "(3) BiCG (Bi-conjugate gradient)\n";
	cout << "(4) BiCG-STAB (Bi-conjugate gradient stabilised)\n";
	cout << "(5) GMRES (generalised minimum residual)\n";
	cout << "(6) GS (Gauss-Seidel)\n";
	cout << ">> ";
	cin >> cmd;
	switch (cmd) {
	case 1: solvertp = LSOLVER_DIRECT;
                break;
	case 2: solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_CG;
                break;
	case 3: solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_BICG;
	        break;
	case 4: solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_BICGSTAB;
	        break;
	case 5: solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_GMRES;
	        break;
	case 6: solvertp = LSOLVER_ITERATIVE, method = ITMETHOD_GAUSSSEIDEL;
	        break;
	}
    }

    if (solvertp == LSOLVER_ITERATIVE) {
        FGenericSparseMatrix::GlobalSelectIterativeSolver (method);
        RGenericSparseMatrix::GlobalSelectIterativeSolver (method);
	SCGenericSparseMatrix::GlobalSelectIterativeSolver_complex (method);
	CGenericSparseMatrix::GlobalSelectIterativeSolver_complex (method);

	if (!pp.GetReal ("LINSOLVER_TOL", iterative_tol)) {
	    do {
	        cout << "\nLinear solver stopping criterion:\n";
		cout << ">> ";
		cin >> iterative_tol;
	    } while (iterative_tol <= 0.0);
	}

	if (!pp.GetInt ("LINSOLVER_MAXIT", iterative_maxit)) {
	    cout << "\nLinear solver iteration limit (0=auto):\n>> ";
	    cin >> iterative_maxit;
	}

	precontp = (PreconType)-1; // undefined
	if (pp.GetString ("LINSOLVER_PRECON", cbuf)) {
	    if (!strcasecmp (cbuf, "NONE")) 
	        precontp = PRECON_NULL;
	    else if (!strcasecmp (cbuf, "DIAG"))
	        precontp = PRECON_DIAG;
#ifdef USE_CUDA_FLOAT
	    else if (!strcasecmp (cbuf, "CUSP_AINV"))
	        precontp = PRECON_CUSP_AINV;
	    else if (!strcasecmp (cbuf, "CUSP_SMOOTHED_AGGREGATION"))
	        precontp = PRECON_CUSP_SMOOTHED_AGGREGATION;
#else
	    else if (!strcasecmp (cbuf, "ICH"))
	        precontp = PRECON_ICH;
	    else if (!strcasecmp (cbuf, "DILU"))
	        precontp = PRECON_DILU;
#endif
	}
	while (precontp == (PreconType)-1) {
	    int cmd;
	    cout << "\nSelect preconditioner for linear solver:\n";
	    cout << "(0) None\n";
	    cout << "(1) Diagonal\n";
#ifdef USE_CUDA_FLOAT
	    cout << "(2) Approximate inverse\n";
	    cout << "(3) Smoothed aggregate\n";
#else
	    cout << "(2) Incomplete Choleski\n";
	    cout << "(3) Diagonal incomplete LU\n";
#endif
	    cin >> cmd;
	    switch (cmd) {
	    case 0: precontp = PRECON_NULL; break;
	    case 1: precontp = PRECON_DIAG; break;
#ifdef USE_CUDA_FLOAT
	    case 2: precontp = PRECON_CUSP_AINV; break;
	    case 3: precontp = PRECON_CUSP_SMOOTHED_AGGREGATION; break;
#else
	    case 2: precontp = PRECON_ICH; break;
	    case 3: precontp = PRECON_DILU; break;
#endif
	    }
	}
#ifdef USE_CUDA_FLOAT
	CuspPreconType cpt;
	switch (precontp) {
	case PRECON_DIAG: cpt = CUSP_PRECON_DIAGONAL; break;
	case PRECON_CUSP_AINV: cpt = CUSP_PRECON_AINV; break;
	case PRECON_CUSP_SMOOTHED_AGGREGATION:
	    cpt = CUSP_PRECON_SMOOTHED_AGGREGATION; break;
	default: cpt = CUSP_PRECON_IDENTITY; break;
	}
	cuda_SetCuspPrecon (cpt);
#endif
    }
}

template<class T>
void TFwdSolver<T>::WriteParams (ParamParser &pp)
{
    switch (dscale) {
    case DATA_LIN:
	pp.PutString ("DATASCALE", "LIN");
	break;
    case DATA_LOG:
	pp.PutString ("DATASCALE", "LOG");
	break;
    }

    if (solvertp == LSOLVER_DIRECT) {

        pp.PutString ("LINSOLVER", "DIRECT");

    } else if (solvertp == LSOLVER_ITERATIVE) {

	switch (method) {
	case ITMETHOD_CG:
	    pp.PutString ("LINSOLVER", "CG"); break;
	case ITMETHOD_BICG:
	    pp.PutString ("LINSOLVER", "BICG"); break;
	case ITMETHOD_BICGSTAB:
	    pp.PutString ("LINSOLVER", "BICGSTAB"); break;
	case ITMETHOD_GMRES:
	    pp.PutString ("LINSOLVER", "GMRES"); break;
	case ITMETHOD_GAUSSSEIDEL:
	    pp.PutString ("LINSOLVER", "GAUSSSEIDEL");break;
	}
	pp.PutReal ("LINSOLVER_TOL", iterative_tol);
	pp.PutInt ("LINSOLVER_MAXIT", iterative_maxit);
	switch (precontp) {
	case PRECON_NULL:
	    pp.PutString ("LINSOLVER_PRECON", "NONE"); break;
	case PRECON_DIAG:
	    pp.PutString ("LINSOLVER_PRECON", "DIAG"); break;
	case PRECON_ICH:
	    pp.PutString ("LINSOLVER_PRECON", "ICH"); break;
	case PRECON_DILU:
	    pp.PutString ("LINSOLVER_PRECON", "DILU"); break;
	case PRECON_CUSP_AINV:
	    pp.PutString ("LINSOLVER_PRECON", "CUSP_AINV"); break;
	case PRECON_CUSP_SMOOTHED_AGGREGATION:
	    pp.PutString ("LINSOLVER_PRECON", "CUSP_SMOOTHED_AGGREGATION");
	    break;
	}
    }
}

template<class T>
void TFwdSolver<T>::SetLinSolver (const char *solver, double tol)
{
    if (!strcasecmp (solver, "DIRECT")) {
	solvertp = LSOLVER_DIRECT;
    } else {
	solvertp = LSOLVER_ITERATIVE;
	iterative_tol = tol;
	if (!strcasecmp (solver, "CG"))
	    method = ITMETHOD_CG;
	else if (!strcasecmp (solver, "BICG"))
	    method = ITMETHOD_BICG;
	else if (!strcasecmp (solver, "BICGSTAB"))
	    method = ITMETHOD_BICGSTAB;
	else if (!strcasecmp (solver, "GMRES"))
	    method = ITMETHOD_GMRES;
	else
	    method = ITMETHOD_GAUSSSEIDEL;
	TGenericSparseMatrix<T>::GlobalSelectIterativeSolver (method);
	CGenericSparseMatrix::GlobalSelectIterativeSolver_complex (method);
    }
}


// =========================================================================
// =========================================================================
// MPI-specific methods

#ifdef MPI_FWDSOLVER

// =========================================================================
// Set MPI parameters

template<class T>
void TFwdSolver<T>::Setup_MPI ()
{
    MPI_Comm_rank (MPI_COMM_WORLD, &rnk);
    MPI_Comm_size (MPI_COMM_WORLD, &sze);
    mpitp = TMPI<T>::MPIType();
    nQ = 0;
}

// =========================================================================
// Deallocate MPI data structures.

template<class T>
void TFwdSolver<T>::Cleanup_MPI ()
{
    if (nQ) {
	delete []Q0;
	delete []Q1;
	nQ = 0;
    }
}

// =========================================================================
// Distribute sources over processes

template<class T>
void TFwdSolver<T>::DistributeSources_MPI (int nq) const
{
    if (nq != nQ) { // re-allocate
	if (nQ) {
	    delete []Q0;
	    delete []Q1;
	}
	if (nQ = nq) {
	    Q0 = new int[nQ];
	    Q1 = new int[nQ];
	}
    }

    if (nq) {
	int i;
	for (i = 0; i < sze; i++) {
	    Q0[i] = (i*nq)/sze;
	    Q1[i] = ((i+1)*nq)/sze;
	}
    }
}

// =========================================================================
// Nonblocking distributed calculation of fields:
// Does not perform synchronisation of results

template<class T>
void TFwdSolver<T>::CalcFields_proc (const TCompRowMatrix<T> &qvec,
    TVector<T> *phi) const
{
    int nq = qvec.nRows();
    if (nq != nQ) DistributeSources_MPI (nq);

    int q, q0, q1;
    
    q0 = Q0[rnk];
    q1 = Q1[rnk];

    for (q = q0; q < q1; q++)
	CalcField (qvec.Row(q), phi[q]);
}

// =========================================================================
// Nonblocking distributed calculation of projections:
// Does not perform synchronisation of results
// Each process only updates the part of the projection vector it is
// responsible for. Only the relevant fields are required by each process,
// so they can be provided by CalcFields_proc

template<class T>
void TFwdSolver<T>::ProjectAll_proc (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl, TVector<T> &proj)
{
	xASSERT (meshptr, "Mesh reference not defined."); 

	// Note: we don't really know the number of sources here, so we use
    // the value stored in the mesh
    int nq = meshptr->nQ;
    if (nq != nQ) DistributeSources_MPI (nq);

    int q, q0, q1;
    
    q0 = Q0[rnk];
    q1 = Q1[rnk];

    for (q = q0; q < q1; q++) {
	TVector<T> projq (proj, meshptr->Qofs[q], meshptr->nQMref[q]);
	projq = ProjectSingle (q, mvec, phi[q], scl);
    }
}

#endif // MPI_FWDSOLVER


// =========================================================================
// ========================================================================
// Nonmember functions

template<class T>
TVector<T> ProjectSingle (const QMMesh *mesh, int q,
    const TCompRowMatrix<T> &mvec, const TVector<T> &phi, DataScale dscale)
{
    // generate a projection from a field

    int i, m;
    TVector<T> proj(mesh->nQMref[q]);

    for (i = 0; i < mesh->nQMref[q]; i++) {
        m = mesh->QMref[q][i];
        proj[i] = dot (phi, mvec.Row(m));
    }
    return (dscale == DATA_LIN ? proj : log(proj));
}

// =========================================================================

template<class T>
TVector<T> ProjectAll (const QMMesh *mesh, const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale dscale)
{
    int i, len, ofs = 0;
    TVector<T> proj(mesh->nQM);

    for (i = 0; i < mesh->nQ; i++) {
	len = mesh->nQMref[i];
	TVector<T> proj_i (proj, ofs, len);
	proj_i = ProjectSingle (mesh, i, mvec, phi[i], dscale);
	ofs += len;
    }
    return proj;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

//template class STOASTLIB SuperLU_data<float>;
//template class STOASTLIB SuperLU_data<double>;
//template class STOASTLIB SuperLU_data<toast::complex>;
//template class STOASTLIB SuperLU_data<scomplex>;

template class STOASTLIB TFwdSolver<float>;
template class STOASTLIB TFwdSolver<double>;
template class STOASTLIB TFwdSolver<std::complex<double> >;
template class STOASTLIB TFwdSolver<std::complex<float> >;

template STOASTLIB FVector ProjectSingle (const QMMesh *mesh, int q,
    const FCompRowMatrix &mvec, const FVector &phi, DataScale dscale);
template STOASTLIB RVector ProjectSingle (const QMMesh *mesh, int q,
    const RCompRowMatrix &mvec, const RVector &phi, DataScale dscale);
template STOASTLIB CVector ProjectSingle (const QMMesh *mesh, int q,
    const CCompRowMatrix &mvec, const CVector &phi, DataScale dscale);
template STOASTLIB SCVector ProjectSingle (const QMMesh *mesh, int q,
    const SCCompRowMatrix &mvec, const SCVector &phi, DataScale dscale);

template STOASTLIB FVector ProjectAll (const QMMesh *mesh,
    const FCompRowMatrix &mvec, const FVector *phi, DataScale dscale);
template STOASTLIB RVector ProjectAll (const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *phi, DataScale dscale);
template STOASTLIB SCVector ProjectAll (const QMMesh *mesh,
    const SCCompRowMatrix &mvec, const SCVector *phi, DataScale dscale);
template STOASTLIB CVector ProjectAll (const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *phi, DataScale dscale);

#endif // NEED_EXPLICIT_INSTANTIATION
