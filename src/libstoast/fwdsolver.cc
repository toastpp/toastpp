#define __FWDSOLVER_CC
#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#ifdef TOAST_MPI
#include "toast_mpi.h"
#endif

using namespace std;
using namespace toast;

// =========================================================================

SuperLU_data::SuperLU_data ()
{
    fact = 'N';
    refact = equed = 'N';
    perm_c = 0;
    perm_r = 0;
    etree  = 0;
    R      = 0;
    C      = 0;
    allocated = false;
    used_LU = false;
}

SuperLU_data::~SuperLU_data ()
{
    Deallocate();
    if (used_LU) {
	Destroy_SuperNode_Matrix (&L);
	Destroy_CompCol_Matrix (&U);
    }
}

void SuperLU_data::Deallocate ()
{
    if (allocated) {
	delete []perm_c;
	delete []perm_r;
	delete []etree;
	delete []R;
	delete []C;
	Destroy_SuperMatrix_Store (&smFW);
	allocated = false;
    }
}

void SuperLU_data::Setup (int n, CCompRowMatrix *CF)
{
    Deallocate();
    doublecomplex *cdat = (doublecomplex*)CF->ValPtr();
    zCreate_CompCol_Matrix (&smFW, CF->nRows(), CF->nCols(),
	    CF->nVal(), cdat, CF->colidx, CF->rowptr, NR, dtypeZ, GE);
    perm_c = new int[n];
    perm_r = new int[n];
    etree  = new int[n];
    R      = new double[n];
    C      = new double[n];
    get_perm_c (0, &smFW, perm_c);
    fact = refact = 'N';
    allocated = true;
}

void SuperLU_data::Solve (SuperMatrix *B, SuperMatrix *X)
{
    int info;
    mem_usage_t mem_usage;
    char equed = 'N';
    char trans = 'N';
    double R = 0.0;
    double C = 0.0;
    double ferr, berr;
    double recip_pivot_growth, rcond;

    zgssvx (&fact, &trans, &refact, &smFW, 0, perm_c, perm_r,
	    etree, &equed, &R, &C, &L, &U, 0, 0, B, X,
	    &recip_pivot_growth, &rcond, &ferr, &berr, &mem_usage, &info);

    fact = 'F';
    refact = 'Y';
    used_LU = true;
}

// =========================================================================

template<class T>
TFwdSolver<T>::TFwdSolver (LSOLVER linsolver, double tol)
{
    dscale   = DATA_LIN;
    solvertp = linsolver;
    iterative_tol = tol;
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;
    pphi = 0;

    SETUP_MPI();
}

template<class T>
TFwdSolver<T>::TFwdSolver (char *solver, double tol)
{
    dscale   = DATA_LIN;
    SetLinSolver (solver, tol);
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;
    pphi = 0;

    SETUP_MPI();
}

template<class T>
TFwdSolver<T>::TFwdSolver (ParamParser &pp)
{
    dscale   = DATA_LIN;
    ReadParams (pp);
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;
    pphi = 0;

    SETUP_MPI();
}

template<class T>
TFwdSolver<T>::~TFwdSolver ()
{
    if (F)      delete F;
    if (FL)     delete FL;
    if (Fd)     delete Fd;
    if (B)      delete B;
    if (precon) delete precon;
    if (pphi)   delete []pphi;

    CLEANUP_MPI();
}

template<class T>
void TFwdSolver<T>::SetDataScaling (DataScale scl)
{
    xASSERT(scl != DATA_DEFAULT, Invalid input argument);
    dscale = scl;
}

template<class T>
DataScale TFwdSolver<T>::GetDataScaling () const
{
    return dscale;
}

template<>
void TFwdSolver<double>::Allocate (const QMMesh &mesh)
{
    int *rowptr, *colidx, nzero;
    int n = mesh.nlen();
    int nq = mesh.nQ;

    meshptr = &mesh;

    // allocate system matrix
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
    F = new RCompRowMatrix (n, n, rowptr, colidx);
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
	precon = new RPrecon_Diag;
    }
}

template<>
void TFwdSolver<complex>::Allocate (const QMMesh &mesh)
{
    int *rowptr, *colidx, nzero;
    int n = mesh.nlen();
    int nq = mesh.nQ;

    meshptr = &mesh;

    // allocate system matrix
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    if (F) delete F;
    F = new CCompRowMatrix (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
	lu_data.Setup(meshptr->nlen(), F);
    } else {
	if (precon) delete precon;
	precon = new CPrecon_Diag;
    }
}

template<>
void TFwdSolver<double>::AssembleSystemMatrix (const Solution &sol,
    double omega, bool elbasis)
{
    // real version
    xASSERT(omega==0, Nonzero omega parameter not allowed here);

    RVector prm(meshptr->nlen());

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
void TFwdSolver<complex>::AssembleSystemMatrix (const Solution &sol,
    double omega, bool elbasis)
{
    // complex version
    RVector prm(meshptr->nlen());

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
void TFwdSolver<T>::AssembleMassMatrix (const Mesh *mesh)
{
    if (!mesh) mesh = meshptr;
    xASSERT(mesh, No mesh information available);

    if (!B) { // allocate on the fly
	int *rowptr, *colidx, nzero, n = mesh->nlen();
	mesh->SparseRowStructure (rowptr, colidx, nzero);
	B = new RCompRowMatrix (n, n, rowptr, colidx);
    } else {
	B->Zero();
    }

    AddToSysMatrix (*mesh, *B, 0, ASSEMBLE_FF);
}

template<>
void TFwdSolver<double>::Reset (const Solution &sol, double omega)
{
    // real version
    AssembleSystemMatrix (sol, omega);
    if (solvertp == LSOLVER_DIRECT)
	CholeskyFactorize (*F, *FL, *Fd, true);
    else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<complex>::Reset (const Solution &sol, double omega)
{
    // complex version
    AssembleSystemMatrix (sol, omega);
    if (solvertp == LSOLVER_DIRECT)
	lu_data.fact = 'N';
    else
	precon->Reset (F);
    if (B) AssembleMassMatrix();
}

template<>
void TFwdSolver<double>::CalcField (const TVector<double> &qvec,
    TVector<double> &phi) const
{
    // calculate the (real) photon density field for a given (real)
    // source distribution. Use only if data type is INTENSITY

    if (solvertp == LSOLVER_DIRECT) {
	CholeskySolve (*FL, *Fd, qvec, phi);
    } else {
	double tol = iterative_tol;
	IterativeSolve (*F, qvec, phi, tol, precon);
    }
}

template<>
void TFwdSolver<complex>::CalcField (const TVector<complex> &qvec,
    TVector<complex> &cphi) const
{
    // calculate the complex field for a given source distribution

#ifdef DO_PROFILE
    times (&tm);
    clock_t time0 = tm.tms_utime;
#endif
    if (solvertp == LSOLVER_DIRECT) {
        SuperMatrix B, X;
	int n = meshptr->nlen();

	doublecomplex *rhsbuf = (doublecomplex*)qvec.data_buffer();
	doublecomplex *xbuf   = (doublecomplex*)cphi.data_buffer();
	zCreate_Dense_Matrix (&B, n, 1, rhsbuf, n, DN, dtypeZ, GE);
	zCreate_Dense_Matrix (&X, n, 1, xbuf, n, DN, dtypeZ, GE);

	lu_data.Solve (&B, &X);

	Destroy_SuperMatrix_Store (&B);
	Destroy_SuperMatrix_Store (&X);
    } else {
        double tol = iterative_tol;
	IterativeSolve (*F, qvec, cphi, tol, precon);
    }
#ifdef DO_PROFILE
    times (&tm);
    solver_time += (double)(tm.tms_utime-time0)/(double)HZ;
#endif
    //cerr << "Solve time = " << solver_time << endl;
}

template<class T>
void TFwdSolver<T>::CalcFields (const TCompRowMatrix<T> &qvec,
    TVector<T> *phi) const
{
    // calculate the fields for all sources

#ifndef TOAST_MPI

    int i, nq = qvec.nRows();
    LOGOUT1_INIT_PROGRESSBAR ("CalcFields", 50, nq);
    for (i = 0; i < nq; i++) {
	CalcField (qvec.Row(i), phi[i]);
	LOGOUT1_PROGRESS(i);
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
void TFwdSolver<T>::CalcFields (const QMMesh &mesh, int nq, int q0, int q1,
    const TCompRowMatrix<T> &qvec, TVector<T> *phi) const
{
    // calculate fields for sources q0 <= q < q1

    for (int i = q0; i < q1; i++) {
	CalcField (qvec.Row(i), phi[i]);
    }
}

// =========================================================================

template<class T>
TVector<T> TFwdSolver<T>::ProjectSingle (int q, const TCompRowMatrix<T> &mvec,
    const TVector<T> &phi, DataScale scl) const
{
    if (scl == DATA_DEFAULT) scl = dscale;
    return ::ProjectSingle (meshptr, q, mvec, phi, scl);
}

// =========================================================================

template<class T>
TVector<T> TFwdSolver<T>::ProjectAll (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl)
{
    if (scl == DATA_DEFAULT) scl = dscale;

#ifndef TOAST_MPI

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

#endif // TOAST_MPI
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

#ifndef TOAST_MPI
    CalcFields (qvec, pphi);
#else
    CalcFields_proc (qvec, pphi);
#endif // TOAST_MPI

    return ProjectAll (mvec, pphi, scl);
}

// =========================================================================

template<class T>
RVector TFwdSolver<T>::ProjectAll_real (const TCompRowMatrix<T> &mvec,
    const TVector<T> *phi, DataScale scl)
{
    return UnfoldComplex (ProjectAll (mvec, phi, scl));
}

// =========================================================================

template<class T>
RVector TFwdSolver<T>::ProjectAll_real (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
    DataScale scl)
{
    return UnfoldComplex (ProjectAll (qvec, mvec, sol, omega, scl));
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
STOASTLIB RVector TFwdSolver<complex>::UnfoldComplex (const CVector &vec)
   const
{
    int n = vec.Dim();
    RVector rvec(n*2);
    RVector rvec_r(rvec,0,n); rvec_r = Re(vec);
    RVector rvec_i(rvec,n,n); rvec_i = Im(vec);
    return rvec;
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
	CGenericSparseMatrix::GlobalSelectIterativeSolver_complex (method);

	if (!pp.GetReal ("LINSOLVER_TOL", iterative_tol)) {
	    do {
	        cout << "\nLinear solver stopping criterion:\n";
		cout << ">> ";
		cin >> iterative_tol;
	    } while (iterative_tol <= 0.0);
	}
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
    }
}

template<class T>
void TFwdSolver<T>::SetLinSolver (char *solver, double tol)
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
	CGenericSparseMatrix::GlobalSelectIterativeSolver_complex (method);
    }
}


// =========================================================================
// =========================================================================
// MPI-specific methods

#ifdef TOAST_MPI

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

#endif // TOAST_MPI


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

// ========================================================================
// OBSOLETE

void Project_cplx (const QMMesh &mesh, int q, const CVector &phi,
    CVector &proj)
{
    int i, m, el, nv, dim, nd, in;
    complex dphi;
    double c2a;

    for (i = 0; i < mesh.nQMref[q]; i++) {
	m = mesh.QMref[q][i];
	el = mesh.Mel[m];
	nv = mesh.elist[el]->nNode();
	dim = mesh.elist[el]->Dimension();
	dphi = 0.0;
	
	RVector fun = mesh.elist[el]->GlobalShapeF (mesh.nlist, mesh.M[m]);
	for (in = 0; in < nv; in++) {
	    nd = mesh.elist[el]->Node[in];
	    c2a = mesh.plist[nd].C2A();   // reflection parameter
	    dphi += phi[nd]*fun[in]*c2a;
	}
	proj[i] = dphi;
    }
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class STOASTLIB TFwdSolver<double>;
template class STOASTLIB TFwdSolver<complex>;

template STOASTLIB RVector ProjectSingle (const QMMesh *mesh, int q,
    const RCompRowMatrix &mvec, const RVector &phi, DataScale dscale);
template STOASTLIB CVector ProjectSingle (const QMMesh *mesh, int q,
    const CCompRowMatrix &mvec, const CVector &phi, DataScale dscale);

template STOASTLIB RVector ProjectAll (const QMMesh *mesh,
    const RCompRowMatrix &mvec, const RVector *phi, DataScale dscale);
template STOASTLIB CVector ProjectAll (const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *phi, DataScale dscale);

#endif // NEED_EXPLICIT_INSTANTIATION
