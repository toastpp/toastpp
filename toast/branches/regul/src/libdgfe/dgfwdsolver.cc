#define __DGFWDSOLVER_CC
#include "dgfwdsolver.h"
#include "nonconformingMesh.h"
#ifdef TOAST_MPI
#include "toast_mpi.h"
#endif
using namespace std;

template<class T>
TDGFwdSolver<T>::TDGFwdSolver (LSOLVER linsolver, double tol)
{
    solvertp = linsolver;
    iterative_tol = tol;
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;

    SETUP_MPI();
}
template<>
TDGFwdSolver<std::complex<double> >::TDGFwdSolver (LSOLVER linsolver,
    double tol)
{
    solvertp = linsolver;
    iterative_tol = tol;
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;

    SETUP_MPI();
}


template<class T>
TDGFwdSolver<T>::TDGFwdSolver (char *solver, double tol)
{
    SetLinSolver (solver, tol);
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;

    SETUP_MPI();
}

template<class T>
TDGFwdSolver<T>::TDGFwdSolver (ParamParser &pp)
{
    ReadParams (pp);
    F  = 0;
    FL = 0;
    Fd = 0;
    B  = 0;
    precon = 0;
    meshptr = 0;

    SETUP_MPI();
}

template<class T>
TDGFwdSolver<T>::~TDGFwdSolver ()
{
    if (F)      delete F;
    if (FL)     delete FL;
    if (Fd)     delete Fd;
    if (B)      delete B;
    if (precon) delete precon;

    CLEANUP_MPI();
}
template<>
TDGFwdSolver<std::complex<double> >::~TDGFwdSolver ()
{
    if (F)      delete F;
    if (FL)     delete FL;
    if (Fd)     delete Fd;
    if (B)      delete B;
    if (precon) delete precon;

    CLEANUP_MPI();
}

// =========================================================================
template<> 
void TDGFwdSolver<double>::Allocate (NonconformingMesh &mesh)
{
    int *rowptr, *colidx, nzero;
    int n = mesh.nlen();
    //int nq = mesh.nQ;

    meshptr = &mesh;
    int dim = meshptr->elist[0]->nNode(); // Hack: assumes the entire mesh is composed of these type of elements

    //cout<<"inside Allocate"<<dim<<endl;
    // allocate system matrix
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    //cout<<"rowptr and colidx have been malloced"<<endl;
    if (F) delete F;
    F = new RCompRowMatrix (dim*n, dim*n, rowptr, colidx);
    //cout<<"memory allocated to the system matrix"<<endl;
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
	F->SymbolicCholeskyFactorize (rowptr, colidx);
	if (FL) delete FL;
	FL = new RCompRowMatrix (dim*n, dim*n, rowptr, colidx);
	if (Fd) delete Fd;
	Fd = new RVector (dim*n);
	delete []rowptr;
	delete []colidx;
    } else {
	if (precon) delete precon;
	precon = new RPrecon_Diag;
    }
}

template<> 
void TDGFwdSolver<std::complex<double> >::Allocate (NonconformingMesh &mesh)
{
    int *rowptr, *colidx, nzero;
   
    meshptr = &mesh;
    int dim = meshptr->elist[0]->nNode(); // Hack: assumes the entire mesh is composed of these type of elements

    // allocate system matrix
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    
    if (F) delete F;
    F = new CCompRowMatrix (dim*meshptr->elen(), dim*meshptr->elen(), rowptr, colidx);
    
    delete []rowptr;
    delete []colidx;

    // allocate factorisations and preconditioners
    if (solvertp == LSOLVER_DIRECT) {
#ifdef ENABLE_DIRECTSOLVER
	lu_data.Setup(dim*(meshptr->elen()), F);
#else
	xERROR("Direct solver not supported");
#endif
    } else {
	if (precon) delete precon;
	precon = new CPrecon_Diag;
    }
}

template<> 
void TDGFwdSolver<double>::AssembleSystemMatrix (const Solution &sol,
    double omega, bool elbasis)
{
    // real version
    xASSERT(omega==0, "Nonzero omega parameter not allowed here");

    RVector prm(meshptr->nlen());

    F->Zero();
    prm = sol.GetParam (OT_CMUA);
    DGAddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
    prm = sol.GetParam (OT_CKAPPA);
    DGAddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
    prm = sol.GetParam (OT_C2A);
    DGAddToSysMatrix (*meshptr, *F, &prm,
		    elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
}
template<> 
void TDGFwdSolver<std::complex<double> >::AssembleSystemMatrix (
    const Solution &sol,
    double omega, bool elbasis)
{
    dASSERT(!meshptr->is_set_up, Setup() needs to be called on the mesh object before performing this operation);	
    RVector prm1(meshptr->elen());
    RVector prm2(meshptr->elen());
    
    F->Zero();
    //cout << "Assembly started ..." << endl;
 
    prm1 = sol.GetParam (OT_CKAPPA);
    DGAddToSysMatrix (*meshptr, *F, &prm1, ASSEMBLE_PDD_EL);
    //cout << "Assembled PDD_EL" << endl;

    prm1 = sol.GetParam (OT_CMUA);
    DGAddToSysMatrix (*meshptr, *F, &prm1, ASSEMBLE_PFF_EL);
    //cout << "Assembled PFF_EL" << endl;
 
    DGAddToSysMatrix (*meshptr, *F, omega, ASSEMBLE_iCFF);
    //cout << "Assembled the element contributions" << endl;
   
    prm1 = sol.GetParam (OT_CKAPPA);
    prm2 = sol.GetParam(OT_N);
    DGAddToSysMatrixInteriorEdgeCont(*meshptr, *F, &prm1, &prm2);
    //cout << "Assembled the interior edge contributions" << endl;

    prm1 = sol.GetParam (OT_C2A);
    DGAddToSysMatrixBoundaryEdgeCont(*meshptr, *F, &prm1);
   //DGAddToSysMatrix(*meshptr, *F, &prm1, ASSEMBLE_BNDPFF_EL);	
    //cout << "Assembled the boundary edge contributions" << endl;

    
    }

template<class T>
void TDGFwdSolver<T>::ReadParams (ParamParser &pp)
{
    char cbuf[256];
    solvertp = LSOLVER_UNDEFINED;
    
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
template<>
void TDGFwdSolver<double>::CalcField (const TVector<double> &qvec, TVector<double> &phi) const
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
void TDGFwdSolver<std::complex<double> >::CalcField (
    const TVector<std::complex<double> > &qvec,
    TVector<std::complex<double> > &cphi) const
{
    // calculate the complex field for a given source distribution

#ifdef DO_PROFILE
    times (&tm);
    clock_t time0 = tm.tms_utime;
#endif
    if (solvertp == LSOLVER_DIRECT) {

#ifdef ENABLE_DIRECTSOLVER
	cout<<"entered DIRECT solver block ..." <<endl;
        SuperMatrix B, X;
	int nNode = meshptr->elist[0]->nNode();
	int n = nNode*meshptr->elen();

	doublecomplex *rhsbuf = (doublecomplex*)qvec.data_buffer();
	doublecomplex *xbuf   = (doublecomplex*)cphi.data_buffer();
	toast_zCreate_Dense_Matrix (&B, n, 1, rhsbuf, n, SLU_DN, SLU_Z,SLU_GE);
	toast_zCreate_Dense_Matrix (&X, n, 1, xbuf, n, SLU_DN, SLU_Z, SLU_GE);

	lu_data.Solve (&B, &X);

	toast_Destroy_SuperMatrix_Store (&B);
	toast_Destroy_SuperMatrix_Store (&X);
#else
	xERROR("Direct solver not supported");
#endif

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

template<>
void TDGFwdSolver<std::complex<double> >::SetLinSolver (char *solver,
    double tol)
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
#ifdef NEED_EXPLICIT_INSTANTIATION

template class STOASTLIB TFwdSolver<double>;
template class STOASTLIB TFwdSolver<toast::complex>;


#endif // NEED_EXPLICIT_INSTANTIATION


