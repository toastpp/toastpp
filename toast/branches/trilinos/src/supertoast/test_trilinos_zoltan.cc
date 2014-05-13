#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <zoltan.h>
#include <iostream>
#include "stoastlib.h"
#include "mesh_zoltan.h"
#include "fwdsolver.h"
#include "source.h"
#include "timing.h"
#include "crmatrix_trilinos.h"
#include "test_trilinos_zoltan_class.h"

#define MAXREGION 100

// ============================================================================
// Local toast prototypes

void SelectMesh (ParamParser &pp, char *meshname, zoltanMesh &mesh);
void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol);
void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
bool ReadNim (char *nimname, RVector &img);
unsigned int simple_hash (unsigned int *key, unsigned int n);

void OutputPartition (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    int myRank, int numProcs, zoltanMesh::GraphData &myGraph,
    const char *fname);

// Solve the linear system(s) AX=B, using GMRES with the preconditioner M.
//
// B: right-hand side(s) of the linear system(s) AX=B.
// X: On input: initial guess(es) for solving AX=B.  On output: solution vector(s).
// A: the sparse matrix (or operator; it need not be a sparse matrix).
// M: if not null, the (right) preconditioner for A.
//
// In this example, MV is a specialization of Tpetra::MultiVector, 
// and OP is a specialization of Tpetra::Operator (the parent class 
// of Tpetra::CrsMatrix, Ifpack2::Preconditioner, and other classes 
// as well).


template<class MV, class OP>
void
solve (std::ostream& out, MV& X, const MV& B, const OP& A, Teuchos::RCP<OP> M) 
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP; 
  using Teuchos::rcp;
  using Teuchos::rcpFromRef; // Make a "weak" RCP from a reference.
  typedef typename MV::scalar_type scalar_type;

  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = parameterList();

  // Set some GMRES parameters.
  //
  // "Num Blocks" = Maximum number of Krylov vectors to store.  This
  // is also the restart length.  "Block" here refers to the ability
  // of this particular solver (and many other Belos solvers) to solve
  // multiple linear systems at a time, even though we may only be
  // solving one linear system in this example.
  //
  // "Maximum Iterations": Maximum total number of iterations,
  // including restarts.
  //
  // "Convergence Tolerance": By default, this is the relative
  // residual 2-norm, although you can change the meaning of the
  // convergence tolerance using other parameters.
  solverParams->set ("Num Blocks", 40);
  solverParams->set ("Maximum Iterations", 400);
  solverParams->set ("Convergence Tolerance", 1.0e-8);

  // Create the GMRES solver using a "factory" and 
  // the list of solver parameters created above.
  Belos::SolverFactory<scalar_type, MV, OP> factory;
  RCP<Belos::SolverManager<scalar_type, MV, OP> > solver = 
    factory.create ("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  typedef Belos::LinearProblem<scalar_type, MV, OP> problem_type;
  RCP<problem_type> problem = 
    rcp (new problem_type (rcpFromRef (A), rcpFromRef (X), rcpFromRef (B)));
  // You don't have to call this if you don't have a preconditioner.
  // If M is null, then Belos won't use a (right) preconditioner.
  //problem->setRightPrec (M);
  //////////////////////////problem->setLeftPrec (M);
  // Tell the LinearProblem to make itself ready to solve.
  problem->setProblem ();

  // Tell the solver what problem you want to solve.
  solver->setProblem (problem);

  // Attempt to solve the linear system.  result == Belos::Converged 
  // means that it was solved to the desired tolerance.  This call 
  // overwrites X with the computed approximate solution.
  Belos::ReturnType result = solver->solve();

  // Ask the solver how many iterations the last solve() took.
  const int numIters = solver->getNumIters();

  if (result == Belos::Converged) {
    out << "The Belos solve took " << numIters << " iteration(s) to reach "
      "a relative residual tolerance of " << 1.0e-8 << "." << std::endl;
  } else {
    out << "The Belos solve took " << numIters << " iteration(s), but did not reach "
      "a relative residual tolerance of " << 1.0e-8 << "." << std::endl;
  }
}

// Get the Ifpack2 preconditioner type and its parameter list.
// You may modify this function to change the preconditioner type
// and its parameters.
//
// The first output argument is the preconditioner name.  In this
// case, it's "ILUT", for Saad's ILUT incomplete factorization
// preconditioner.
//
// The second output argument is the parameter list for the
// preconditioner.  Give it to the preconditioner's setParameters()
// method.  The parameter list this function returns tells ILUT to
// use fill level 2, drop tolerance 0, and absolute threshold 0.1.
//
// Note that with Ifpack2, the type of preconditioner is separate
// from the ParameterList for that preconditioner.
void 
getPrecondTypeAndParameters (std::string& precondType, Teuchos::ParameterList& pl)
{
  using Teuchos::ParameterList;

  // The name of the type of preconditioner to use.
  precondType = "ILUT";

  // Ifpack2 expects arguments of type 'double' here, regardless of
  // the scalar or magnitude types of the entries of the sparse
  // matrix.
  const double fillLevel = 2.0; //2.0;
  const double dropTol = 0.0;
  const double absThreshold = 0.1;

  pl.set ("fact: ilut level-of-fill", fillLevel);
  pl.set ("fact: drop tolerance", dropTol);
  //pl.set ("fact: absolute threshold", absThreshold);
}

// This function encapsulates creation of an Ifpack2 preconditioner
// from a Tpetra::CrsMatrix.  It returns the preconditioner as a
// Tpetra::Operator, which is the parent class of
// Ifpack2::Preconditioner.
//
// The template parameter TpetraMatrixType must be a specialization of
// Tpetra::CrsMatrix.  We template this function on the matrix type,
// rather than on the five template arguments of Tpetra::CrsMatrix,
// because CrsMatrix has some nice typedefs that let us retrieve those
// template arguments.  It's easier to template on one thing than on
// five things!  Recall that if T is a template parameter, and if T is
// a class with a typedef inside T::U, then you have to use "typename"
// to get at the typedef U.
//
// You don't have to use this function to make an Ifpack2
// preconditioner.  We just find it easier to read the code example if
// we wrap up the preconditioner creation in its own little function.
template<class TpetraMatrixType>
Teuchos::RCP<Tpetra::Operator<typename TpetraMatrixType::scalar_type,
                              typename TpetraMatrixType::local_ordinal_type,
                              typename TpetraMatrixType::global_ordinal_type,
                              typename TpetraMatrixType::node_type> >
createPreconditioner (const Teuchos::RCP<const TpetraMatrixType>& A,
                      const std::string& precondType,
                      const Teuchos::ParameterList& plist,
                      std::ostream& out,
                      std::ostream& err)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using std::endl;

  // Fetch the typedefs defined by Tpetra::CrsMatrix.
  typedef typename TpetraMatrixType::scalar_type scalar_type;
  typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraMatrixType::node_type node_type;

  // Ifpack2's generic Preconditioner interface implements
  // Tpetra::Operator.  A Tpetra::Operator is an abstraction of a
  // function mapping a (Multi)Vector to a (Multi)Vector, with the
  // option of applying the transpose or conjugate transpose of the
  // operator.  Tpetra::CrsMatrix implements Operator as well.
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, 
                           global_ordinal_type, node_type> op_type;

  // These are just some convenience typedefs.
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
  // creates a Preconditioner object, but users of iterative methods
  // want a Tpetra::Operator.  That's why create() returns an Operator
  // instead of a Preconditioner.
  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, 
                                  global_ordinal_type, node_type> prec_type;

  // Create timers to show how long it takes for Ifpack2 to do various operations.
  RCP<Time> initTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::initialize");
  RCP<Time> computeTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::compute");
  RCP<Time> condestTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::condest");

  err << "Creating ILUT preconditioner" << endl 
      << "-- Configuring" << endl;
  //
  // Create the preconditioner and set parameters.
  //
  // This doesn't actually _compute_ the preconditioner.
  // It just sets up the specific type of preconditioner and
  // its associated parameters (which depend on the type).
  //
  RCP<prec_type> prec;
  Ifpack2::Factory factory;
  // Set up the preconditioner of the given type.
  prec = factory.create (precondType, A);
  prec->setParameters (plist);

  err << "-- Initializing" << endl;
  {
    TimeMonitor mon (*initTimer);
    prec->initialize();
  }

  // THIS ACTUALLY COMPUTES THE PRECONDITIONER
  // (e.g., does the incomplete factorization).
  err << "-- Computing" << endl;
  {
    TimeMonitor mon (*computeTimer);
    prec->compute();
  }

  if (precondType != "RELAXATION") {
    err << "-- Estimating condition number" << endl;
    magnitude_type condest = STM::one();
    {
      TimeMonitor mon (*condestTimer);
      condest = prec->computeCondEst (Ifpack2::Cheap);
    }
    out << endl << "Ifpack2 preconditioner's estimated condition number: " << condest << endl;
  }
  return prec;
}

template<class scalar_type>
void assemble_clbk (zoltanMesh *mesh, void *data,
     TCompRowMatrix<scalar_type> &localMatrix)
{
    Solution *sol = (Solution*)data;
    RVector prm;
    prm = sol->GetParam(OT_CMUA);
    mesh->AddToSysMatrix (localMatrix, &prm, ASSEMBLE_PFF);
    prm = sol->GetParam(OT_CKAPPA);
    mesh->AddToSysMatrix (localMatrix, &prm, ASSEMBLE_PDD);
    prm = sol->GetParam(OT_C2A);
    mesh->AddToSysMatrix (localMatrix, &prm, ASSEMBLE_BNDPFF);
}

template<class TpetraMatrixType>
TCompRowMatrixTrilinos<typename TpetraMatrixType::scalar_type> createToastMatrix (
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<typename TpetraMatrixType::node_type>& node,
    int m, int n, zoltanMesh *mesh, Solution *sol)
{
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using Teuchos::tuple;

    typedef TpetraMatrixType matrix_type;

    int i, j;
    int *rowptr, *colidx, nzero;
    int nlen = mesh->nlen();
    mesh->SparseRowStructure(rowptr, colidx, nzero);
    zoltanMesh::GraphData *graph = mesh->GetGraphData();

    // Fetch the timer for sparse matrix creation.
    RCP<Time> timer = TimeMonitor::lookupCounter ("Sparse matrix creation");
    if (timer.is_null())
	timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Time the whole scope of this routine, not counting timer lookup.
    TimeMonitor monitor (*timer);

    // Fetch typedefs from the Tpetra::CrsMatrix.
    typedef typename TpetraMatrixType::scalar_type scalar_type;
    typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename TpetraMatrixType::node_type node_type;

    // create an ArrayView from the graph data
    int *ndtp, nnd = 0;
    mesh->GetNodeTypeList (&ndtp);
    for (i = 0; i < nlen; i++)
	if (ndtp[i] == 2) nnd++;

    std::vector<int> myNode(nnd);
    for (i = j = 0; i < nlen; i++)
	     if (ndtp[i] == 2) myNode[j++] = i;
    ArrayView<int> myNodeView(myNode);

    // The type of the Tpetra::Map that describes how the matrix is distributed.
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    // The global number of rows in the matrix A to create.
    //const Tpetra::global_size_t numGlobalElements = (Tpetra::global_size_t)m;
    const Tpetra::global_size_t numGlobalElements = (const Tpetra::global_size_t)m;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    const global_ordinal_type indexBase = 0;
    RCP<const map_type > map = 
	rcp (new map_type (numGlobalElements, myNodeView, indexBase, comm,
			   node));


    TCompRowMatrixTrilinos<scalar_type> A(nlen,nlen,rowptr,colidx,map);
    A.Assemble (comm, node, mesh, (void*)sol, assemble_clbk<scalar_type>);

#ifdef UNDEF
    // Create a Tpetra::Matrix using the Map, with a static allocation
    // dictated by NumNz.
    RCP<matrix_type> A = rcp (new matrix_type (map, NumNz,
					       Tpetra::StaticProfile));
  

    // assemble local matrix into Tpetra matrix
    double *Alocal_val = Alocal.ValPtr();
    for (i = 0; i < numMyElements; i++) {
	int rp = myRowptr[i];
	ArrayView<int> cols (myColidx+rp, NumNz[i]);
	ArrayView<double> vals (Alocal_val+rp, NumNz[i]);
	A->insertGlobalValues (myGlobalElements[i], cols, vals);
    }

    // We are done with NumNZ; free it.
    NumNz = Teuchos::null;

    delete []rowptr;
    delete []colidx;
    delete []myRowptr;
    delete []myColidx;

    // Finish up the matrix.
    A->fillComplete ();
#endif
    return A;
}

#ifdef UNDEF
// Create a CrsMatrix from a mesh graph
template<class TpetraMatrixType>
Teuchos::RCP<const TpetraMatrixType>
createToastMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<typename TpetraMatrixType::node_type>& node,
    int m, int n, zoltanMesh *mesh, Solution *sol)
{
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using Teuchos::tuple;

    typedef TpetraMatrixType matrix_type;

    int i, j;
    int *rowptr, *colidx, nzero;
    mesh->SparseRowStructure(rowptr, colidx, nzero);
    zoltanMesh::GraphData *graph = mesh->GetGraphData();

    // Fetch the timer for sparse matrix creation.
    RCP<Time> timer = TimeMonitor::lookupCounter ("Sparse matrix creation");
    if (timer.is_null())
	timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Time the whole scope of this routine, not counting timer lookup.
    TimeMonitor monitor (*timer);

    // Fetch typedefs from the Tpetra::CrsMatrix.
    typedef typename TpetraMatrixType::scalar_type scalar_type;
    typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename TpetraMatrixType::node_type node_type;

    // create an ArrayView from the graph data
    int *ndtp, nnd = 0;
    mesh->GetNodeTypeList (&ndtp);
    for (i = 0; i < mesh->nlen(); i++)
	if (ndtp[i] == 2) nnd++;

    std::vector<int> myNode(nnd);
    for (i = j = 0; i < mesh->nlen(); i++)
	     if (ndtp[i] == 2) myNode[j++] = i;
    ArrayView<int> myNodeView(myNode);

    // The type of the Tpetra::Map that describes how the matrix is distributed.
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    // The global number of rows in the matrix A to create.
    //const Tpetra::global_size_t numGlobalElements = (Tpetra::global_size_t)m;
    const Tpetra::global_size_t numGlobalElements = (const Tpetra::global_size_t)m;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    const global_ordinal_type indexBase = 0;
    RCP<const map_type > map = 
	rcp (new map_type (numGlobalElements, myNodeView, indexBase, comm,
			   node));


    // Get update list and the number of equations that this MPI process
    // owns.
    const size_t numMyElements = map->getNodeNumElements();
    ArrayView<const global_ordinal_type> myGlobalElements =
	map->getNodeElementList();

    // NumNz[i] will be the number of nonzero entries for the i-th
    // global equation on this MPI process.
    ArrayRCP<size_t> NumNz = arcp<size_t> (numMyElements);
    int numNzero = 0;

    // Assign number of elements for each row we are owning
    for (size_t i = 0; i < numMyElements; ++i) {
	global_ordinal_type r = myGlobalElements[i];
	NumNz[i] = (size_t)(rowptr[r+1]-rowptr[r]);
	numNzero += NumNz[i];
    }

    // Create a local matrix for the assembly
    int *myRowptr = new int[numMyElements+1];
    int *myColidx = new int[numNzero];
    myRowptr[0] = 0;
    for (i = 0; i < numMyElements; i++) {
	myRowptr[i+1] = myRowptr[i] + NumNz[i];
	for (j = 0; j < NumNz[i]; j++)
	    myColidx[myRowptr[i]+j] = colidx[rowptr[myGlobalElements[i]]+j];
    }
    RCompRowMatrix Alocal(numMyElements, mesh->nlen(), myRowptr, myColidx);

    // Element-wise assembly
    RVector prm;
    prm = sol->GetParam(OT_CMUA);
    mesh->AddToSysMatrix (Alocal, &prm, ASSEMBLE_PFF);
    prm = sol->GetParam(OT_CKAPPA);
    mesh->AddToSysMatrix (Alocal, &prm, ASSEMBLE_PDD);
    prm = sol->GetParam(OT_C2A);
    mesh->AddToSysMatrix (Alocal, &prm, ASSEMBLE_BNDPFF);

    // Create a Tpetra::Matrix using the Map, with a static allocation
    // dictated by NumNz.
    RCP<matrix_type> A = rcp (new matrix_type (map, NumNz,
					       Tpetra::StaticProfile));
  

    // assemble local matrix into Tpetra matrix
    double *Alocal_val = Alocal.ValPtr();
    for (i = 0; i < numMyElements; i++) {
	int rp = myRowptr[i];
	ArrayView<int> cols (myColidx+rp, NumNz[i]);
	ArrayView<double> vals (Alocal_val+rp, NumNz[i]);
	A->insertGlobalValues (myGlobalElements[i], cols, vals);
    }

    // We are done with NumNZ; free it.
    NumNz = Teuchos::null;

    delete []rowptr;
    delete []colidx;
    delete []myRowptr;
    delete []myColidx;

    // Finish up the matrix.
    A->fillComplete ();
    return A;
}
#endif

// ============================================================================

int main (int argc, char *argv[]) 
{
    ParamParser pp;
    char meshname[256];

    // load the mesh
    pp.Open (argv[1]);
    pp.LogOpen ("fwdfem.out");

    using std::endl;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;

    Teuchos::oblackholestream blackHole;
    Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
    RCP<const Teuchos::Comm<int> > comm = 
	Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Set up Tpetra typedefs.
    typedef double scalar_type;
    typedef int local_ordinal_type;
    typedef int global_ordinal_type;
    //typedef long global_ordinal_type;
    typedef Kokkos::DefaultNode::DefaultNodeType node_type;

    // Create the Kokkos Node instance.  
    // Tpetra objects use this object for intranode parallel operations.
    RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode ();

    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    std::ostream& out = (myRank == 0) ? std::cout : blackHole;
    std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

    // Make a timer for sparse matrix creation.
    //
    // If you are using Trilinos 10.6 instead of the development branch
    // (10.7), just delete this line of code, and make the other change
    // mentioned above.
    RCP<Time> sparseMatrixCreationTimer = 
    	TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Run the whole example: create the sparse matrix, and compute the
    // preconditioner.
    // Print out the Tpetra software version information.
    out << Tpetra::version() << endl << endl;

    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
			      global_ordinal_type, node_type> matrix_type;
    typedef Tpetra::Operator<scalar_type, local_ordinal_type,
			     global_ordinal_type, node_type> op_type;
    typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
				global_ordinal_type, node_type> vec_type;


    float version;
    struct Zoltan_Struct *zz;
    Zoltan_Initialize(argc, argv, &version);
    zz = Zoltan_Create(MPI_COMM_WORLD);

    // Load the mesh
    zoltanMesh zmesh(comm);
    SelectMesh (pp, meshname, zmesh);
    int nq = zmesh.nQ;

    int *nndel, **ndel;
    zmesh.NodeToElementMap (&nndel, &ndel);

    // source and measurement operator parameters
    SourceMode srctp;
    int i, qprof, mprof;
    double qwidth, mwidth;
    SelectSourceProfile (pp, qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    RVector qvec = QVec_Gaussian (zmesh, zmesh.Q[0], qwidth, srctp);

    walltic();

    // Partition the mesh over processes
    zoltanMesh::GraphData &myGraph = *zmesh.GetGraphData();

    // Save mesh partition before Zoltan
#ifdef DBG_OUTPUT
    OutputPartition (comm, myRank, numProcs, myGraph, "pre_partition.dat");
#endif

    // Perform the repartitioning
    zmesh.Partition (zz);

    // Migrate vertices to new partitions
    zmesh.Migrate (zz);

    // Save the partition after Zoltan
#ifdef DBG_OUTPUT
    OutputPartition (comm, myRank, numProcs, myGraph, "post_partition.dat");
#endif

    // Load nodal parameters
    int nnd = zmesh.nlen();
    Solution sol (OT_NPARAM, nnd);
    SelectInitialParams (pp, zmesh, sol);

    TCompRowMatrixTrilinos<scalar_type> Atoast = createToastMatrix<matrix_type> (
        comm, node, nnd, nnd, &zmesh, &sol);

    RCP<const matrix_type> A = Atoast.GetMatrix();

    // Create and populate distributed system matrix
    //RCP<const matrix_type> A = createToastMatrix<matrix_type> (comm, node,
    //    nnd, nnd, &zmesh, &sol);

    double t_assemble = walltoc();
    walltic();

    // Get the preconditioner type and its parameter list.  Modify the
    // definition of that function if you want to use a different
    // preconditioner or parameters.
    std::string precondType;
    ParameterList plist;
    getPrecondTypeAndParameters (precondType, plist);

    // Compute the preconditioner using the matrix A.
    // The matrix A itself is not modified.
    RCP<op_type> M = createPreconditioner<matrix_type> (A, precondType, plist,
							out, err);

    // Create vectors ("multivectors" may store one or more vectors).
    RCP<vec_type> X = rcp (new vec_type (A->getDomainMap (), nq));
        // Set to zeros by default.

    // Build RHS vectors from source vectors
    RCP<vec_type> B = rcp (new vec_type (A->getRangeMap(), nq));
    ArrayView<const global_ordinal_type> mapB = 
	A->getRangeMap()->getNodeElementList();
    for (size_t q = 0; q < nq; q++) {
	RVector qvec = QVec_Gaussian (zmesh, zmesh.Q[q], qwidth, srctp);
	for (size_t i = 0; i < B->getLocalLength(); i++) {
	    size_t j = mapB[i];
	    B->replaceLocalValue(i,q,qvec[j]);
	}
    }

    // Solve the linear system using Belos.
    solve<vec_type, op_type> (out, *X, *B, *A, M);

    double t_solve = walltoc();

    // Summarize global performance timing results, for all timers
    // created using TimeMonitor::getNewCounter().
    TimeMonitor::summarize (out);

    // Retrieve result
    RVector xlocal(nnd), x(nnd);
    ArrayRCP<const double> xval = X->get1dView();
    const size_t localLength = X->getLocalLength();
    ArrayView<const global_ordinal_type> myGlobalElements =
	A->getDomainMap()->getNodeElementList();
    for (size_t i = 0; i < localLength; i++) {
	size_t j = myGlobalElements[i];
	xlocal[j] = xval[i];
    }
    Teuchos::reduceAll<int,double> (*comm, Teuchos::REDUCE_SUM,
    	    nnd, xlocal.data_buffer(), x.data_buffer());

    if (myRank == 0) {
	std::ofstream ofs1("dbg1.dat");
	ofs1 << x << std::endl;
	ofs1.close();
	std::cerr << "Wallclock timings: assemble = " << t_assemble << ", solve = " << t_solve << std::endl;
    }

    Zoltan_Destroy (&zz);
    return 0;
}



// ============================================================================

void SelectMesh (ParamParser &pp, char *meshname, zoltanMesh &mesh)
{
    char qmname[256];

    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nMesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();

    if (!pp.GetString ("QMFILE", qmname)) {
        cout << "\nQM file name:\n>> ";
	cin >> qmname;
    }
    ifstream qmf (qmname);
    mesh.LoadQM (qmf);

    // write back
    pp.PutString ("MESHFILE", meshname);
    pp.PutString ("QMFILE", qmname);
}

// ============================================================================

int ScanRegions (const Mesh &mesh, int *nregnode)
{
    int i, reg, nreg;
    for (i = 0; i < MAXREGION; i++) nregnode[i] = 0;
    for (i = 0; i < mesh.nlen(); i++) {
	reg = mesh.nlist[i].Region();
	if (reg >= 0 && reg < MAXREGION) nregnode[reg]++;
    }
    for (nreg = i = 0; i < MAXREGION; i++)
	if (nregnode[i]) nreg++;
    return nreg;
}

// ============================================================================

void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    RVector param[3];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];
    const char *resetstr[3] = {"RESET_MUA", "RESET_MUS", "RESET_N"};
    const ParameterType prmtp[3] = {PRM_MUA, PRM_MUS, PRM_N};
    for (p = 0; p < 3; p++) {

	param[p].New(mesh.nlen());
	if (pp.GetString (resetstr[p], cbuf)) {
	    pp.PutString (resetstr[p], cbuf);
	    if (!strcasecmp (cbuf, "MESH")) {
		param[p] = mesh.plist.Param(prmtp[p]);
	    } else if (!strncasecmp (cbuf, "HOMOG", 5)) {
		sscanf (cbuf+5, "%lf", &prm);
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    sscanf (valstr, "%lf", reg_prm+n);
		    valstr = strtok (NULL, " \t");
		}
		nreg = ScanRegions (mesh, nregnode);
		for (i = k = 0; k < n && i < MAXREGION; i++) {
		    if (nregnode[i]) {
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = reg_prm[k];
			k++;
		    }
		}	     
	    } else if (!strncasecmp (cbuf, "NIM", 3)) {
		ReadNim (cbuf+4, param[p]);
	    }
	} else {
	    cout << "\nSelect initial distribution for " << resetstr[p]
		 << endl;
	    cout << "(1) Use values stored in mesh\n";
	    cout << "(2) Global homogeneous\n";
	    cout << "(3) Homogeneous in regions\n";
	    cout << "(4) Nodal image file (NIM)\n";
	    cout << "[1|2|3|4] >> ";
	    cin >> resettp;
	    switch (resettp) {
	    case 1:
		param[p] = mesh.plist.Param(prmtp[p]);
		strcpy (cbuf, "MESH");
		break;
	    case 2:
		cout << "\nGlobal value:\n>> ";
		cin >> prm;
		param[p] = prm;
		sprintf (cbuf, "HOMOG %f", prm);
		break;
	    case 3:
		nreg = ScanRegions (mesh, nregnode);
		strcpy (cbuf, "REGION_HOMOG");
		cout << "\nFound " << nreg << " regions\n";
		for (i = 0; i < MAXREGION; i++) {
		    if (nregnode[i]) {
			cout << "Value for region " << i << " (" << nregnode[i]
			     << " nodes):\n>> ";
			cin >> prm;
			sprintf (cbuf+strlen(cbuf), " %f", prm);
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = prm;
		    }
		}
		break;
	    case 4:
		cout << "\nNIM file name:\n>> ";
		strcpy (cbuf, "NIM ");
		cin >> cbuf+4;
		ReadNim (cbuf+4, param[p]);
		break;
	    }
	    pp.PutString (resetstr[p], cbuf);
	}
    }
    msol.SetParam (OT_CMUA,   param[0]*c0/param[2]);
    msol.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    msol.SetParam (OT_N, param[2]);
    for (i = 0; i < param[OT_C2A].Dim(); i++)
	param[OT_C2A][i] = c0/(2*param[2][i]*A_Keijzer(param[OT_C2A][i]));
    msol.SetParam (OT_C2A, param[OT_C2A]);
}

// ============================================================================

void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp)
{
    char cbuf[256];
    int cmd;

    bool typeok = false;
    if (pp.GetString ("SOURCETYPE", cbuf)) {
	if (!strcasecmp (cbuf, "NEUMANN")) {
	    srctp = SRCMODE_NEUMANN;
	    typeok = true;
	} else if (!strcasecmp (cbuf, "ISOTROPIC")) {
	    srctp = SRCMODE_ISOTROPIC;
	    typeok = true;
	}
    }
    while (!typeok) {
	cout << "\nSource type:\n";
	cout << "(1) Neumann boundary source\n";
	cout << "(2) Isotropic point source\n";
	cout << "[1|2] >> ";
	cin  >> cmd;
	switch (cmd) {
	    case 1: srctp = SRCMODE_NEUMANN;   typeok = true; break;
	    case 2: srctp = SRCMODE_ISOTROPIC; typeok = true; break;
	}
    }
    pp.PutString ("SOURCETYPE",
        srctp == SRCMODE_NEUMANN ? "NEUMANN" : "ISOTROPIC");

    qtype = -1;
    if (pp.GetString ("SOURCEPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    qtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    qtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    qtype = 2;
	}
    }
    while (qtype < 0) {
	cout << "\nSource profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> qtype;
	qtype -= 1;
    }
    if (qtype > 0 && !pp.GetReal ("SOURCEWIDTH", qwidth)) {
	switch (qtype) {
	case 1:
	    cout << "\nSource 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nSource support radius [mm]:\n>> ";
	    break;
	}
	cin >> qwidth;
    }
    switch (qtype) {
    case 0:
	pp.PutString ("SOURCEPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("SOURCEPROFILE", "GAUSSIAN");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    case 2:
	pp.PutString ("SOURCEPROFILE", "COSINE");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    }
}

// ============================================================================

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth)
{
    char cbuf[256];
    mtype = -1;
    if (pp.GetString ("MEASUREMENTPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    mtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    mtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    mtype = 2;
	}
    }
    while (mtype < 0) {
	cout << "\nMeasurement profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> mtype;
	mtype -= 1;
    }
    if (mtype > 0 && !pp.GetReal ("MEASUREMENTWIDTH", mwidth)) {
	switch (mtype) {
	case 1:
	    cout << "\nMeasurement 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nMeasurement support radius [mm]:\n>> ";
	    break;
	}
	cin >> mwidth;
    }
    switch (mtype) {
    case 0:
	pp.PutString ("MEASUREMENTPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("MEASUREMENTPROFILE", "GAUSSIAN");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    case 2:
	pp.PutString ("MEASUREMENTPROFILE", "COSINE");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    }
}

// ============================================================================

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}

// ============================================================================

template<class TpetraMatrixType>
void AddToSysMatrix_Trilinos (const Mesh &mesh, TpetraMatrixType &M,
    const RVector *coeff, int mode)
{
    using Teuchos::RCP;
    using Teuchos::ArrayView;

    typedef TpetraMatrixType matrix_type;
    typedef typename TpetraMatrixType::node_type node_type;
    typedef Tpetra::Map<int, int, node_type> map_type;

    RCP<map_type> map = M->getRangeMap();
    const size_t numMyElements = map->getNodeNumElements();
    ArrayView<const int> myRows = map->getNodeElementList();

    for (size_t i = 0; i < numMyElements; i++) {
	int r = myRows[i];
    }
}


// ============================================================================

void OutputPartition (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    int myRank, int numProcs, zoltanMesh::GraphData &myGraph,
    const char *fname)
{
    int i, rc;
    int *graph_proc;
    ZOLTAN_ID_TYPE *graph_gid;
    int *numVertices;
    int *displace;
    int totVertices;

    int *parts = new int[myGraph.numMyVertices];
    rc = Zoltan_DD_Find (myGraph.dd, myGraph.vertexGID, NULL, NULL, parts,
			 myGraph.numMyVertices, NULL);

    numVertices = new int[numProcs];
    displace = new int[numProcs];
    Teuchos::gather<int,int>(&myGraph.numMyVertices, 1, numVertices, 1, 0,
			     *comm);
    //Teuchos::gatherAll<int,int> (*comm, 1, &myGraph.numMyVertices,
    //				 1, numVertices);

    //MPI_Gather (&myGraph.numMyVertices, 1, MPI_INT, numVertices, 1, MPI_INT,
    //		0, MPI_COMM_WORLD);

    if (myRank == 0) {
        for (i = totVertices = 0; i < numProcs; i++)
	    totVertices += numVertices[i];
	for (i = 1, displace[0] = 0; i < numProcs; i++)
	    displace[i] = displace[i-1] + numVertices[i-1];

	graph_proc = new int[totVertices];
	graph_gid  = new ZOLTAN_ID_TYPE[totVertices];
    }

    MPI_Gatherv (parts, myGraph.numMyVertices, MPI_INT, graph_proc, numVertices,
		 displace, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv (myGraph.vertexGID, myGraph.numMyVertices, ZOLTAN_ID_MPI_TYPE,
		 graph_gid, numVertices, displace, ZOLTAN_ID_MPI_TYPE, 0,
		 MPI_COMM_WORLD);

    if (myRank == 0) {
	ofstream ofs (fname);
	for (i = 0; i < totVertices; i++)
	    ofs << graph_gid[i] << ' ' << graph_proc[i] << endl;
	delete []numVertices;
	delete []displace;
	delete []graph_proc;
	delete []graph_gid;
    }

    delete []parts;
}

#ifdef UNDEF
template<class MT> class TCompRowMatrixTrilinosTest
{
public:
    typedef Kokkos::DefaultNode::DefaultNodeType node_type;
    typedef Tpetra::CrsMatrix<double, int, int, node_type> matrix_type;
    typedef Tpetra::Map<int, int, node_type> map_type;
  
    TCompRowMatrixTrilinosTest() {}

    void tmp (const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
	      const Teuchos::RCP<typename matrix_type::node_type>& node);

private:
    Teuchos::RCP<const map_type> map;
};

template<class MT>
void TCompRowMatrixTrilinosTest<MT>::tmp (
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<typename matrix_type::node_type>& node)
{
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::arcp;

    const Tpetra::global_size_t numGlobalElements = (const Tpetra::global_size_t)100;
    const size_t numMyElements = 10;
    const int indexBase = 0;

    std::vector<int> myNode(numGlobalElements);
    ArrayView<int> myNodeView(myNode);

    ArrayRCP<size_t> NumNz = arcp<size_t> (numMyElements);

    RCP<matrix_type> A = rcp (new matrix_type (map, NumNz,
					       Tpetra::StaticProfile));
}
#endif
