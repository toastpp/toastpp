#include "stoastlib.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "zoltan.h"
#include "timing.h"

using namespace std;

// =========================================================================
// local prototypes

void SelectMesh (ParamParser &pp, char *meshname, MeshMPI &mesh);

// Auxiliary functions
void OutputPartition (int myRank, int numProcs, MeshMPI::GraphData &myGraph,
    const char *fname);

// =========================================================================
// MAIN 

int main (int argc, char *argv[])
{
    ParamParser pp;
    if (argc > 1) pp.Open (argv[1]);

    double qwidth = 2.0;
    SourceMode srctp = SRCMODE_NEUMANN;
    double t, tol = 1e-12;
    int maxit = 1000;

    int i, j, rc, nlen;
    int myRank, numProcs;
    struct Zoltan_Struct *zz;
    MPI_Status status;
    float version;
    int *nndel, **ndel;

    // Initialise MPI and Zoltan
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    rc = Zoltan_Initialize(argc, argv, &version);
    if (rc != ZOLTAN_OK) {
        cerr << "Failed to initialise Zoltan." << endl;
	MPI_Finalize();
	exit(1);
    }
    
    // Create the Zoltan environment    
    zz = Zoltan_Create(MPI_COMM_WORLD);

    // Load the mesh
    char meshname[256];
    MeshMPI mesh;
    SelectMesh (pp, meshname, mesh);
    nlen = mesh.nlen();
    mesh.NodeToElementMap (&nndel, &ndel);

    // Define the RHS
    int nq = mesh.nQ;
    CCompRowMatrix qvec(nq,nlen);
    for (i = 0; i < nq; i++) {
        CVector q(nlen);
	SetReal (q, QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp));
	qvec.SetRow (i, q);
    }

    // Partition the mesh over processes
    MeshMPI::GraphData &myGraph = *mesh.GetGraphData();

    // Save mesh partition before Zoltan
    OutputPartition (myRank, numProcs, myGraph, "pre_partition.dat");

    // Perform the repartitioning
    mesh.Partition (zz);

    // Migrate vertices to new partitions
    mesh.Migrate (zz);

    // Save the partition after Zoltan
    OutputPartition (myRank, numProcs, myGraph, "post_partition.dat");

    // Part-assemble by performing integrals over owned elements
    int *rowptr, *colidx, nzero;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    int *mynode = new int[myGraph.numMyVertices];
    for (i = 0; i < myGraph.numMyVertices; i++)
        mynode[i] = myGraph.vertexGID[i]-1;
    CCompRowMatrixMPI smat(nlen, nlen, rowptr, colidx,
	myGraph.numMyVertices, mynode);

    // assemble the partial system matrix
    if (myRank == 0)
        tic();

    RVector cmua(nlen);
    RVector ckap(nlen);
    RVector c2a(nlen);
    RVector omega(nlen);
    double c0 = 0.3;
    double ref = 1.4;
    double freq = 100;
    double c = c0/ref;
    double mua = 0.025;
    double mus = 2.0;
    double kap = 1.0/(3.0*(mua+mus));
    cmua = c*mua;
    ckap = c*kap;
    c2a  = Parameter::C2A(REFLECTION_KEIJZER, ref);
    omega = freq * 2.0*Pi*1e-6;
    mesh.AddToSysMatrix (smat, &cmua, ASSEMBLE_PFF);
    mesh.AddToSysMatrix (smat, &ckap, ASSEMBLE_PDD);
    mesh.AddToSysMatrix (smat, &c2a, ASSEMBLE_BNDPFF);
    mesh.AddToSysMatrix (smat, &omega, ASSEMBLE_iPFF);

    if (myRank == 0) {
        t = toc();
	cerr << "assembly time: " << t << endl;
    }

    // Lets try a distributed matrix-vector product on the system matrix
    CVector x(nlen);
    for (i = 0; i < nlen; i++)
        x[i] = toast::complex(1,0);
    
    CVector *dphi = new CVector[nq];
    for (i = 0; i < nq; i++)
        dphi[i].New (nlen);

    if (myRank == 0)
        tic();

    CPrecon_Null precon;
    for (i = 0; i < nq; i++) {
        GMRES (smat, qvec.Row(i), dphi[i], tol, &precon, 20, maxit);
    }
    
    if (myRank == 0) {
        t = toc();
	cerr << "solver time: " << t << endl;
    }

    if (myRank == 0) {
        ofstream ofs ("dbg.dat");
	ofs << dphi[0] << endl;
    }

    Zoltan_Destroy (&zz);
    MPI_Finalize();

    return 0;
}                                                                              

// ============================================================================

void SelectMesh (ParamParser &pp, char *meshname, MeshMPI &mesh)
{
    char qmname[256];
    
    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nMesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup ();
    
    if (!pp.GetString ("QMFILE", qmname)) {
        cout << "\nQM file name:\n";
	cin >> qmname;
    }
    ifstream qmf (qmname);
    mesh.LoadQM (qmf);
}

// ============================================================================

void OutputPartition (int myRank, int numProcs, MeshMPI::GraphData &myGraph,
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

    if (myRank == 0) {
        numVertices = new int[numProcs];
	displace = new int[numProcs];
    }

    MPI_Gather (&myGraph.numMyVertices, 1, MPI_INT, numVertices, 1, MPI_INT,
		0, MPI_COMM_WORLD);

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

