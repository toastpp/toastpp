#include "stoastlib.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "zoltan.h"

using namespace std;

// =========================================================================

struct GRAPH_DATA {
    int numMyVertices;         // total vertices in my partition
    int numAllNbors;           // total number of neighbours of my vertices
    ZOLTAN_ID_TYPE *vertexGID; // global ID for each of my vertices (1-based)
    int *nborIndex;            // nborIndex[i] is start of neighbours for vtx i
    ZOLTAN_ID_TYPE *nborGID;   // nborGID[nborIndex[i]] is first neighbour of i
    int *nborPart;             // process owning each nbor in nborGID
    int vertex_capacity;       // size of vertexGID buffer
    int nbor_capacity;         // size of nborGID & nborPart buffers
    struct Zoltan_DD_Struct *dd;
};


// =========================================================================
// local prototypes

// Application defined query functions for partitioning
static int get_number_of_vertices(void *data, int *ierr);
static void get_vertex_list(void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr);
static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int *numEdges, int *ierr);
static void get_edge_list(void *data, int sizeGID, int sizeLID,
    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int *num_edges, ZOLTAN_ID_PTR nborGID, int *nborPart,
    int wgt_dim, float *ewgts, int *ierr);

// Application defined query functions for migrating
static void get_message_sizes(void *data, int gidSize, int lidSize, int num_ids,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *sizes, int *ierr);
static void pack_object_messages(void *data, int gidSize, int lidSize,
    int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, int *dests,
    int *sizes, int *idx, char *buf, int *ierr);
static void mid_migrate(void *data, int gidSize, int lidSize,
    int numImport, ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID,
    int *importProc, int *importPart, int numExport,
    ZOLTAN_ID_PTR exportGlobalID, ZOLTAN_ID_PTR exportLocalID, int *exportProc,
    int *exportPart, int *ierr);
static void unpack_object_messages(void *data, int gidSize, int num_ids,
    ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf, int *ierr);

// Auxiliary functions
void PartitionMesh (int myRank, int numProcs, const QMMesh &mesh,
    GRAPH_DATA *graph);
void BuildNodeToElementMap (const QMMesh &mesh, int **nndel, int ***ndel);
void ComputeElementList (int myRank, int *nndel, int **ndel,
    ZOLTAN_ID_TYPE *gid, int ngid, int **_ellist, int *_nel);
void OutputPartition (int myRank, int numProcs, GRAPH_DATA &myGraph,
    const char *fname);
unsigned int simple_hash (unsigned int *key, unsigned int n);

void AddToElMatrix (const Mesh &mesh, int el, CGenericSparseMatrix &M,
		    const RVector *coeff, int mode, int *nodeType);

int rnk; // temporaray

// =========================================================================
// MAIN 

int main (int argc, char *argv[])
{
    //const char *fname = "PedroToastFormat_shifted.msh";

    const char *fname = "circle25_32.msh";
    const char *qmname = "circle25_32x32.qm";
    double qwidth = 2.0;
    SourceMode srctp = SRCMODE_NEUMANN;
    double tol = 1e-12;
    int maxit = 1000;

    int i, j, rc, nlen;
    int myRank, numProcs;
    struct Zoltan_Struct *zz;
    struct Zoltan_DD_Struct *dd;
    MPI_Status status;
    float version;
    GRAPH_DATA myGraph;
    int *parts = NULL;
    ZOLTAN_ID_PTR lids = NULL;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    int start_gid, num_nbors;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids,
      exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int gid_length = 1;   // our global IDs consist of 1 integer
    int lid_length = 1;   // our local IDs consist of 1 integer
    int *nndel, **ndel;

    // Initialise MPI and Zoltan
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    rnk = myRank; // temporary

    rc = Zoltan_Initialize(argc, argv, &version);
    if (rc != ZOLTAN_OK) {
        cerr << "Failed to initialise Zoltan." << endl;
	MPI_Finalize();
	exit(1);
    }
    
    // Load the mesh
    QMMesh mesh;
    ifstream ifs(fname);
    ifs >> mesh;
    mesh.Setup();
    ifstream qmf(qmname);
    mesh.LoadQM(qmf);
    nlen = mesh.nlen();
    BuildNodeToElementMap (mesh, &nndel, &ndel);

    // Define the RHS
    int nq = mesh.nQ;
    CCompRowMatrix qvec(nq,nlen);
    for (i = 0; i < nq; i++) {
        CVector q(nlen);
	SetReal (q, QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp));
	qvec.SetRow (i, q);
    }

    // Partition the mesh over processes
    PartitionMesh (myRank, numProcs, mesh, &myGraph);

    // Create a distributed data directory which maps vertex global IDs to
    // their current partition number
    rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD,
        gid_length,            // length of a global ID
	lid_length,            // length of a local ID
	0,                     // length of user data
	myGraph.numMyVertices, // hash table size
	0);                    // debug level

    parts = new int[myGraph.numMyVertices];
    lids  = new ZOLTAN_ID_TYPE[myGraph.numMyVertices];

    for (i = 0; i < myGraph.numMyVertices; i++) {
        parts[i] = myRank;           // part number of this vertex
	lids[i] = (ZOLTAN_ID_TYPE)i; // local ID on my process for this vertex
    }

    rc = Zoltan_DD_Update(dd,
        myGraph.vertexGID,
	lids,
	NULL,
	parts,
	myGraph.numMyVertices);

    myGraph.dd = dd;

    // Create the Zoltan environment    
    zz = Zoltan_Create(MPI_COMM_WORLD);

    // General parameters
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

    // Graph parameters
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2"); 
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

    // Register query functions for partitioning
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &myGraph);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &myGraph);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &myGraph);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &myGraph);

    // Register query functions for migration
    Zoltan_Set_Obj_Size_Multi_Fn(zz, get_message_sizes, &myGraph);
    Zoltan_Set_Pack_Obj_Multi_Fn(zz, pack_object_messages, &myGraph);
    Zoltan_Set_Unpack_Obj_Multi_Fn(zz, unpack_object_messages, &myGraph);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, mid_migrate, &myGraph);
    
    // Save mesh partition before Zoltan
    OutputPartition (myRank, numProcs, myGraph, "pre_partition.dat");

    // Perform partition
    rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output)
	&changes,        // 1 if partitioning was changed, 0 otherwise
	&numGidEntries,  // Number of integers used for a global ID
	&numLidEntries,  // Number of integers used for a local ID
	&numImport,      // Number of vertices to be sent to me
	&importGlobalGids,  // Global IDs of vertices to be sent to me
	&importLocalGids,   // Local IDs of vertices to be sent to me
	&importProcs,    // Process rank for source of each incoming vertex
	&importToPart,   // New partition for each incoming vertex
	&numExport,      // Number of vertices I must send to other processes
        &exportGlobalGids,  // Global IDs of the vertices I must send
	&exportLocalGids,   // Local IDs of the vertices I must send
	&exportProcs,    // Process to which I send each of the vertices
	&exportToPart);  // Partition to which each vertex will belong

    if (rc != ZOLTAN_OK) {
        printf("sorry...\n");
	MPI_Finalize();
	Zoltan_Destroy(&zz);
	exit(0);
    }

    // Update the data directory with the new partition numbers
    for (i = 0; i < numExport; i++)
        parts[exportLocalGids[i]] = exportToPart[i];
    
    rc = Zoltan_DD_Update(dd, myGraph.vertexGID, lids, NULL, parts,
        myGraph.numMyVertices);

    // Migrate vertices to new partitions
    rc = Zoltan_Migrate(zz,
        numImport, importGlobalGids, importLocalGids,
	importProcs, importToPart,
	numExport, exportGlobalGids, exportLocalGids,
	exportProcs, exportToPart);

    // Use the data dictionary to find neighbours' partitions
    start_gid = myGraph.numMyVertices - numImport;
    num_nbors = myGraph.nborIndex[myGraph.numMyVertices] -
        myGraph.nborIndex[start_gid];

    rc = Zoltan_DD_Find(dd,
        (ZOLTAN_ID_PTR)(myGraph.nborGID + start_gid), NULL, NULL, 
        myGraph.nborPart + start_gid, num_nbors, NULL);

    int nel, *ellist;
    ComputeElementList (myRank, nndel, ndel, myGraph.vertexGID,
        myGraph.numMyVertices, &ellist, &nel);

    // Save the partition after Zoltan
    OutputPartition (myRank, numProcs, myGraph, "post_partition.dat");

    // Mark node types (2=my node, 1=ghost node, 0=not my node)
    int *nodeType = new int[nlen];
    for (i = 0; i < nlen; i++)
        nodeType[i] = 0;
    for (i = 0; i < nel; i++) {
        Element *pel = mesh.elist[ellist[i]];
	for (j = 0; j < pel->nNode(); j++)
	  nodeType[pel->Node[j]] = 1;
    }
    for (i = 0; i < myGraph.numMyVertices; i++)
        nodeType[myGraph.vertexGID[i]-1] = 2;




   // Part-assemble by performing integrals over owned elements
    int *rowptr, *colidx, nzero;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    int *mynode = new int[myGraph.numMyVertices];
    for (i = 0; i < myGraph.numMyVertices; i++)
        mynode[i] = myGraph.vertexGID[i]-1;
     CCompRowMatrixMPI smat(nlen, nlen, rowptr, colidx,
	myGraph.numMyVertices, mynode);

    // assemble the partial system matrix
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
    for (i = 0; i < nel; i++) {
        int el = ellist[i];
	AddToElMatrix (mesh, el, smat, &cmua, ASSEMBLE_PFF, nodeType);
	AddToElMatrix (mesh, el, smat, &ckap, ASSEMBLE_PDD, nodeType);
	AddToElMatrix (mesh, el, smat, &c2a, ASSEMBLE_BNDPFF, nodeType);
	AddToElMatrix (mesh, el, smat, &omega, ASSEMBLE_iPFF, nodeType);
    }

    // Lets try a distributed matrix-vector product on the system matrix
    CVector x(nlen);
    for (i = 0; i < nlen; i++)
        x[i] = toast::complex(1,0);
    
    CVector *dphi = new CVector[nq];
    for (i = 0; i < nq; i++)
        dphi[i].New (nlen);

    CPrecon_Null precon;
    for (i = 0; i < nq; i++) {
        GMRES (smat, qvec.Row(i), dphi[i], tol, &precon, 20, maxit);
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
// Application defined query functions

static int get_number_of_vertices(void *data, int *ierr)
{
    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return graph->numMyVertices;
}


// ============================================================================

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
			    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			    int wgt_dim, float *obj_wgts, int *ierr)
{
    int i;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */
    for (i=0; i<graph->numMyVertices; i++) {
        globalID[i] = graph->vertexGID[i];
	localID[i] = (ZOLTAN_ID_TYPE)i;
    }
}


// ============================================================================

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
			       int num_obj,
			       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			       int *numEdges, int *ierr)
{
    int i;
    ZOLTAN_ID_TYPE idx;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyVertices)){
        *ierr = ZOLTAN_FATAL;
	return;
    }

    for (i=0;  i < num_obj ; i++) {
        idx = localID[i];
	numEdges[i] = graph->nborIndex[idx+1] - graph->nborIndex[idx];
    }

    *ierr = ZOLTAN_OK;
    return;
}


// ============================================================================

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborPart,
        int wgt_dim, float *ewgts, int *ierr)
{
    int i, j, from, to;
    int *nextProc;
    ZOLTAN_ID_TYPE *nextNbor;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    if ( (sizeGID != 1) || (sizeLID != 1) || 
	 (num_obj != graph->numMyVertices)||
	 (wgt_dim != 0)) {
        *ierr = ZOLTAN_FATAL;
	return;
    }

    nextNbor = nborGID;
    nextProc = nborPart;

    for (i=0; i < num_obj; i++) {
        /*
	 * In this example, we are not setting edge weights.  Zoltan will
	 * set each edge to weight 1.0.
	 */
        to = graph->nborIndex[localID[i]+1];
	from = graph->nborIndex[localID[i]];
	if ((to - from) != num_edges[i]) {
	    *ierr = ZOLTAN_FATAL;
	    return;
	}

	for (j=from; j < to; j++) {
	    *nextNbor++ = graph->nborGID[j];
	    *nextProc++ = graph->nborPart[j];
	}
    }
    return;
}


// ============================================================================

static void get_message_sizes(void *data, int gidSize, int lidSize, int num_ids,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *sizes, int *ierr)
{
    GRAPH_DATA *graph;
    int i, len;

    graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    for (i = 0; i < num_ids; i++) {
        len = graph->nborIndex[localID[i]+1] - graph->nborIndex[localID[i]];
	sizes[i] = len * sizeof(ZOLTAN_ID_TYPE);
    }
}


// ============================================================================

static void pack_object_messages(void *data, int gidSize, int lidSize,
    int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, 
    int *dests, int *sizes, int *idx, char *buf, int *ierr)
{
    int i, j, num_nbors;
    ZOLTAN_ID_TYPE *nbors=NULL, *ibuf=NULL;
    ZOLTAN_ID_TYPE lid;
    GRAPH_DATA *graph;
    *ierr = ZOLTAN_OK;
    graph = (GRAPH_DATA *)data;

    // For each exported vertex, write its neighbor global IDs to the
    //supplied buffer

    for (i = 0; i < num_ids; i++){
        lid = localIDs[i];
	nbors = graph->nborGID + graph->nborIndex[lid];
	num_nbors = graph->nborIndex[lid+1] - graph->nborIndex[lid];

	ibuf = (ZOLTAN_ID_TYPE *)(buf + idx[i]);

	for (j=0; j < num_nbors; j++){
	    ibuf[j] = *nbors++;
	}
    }
}

// ============================================================================

static void mid_migrate(void *data, int gidSize, int lidSize, int numImport,
    ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID, int *importProc,
    int *importPart, int numExport, ZOLTAN_ID_PTR exportGlobalID,
    ZOLTAN_ID_PTR exportLocalID, int *exportProc, int *exportPart, int *ierr)
{
    GRAPH_DATA *graph; 
    int i, len, next_vertex, next_nbor;
    int *exports;

    *ierr = ZOLTAN_OK;
    graph = (GRAPH_DATA *)data;

    // The exported vertices have been packed. 
    // Remove them from our local graph.

    exports = (int *)calloc(sizeof(int) , graph->numMyVertices);
    for (i = 0; i <numExport; i++) {
        exports[exportLocalID[i]] = 1; 
    }

    next_vertex = 0;
    next_nbor = 0;

    graph->nborIndex[0] = 0;

    for (i = 0; i < graph->numMyVertices; i++) {
        if (exports[i] == 0) {
	    len = graph->nborIndex[i+1] - graph->nborIndex[i];

	    if (i > next_vertex) {
	        graph->vertexGID[next_vertex] = graph->vertexGID[i];
		if (len > 0) {
		    memcpy(graph->nborGID + next_nbor,
                        graph->nborGID + graph->nborIndex[i],
                        sizeof(ZOLTAN_ID_TYPE) * len);
		    memcpy(graph->nborPart + next_nbor,
                        graph->nborPart + graph->nborIndex[i],
                        sizeof(int) * len);
		}
		graph->nborIndex[next_vertex+1] =
                    graph->nborIndex[next_vertex] + len;
	    }
	    next_nbor += len;
	    next_vertex++;
	}
    }

    free(exports);
    graph->numMyVertices = next_vertex;
}

// ============================================================================

static void unpack_object_messages(void *data, int gidSize, int num_ids,
    ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf, int *ierr)
{
    int i, len, num_nbors, num_vertex, next_vertex, next_nbor;
    ZOLTAN_ID_TYPE *ibuf=NULL;
    GRAPH_DATA *graph;
    *ierr = ZOLTAN_OK;
    graph = (GRAPH_DATA *)data;

    // Add incoming vertices to local graph
    next_vertex = graph->numMyVertices;
    next_nbor = graph->nborIndex[next_vertex];

    num_nbors = 0;
    for (i = 0; i < num_ids; i++) {
        num_nbors += size[i];
    }
    num_nbors /= sizeof(int);        /* number of incoming neighbors */

    num_nbors += next_nbor;          /* plus number of existing neighbors */

    num_vertex = next_vertex + num_ids;

    if (num_vertex > graph->vertex_capacity) {
        graph->vertexGID = (ZOLTAN_ID_TYPE *)realloc(graph->vertexGID,
            sizeof(ZOLTAN_ID_TYPE) * num_vertex);
	graph->nborIndex = (int *)realloc(graph->nborIndex,
            sizeof(int) * (num_vertex+1));
	graph->vertex_capacity = num_vertex; 
    }

    if (num_nbors > graph->nbor_capacity) {
        graph->nborGID = (ZOLTAN_ID_TYPE *)realloc(graph->nborGID,
            sizeof(ZOLTAN_ID_TYPE) * num_nbors);
	graph->nborPart = (int *)realloc(graph->nborPart,
            sizeof(int) * num_nbors);
	graph->nbor_capacity = num_nbors;
    }

    for (i = 0; i < num_ids; i++) {
        graph->vertexGID[next_vertex] = globalIDs[i];
	len = size[i] / sizeof(int);

	if (len > 0) {
	    ibuf = (ZOLTAN_ID_TYPE *)(buf + idx[i]);
	    memcpy(graph->nborGID + next_nbor, ibuf,
                len * sizeof(ZOLTAN_ID_TYPE));
	}
	graph->nborIndex[next_vertex+1] = graph->nborIndex[next_vertex] + len;
	next_vertex++;
	next_nbor += len;
    }

    graph->numMyVertices += num_ids;
}

// ============================================================================
// Partition a mesh by splitting its graph into sub-graphs and distributing
// them over processes. The sub-graph data are contained in the GRAPH_DATA
// struct for each process.

void PartitionMesh (int myRank, int numProcs, const QMMesh &mesh,
		    GRAPH_DATA *graph)
{
    int numGlobalVertices;
    int numGlobalNeighbours;
    int *nnbrs;
    int **nbrs;
    int *idx;
    int i, j, num, nnbors, procID;
    int send_count[2];
    unsigned int id;
    GRAPH_DATA *send_graph;
    int ack = 0, ack_tag = 5, count_tag = 10, id_tag = 15;
    MPI_Status status;

    if (myRank == 0) {
        numGlobalVertices = mesh.nlen();

	// Set up the neighbour graphs
	mesh.NodeNeighbourList (&nnbrs, &nbrs);
	for (i = numGlobalNeighbours = 0; i < numGlobalVertices; i++)
	    numGlobalNeighbours += nnbrs[i];

	// Allocate arrays for entire graph
	graph->vertexGID = new ZOLTAN_ID_TYPE[numGlobalVertices];
	graph->nborIndex = new int[numGlobalVertices+1];
	graph->nborGID   = new ZOLTAN_ID_TYPE[numGlobalNeighbours];
	graph->nborPart  = new int[numGlobalNeighbours];

	graph->nborIndex[0] = 0;

	for (i = 0; i < numGlobalVertices; i++) {
	  graph->vertexGID[i] = (ZOLTAN_ID_TYPE)(i+1);
	    for (j = 0; j < nnbrs[i]; j++) {
	      graph->nborGID[graph->nborIndex[i]+j] =
		  (ZOLTAN_ID_TYPE)(nbrs[i][j]+1);
	    }
	    graph->nborIndex[i+1] = graph->nborIndex[i] + nnbrs[i];
	}
	
	// Assign each vertex to a process using a hash function
	for (i = 0; i < numGlobalNeighbours; i++) {
	    id = (unsigned int)graph->nborGID[i];
	    graph->nborPart[i] = simple_hash(&id, numProcs);
	}

	// Create a subgraph for each process
	send_graph = new GRAPH_DATA[numProcs];
	memset (send_graph, 0, sizeof(GRAPH_DATA)*numProcs);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    procID = simple_hash(&id, numProcs);
	    send_graph[procID].numMyVertices++;
	}

	for (i = 0; i < numProcs; i++) {
	    num = send_graph[i].numMyVertices;
	    send_graph[i].vertexGID = new ZOLTAN_ID_TYPE[num];
	    send_graph[i].nborIndex = new int[num+1];
	    memset (send_graph[i].nborIndex, 0, sizeof(int)*(num+1));
	}

	idx = new int[numProcs];
	memset (idx, 0, sizeof(int)*numProcs);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
	    procID = simple_hash(&id, numProcs);
	    j = idx[procID];
	    send_graph[procID].vertexGID[j] = (ZOLTAN_ID_TYPE)id;
	    send_graph[procID].nborIndex[j+1] =
	        send_graph[procID].nborIndex[j] + nnbors;
	    idx[procID] = j+1;
	}

	for (i = 0; i < numProcs; i++) {
	    num = send_graph[i].nborIndex[send_graph[i].numMyVertices];
	    send_graph[i].nborGID = new ZOLTAN_ID_TYPE[num];
	    send_graph[i].nborPart = new int[num];
	    send_graph[i].numAllNbors = num;
	}

	memset (idx, 0, sizeof(int)*numProcs);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
	    procID = simple_hash(&id, numProcs);
	    j = idx[procID];
	    if (nnbors > 0) {
	        memcpy(send_graph[procID].nborGID+j,
		       graph->nborGID+graph->nborIndex[i],
		       nnbors*sizeof(ZOLTAN_ID_TYPE));
		memcpy(send_graph[procID].nborPart+j,
		       graph->nborPart+graph->nborIndex[i],
		       nnbors*sizeof(int));
		idx[procID] = j + nnbors;
	    }
	}
	delete idx;

	// Process zero subgraph
	delete graph->vertexGID;
	delete graph->nborIndex;
	delete graph->nborGID;
	delete graph->nborPart;

	*graph = send_graph[0];

	// Send other subgraphs to their processes
	for (i = 1; i < numProcs; i++) {
	    send_count[0] = send_graph[i].numMyVertices;
	    send_count[1] = send_graph[i].numAllNbors;

	    MPI_Send(send_count, 2, MPI_INT, i, count_tag, MPI_COMM_WORLD);
	    MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);
	    if (send_count[0] > 0) {
	        MPI_Send(send_graph[i].vertexGID, send_count[0],
			 ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
		delete send_graph[i].vertexGID;

		MPI_Send(send_graph[i].nborIndex, send_count[0]+1, MPI_INT, i,
			 id_tag+1, MPI_COMM_WORLD);
		delete send_graph[i].nborIndex;

		if (send_count[1] > 0) {
		    MPI_Send(send_graph[i].nborGID, send_count[1],
			     ZOLTAN_ID_MPI_TYPE, i, id_tag+2, MPI_COMM_WORLD);
		    delete send_graph[i].nborGID;

		    MPI_Send(send_graph[i].nborPart, send_count[1], MPI_INT, i,
			     id_tag+3, MPI_COMM_WORLD);
		    delete send_graph[i].nborPart;
		}
	    }
	}
	delete send_graph;

	// Signal all procs it is ok to go on
	ack = 0;
	for (i = 1; i < numProcs; i++) {
	    MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
    } else { // myRank > 0
        // pick up the graph for my process
        MPI_Recv(send_count, 2, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);
	if (send_count[0] < 0) {
	    MPI_Finalize();
	    exit(1);
	}

	ack = 0;
	graph->numMyVertices = send_count[0];
	graph->numAllNbors = send_count[1];

	if (send_count[0] > 0) {
	    graph->vertexGID = new ZOLTAN_ID_TYPE[send_count[0]];
	    graph->nborIndex = new int[send_count[0]+1];
	    if (send_count[1] > 0) {
	        graph->nborGID = new ZOLTAN_ID_TYPE[send_count[1]];
		graph->nborPart = new int[send_count[1]];
	    }
	}

	MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

	if (send_count[0] > 0) {
	    MPI_Recv(graph->vertexGID, send_count[0], ZOLTAN_ID_MPI_TYPE, 0,
		     id_tag, MPI_COMM_WORLD, &status);
	    MPI_Recv(graph->nborIndex, send_count[0]+1, MPI_INT, 0, id_tag+1,
		     MPI_COMM_WORLD, &status);

	    if (send_count[1] > 0) {
	        MPI_Recv(graph->nborGID, send_count[1], ZOLTAN_ID_MPI_TYPE, 0,
			 id_tag+2, MPI_COMM_WORLD, &status);
		MPI_Recv(graph->nborPart, send_count[1], MPI_INT, 0, id_tag+3,
			 MPI_COMM_WORLD, &status);
	    }
	}

	// ok to go on?
	MPI_Recv (&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (ack < 0) {
	    MPI_Finalize();
	    exit(1);
	}
    }
}

// ============================================================================

void update_element_list (int myRank, const QMMesh &mesh, GRAPH_DATA &myGraph)
{
    // update the list of elements required by this process from the list
    // of GIDs

  
}

// ============================================================================

void BuildNodeToElementMap (const QMMesh &mesh, int **_nndel, int ***_ndel)
{
    int i, j, n;
    int *nndel = new int[mesh.nlen()];  // number of elements for each node
    int **ndel = new int*[mesh.nlen()]; // list of element indices for each node

    for (i = 0; i < mesh.nlen(); i++)
        nndel[i] = 0;
    for (i = 0; i < mesh.elen(); i++) {
        Element *pel = mesh.elist[i];
	for (j = 0; j < pel->nNode(); j++)
	    nndel[pel->Node[j]]++;
    }
    for (i = 0; i < mesh.nlen(); i++) {
        ndel[i] = new int[nndel[i]];
	nndel[i] = 0;
    }
    for (i = 0; i < mesh.elen(); i++) {
        Element *pel = mesh.elist[i];
	for (j = 0; j < pel->nNode(); j++) {
	    n = pel->Node[j];
	    ndel[n][nndel[n]++] = i;
	}
    }
    *_nndel = nndel;
    *_ndel = ndel;
}

// ============================================================================

int qcomp (const void *arg1, const void *arg2)
{
    int i1 = *(int*)arg1;
    int i2 = *(int*)arg2;
    return (i1 < i2 ? -1 : i1 > i2 ? 1 : 0);
}

void ComputeElementList (int myRank, int *nndel, int **ndel,
    ZOLTAN_ID_TYPE *gid, int ngid, int **_ellist, int *_nel)
{
    int i, j;
    int nel = 0;
    for (i = 0; i < ngid; i++)
        nel += nndel[gid[i]-1];
    int *ellist = new int[nel];
    nel = 0;
    for (i = 0; i < ngid; i++) {
        for (j = 0; j < nndel[gid[i]-1]; j++)
	    ellist[nel++] = ndel[gid[i]-1][j];
    }

    // sort the list
    qsort (ellist, nel, sizeof(int), qcomp);

    // shrink to unique values
    for (j = 0, i = 1; i < nel; i++) {
        if (ellist[i] != ellist[j]) {
	    if (i > j+1) ellist[j+1] = ellist[i];
	    j++;
	}
    }
    nel = j+1;

    *_ellist = ellist;
    *_nel = nel;
}

// ============================================================================

void OutputPartition (int myRank, int numProcs, GRAPH_DATA &myGraph,
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

// ============================================================================

unsigned int simple_hash (unsigned int *key, unsigned int n)
{
    unsigned int h, rest, *p, bytes, num_bytes;
    char *byteptr;

    num_bytes = (unsigned int) sizeof(int);

    /* First hash the int-sized portions of the key */
    h = 0;
    for (p = (unsigned int *)key, bytes=num_bytes;
	 bytes >= (unsigned int) sizeof(int);
	 bytes-=sizeof(int), p++){
      h = (h*2654435761U) ^ (*p);
    }

    /* Then take care of the remaining bytes, if any */
    rest = 0;
    for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
      rest = (rest<<8) | (*byteptr);
    }

    /* Merge the two parts */
    if (rest)
      h = (h*2654435761U) ^ rest;

    /* Return h mod n */
    return (h%n);
}

void AddToElMatrix (const Mesh &mesh, int el, CGenericSparseMatrix &M,
		    const RVector *coeff, int mode, int *nodeType)
{
    int i, j, is, js, nnode;
    toast::complex entry;

    nnode = mesh.elist[el]->nNode();
    for (i = 0; i < nnode; i++) {
	is = mesh.elist[el]->Node[i];
	if (nodeType[is] != 2) continue;
	for (j = 0; j < nnode; j++) {
	    js = mesh.elist[el]->Node[j];
	    switch (mode) {
	    case ASSEMBLE_FF:
		entry.re = mesh.elist[el]->IntFF (i, j);
		break;
	    case ASSEMBLE_DD:
		entry.re = mesh.elist[el]->IntDD (i, j);
		break;
	    case ASSEMBLE_PFF:
		entry.re = mesh.elist[el]->IntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PDD:
		entry.re = mesh.elist[el]->IntPDD (i, j, *coeff);
		break;
	    case ASSEMBLE_BNDPFF:
		entry.re = mesh.elist[el]->BndIntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PFF_EL:
		entry.re = mesh.elist[el]->IntFF (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_PDD_EL:
		entry.re = mesh.elist[el]->IntDD (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_BNDPFF_EL:
		entry.re = mesh.elist[el]->BndIntFF (i, j) * (*coeff)[el];
		break;

	    case ASSEMBLE_iPFF:
	        entry.im = mesh.elist[el]->IntPFF (i, j, *coeff);
		break;
	    }
	    M.Add (is, js, entry);
	}
    }
}

