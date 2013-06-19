#define FELIB_IMPLEMENTATION

#include "mathlib.h"
#include "felib.h"
#include "toast_mpi.h"
#include "mesh_mpi.h"

// ============================================================================
// Local prototypes

unsigned int simple_hash (unsigned int *key, unsigned int n);


// ============================================================================
// class MeshMPI

MeshMPI::MeshMPI (): QMMesh ()
{
    rank = TMPI::Rank();
    nproc = TMPI::Size();
    nndel = NULL;
    ndel = NULL;
    parts = NULL;
    lids = NULL;
    nodeType = NULL;
    nel = 0;
}

MeshMPI::~MeshMPI ()
{
    ClearNodeToElementMap ();
    if (parts) delete []parts;
    if (lids) delete []lids;
    if (nel) delete []ellist;
    if (nodeType) delete []nodeType;
}

void MeshMPI::Setup (bool mark_boundary)
{
    QMMesh::Setup (mark_boundary);
    ClearNodeToElementMap ();
    ComputeNodeToElementMap (&nndel, &ndel);
    ComputePartition (&myGraph);
    CreateZoltanDictionary ();
}

void MeshMPI::NodeToElementMap (int **_nndel, int ***_ndel) const
{
    *_nndel = nndel;
    *_ndel = ndel;
}

void MeshMPI::ComputeNodeToElementMap (int **_nndel, int ***_ndel) const
{
    int i, j, n;
    int *nndel = new int[nlen()];  // number of elements for each node
    int **ndel = new int*[nlen()]; // list of element indices for each node

    for (i = 0; i < nlen(); i++)
        nndel[i] = 0;
    for (i = 0; i < elen(); i++) {
        Element *pel = elist[i];
	for (j = 0; j < pel->nNode(); j++)
	    nndel[pel->Node[j]]++;
    }
    for (i = 0; i < nlen(); i++) {
        ndel[i] = new int[nndel[i]];
	nndel[i] = 0;
    }
    for (i = 0; i < elen(); i++) {
        Element *pel = elist[i];
	for (j = 0; j < pel->nNode(); j++) {
	    n = pel->Node[j];
	    ndel[n][nndel[n]++] = i;
	}
    }
    *_nndel = nndel;
    *_ndel = ndel;
}

void MeshMPI::ClearNodeToElementMap ()
{
    if (nndel) {
        delete []nndel;
	nndel = NULL;
    }
    if (ndel) {
        for (int i = 0; i < nlen(); i++)
	    delete []ndel[i];
	delete []ndel;
	ndel = NULL;
    }
}

// ============================================================================

int qcomp (const void *arg1, const void *arg2)
{
    int i1 = *(int*)arg1;
    int i2 = *(int*)arg2;
    return (i1 < i2 ? -1 : i1 > i2 ? 1 : 0);
}

void MeshMPI::ComputeElementList (ZOLTAN_ID_TYPE *gid, int ngid,
    int **_ellist, int *_nel) const
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

void MeshMPI::ComputeNodeTypeList (int **ntype) const
{
    int i, j, n = nlen();

    int *ntp = new int[n];
    for (i = 0; i < n; i++)
        ntp[i] = 0;
    for (i = 0; i < nel; i++) {
        Element *pel = elist[ellist[i]];
	for (j = 0; j < pel->nNode(); j++)
	  ntp[pel->Node[j]] = 1;
    }
    for (i = 0; i < myGraph.numMyVertices; i++)
        ntp[myGraph.vertexGID[i]-1] = 2;

    *ntype = ntp;
}

// ============================================================================
// Partition a mesh by splitting its graph into sub-graphs and distributing
// them over processes. The sub-graph data are contained in the GraphData
// struct for each process.

void MeshMPI::ComputePartition (GraphData *graph)
{
    int numGlobalVertices;
    int numGlobalNeighbours;
    int *nnbrs;
    int **nbrs;
    int *idx;
    int i, j, num, nnbors, procID;
    int send_count[2];
    unsigned int id;
    GraphData *send_graph;
    int ack = 0, ack_tag = 5, count_tag = 10, id_tag = 15;
    MPI_Status status;

    if (rank == 0) {
        numGlobalVertices = nlen();

	// Set up the neighbour graphs
	NodeNeighbourList (&nnbrs, &nbrs);
	for (i = numGlobalNeighbours = 0; i < numGlobalVertices; i++)
	    numGlobalNeighbours += nnbrs[i];

	// Allocate arrays for entire graph
	graph->vertexGID = new ZOLTAN_ID_TYPE[numGlobalVertices];
	graph->nborIndex = new int[numGlobalVertices+1];
	graph->nborGID   = new ZOLTAN_ID_TYPE[numGlobalNeighbours];
	graph->nborPart  = new int[numGlobalNeighbours];

	graph->vertex_capacity = numGlobalVertices;
	graph->nbor_capacity = numGlobalNeighbours;

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
	    graph->nborPart[i] = simple_hash(&id, nproc);
	}

	// Create a subgraph for each process
	send_graph = new GraphData[nproc];
	memset (send_graph, 0, sizeof(GraphData)*nproc);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    procID = simple_hash(&id, nproc);
	    send_graph[procID].numMyVertices++;
	}

	for (i = 0; i < nproc; i++) {
	    num = send_graph[i].numMyVertices;
	    send_graph[i].vertexGID = new ZOLTAN_ID_TYPE[num];
	    send_graph[i].nborIndex = new int[num+1];
	    memset (send_graph[i].nborIndex, 0, sizeof(int)*(num+1));
	}

	idx = new int[nproc];
	memset (idx, 0, sizeof(int)*nproc);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
	    procID = simple_hash(&id, nproc);
	    j = idx[procID];
	    send_graph[procID].vertexGID[j] = (ZOLTAN_ID_TYPE)id;
	    send_graph[procID].nborIndex[j+1] =
	        send_graph[procID].nborIndex[j] + nnbors;
	    idx[procID] = j+1;
	}

	for (i = 0; i < nproc; i++) {
	    num = send_graph[i].nborIndex[send_graph[i].numMyVertices];
	    send_graph[i].nborGID = new ZOLTAN_ID_TYPE[num];
	    send_graph[i].nborPart = new int[num];
	    send_graph[i].numAllNbors = num;
	}

	memset (idx, 0, sizeof(int)*nproc);

	for (i = 0; i < numGlobalVertices; i++) {
	    id = (unsigned int)graph->vertexGID[i];
	    nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
	    procID = simple_hash(&id, nproc);
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

	// MS: is the following correct? (was missing in original example)
	graph->vertex_capacity = graph->numMyVertices;
	graph->nbor_capacity = graph->nborIndex[graph->numMyVertices];

	// Send other subgraphs to their processes
	for (i = 1; i < nproc; i++) {
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
	for (i = 1; i < nproc; i++) {
	    MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
    } else { // rank > 0
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
	graph->vertex_capacity = send_count[0];
	graph->nbor_capacity = send_count[1];
    }
}

// ============================================================================

void MeshMPI::Partition (struct Zoltan_Struct *zz)
{
  int i, rc;

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
    Zoltan_Set_Num_Obj_Fn(zz, Zoltan_get_number_of_vertices, &myGraph);
    Zoltan_Set_Obj_List_Fn(zz, Zoltan_get_vertex_list, &myGraph);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, Zoltan_get_num_edges_list, &myGraph);
    Zoltan_Set_Edge_List_Multi_Fn(zz, Zoltan_get_edge_list, &myGraph);

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
    
    rc = Zoltan_DD_Update(myGraph.dd, myGraph.vertexGID, lids, NULL, parts,
        myGraph.numMyVertices);
}

// ============================================================================

void MeshMPI::Migrate (struct Zoltan_Struct *zz)
{
    int rc;
    int start_gid, num_nbors;

    // Register query functions for migration
    Zoltan_Set_Obj_Size_Multi_Fn(zz, Zoltan_get_message_sizes, &myGraph);
    Zoltan_Set_Pack_Obj_Multi_Fn(zz, Zoltan_pack_object_messages, &myGraph);
    Zoltan_Set_Unpack_Obj_Multi_Fn(zz, Zoltan_unpack_object_messages, &myGraph);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, Zoltan_mid_migrate, &myGraph);
    
    rc = Zoltan_Migrate (zz,
        numImport, importGlobalGids, importLocalGids,
	importProcs, importToPart,
	numExport, exportGlobalGids, exportLocalGids,
	exportProcs, exportToPart);

    // Use the data dictionary to find neighbours' partitions
    start_gid = myGraph.numMyVertices - numImport;
    num_nbors = myGraph.nborIndex[myGraph.numMyVertices] -
        myGraph.nborIndex[start_gid];

    rc = Zoltan_DD_Find(myGraph.dd,
        (ZOLTAN_ID_PTR)(myGraph.nborGID + start_gid), NULL, NULL, 
        myGraph.nborPart + start_gid, num_nbors, NULL);

    // Update the list of owned elements
    if (nel) delete []ellist;
    ComputeElementList (myGraph.vertexGID, myGraph.numMyVertices,
        &ellist, &nel);

    // Update node type list
    if (nodeType) delete []nodeType;
    ComputeNodeTypeList (&nodeType);
}

// ============================================================================

void MeshMPI::AddToElMatrix (int el, CCompRowMatrixMPI &M, RVector *coeff,
    int mode) const
{
    int i, j, is, js, nnode;
    toast::complex entry;

    nnode = elist[el]->nNode();
    for (i = 0; i < nnode; i++) {
	is = elist[el]->Node[i];
	if (nodeType[is] != 2) continue;
	for (j = 0; j < nnode; j++) {
	    js = elist[el]->Node[j];
	    switch (mode) {
	    case ASSEMBLE_FF:
		entry.re = elist[el]->IntFF (i, j);
		break;
	    case ASSEMBLE_DD:
		entry.re = elist[el]->IntDD (i, j);
		break;
	    case ASSEMBLE_PFF:
		entry.re = elist[el]->IntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PDD:
		entry.re = elist[el]->IntPDD (i, j, *coeff);
		break;
	    case ASSEMBLE_BNDPFF:
		entry.re = elist[el]->BndIntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PFF_EL:
		entry.re = elist[el]->IntFF (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_PDD_EL:
		entry.re = elist[el]->IntDD (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_BNDPFF_EL:
		entry.re = elist[el]->BndIntFF (i, j) * (*coeff)[el];
		break;

	    case ASSEMBLE_iPFF:
	        entry.im = elist[el]->IntPFF (i, j, *coeff);
		break;
	    }
	    M.Add (is, js, entry);
	}
    }

}

void MeshMPI::AddToSysMatrix (CCompRowMatrixMPI &M, RVector *coeff, int mode)
    const
{
    int i;
    for (i = 0; i < nel; i++) {
        int el = ellist[i];
	AddToElMatrix (el, M, coeff, mode);
    }
}

// ============================================================================

void MeshMPI::CreateZoltanDictionary ()
{
    int gid_length = 1;   // our global IDs consist of 1 integer
    int lid_length = 1;   // our local IDs consist of 1 integer
    int i, rc;
    struct Zoltan_DD_Struct *dd;

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
        parts[i] = rank;             // part number of this vertex
	lids[i] = (ZOLTAN_ID_TYPE)i; // local ID on my process for this vertex
    }

    rc = Zoltan_DD_Update(dd,
        myGraph.vertexGID,
	lids,
	NULL,
	parts,
	myGraph.numMyVertices);

    myGraph.dd = dd;
}

// ============================================================================
// Zoltan access functions

int MeshMPI::Zoltan_get_number_of_vertices (void *data, int *ierr)
{
    GraphData *graph = (GraphData*)data;
    *ierr = ZOLTAN_OK;
    return graph->numMyVertices;
}

void MeshMPI::Zoltan_get_vertex_list (void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
    float *obj_wgts, int *ierr)
{
    int i;

    GraphData *graph = (GraphData*)data;
    *ierr = ZOLTAN_OK;

    /* For now, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */
    for (i = 0; i < graph->numMyVertices; i++) {
        globalID[i] = graph->vertexGID[i];
	localID[i] = (ZOLTAN_ID_TYPE)i;
    }
}

void MeshMPI::Zoltan_get_num_edges_list (void *data, int sizeGID, int sizeLID,
    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int *numEdges, int *ierr)
{
    int i;
    ZOLTAN_ID_TYPE idx;

    GraphData *graph = (GraphData*)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyVertices)){
        *ierr = ZOLTAN_FATAL;
	return;
    }

    for (i = 0;  i < num_obj ; i++) {
        idx = localID[i];
	numEdges[i] = graph->nborIndex[idx+1] - graph->nborIndex[idx];
    }

    *ierr = ZOLTAN_OK;
    return;
}

void MeshMPI::Zoltan_get_edge_list(void *data, int sizeGID, int sizeLID,
    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int *num_edges,
    ZOLTAN_ID_PTR nborGID, int *nborPart,
    int wgt_dim, float *ewgts, int *ierr)
{
    int i, j, from, to;
    int *nextProc;
    ZOLTAN_ID_TYPE *nextNbor;

    GraphData *graph = (GraphData*)data;
    *ierr = ZOLTAN_OK;

    if ( (sizeGID != 1) || (sizeLID != 1) || 
	 (num_obj != graph->numMyVertices)||
	 (wgt_dim != 0)) {
        *ierr = ZOLTAN_FATAL;
	return;
    }

    nextNbor = nborGID;
    nextProc = nborPart;

    for (i = 0; i < num_obj; i++) {
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

void MeshMPI::Zoltan_get_message_sizes(void *data, int gidSize, int lidSize,
    int num_ids, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *sizes,
    int *ierr)
{
    GraphData *graph;
    int i, len;

    graph = (GraphData*)data;
    *ierr = ZOLTAN_OK;

    for (i = 0; i < num_ids; i++) {
        len = graph->nborIndex[localID[i]+1] - graph->nborIndex[localID[i]];
	sizes[i] = len * sizeof(ZOLTAN_ID_TYPE);
    }
}

void MeshMPI::Zoltan_pack_object_messages (void *data, int gidSize, int lidSize,
    int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, 
    int *dests, int *sizes, int *idx, char *buf, int *ierr)
{
    int i, j, num_nbors;
    ZOLTAN_ID_TYPE *nbors=NULL, *ibuf=NULL;
    ZOLTAN_ID_TYPE lid;
    GraphData *graph;
    *ierr = ZOLTAN_OK;
    graph = (GraphData*)data;

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

void MeshMPI::Zoltan_unpack_object_messages (void *data, int gidSize,
    int num_ids, ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf,
    int *ierr)
{
    int i, len, num_nbors, num_vertex, next_vertex, next_nbor;
    ZOLTAN_ID_TYPE *ibuf=NULL;
    GraphData *graph;
    *ierr = ZOLTAN_OK;
    graph = (GraphData*)data;

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

void MeshMPI::Zoltan_mid_migrate (void *data, int gidSize, int lidSize,
    int numImport, ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID,
    int *importProc, int *importPart, int numExport,
    ZOLTAN_ID_PTR exportGlobalID, ZOLTAN_ID_PTR exportLocalID, int *exportProc,
    int *exportPart, int *ierr)
{
    GraphData *graph; 
    int i, len, next_vertex, next_nbor;
    int *exports;

    *ierr = ZOLTAN_OK;
    graph = (GraphData*)data;

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
// Local functions

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
