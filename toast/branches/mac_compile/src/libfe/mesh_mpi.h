// -*-C++-*-
// ==========================================================================
// Module libfe
// File mesh.h
// Declaration of class MeshMPI (distributed mesh)
// Uses Zoltan library for mesh partitioning
// ==========================================================================

#ifndef __MESH_MPI_H
#define __MESH_MPI_H

#include "mesh.h"
#include "zoltan.h"

// ==========================================================================
// Prototypes

// ==========================================================================
// class Mesh
// ==========================================================================
/**
 * \brief Finite-element mesh management
 *
 * Defines node coordinates (via NodeList) and element connectivity (via
 * ElementList) and provides general mesh-related methods.
 * This version also defines node-based parameter coefficient sets for
 * optical tomography (via ParameterList), but this should be removed.
 */

class FELIB MeshMPI: public QMMesh {
public:
    struct GraphData {
      int numMyVertices;        // total vertices in my partition
      int numAllNbors;          // total number of neighbours of my vertices
      ZOLTAN_ID_TYPE *vertexGID;// global IDs for my vertices (1-based)
      int *nborIndex;           // nborIndex[i] is start of neighbours for vtx i
      ZOLTAN_ID_TYPE *nborGID;  // nborGID[nborIndex[i]] is first neighbour of i
      int *nborPart;            // process owning each nbor in nborGID
      int vertex_capacity;      // size of vertexGID buffer
      int nbor_capacity;        // size of nborGID & nborPart buffers
      struct Zoltan_DD_Struct *dd;
    };

    MeshMPI ();
    ~MeshMPI ();
    void Setup (bool mark_boundary=true);
    void NodeToElementMap (int **nndel, int ***ndel) const;
    GraphData *GetGraphData() { return &myGraph; }
    void Partition (struct Zoltan_Struct *zz);
    void Migrate (struct Zoltan_Struct *zz);
    Zoltan_DD_Struct *ZoltanDictionary () { return myGraph.dd; }
    int GetOwnedElements (int **_ellist) const
    { *_ellist = ellist; return nel; }
    int GetNodeTypeList (int **_ntp) const { *_ntp = nodeType; }

    void AddToElMatrix (int el, CCompRowMatrixMPI &M, RVector *coeff,
        int mode) const;
    void AddToSysMatrix (CCompRowMatrixMPI &M, RVector *coeff, int mode) const;

    // temporary
    int *Parts() const { return parts; }
    ZOLTAN_ID_PTR Lids() const { return lids; }

protected:
    void ComputeNodeToElementMap (int **nndel, int ***ndel) const;
    void ClearNodeToElementMap ();
    void ComputeElementList (ZOLTAN_ID_TYPE *gid, int ngid,
	int **_ellist, int *_nel) const;
    void ComputeNodeTypeList (int **ntype) const;
    void ComputePartition (GraphData *graph);
    void CreateZoltanDictionary ();

private:
    // Zoltan access functions
    static int Zoltan_get_number_of_vertices (void *data, int *ierr);
    static void Zoltan_get_vertex_list (void *data, int sizeGID, int sizeLID,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
        float *obj_wgts, int *ierr);
    static void Zoltan_get_num_edges_list (void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *numEdges, int *ierr);
    static void Zoltan_get_edge_list (void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges, ZOLTAN_ID_PTR nborGID, int *nborPart,
        int wgt_dim, float *ewgts, int *ierr);
    static void Zoltan_get_message_sizes (void *data, int gidSize, int lidSize,
        int num_ids, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *sizes,
	int *ierr);
    static void Zoltan_pack_object_messages (void *data, int gidSize,
        int lidSize, int num_ids, ZOLTAN_ID_PTR globalIDs,
        ZOLTAN_ID_PTR localIDs, int *dests, int *sizes, int *idx, char *buf,
        int *ierr);
    static void Zoltan_unpack_object_messages (void *data, int gidSize,
        int num_ids, ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf,
	int *ierr);
    static void Zoltan_mid_migrate (void *data, int gidSize, int lidSize,
        int numImport, ZOLTAN_ID_PTR importGlobalID,
        ZOLTAN_ID_PTR importLocalID, int *importProc, int *importPart,
        int numExport, ZOLTAN_ID_PTR exportGlobalID,
        ZOLTAN_ID_PTR exportLocalID, int *exportProc, int *exportPart,
	int *ierr);

    int rank;
    int nproc;
    int *nndel;        // number of elements for each node
    int **ndel;        // list of element indices for each node
    int nel;           // number of elements processed by this process
    int *ellist;       // list of elements processed by this process
    int *nodeType;     // 2=owned, 1=ghost, 0=not my node
    GraphData myGraph; // mesh partition data for this process
    int *parts;
    ZOLTAN_ID_PTR lids;

    // the following data are returned by the Zoltan partitioner.
    // Maybe store in local method if not required globally
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids,
      exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
};

#endif // !__MESH_MPI_H
