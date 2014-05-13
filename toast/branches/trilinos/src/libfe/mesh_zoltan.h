// -*-C++-*-
// ==========================================================================
// Module libfe
// File mesh.h
// Declaration of class zoltanMesh (distributed mesh)
// Uses Zoltan library for mesh partitioning
// ==========================================================================

#ifndef __MESH_ZOLTAN_H
#define __MESH_ZOLTAN_H

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <iostream>
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

class FELIB zoltanMesh: public QMMesh {
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

    zoltanMesh ();
    zoltanMesh (const Teuchos::RCP<const Teuchos::Comm<int> >& comm);
    ~zoltanMesh ();
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

    void AddToElMatrix (int el, RCompRowMatrix &Mlocal,
        RVector *coeff, int mode) const
    {
	int i, j, is, js, isl, nnode;
	double entry;
	
	nnode = elist[el]->nNode();
	for (i = 0; i < nnode; i++) {
	    is = elist[el]->Node[i];
	    if (nodeType[is] != 2) continue; // node not owned
	    isl = nodeGlobalToLocal[is];
	    for (j = 0; j < nnode; j++) {
		js = elist[el]->Node[j];
		switch (mode) {
		case ASSEMBLE_FF:
		    entry = elist[el]->IntFF (i, j);
		    break;
		case ASSEMBLE_DD:
		    entry = elist[el]->IntDD (i, j);
		    break;
		case ASSEMBLE_PFF:
		    entry = elist[el]->IntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_PDD:
		    entry = elist[el]->IntPDD (i, j, *coeff);
		    break;
		case ASSEMBLE_BNDPFF:
		    entry = elist[el]->BndIntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_PFF_EL:
		    entry = elist[el]->IntFF (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_PDD_EL:
		    entry = elist[el]->IntDD (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_BNDPFF_EL:
		    entry = elist[el]->BndIntFF (i, j) * (*coeff)[el];
		    break;
		}
		Mlocal.Add (isl, js, entry);
	    }
	}
    }

    void AddToSysMatrix (RCompRowMatrix &Mlocal,
        RVector *coeff, int mode) const
    {
	int i, j, k, ncol, nnode, is, js, isl;
	double entry;

	for (i = 0; i < nel; i++) {
	    int el = ellist[i];
	    AddToElMatrix (el, Mlocal, coeff, mode);
	}
    }

    // temporary
    int *Parts() const { return parts; }
    ZOLTAN_ID_PTR Lids() const { return lids; }

protected:
    void ComputeNodeToElementMap (int **nndel, int ***ndel) const;
    void ClearNodeToElementMap ();
    void ComputeElementList (ZOLTAN_ID_TYPE *gid, int ngid,
	int **_ellist, int *_nel) const;
    void ComputeNodeTypeList (int **ntype) const;
    void ComputeNodeMaps (int **glob2loc) const;
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
    int *nodeGlobalToLocal; // maps from global to local node indices
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

#endif // !__MESH_ZOLTAN_H
