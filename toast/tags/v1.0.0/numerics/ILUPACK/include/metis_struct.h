#ifndef _METIS_STRUCT_H_
#define _METIS_STRUCT_H_
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h,v 1.1 1998/11/27 17:59:31 karypis Exp $
 */

/* Undefine the following #define in order to use short integer as the idxtype */
#include "long_integer.h"
#include "metis_idxtype.h"

/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;


/*************************************************************************
* The following data structure will hold a node of a doubly-linked list.
**************************************************************************/
struct ListNodeType {
  integer id;                       	/* The id value of the node */
  struct ListNodeType *prev, *next;     /* It's a doubly-linked list */
};

typedef struct ListNodeType ListNodeType;



/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct PQueueType {
  integer type;                     /* The type of the representation used */
  integer nnodes;
  integer maxnodes;
  integer mustfree;

  /* Linear array version of the data structures */
  integer pgainspan, ngainspan;     /* plus and negative gain span */
  integer maxgain;
  ListNodeType *nodes;
  ListNodeType **buckets;

  /* Heap version of the data structure */
  KeyValueType *heap;
  idxtype *locator;
};

typedef struct PQueueType PQueueType;


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edegreedef {
  idxtype pid;
  idxtype ed;
};
typedef struct edegreedef EDegreeType;


/*************************************************************************
* The following data structure stores an edge for vol
**************************************************************************/
struct vedegreedef {
  idxtype pid;
  idxtype ed, ned;
  idxtype gv;
};
typedef struct vedegreedef VEDegreeType;


/*************************************************************************
* This data structure holds various working space data
**************************************************************************/
struct workspacedef {
  idxtype *core;			/* Where pairs, indices, and degrees are coming from */
  integer maxcore, ccore;

  EDegreeType *edegrees;
  VEDegreeType *vedegrees;
  integer cdegree;

  idxtype *auxcore;			/* This points to the memory of the edegrees */

  idxtype *pmat;			/* An array of k^2 used for eliminating domain 
                                           connectivity in k-way refinement */
};

typedef struct workspacedef WorkSpaceType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 integer id, ed;            	/* ID/ED of nodes */
 integer ndegrees;          	/* The number of different ext-degrees */
 EDegreeType *edegrees;     	/* List of edges */
};

typedef struct rinfodef RInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* vol-based partition
**************************************************************************/
struct vrinfodef {
 integer id, ed, nid;            	/* ID/ED of nodes */
 integer gv;            		/* IV/EV of nodes */
 integer ndegrees;          	/* The number of different ext-degrees */
 VEDegreeType *edegrees;     	/* List of edges */
};

typedef struct vrinfodef VRInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 idxtype edegrees[2];  
};

typedef struct nrinfodef NRInfoType;


/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  idxtype *gdata, *rdata;	/* Memory pools for graph and refinement data.
                                   This is where memory is allocated and used
                                   the rest of the fields in this structure */

  integer nvtxs, nedges;		/* The # of vertices and edges in the graph */
  idxtype *xadj;		/* Pointers to the locally stored vertices */
  idxtype *vwgt;		/* Vertex weights */
  idxtype *vsize;		/* Vertex sizes for min-volume formulation */
  idxtype *adjncy;		/* Array that stores the adjacency lists of nvtxs */
  idxtype *adjwgt;		/* Array that stores the weights of the adjacency lists */

  idxtype *adjwgtsum;		/* The sum of the adjacency weight of each vertex */

  idxtype *label;

  idxtype *cmap;

  /* Partition parameters */
  integer mincut, minvol;
  idxtype *where, *pwgts;
  integer nbnd;
  idxtype *bndptr, *bndind;

  /* Bisection refinement parameters */
  idxtype *id, *ed;

  /* K-way refinement parameters */
  RInfoType *rinfo;

  /* K-way volume refinement parameters */
  VRInfoType *vrinfo;

  /* Node refinement information */
  NRInfoType *nrinfo;


  /* Additional info needed by the MOC routines */
  integer ncon;			/* The # of constrains */ 
  float *nvwgt;			/* Normalized vertex weights */
  float *npwgts;		/* The normalized partition weights */

  struct graphdef *coarser, *finer;
};

typedef struct graphdef GraphType;



/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/*************************************************************************
* The following structure stores information used by Metis
**************************************************************************/
struct controldef {
  integer CoarsenTo;		/* The # of vertices in the coarsest graph */
  integer dbglvl;			/* Controls the debuging output of the program */
  integer CType;			/* The type of coarsening */
  integer IType;			/* The type of initial partitioning */
  integer RType;			/* The type of refinement */
  integer maxvwgt;			/* The maximum allowed weight for a vertex */
  float nmaxvwgt;		/* The maximum allowed weight for a vertex for each constrain */
  integer optype;			/* Type of operation */
  integer pfactor;			/* .1*prunning factor */
  integer nseps;			/* The number of separators to be found during multiple bisections */
  integer oflags;

  WorkSpaceType wspace;		/* Work Space Informations */

  /* Various Timers */
  timer TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, UncoarsenTmr, 
        SepTmr, RefTmr, ProjectTmr, SplitTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;

};

typedef struct controldef CtrlType;


/*************************************************************************
* The following data structure stores max-partition weight info for 
* Vertical MOC k-way refinement
**************************************************************************/
struct vpwgtdef {
  float max[2][MAXNCON];
  integer imax[2][MAXNCON];
};

typedef struct vpwgtdef VPInfoType;


#endif
