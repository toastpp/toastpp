#ifndef _METIS_PROTO_H_
#define _METIS_PROTO_H_
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.1 1998/11/27 17:59:28 karypis Exp $
 *
 */

#include "long_integer.h"





/* balance.c */
void Balance2Way(CtrlType *, GraphType *, integer *, float);
void Bnd2WayBalance(CtrlType *, GraphType *, integer *);
void General2WayBalance(CtrlType *, GraphType *, integer *);

/* bucketsort.c */
void BucketSortKeysInc(integer, integer, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void CreateCoarseGraph(CtrlType *, GraphType *, integer, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, integer, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, integer, idxtype *, idxtype *);
GraphType *SetUpCoarseGraph(GraphType *, integer, integer);
void ReAdjustMemory(GraphType *, GraphType *, integer);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* compress.c */
void CompressGraph(CtrlType *, GraphType *, integer, idxtype *, idxtype *, idxtype *, idxtype *);
void PruneGraph(CtrlType *, GraphType *, integer, idxtype *, idxtype *, idxtype *, float);

/* debug.c */
integer ComputeCut(GraphType *, idxtype *);
integer CheckBnd(GraphType *);
integer CheckBnd2(GraphType *);
integer CheckNodeBnd(GraphType *, integer);
integer CheckRInfo(RInfoType *);
integer CheckNodePartitionParams(GraphType *);
integer IsSeparable(GraphType *);

/* estmem.c */
void METIS_EstimateMemory(integer *, idxtype *, idxtype *, integer *, integer *, integer *);
void EstimateCFraction(integer, idxtype *, idxtype *, float *, float *);
integer ComputeCoarseGraphSize(integer, idxtype *, idxtype *, integer, idxtype *, idxtype *, idxtype *);

/* fm.c */
void FM_2WayEdgeRefine(CtrlType *, GraphType *, integer *, integer);

/* fortran.c */
void Change2CNumbering(integer, idxtype *, idxtype *);
void Change2FNumbering(integer, idxtype *, idxtype *, idxtype *);
void Change2FNumbering2(integer, idxtype *, idxtype *);
void Change2FNumberingOrder(integer, idxtype *, idxtype *, idxtype *, idxtype *);
void ChangeMesh2CNumbering(integer, idxtype *);
void ChangeMesh2FNumbering(integer, idxtype *, integer, idxtype *, idxtype *);
void ChangeMesh2FNumbering2(integer, idxtype *, integer, integer, idxtype *, idxtype *);

/* frename.c */
void METIS_PARTGRAPHRECURSIVE(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphrecursive(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphrecursive_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphrecursive__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphrecursive(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphrecursive_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphrecursive__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void METIS_PARTGRAPHKWAY(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphkway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphkway_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void metis_partgraphkway__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void METIS_WPARTGRAPHKWAY(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphkway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphkway_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void metis_wpartgraphkway__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void METIS_EDGEND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_edgend(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_edgend_(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_edgend__(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_NODEND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodend(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodend_(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodend__(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_NODEWND(integer *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodewnd(integer *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodewnd_(integer *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void metis_nodewnd__(integer *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshnodal(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshnodal_(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshnodal__(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshdual(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshdual_(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void metis_partmeshdual__(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void METIS_MESHTONODAL(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtonodal(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtonodal_(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtonodal__(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void METIS_MESHTODUAL(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtodual(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtodual_(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void metis_meshtodual__(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(integer *, idxtype *, idxtype *, integer *, integer *, integer *);
void metis_estimatememory(integer *, idxtype *, idxtype *, integer *, integer *, integer *);
void metis_estimatememory_(integer *, idxtype *, idxtype *, integer *, integer *, integer *);
void metis_estimatememory__(integer *, idxtype *, idxtype *, integer *, integer *, integer *);
void METIS_MCPARTGRAPHRECURSIVE(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_mcpartgraphrecursive(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_mcpartgraphrecursive_(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_mcpartgraphrecursive__(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void METIS_MCPARTGRAPHKWAY(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_mcpartgraphkway(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_mcpartgraphkway_(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_mcpartgraphkway__(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void METIS_PARTGRAPHVKWAY(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_partgraphvkway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_partgraphvkway_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void metis_partgraphvkway__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void METIS_WPARTGRAPHVKWAY(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_wpartgraphvkway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_wpartgraphvkway_(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void metis_wpartgraphvkway__(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);

/* graph.c */
void SetUpGraph(GraphType *, integer, integer, integer, idxtype *, idxtype *, idxtype *, idxtype *, integer);
void SetUpGraphKway(GraphType *, integer, idxtype *, idxtype *);
void SetUpGraph2(GraphType *, integer, integer, idxtype *, idxtype *, float *, idxtype *);
void VolSetUpGraph(GraphType *, integer, integer, integer, idxtype *, idxtype *, idxtype *, idxtype *, integer);
void RandomizeGraph(GraphType *);
integer IsConnectedSubdomain(CtrlType *, GraphType *, integer, integer);
integer IsConnected(CtrlType *, GraphType *, integer);
integer IsConnected2(GraphType *, integer);
integer FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);

/* initpart.c */
void Init2WayPartition(CtrlType *, GraphType *, integer *, float);
void InitSeparator(CtrlType *, GraphType *, float);
void GrowBisection(CtrlType *, GraphType *, integer *, float);
void GrowBisectionNode(CtrlType *, GraphType *, float);
void RandomBisection(CtrlType *, GraphType *, integer *, float);

/* kmetis.c */
void METIS_PartGraphKway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void METIS_WPartGraphKway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
integer MlevelKWayPartitioning(CtrlType *, GraphType *, integer, idxtype *, float *, float);

/* kvmetis.c */
void METIS_PartGraphVKway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void METIS_WPartGraphVKway(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
integer MlevelVolKWayPartitioning(CtrlType *, GraphType *, integer, idxtype *, float *, float);

/* kwayfm.c */
void Random_KWayEdgeRefine(CtrlType *, GraphType *, integer, float *, float, integer, integer);
void Greedy_KWayEdgeRefine(CtrlType *, GraphType *, integer, float *, float, integer);
void Greedy_KWayEdgeBalance(CtrlType *, GraphType *, integer, float *, float, integer);

/* kwayrefine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, integer, float *, float);
void AllocateKWayPartitionMemory(CtrlType *, GraphType *, integer);
void ComputeKWayPartitionParams(CtrlType *, GraphType *, integer);
void ProjectKWayPartition(CtrlType *, GraphType *, integer);
integer IsBalanced(idxtype *, integer, float *, float);
void ComputeKWayBoundary(CtrlType *, GraphType *, integer);
void ComputeKWayBalanceBoundary(CtrlType *, GraphType *, integer);

/* kwayvolfm.c */
void Random_KWayVolRefine(CtrlType *, GraphType *, integer, float *, float, integer, integer);
void Random_KWayVolRefineMConn(CtrlType *, GraphType *, integer, float *, float, integer, integer);
void Greedy_KWayVolBalance(CtrlType *, GraphType *, integer, float *, float, integer);
void Greedy_KWayVolBalanceMConn(CtrlType *, GraphType *, integer, float *, float, integer);
void KWayVolUpdate(CtrlType *, GraphType *, integer, integer, integer, idxtype *, idxtype *, idxtype *);
void ComputeKWayVolume(GraphType *, integer, idxtype *, idxtype *, idxtype *);
integer ComputeVolume(GraphType *, idxtype *);
void CheckVolKWayPartitionParams(CtrlType *, GraphType *, integer);
void ComputeVolSubDomainGraph(GraphType *, integer, idxtype *, idxtype *);
void EliminateVolSubDomainEdges(CtrlType *, GraphType *, integer, float *);
void EliminateVolComponents(CtrlType *, GraphType *, integer, float *, float);

/* kwayvolrefine.c */
void RefineVolKWay(CtrlType *, GraphType *, GraphType *, integer, float *, float);
void AllocateVolKWayPartitionMemory(CtrlType *, GraphType *, integer);
void ComputeVolKWayPartitionParams(CtrlType *, GraphType *, integer);
void ComputeKWayVolGains(CtrlType *, GraphType *, integer);
void ProjectVolKWayPartition(CtrlType *, GraphType *, integer);
void ComputeVolKWayBoundary(CtrlType *, GraphType *, integer);
void ComputeVolKWayBalanceBoundary(CtrlType *, GraphType *, integer);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);

/* mbalance.c */
void MocBalance2Way(CtrlType *, GraphType *, float *, float);
void MocGeneral2WayBalance(CtrlType *, GraphType *, float *, float);

/* mbalance2.c */
void MocBalance2Way2(CtrlType *, GraphType *, float *, float *);
void MocGeneral2WayBalance2(CtrlType *, GraphType *, float *, float *);
void SelectQueue3(integer, float *, float *, integer *, integer *, PQueueType [MAXNCON][2], float *);

/* mcoarsen.c */
GraphType *MCCoarsen2Way(CtrlType *, GraphType *);

/* memory.c */
void AllocateWorkSpace(CtrlType *, GraphType *, integer);
void FreeWorkSpace(CtrlType *, GraphType *);
integer WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, integer);
void idxwspacefree(CtrlType *, integer);
float *fwspacemalloc(CtrlType *, integer);
void fwspacefree(CtrlType *, integer);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* mesh.c */
void METIS_MeshToDual(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void METIS_MeshToNodal(integer *, integer *, idxtype *, integer *, integer *, idxtype *, idxtype *);
void GENDUALMETIS(integer, integer, integer, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(integer, integer, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(integer, integer, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(integer, integer, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(integer, integer, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */
void METIS_PartMeshNodal(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);
void METIS_PartMeshDual(integer *, integer *, idxtype *, integer *, integer *, integer *, integer *, idxtype *, idxtype *);

/* mfm.c */
void MocFM_2WayEdgeRefine(CtrlType *, GraphType *, float *, integer);
void SelectQueue(integer, float *, float *, integer *, integer *, PQueueType [MAXNCON][2]);
integer BetterBalance(integer, float *, float *, float *);
float Compute2WayHLoadImbalance(integer, float *, float *);
void Compute2WayHLoadImbalanceVec(integer, float *, float *, float *);

/* mfm2.c */
void MocFM_2WayEdgeRefine2(CtrlType *, GraphType *, float *, float *, integer);
void SelectQueue2(integer, float *, float *, integer *, integer *, PQueueType [MAXNCON][2], float *);
integer IsBetter2wayBalance(integer, float *, float *, float *);

/* mincover.o */
void MinCover(idxtype *, idxtype *, integer, integer, idxtype *, integer *);
integer MinCover_Augment(idxtype *, idxtype *, integer, idxtype *, idxtype *, idxtype *, integer);
void MinCover_Decompose(idxtype *, idxtype *, integer, integer, idxtype *, idxtype *, integer *);
void MinCover_ColDFS(idxtype *, idxtype *, integer, idxtype *, idxtype *, integer);
void MinCover_RowDFS(idxtype *, idxtype *, integer, idxtype *, idxtype *, integer);

/* minitpart.c */
void MocInit2WayPartition(CtrlType *, GraphType *, float *, float);
void MocGrowBisection(CtrlType *, GraphType *, float *, float);
void MocRandomBisection(CtrlType *, GraphType *, float *, float);
void MocInit2WayBalance(CtrlType *, GraphType *, float *);
integer SelectQueueoneWay(integer, float *, float *, integer, PQueueType [MAXNCON][2]);

/* minitpart2.c */
void MocInit2WayPartition2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisection2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisectionNew2(CtrlType *, GraphType *, float *, float *);
void MocInit2WayBalance2(CtrlType *, GraphType *, float *, float *);
integer SelectQueueOneWay2(integer, float *, PQueueType [MAXNCON][2], float *);

/* mkmetis.c */
void METIS_mCPartGraphKway(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
integer MCMlevelKWayPartitioning(CtrlType *, GraphType *, integer, idxtype *, float *);

/* mkwayfmh.c */
void MCRandom_KWayEdgeRefineHorizontal(CtrlType *, GraphType *, integer, float *, integer);
void MCGreedy_KWayEdgeBalanceHorizontal(CtrlType *, GraphType *, integer, float *, integer);
integer AreAllHVwgtsBelow(integer, float, float *, float, float *, float *);
integer AreAllHVwgtsAbove(integer, float, float *, float, float *, float *);
void ComputeHKWayLoadImbalance(integer, integer, float *, float *);
integer MocIsHBalanced(integer, integer, float *, float *);
integer IsHBalanceBetterFT(integer, integer, float *, float *, float *, float *);
integer IsHBalanceBetterTT(integer, integer, float *, float *, float *, float *);

/* mkwayrefine.c */
void MocRefineKWayHorizontal(CtrlType *, GraphType *, GraphType *, integer, float *);
void MocAllocateKWayPartitionMemory(CtrlType *, GraphType *, integer);
void MocComputeKWayPartitionParams(CtrlType *, GraphType *, integer);
void MocProjectKWayPartition(CtrlType *, GraphType *, integer);
void MocComputeKWayBalanceBoundary(CtrlType *, GraphType *, integer);

/* mmatch.c */
void MCMatch_RM(CtrlType *, GraphType *);
void MCMatch_HEM(CtrlType *, GraphType *);
void MCMatch_SHEM(CtrlType *, GraphType *);
void MCMatch_SHEBM(CtrlType *, GraphType *, integer);
void MCMatch_SBHEM(CtrlType *, GraphType *, integer);
float BetterVBalance(integer, integer, float *, float *, float *);
integer AreAllVwgtsBelowFast(integer, float *, float *, float);

/* mmd.c */
void genmmd(integer, idxtype *, idxtype *, idxtype *, idxtype *, integer , idxtype *, idxtype *, idxtype *, idxtype *, integer, integer *);
void mmdelm(integer, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, integer, integer);
integer  mmdint(integer, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(integer, idxtype *, idxtype *, idxtype *);
void mmdupd(integer, integer, idxtype *, idxtype *, integer, integer *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, integer, integer *tag);

/* mpmetis.c */
void METIS_mCPartGraphRecursive(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *);
void METIS_mCHPartGraphRecursive(integer *, integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *);
void METIS_mCPartGraphRecursiveInternal(integer *, integer *, idxtype *, idxtype *, float *, idxtype *, integer *, integer *, integer *, idxtype *);
void METIS_mCHPartGraphRecursiveInternal(integer *, integer *, idxtype *, idxtype *, float *, idxtype *, integer *, float *, integer *, integer *, idxtype *);
integer MCMlevelRecursiveBisection(CtrlType *, GraphType *, integer, idxtype *, float, integer);
integer MCHMlevelRecursiveBisection(CtrlType *, GraphType *, integer, idxtype *, float *, integer);
void MCMlevelEdgeBisection(CtrlType *, GraphType *, float *, float);
void MCHMlevelEdgeBisection(CtrlType *, GraphType *, float *, float *);

/* mrefine.c */
void MocRefine2Way(CtrlType *, GraphType *, GraphType *, float *, float);
void MocAllocate2WayPartitionMemory(CtrlType *, GraphType *);
void MocCompute2WayPartitionParams(CtrlType *, GraphType *);
void MocProject2WayPartition(CtrlType *, GraphType *);

/* mrefine2.c */
void MocRefine2Way2(CtrlType *, GraphType *, GraphType *, float *, float *);

/* mutil.c */
integer AreAllVwgtsBelow(integer, float, float *, float, float *, float);
integer AreAnyVwgtsBelow(integer, float, float *, float, float *, float);
integer AreAllVwgtsAbove(integer, float, float *, float, float *, float);
float ComputeLoadImbalance(integer, integer, float *, float *);
integer AreAllBelow(integer, float *, float *);

/* myqsort.c */
void iidxsort(integer, idxtype *);
void iintsort(integer, integer *);
void ikeysort(integer, KeyValueType *);
void ikeyvalsort(integer, KeyValueType *);

/* ometis.c */
void METIS_EdgeND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_NodeND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_NodeWND(integer *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, integer);
void MlevelNestedDissectionCC(CtrlType *, GraphType *, idxtype *, float, integer);
void MlevelNodeBisectionMultiple(CtrlType *, GraphType *, integer *, float);
void MlevelNodeBisection(CtrlType *, GraphType *, integer *, float);
void SplitGraphOrder(CtrlType *, GraphType *, GraphType *, GraphType *);
void MMDOrder(CtrlType *, GraphType *, idxtype *, integer);
integer SplitGraphOrderCC(CtrlType *, GraphType *, GraphType *, integer, idxtype *, idxtype *);

/* parmetis.c */
void METIS_PartGraphKway2(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void METIS_WPartGraphKway2(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
void METIS_NodeNDP(integer, idxtype *, idxtype *, integer, integer *, idxtype *, idxtype *, idxtype *);
void MlevelNestedDissectionP(CtrlType *, GraphType *, idxtype *, integer, integer, integer, idxtype *);
void METIS_NodeComputeSeparator(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *); 
void METIS_EdgeComputeSeparator(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, idxtype *); 

/* pmetis.c */
void METIS_PartGraphRecursive(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, integer *, integer *, idxtype *); 
void METIS_WPartGraphRecursive(integer *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, integer *, integer *, float *, integer *, integer *, idxtype *); 
integer MlevelRecursiveBisection(CtrlType *, GraphType *, integer, idxtype *, float *, float, integer);
void MlevelEdgeBisection(CtrlType *, GraphType *, integer *, float);
void SplitGraphPart(CtrlType *, GraphType *, GraphType *, GraphType *);
void SetUpSplitGraph(GraphType *, GraphType *, integer, integer);

/* pqueue.c */
void PQueueInit(CtrlType *ctrl, PQueueType *, integer, integer);
void PQueueReset(PQueueType *);
void PQueueFree(CtrlType *ctrl, PQueueType *);
integer PQueueGetSize(PQueueType *);
integer PQueueInsert(PQueueType *, integer, integer);
integer PQueueDelete(PQueueType *, integer, integer);
integer PQueueUpdate(PQueueType *, integer, integer, integer);
void PQueueUpdateUp(PQueueType *, integer, integer, integer);
integer PQueueGetMax(PQueueType *);
integer PQueueSeeMax(PQueueType *);
integer PQueueGetKey(PQueueType *);
integer CheckHeap(PQueueType *);

/* refine.c */
void Refine2Way(CtrlType *, GraphType *, GraphType *, integer *, float ubfactor);
void Allocate2WayPartitionMemory(CtrlType *, GraphType *);
void Compute2WayPartitionParams(CtrlType *, GraphType *);
void Project2WayPartition(CtrlType *, GraphType *);

/* separator.c */
void ConstructSeparator(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator0(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator(CtrlType *, GraphType *, float);

/* sfm.c */
void FM_2WayNodeRefine(CtrlType *, GraphType *, float, integer);
void FM_2WayNodeRefineEqWgt(CtrlType *, GraphType *, integer);
void FM_2WayNodeRefine_OneSided(CtrlType *, GraphType *, float, integer);
void FM_2WayNodeBalance(CtrlType *, GraphType *, float);
integer ComputeMaxNodeGain(integer, idxtype *, idxtype *, idxtype *);

/* srefine.c */
void Refine2WayNode(CtrlType *, GraphType *, GraphType *, float);
void Allocate2WayNodePartitionMemory(CtrlType *, GraphType *);
void Compute2WayNodePartitionParams(CtrlType *, GraphType *);
void Project2WayNodePartition(CtrlType *, GraphType *);

/* stat.c */
void ComputePartitionInfo(GraphType *, integer, idxtype *);
void ComputePartitionInfoBipartite(GraphType *, integer, idxtype *);
void ComputePartitionBalance(GraphType *, integer, idxtype *, float *);
float ComputeElementBalance(integer, integer, idxtype *);

/* subdomains.c */
void Random_KWayEdgeRefineMConn(CtrlType *, GraphType *, integer, float *, float, integer, integer);
void Greedy_KWayEdgeBalanceMConn(CtrlType *, GraphType *, integer, float *, float, integer);
void PrintSubDomainGraph(GraphType *, integer, idxtype *);
void ComputeSubDomainGraph(GraphType *, integer, idxtype *, idxtype *);
void EliminateSubDomainEdges(CtrlType *, GraphType *, integer, float *);
void MoveGroupMConn(CtrlType *, GraphType *, idxtype *, idxtype *, integer, integer, integer, idxtype *);
void EliminateComponents(CtrlType *, GraphType *, integer, float *, float);
void MoveGroup(CtrlType *, GraphType *, integer, integer, integer, idxtype *, idxtype *);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* util.c */
void errexit(char *,...);
#ifndef DMALLOC
integer *imalloc(integer, char *);
idxtype *idxmalloc(integer, char *);
float *fmalloc(integer, char *);
integer *ismalloc(integer, integer, char *);
idxtype *idxsmalloc(integer, idxtype, char *);
void *GKmalloc(size_t, char *);
#endif
/*void GKfree(void **,...); */
integer *iset(integer n, integer val, integer *x);
idxtype *idxset(integer n, idxtype val, idxtype *x);
float   *sset(integer n, float val, float *x);
integer iamax(integer, integer *);
integer idxamax(integer, idxtype *);
integer idxamax_strd(integer, idxtype *, integer);
integer samax(integer, float *);
integer samax2(integer, float *);
integer idxamin(integer, idxtype *);
integer samin(integer, float *);
integer idxsum(integer, idxtype *);
integer idxsum_strd(integer, idxtype *, integer);
void    idxadd(integer, idxtype *, idxtype *);
integer charsum(integer, char *);
integer isum(integer, integer *);
float   ssum(integer, float *);
float   ssum_strd(integer n, float *x, integer);
void    sscale(integer n, float, float *x);
float   snorm2(integer, float *);
float   sdot(integer n, float *, float *);
void    saxpy(integer, float, float *, integer, float *, integer);
void RandomPermute(integer, idxtype *, integer);
double drand48();
void srand48(long);
integer ispow2(integer);
void InitRandom(integer);
integer metis_log2(integer);










/***************************************************************
* Programs Directory
****************************************************************/

/* io.c */
void ReadGraph(GraphType *, char *, integer *);
void WritePartition(char *, idxtype *, integer, integer);
void WriteMeshPartition(char *, integer, integer, idxtype *, integer, idxtype *);
void WritePermutation(char *, idxtype *, integer);
integer CheckGraph(GraphType *);
idxtype *ReadMesh(char *, integer *, integer *, integer *);
void WriteGraph(char *, integer, idxtype *, idxtype *);

/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
integer smbfct(integer, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, integer *, idxtype *, idxtype *, integer *);


/***************************************************************
* Test Directory
****************************************************************/
void Test_PartGraph(integer, idxtype *, idxtype *);
integer VerifyPart(integer, idxtype *, idxtype *, idxtype *, idxtype *, integer, integer, idxtype *);
integer VerifyWPart(integer, idxtype *, idxtype *, idxtype *, idxtype *, integer, float *, integer, idxtype *);
void Test_PartGraphV(integer, idxtype *, idxtype *);
integer VerifyPartV(integer, idxtype *, idxtype *, idxtype *, idxtype *, integer, integer, idxtype *);
integer VerifyWPartV(integer, idxtype *, idxtype *, idxtype *, idxtype *, integer, float *, integer, idxtype *);
void Test_PartGraphmC(integer, idxtype *, idxtype *);
integer VerifyPartmC(integer, integer, idxtype *, idxtype *, idxtype *, idxtype *, integer, float *, integer, idxtype *);
void Test_ND(integer, idxtype *, idxtype *);
integer VerifyND(integer, idxtype *, idxtype *);



#endif
