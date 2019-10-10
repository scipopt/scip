/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce.h
 * @brief  includes various reduction methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_REDUCE_H_
#define APPLICATIONS_STP_SRC_REDUCE_H_

#define EXT_CLOSENODES_MAXN 64


#include "scip/scip.h"
#include "graph.h"

/** path root state */
typedef struct pathroot_state
{
   int pathroot_id;
   int pathroot_nrecomps;
} PRSTATE;

/** distance data */
typedef struct distance_data
{
   //    SCIP_Bool*  pathrootsd_isdirty;
   DHEAP* dheap;
   RANGE* closenodes_range;          /* of size nnodes */
   int* closenodes_indices;          /* of size closenodes_totalsize */
   SCIP_Real* closenodes_distances;  /* of size closenodes_totalsize */
   int* pathroot_blocksizes;         /* of size nedges / 2 */
   int* pathroot_blocksizesmax;      /* of size nedges / 2 */
   PRSTATE** pathroot_blocks;        /* of size nedges / 2 */
   SCIP_Bool* pathroot_isdirty;      /* of size nnodes */
   int* pathroot_nrecomps;           /* of size nnodes */
   DIJK* dijkdata;
   int closenodes_totalsize;
   int closenodes_maxnumber;
} DISTDATA;


/** reduced cost result data */
typedef struct reduce_costs_data
{
   const SCIP_Real* const redEdgeCost;           /**< reduced costs */
   const SCIP_Real* const rootToNodeDist;        /**< shortest path distances from root  */
   const PATH* const nodeTo3TermsPaths;          /**< paths to three nearest terminals */
   const SCIP_Real cutoff;                       /**< reduced cost cutoff value or -1.0 if not used */
   const int redCostRoot;                        /**< graph root for reduced cost calculation */
} REDCOST;


/* reduce.c
 */
extern SCIP_RETCODE level0(SCIP*, GRAPH*);
extern SCIP_RETCODE level0infeas(SCIP*, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE level0RpcRmw(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE level0RpcRmwInfeas(SCIP*, GRAPH*, SCIP_Real*, SCIP_Bool*);
extern SCIP_RETCODE reduceStp(SCIP*, GRAPH*, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reducePc(SCIP*, const int*, GRAPH*, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopStp(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Real, SCIP_Bool, SCIP_Bool, SCIP_Bool, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopPc(SCIP*, const int*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Bool, SCIP_Bool, SCIP_Bool, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopMw(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, STP_Bool, STP_Bool, STP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE reduce(SCIP*, GRAPH*, SCIP_Real*, int, int, SCIP_Bool);

/* reduce_alt.c
 */
extern void    reduce_ans(SCIP*, GRAPH*, int*, int*);
extern void    reduce_ansAdv(SCIP*, GRAPH*, int*, int*, SCIP_Bool);
extern void    reduce_ansAdv2(SCIP*, GRAPH*, int*, int*);
extern void    reduce_nnp(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_sdsp(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, int*);
extern SCIP_RETCODE    reduce_sdStar(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdStarPc(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalk(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalk_csr(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkTriangle(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt2(SCIP*, int, const int*, GRAPH*, int*,  SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdspSap(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_sd(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, SCIP_Bool, int*);
extern SCIP_RETCODE    reduce_sdPc(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_getSd(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_getSdPcMw(SCIP*, const GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE    reduce_nts(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_bd34(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, SCIP_Real*);
extern SCIP_RETCODE    reduce_bd34WithSd(SCIP*, GRAPH*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nv(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nvAdv(SCIP*, const int*, GRAPH*, PATH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_sl(SCIP*, const int*, GRAPH*, PATH*, double*, int*, int*, int*, int*, STP_Bool*, int*, int*);
extern SCIP_RETCODE    reduce_ledge(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_cnsAdv(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_npv(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_chain2(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);

/* reduce_bnd.c
 */
extern SCIP_RETCODE    reduce_da(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, STP_Bool*, int*, int, SCIP_RANDNUMGEN*, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_daSlackPrune(SCIP*, SCIP_VAR**, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    reduce_daPcMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, int*, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_RANDNUMGEN*, SCIP_Real, SCIP_Bool);
extern SCIP_RETCODE    reduce_bound(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundMw(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundPrune(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, const int*, const int*, int*, int);
extern SCIP_RETCODE    reduce_boundHop(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopR(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopRc(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);


/* reduce_ext.c
 */
extern SCIP_RETCODE    reduce_deleteConflictEdges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_extendedCheck3Tree(SCIP*, const GRAPH*, int, const SCIP_Real*, const SCIP_Real*, const PATH*, const int*, SCIP_Real, const int*, int, int*, SCIP_Real*, SCIP_Bool*, unsigned int*, int*, SCIP_Bool*);
extern int             reduce_extendedEdge(SCIP*, GRAPH*, const PATH*, const SCIP_Real*, const SCIP_Real*, const int*, SCIP_Real, int, int*, STP_Bool*, SCIP_Bool);


/* reduce_ext2.c
 */
extern SCIP_RETCODE    reduce_extendedEdge2(SCIP*, GRAPH*, const PATH*, const SCIP_Real*, const SCIP_Real*, const int*, SCIP_Real, int, SCIP_Bool,  STP_Bool*, int*);
extern SCIP_RETCODE    reduce_extendedCheckArc(SCIP*, const GRAPH*, const REDCOST*, STP_Bool*,  const SCIP_Bool*, int, SCIP_Bool, DISTDATA*, SCIP_Real*, int*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_extendedCheckEdge(SCIP*, const GRAPH*, const REDCOST*, const STP_Bool*,  const SCIP_Bool*, int, SCIP_Bool, DISTDATA*, SCIP_Real*, int*, SCIP_Bool*);


/* reduce_test.c
 */
extern SCIP_RETCODE    testAll(SCIP*);
extern SCIP_RETCODE    pseudoAncestors_test(SCIP*);
extern SCIP_RETCODE    pseudoDel_test(SCIP*);
extern SCIP_RETCODE    reduce_extTest(SCIP*);
extern SCIP_RETCODE    dheap_Test(SCIP*);
extern SCIP_RETCODE    reduce_sdPcMwTest(SCIP*);
extern SCIP_RETCODE    heur_extendPcMwOuterTest(SCIP*);


/* reduce_simple.c
 */
extern SCIP_RETCODE    reduce_simple_deleteMultiedges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_simple_contract0Edges(SCIP*, GRAPH*, SCIP_Bool);
extern SCIP_RETCODE    reduce_simple(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE    reduce_simple_hc(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_mw(SCIP*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_pc(SCIP*, const int*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE    reduce_simple_aritculations(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_fixedConflict(SCIP*, const int*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_simple_sap(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_rpt(SCIP*, GRAPH*, SCIP_Real*, int*);

/* reduce_util.c
 */
extern SCIP_RETCODE    reduce_distDataInit(SCIP*, const GRAPH*, int, SCIP_Bool, DISTDATA*);
extern SCIP_Real       reduce_distDataGetSD(SCIP*, const GRAPH*, int, int, DISTDATA*);
extern void            reduce_distDataFreeMembers(SCIP*, const GRAPH*, DISTDATA*);
extern void            reduce_distDataDeleteEdge(SCIP*, const GRAPH*, int, DISTDATA*);


#endif /* APPLICATIONS_STP_SRC_REDUCE_H_ */
