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

/**@file   extreduce.h
 * @brief  This file implements extended reduction techniques for several Steiner problems.
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_EXTREDUCE_H_
#define APPLICATIONS_STP_SRC_EXTREDUCE_H_

#define EXT_CLOSENODES_MAXN 64
#define STP_EXTTREE_MAXNLEAVES 20

#include "scip/scip.h"
#include "graph.h"
#include "reduce.h"
#include "completegraph.h"


/** permanent extension data */
typedef struct extension_data_permanent
{
   CMST*                 cmst;               /**< complete MST */
   CGRAPH*               cgraph;             /**< complete graph */
   int*                  cgraphEdgebuffer;   /**< complete graph edge buffer */
   STP_Bool*             edgedeleted;        /**< (non-owned!) edge array to mark which directed edge can be removed */
   SCIP_Bool*            isterm;             /**< marks whether node is a terminal (or proper terminal for PC) */
   SCIP_Real*            bottleneckDistNode; /**< needs to be set to -1.0 (size nnodes) */
   SCIP_Real*            pcSdToNode;         /**< needs to be set to -1.0 for all nodes, or NULL if not Pc */
   int*                  tree_deg;           /**< -1 for forbidden nodes (e.g. PC terminals), nnodes for tail, 0 otherwise; in method: position ( > 0) for nodes in tree */
   int                   nnodes;             /**< number of nodes */
} EXTPERMA;


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


/* extreduce_base.c
 */
extern SCIP_RETCODE    extreduce_deleteEdges(SCIP*, const REDCOST*, const int*, GRAPH*, STP_Bool*, int*);
extern SCIP_RETCODE    extreduce_checkArc(SCIP*, const GRAPH*, const REDCOST*, int, SCIP_Bool, DISTDATA*, EXTPERMA*, SCIP_Bool*);
extern SCIP_RETCODE    extreduce_checkEdge(SCIP*, const GRAPH*, const REDCOST*, int, SCIP_Bool, DISTDATA*, EXTPERMA*, SCIP_Bool*);


/* extreduce_util.c
 */
extern SCIP_RETCODE    extreduce_distDataInit(SCIP*, const GRAPH*, int, SCIP_Bool, DISTDATA*);
extern SCIP_Real       extreduce_distDataGetSD(SCIP*, const GRAPH*, int, int, DISTDATA*);
extern void            extreduce_distDataFreeMembers(SCIP*, const GRAPH*, DISTDATA*);
extern void            extreduce_distDataDeleteEdge(SCIP*, const GRAPH*, int, DISTDATA*);
extern SCIP_RETCODE    extreduce_extPermaInit(SCIP*, const GRAPH*, STP_Bool*, EXTPERMA*);
extern SCIP_Bool       extreduce_extPermaIsClean(const GRAPH*, const EXTPERMA*);
extern void            extreduce_extPermaFreeMembers(SCIP*, EXTPERMA*);




#endif /* APPLICATIONS_STP_SRC_EXTREDUCE_H_ */
