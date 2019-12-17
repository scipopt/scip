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

/**@file   completegraph.c
 * @brief  includes complete graph methods, in particular for MST calculation
 * @author Daniel Rehfeldt
 *
 * Complete graph methods, in particular for MST calculation. Assumes a maximum size of the graph.
 * Only sensible for small maximum sizes.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#include "completegraph.h"


/** gets start position in edge cost array */
static inline
int getStartPos(
   int                   nodepos,            /**< node position */
   int                   nnodes_max          /**< maximum number of nodes */
)
{
   assert(nodepos >= 0);
   assert(nodepos < nnodes_max);

   return (nodepos * nnodes_max);
}


/** gets end position in edge cost array */
static inline
int getEndPos(
   int                   nodepos,            /**< node position */
   int                   nnodes_curr,        /**< current number of nodes */
   int                   nnodes_max,         /**< maximum number of nodes */
)
{
   assert(nodepos >= 0);
   assert(nnodes_curr <= nnodes_max);
   assert(nodepos < nnodes_curr);

   return (nodepos * nnodes_max + nnodes_curr);
}



/** initialize complete, undirected graph */
SCIP_RETCODE cgraph_init(
   SCIP*                 scip,               /**< SCIP data structure */
   CGRAPH**              cgraph,             /**< new graph */
   int                   maxnnodes           /**< maximum number of nodes */
   )
{
   CGRAPH* g;

   assert(scip && cgraph);
   assert(maxnnodes > 1);

   SCIP_CALL( SCIPallocMemory(scip, cgraph) );

   g = *cgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->edgecosts), maxnnodes * (maxnnodes - 1)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->nodeids), maxnnodes) );

   g->nnodes_curr = 0;
   g->nnodes_max = maxnnodes;

   return SCIP_OKAY;
}


/** adds node */
void cgraph_node_append(
   CGRAPH*               cgraph,             /**< new graph */
   const SCIP_Real*      adjcosts,           /**< array with edge costs to neigbors of size nnodes_curr (use -1.0 as debug sentinel) */
   int                   nodeid              /**< the node id */
   )
{
   SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int nnodes_curr = cgraph->nnodes_curr;
   const int nnodes_max = cgraph->nnodes_curr;
   const int start_new = getStartPos(nnodes_curr, nnodes_max);

   assert(adjcosts);
   assert(nodeid >= 0);
   assert(nnodes_curr < nnodes_max);

   BMScopyMemoryArray(edgecosts + start_new, adjcosts, nnodes_curr);

#ifndef NDEBUG
   for( int j = 0; j < nnodes_curr; j++ )
      assert(adjcosts[j] >= 0.0);
#endif

   /* adapt all other edges */
   for( int i = 0; i < nnodes_curr; i++ )
   {
      const int end = getEndPos(i, nnodes_curr + 1, nnodes_max);
      edgecosts[end] = adjcosts[i];
   }

   cgraph->nodeids[nnodes_curr] = nodeid;
   cgraph->edgecosts[nnodes_curr] = FARAWAY;

   cgraph->nnodes_curr++;
}


/** deletes node */
void cgraph_node_deleteTop(
   CGRAPH*               cgraph              /**< new graph */
   )
{
   assert(cgraph);
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);
   assert(cgraph->nnodes_curr > 0);

   cgraph->nnodes_curr--;
}


/** deletes node */
void cgraph_node_exchange(
   CGRAPH*               cgraph,             /**< new graph */
   const SCIP_Real*      adjcosts,           /**< array with edge costs to neigbors of size nnodes_curr (FARAWAY at own position) */
   int                   nodepos,            /**< the node position */
   int                   nodeid_old,         /**< node id (for debugging) */
   int                   nodeid_new,         /**< node id of new entry */
   )
{
   assert(cgraph);
   assert(nodeid_new >= 0);
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);


}

/** Get edge cost.
 * Note: quite slow, only use for debugging! */
SCIP_Real cgraph_edge_getCost(
   const CGRAPH*         cgraph,             /**< new graph */
   int                   nodeid_tail,        /**< the node id */
   int                   nodeid_head        /**< the node id */
   )
{
   assert(cgraph);

   return 0.0;
}

/** free complete graph */
void cgraph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   CGRAPH**              cgraph              /**< new graph */
   )
{
   CGRAPH* g;

   assert(scip && cgraph);

   g = *cgraph;

   SCIPfreeMemoryArray(scip, &(g->nodeids));
   SCIPfreeMemoryArray(scip, &(g->edgecosts));

   SCIPfreeMemory(scip, cgraph);
}
