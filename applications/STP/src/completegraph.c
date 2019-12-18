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
#include "portab.h"


#define NODE_ID_UNDEFINED -2
#define EDGECOST_UNDEFINED_VALUE -1.0


/** gets current number of nodes */
static inline
int getNnodesCurr(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   assert(cgraph);
   return cgraph->nnodes_curr;
}


/** gets maximum number of nodes */
static inline
int getNnodesMax(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   assert(cgraph);
   return cgraph->nnodes_max;
}


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
   int                   nnodes_max          /**< maximum number of nodes */
)
{
   assert(nodepos >= 0);
   assert(nnodes_curr <= nnodes_max);
   assert(nodepos < nnodes_curr);

   return (nodepos * nnodes_max + nnodes_curr);
}


/** is the graph valid? */
SCIP_Bool cgraph_valid(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   const int nnodes_curr = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   const SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int* const nodeids = cgraph->nodeids;

   assert(edgecosts && nodeids);
   assert(nnodes_curr <= nnodes_max);
   assert(nnodes_curr >= 0);

   for( int i = getStartPos(nnodes_curr, nnodes_max); i < getEndPos(nnodes_max - 1, nnodes_max - 1, nnodes_max); i++ )
   {
      if( !EQ(EDGECOST_UNDEFINED_VALUE, edgecosts[i]) )
      {
         SCIPdebugMessage("unused edge value is wrong: %f \n", edgecosts[i]);
         return FALSE;
      }
   }

   for( int i = 0; i < nnodes_curr; i++ )
   {
      if( NODE_ID_UNDEFINED == nodeids[i] )
      {
         SCIPdebugMessage("node %d is wrongly set to NODE_ID_UNDEFINED \n", i);
         return FALSE;
      }
   }

   for( int i = nnodes_curr; i < nnodes_max; i++ )
   {
      if( NODE_ID_UNDEFINED != nodeids[i] )
      {
         SCIPdebugMessage("node %d is not set to NODE_ID_UNDEFINED \n", i);
         return FALSE;
      }
   }

   return TRUE;
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

#ifndef NDEBUG
   for( int i = 0; i < maxnnodes * (maxnnodes - 1); i++ )
      g->edgecosts[i] = EDGECOST_UNDEFINED_VALUE;

   for( int i = 0; i < maxnnodes; i++ )
      g->nodeids[i] = NODE_ID_UNDEFINED;
#endif

   assert(cgraph_valid(g));

   return SCIP_OKAY;
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

#ifndef NDEBUG
   assert(cgraph_valid(cgraph));
   assert(adjcosts && edgecosts && cgraph->nodeids);
   assert(nodeid >= 0);
   assert(nnodes_curr < nnodes_max);
   assert(NODE_ID_UNDEFINED == cgraph->nodeids[nnodes_curr]);

   for( int j = 0; j < nnodes_curr; j++ )
      assert(GE(adjcosts[j], 0.0));
#endif

   BMScopyMemoryArray(edgecosts + start_new, adjcosts, nnodes_curr);

   /* adapt all other edges (going to the new node) */
   for( int i = 0; i < nnodes_curr; i++ )
   {
      const int end = getEndPos(i, nnodes_curr + 1, nnodes_max);

      assert(EQ(edgecosts[end], EDGECOST_UNDEFINED_VALUE));

      edgecosts[end] = adjcosts[i];
   }

   edgecosts[nnodes_curr] = FARAWAY;
   cgraph->nodeids[nnodes_curr] = nodeid;

   cgraph->nnodes_curr++;

   assert(cgraph_valid(cgraph));
}


/** deletes node at the top */
void cgraph_node_deleteTop(
   CGRAPH*               cgraph              /**< new graph */
   )
{
   assert(cgraph_valid(cgraph));
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);
   assert(cgraph->nnodes_curr > 0);

   cgraph->nnodes_curr--;

#ifndef NDEBUG
   assert(NODE_ID_UNDEFINED != cgraph->nodeids[cgraph->nnodes_curr]);

   cgraph->nodeids[cgraph->nnodes_curr] = NODE_ID_UNDEFINED;

   for( int i = 0; i < cgraph->nnodes_curr; i++ )
   {
      const int end = getEndPos(i, cgraph->nnodes_curr, cgraph->nnodes_max);

      assert(end == cgraph->nnodes_curr);
      assert(!EQ(cgraph->edgecosts[end], EDGECOST_UNDEFINED_VALUE));

      cgraph->edgecosts[end] = EDGECOST_UNDEFINED_VALUE;
   }
#endif

   assert(cgraph_valid(cgraph));
}


/** deletes node */
void cgraph_node_exchange(
   CGRAPH*               cgraph,             /**< new graph */
   const SCIP_Real*      adjcosts,           /**< array with edge costs to neigbors of size nnodes_curr (FARAWAY at own position) */
   int                   nodepos,            /**< the node position */
   int                   nodeid_old,         /**< node id (for debugging) */
   int                   nodeid_new          /**< node id of new entry */
   )
{
   int todo; // what do we really need? nodeid_old, nodeid_new?

   assert(cgraph);
   assert(cgraph_valid(cgraph));
   assert(nodeid_new >= 0);
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);

}

/** Get edge cost.
 * Note: quite slow, only use for debugging! */
SCIP_Real cgraph_edge_getCost(
   const CGRAPH*         cgraph,             /**< new graph */
   int                   nodeid_tail,        /**< the node id */
   int                   nodeid_head         /**< the node id */
   )
{
   assert(cgraph);
   assert(cgraph_valid(cgraph));

   return 0.0;
}
