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

/**@file   dpborder_core.c
 * @brief  Core of dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 *
 * Contains core methods.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "dpborder.h"
#include "dpborderinterns.h"
#include "stpvector.h"

/*
 * Local methods
 */

#ifdef SCIP_DEBUG
/** prints border nodes */
static
void printBorder(
   const GRAPH*          graph,              /**< graph */
   int                   iteration,          /**< iteration number */
   const DPBORDER*       dpborder            /**< border */
)
{
   const DPBLEVEL* const toplevel = dpborder->borderlevels[iteration];
   const int nbordernodes = toplevel->nbordernodes;

   printf("\n---PRINTING BORDER: \n");
   printf("extnode=%d, nbordernodes=%d \n", toplevel->extnode, nbordernodes);
   printf("bordernodes: \n");

   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int node = dpborder->bordernodes[i];
      assert(dpborder->nodes_isBorder[node]);

      graph_knot_printInfo(graph, node);
   }

   printf("prev-bordernodes: \n");

   for( int i = 0; i < StpVecGetSize(dpborder->prevbordernodes); i++ )
   {
      const int node = dpborder->prevbordernodes[i];
      assert(!dpborder->nodes_isBorder[node] || node == toplevel->extnode);

      graph_knot_printInfo(graph, node);
   }

}
#endif


/** adapts degree of nodes */
static
void updateBorder(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   iteration,          /**< iteration number */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBLEVEL* const toplevel = dpborder->borderlevels[iteration];
   SCIP_Bool* const nodes_isBorder = dpborder->nodes_isBorder;
   int* const nodes_outdegree = dpborder->nodes_outdeg;
   const CSR* const csr = graph->csr_storage;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   const int extnode = toplevel->extnode;
   int nbordernodes = 0;

   assert(iteration >= 1);
   assert(!nodes_isBorder[extnode]);
   assert(dpborder->borderlevels[iteration - 1]->nbordernodes == StpVecGetSize(dpborder->bordernodes));

   if( dpborder->prevbordernodes )
      StpVecClear(dpborder->prevbordernodes);

   for( int e = start_csr[extnode]; e != start_csr[extnode + 1]; e++ )
   {
      const int head = head_csr[e];
      if( nodes_isBorder[head] )
      {
         nodes_outdegree[extnode]--;
         nodes_outdegree[head]--;
         assert(nodes_outdegree[extnode] >= 0);
         assert(nodes_outdegree[head] >= 0);

         if( nodes_outdegree[head] == 0 )
         {
            StpVecPushBack(scip, dpborder->prevbordernodes, head);
            nodes_isBorder[head] = FALSE;
         }
      }
   }

   if( nodes_outdegree[extnode] == 0 )
      StpVecPushBack(scip, dpborder->prevbordernodes, extnode);

   /* remove outdated border nodes */
   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];
      assert(graph_knot_isInRange(graph, bordernode));

      if( nodes_isBorder[bordernode] )
      {
         /* NOTE: might happen for previous border node */
         if( nodes_outdegree[bordernode] == 0 )
         {
            nodes_isBorder[bordernode] = FALSE;
            continue;
         }

         dpborder->bordernodes[nbordernodes++] = dpborder->bordernodes[i];
      }
   }

   for( int i = StpVecGetSize(dpborder->bordernodes) - nbordernodes; i > 0; i-- )
      StpVecPopBack(dpborder->bordernodes);

   if( nodes_outdegree[extnode] != 0 )
   {
      nbordernodes++;
      StpVecPushBack(scip, dpborder->bordernodes, extnode);
      nodes_isBorder[extnode] = TRUE;
   }

   toplevel->nbordernodes = nbordernodes;
   assert(toplevel->nbordernodes == StpVecGetSize(dpborder->bordernodes));
}


/** adds border level */
static
SCIP_RETCODE addLevel(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   iteration,          /**< iteration number */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBLEVEL* level;
   const DPBSEQUENCE* const dpbsequence = dpborder->dpbsequence;
   const int* const nodessequence = dpbsequence->nodessquence;
   const int node = nodessequence[iteration];

   assert(graph_knot_isInRange(graph, node));

   SCIP_CALL( dpborder_dpblevelInit(scip, &level) );
   StpVecPushBack(scip, dpborder->borderlevels, level);

   level->extnode = node;
   level->exnodeIsTerm = Is_term(graph->term[node]);

   if( level->exnodeIsTerm )
      dpborder->ntermsvisited++;

   updateBorder(scip, graph, iteration, dpborder);

#ifdef SCIP_DEBUG
   printBorder(graph, iteration, dpborder);
#endif

   return SCIP_OKAY;
}


/** adds first border level */
static
SCIP_RETCODE addLevelFirst(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBLEVEL* level;
   const DPBSEQUENCE* const dpbsequence = dpborder->dpbsequence;
   const int* const nodessequence = dpbsequence->nodessquence;
   const int root = nodessequence[0];

   assert(Is_term(graph->term[root]));
   assert(!dpborder->nodes_isBorder[root]);

   SCIPdebugMessage("adding initial level for DP-root=%d \n", root);

   dpborder->nodes_isBorder[root] = TRUE;
   StpVecPushBack(scip, dpborder->bordernodes, root);
   StpVecPushBack(scip, dpborder->global_partstarts, 0);
   dpborder->ntermsvisited = 1;

   SCIP_CALL( dpborder_dpblevelInit(scip, &level) );
   StpVecPushBack(scip, dpborder->borderlevels, level);

   level->extnode = root;
   level->exnodeIsTerm = TRUE;
   level->nbordernodes = 1;

   assert(StpVecGetSize(dpborder->bordernodes) == 1);

   return SCIP_OKAY;
}


/** initializes for solve */
static
SCIP_RETCODE initSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   const DPBSEQUENCE* const dpbsequence = dpborder->dpbsequence;
   const uint64_t maxnpartitions = dpbsequence->maxnpartitions;

   assert(dpborder->nterms == graph->terms && dpborder->ntermsvisited == 0);
   assert(!dpborder->bordernodes);

   if( maxnpartitions > BPBORDER_MAXNPARTITIONS || maxnpartitions > INT_MAX )
   {
      SCIPerrorMessage("too many partitions! \n");
      return SCIP_ERROR;
   }

   dpborder->global_partcap = maxnpartitions / 2;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpborder->global_partitions), dpborder->global_partcap) );

   SCIP_CALL( addLevelFirst(scip, graph, dpborder) );

   return SCIP_OKAY;
}


/** computes node ordering and updates if better */
static
SCIP_RETCODE computeOrderingFromNode(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   root,               /**< node to start from */
   DPBSEQUENCE*          dpbsequence         /**< sequence */
)
{
   SCIP_Bool* RESTRICT nodes_isVisited;
   int* RESTRICT nodessquence = dpbsequence->nodessquence;
   const CSR* const csr = graph->csr_storage;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   const int nnodes = graph_get_nNodes(graph);
   int nvisited;
   int nscanned;

   assert(graph_knot_isInRange(graph, root));
   assert(Is_term(graph->term[root]));

   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodes_isVisited, nnodes) );
   assert(!nodes_isVisited[root]);

   nvisited = 1;
   nscanned = 0;
   nodessquence[0] = root;
   nodes_isVisited[root] = TRUE;

   while( nvisited < nnodes )
   {
      const int node = nodessquence[nscanned++];
      assert(nscanned <= nvisited);

      for( int e = start_csr[node]; e != start_csr[node + 1]; e++ )
      {
         const int head = head_csr[e];
         if( !nodes_isVisited[head] )
         {
            nodes_isVisited[head] = TRUE;
            nodessquence[nvisited++] = head;
         }
      }
   }

   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


// compute ordering


/** computes node ordering */
SCIP_RETCODE dpborder_coreComputeOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   // todo try some sources...

   SCIP_CALL( computeOrderingFromNode(scip, graph, graph->source, dpborder->dpbsequence) );

   return SCIP_OKAY;
}


/** solves problem */
SCIP_RETCODE dpborder_coreSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder,           /**< border */
   SCIP_Bool*            wasSolved           /**< was problem solved to optimality? */
)
{
   const int nnodes = graph_get_nNodes(graph);

   assert(scip && graph && dpborder && wasSolved);
   assert(dpborder->dpbsequence);

   *wasSolved = TRUE;

   SCIP_CALL( initSolve(scip, graph, dpborder) );

   /* DP loop */
   for( int s = 1; s < nnodes; s++ )
   {
      SCIP_CALL( addLevel(scip, graph, s, dpborder) );

      // todo update partitions

      if( SCIPisStopped(scip) )
      {
         *wasSolved = FALSE;
         break;
      }
   }


   if( *wasSolved )
   {
      SCIPdebugMessage("OBJ=%f \n", dpborder->global_obj);
      assert(graph->terms == dpborder->ntermsvisited);
   }


   return SCIP_OKAY;
}
