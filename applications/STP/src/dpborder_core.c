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

   printf("---PRINTING BORDER: \n");
   printf("extnode=%d, nbordernodes=%d \n", toplevel->extnode, nbordernodes);
   printf("bordernodes: \n");

   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int node = dpborder->bordernodes[i];
      assert(dpborder->nodes_isBorder[node]);

      printf("char=%d node=%d \n", i, node);

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


/** does reallocation of global array if necessary */
static inline
SCIP_RETCODE partitionTryRealloc(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   newadds,
   DPBORDER*             dpborder             /**< border */
   )
{
   const int size = dpborder->global_partstarts[dpborder->global_npartitions];
   assert(size >= 0);
   assert(newadds > 0);


   if( newadds + size > dpborder->global_partcap )
   {
      dpborder->global_partcap *= 2;
      assert(newadds + dpborder->global_npartitions > size);
      SCIPdebugMessage("reallocating memory (to %d) \n", dpborder->global_partcap);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(dpborder->global_partitions), dpborder->global_partcap) );
   }

   return SCIP_OKAY;
}


/** gets partition range */
static inline
void partitionGetRangeGlobal(
   const DPBORDER*       dpborder,            /**< border */
   int                   globalposition,     /**< position of partition */
   int*                  start,
   int*                  end
   )
{
   const int* const starts = dpborder->global_partstarts;
   assert(globalposition >= 0);
   assert(globalposition < dpborder->global_npartitions);

   *start = starts[globalposition];
   *end = starts[globalposition + 1];

   assert(*start < *end);
}


/** adapts border */
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

   /* NOTE: needs to be called before border is updated! */
   dpborder_buildBorderDists(graph, dpborder);

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

   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];
      assert(graph_knot_isInRange(graph, bordernode));

      /* NOTE: might happen for previous border node */
      if( nodes_isBorder[bordernode] && nodes_outdegree[bordernode] == 0 )
         nodes_isBorder[bordernode] = FALSE;
   }

   dpborder_buildBorderMap(dpborder);

   /* remove outdated border nodes */
   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];

      if( nodes_isBorder[bordernode] )
         dpborder->bordernodes[nbordernodes++] = dpborder->bordernodes[i];
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


/** updates partition from given partition */
static
SCIP_RETCODE updateFromPartition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   globalposition,     /**< position of partition */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBPART partition;
   int* subbuffer;
   STP_Vectype(int) candstarts;
   DPBLEVEL* const toplevel = dpborder_getTopLevel(dpborder);
   DPB_Ptype* global_partitions = dpborder->global_partitions;
   int part_start;
   int part_end;
   int ncands;
   uint32_t powsize;
   const SCIP_Real part_cost = dpborder->global_partcosts[globalposition];
   const SCIP_Bool allTermsAreVisited = (dpborder->nterms == dpborder->ntermsvisited);

   partitionGetRangeGlobal(dpborder, globalposition, &part_start, &part_end);

   // todo check for P \ L
   if( !toplevel->exnodeIsTerm && 0  )
   {

   }

   partition.partchars = &(global_partitions[part_start]);
   partition.partsize = (part_end - part_start);
   partition.delimiter = dpborder_getDelimiter(dpborder, StpVecGetSize(dpborder->borderlevels) - 2);

   candstarts = dpborder_partGetCandstarts(scip, &partition, dpborder);
   assert(StpVecGetSize(candstarts) > 0);

   ncands = StpVecGetSize(candstarts);
   assert(ncands > 0);
   powsize = (uint32_t) pow(2.0, ncands);
   SCIP_CALL( SCIPallocBufferArray(scip, &subbuffer, BPBORDER_MAXBORDERSIZE) );

   /* make sure that partition storage is large enough */
   SCIP_CALL( partitionTryRealloc(scip, powsize * BPBORDER_MAXBORDERSIZE * 2, dpborder) );
   global_partitions = dpborder->global_partitions;
   partition.partchars = &(global_partitions[part_start]);

   SCIPdebugMessage("partition ncands=%d \n", ncands);

   dpborder_partPrint(&partition);

   /* loop over all subsets (also the empty set) */
   for( uint32_t counter = powsize; counter >= 1; counter-- )
   {
      int nsub = 0;
      int globalposition_new;
      SCIP_Real cost_new;
      const uint32_t mask = counter - 1;

      SCIPdebugMessage("new partition:  \n");
      for( uint32_t j = 0; j < (uint32_t) ncands; j++ )
      {
         if( mask & ((uint32_t) 1 << j) )
         {
            assert(nsub < BPBORDER_MAXBORDERSIZE);
            subbuffer[nsub++] = candstarts[j];
            SCIPdebugMessage("...%d \n", candstarts[j]);
         }
      }

      globalposition_new = dpborder_partGetIdxNew(scip, &partition, subbuffer, nsub, dpborder);

      /* non-valid partition? */
      if( globalposition_new == -1 )
      {
         SCIPdebugMessage("partition is invalid... \n");
         continue;
      }

      assert(globalposition_new >= 0);

      // todo compute connection cost
      cost_new = 0.0;
      cost_new += part_cost;

      if( GT(dpborder->global_partcosts[globalposition_new], cost_new) )
      {
         // todo update

         if( allTermsAreVisited )
         {
            // check whether size of parition is 1!
            // todo extra method to get cardinality!
            if( 0 )
            {
               dpborder->global_obj = MIN(dpborder->global_obj, cost_new);
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &subbuffer);
   StpVecFree(scip, candstarts);

   return SCIP_OKAY;
}


/** adds partitions for new border */
static
SCIP_RETCODE addPartitions(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   iteration,          /**< iteration number */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBLEVEL* const toplevel = dpborder->borderlevels[iteration];
   DPBLEVEL* const prevlevel = dpborder->borderlevels[iteration - 1];
   const int global_start = prevlevel->globalstartidx;
   const int global_end = toplevel->globalstartidx;

   SCIPdebugMessage("adding partitions: \n");

   {
      int todo;

      if( iteration == 3 )
      {
         StpVecPushBack(scip, dpborder->global_partcosts, 2.0);
         StpVecPushBack(scip, dpborder->global_partstarts, 4);
         dpborder->global_partitions[0] = 0;
         dpborder->global_partitions[1] = 1;
         dpborder->global_partitions[2] = 3;
         dpborder->global_partitions[3] = 2;


         //StpVecPushBack(scip, dpborder->global_partitions, 0);
         //StpVecPushBack(scip, dpborder->global_partitions, 1);
         //StpVecPushBack(scip, dpborder->global_partitions, 3); // delim
         //StpVecPushBack(scip, dpborder->global_partitions, 2);

         dpborder->global_npartitions++;
         dpborder->borderchardists[1] = 2.0;
         dpborder->borderchardists[2] = 2.0;

         SCIP_CALL( updateFromPartition(scip, 0, graph, dpborder) );
         assert(0);

      }


      return SCIP_OKAY;
   }

   assert(iteration >= 1);
   assert(global_start < global_end);

   /* loop over all valid partitions of previous border */
   for( int i = global_start; i != global_end; i++ )
   {
      SCIP_CALL( updateFromPartition(scip, i, graph, dpborder) );
   }

   return SCIP_OKAY;
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

#ifdef SCIP_DEBUG
   printf("\n----- Starting level %d ----- \n", iteration);
#endif

   SCIP_CALL( dpborder_dpblevelInit(scip, &level) );
   StpVecPushBack(scip, dpborder->borderlevels, level);

   level->extnode = node;
   level->exnodeIsTerm = Is_term(graph->term[node]);
   level->globalstartidx = dpborder->global_npartitions;

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
   level->globalstartidx = 0;

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

   dpborder->global_partcap = MAX(maxnpartitions / 2, BPBORDER_MAXBORDERSIZE);
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
      SCIP_CALL( addPartitions(scip, graph, s, dpborder) );

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
