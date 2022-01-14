/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
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

//#define SCIP_DEBUG
#include "dpborder.h"
#include "dpborderinterns.h"
#include "stpvector.h"
#include "misc_stp.h"

#define DPB_ORDERMULT_PREVS    2
#define DPB_ORDERMULT_TERM     5
#define DPB_ORDERMULT_OUTDEG   1
#define DPB_ORDERMULT_OUTDELTA 1

#define DPB_ORDER_MAXNROOTS 100


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
   int                   newadds,            /**< number of additions */
   DPBORDER*             dpborder            /**< border */
   )
{
   const int size = dpborder->global_partstarts[dpborder->global_npartitions];
   assert(size >= 0);
   assert(newadds > 0);

   if( newadds + size > dpborder->global_partcap )
   {
      while( newadds + size > dpborder->global_partcap )
         dpborder->global_partcap *= DPBORDER_GROWTH_FACTOR;

      SCIPdebugMessage("reallocating memory (to %d) \n", dpborder->global_partcap);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(dpborder->global_partitions), dpborder->global_partcap) );

      hashmap_updateKeyarr(dpborder->global_partitions, &dpborder->hashmap);
   }

   return SCIP_OKAY;
}


/** gets partition range */
static inline
void partitionGetRangeGlobal(
   const DPBORDER*       dpborder,            /**< border */
   int                   globalposition,      /**< position of partition */
   int*                  start,               /**< pointer to range start (OUT) */
   int*                  end                  /**< pointer to range end (OUT) */
   )
{
   const int* const starts = dpborder->global_partstarts;
   assert(globalposition >= 0);
   assert(globalposition < dpborder->global_npartitions);

   *start = starts[globalposition];
   *end = starts[globalposition + 1];

   assert(*start < *end);
}



/** builds map between old and new border char representation */
static
void borderBuildCharMap(
   int                   iteration,          /**< iteration number */
   DPBORDER*             dpborder            /**< border */
)
{
   const int extnode = dpborder_getTopLevel(dpborder)->extnode;
   const SCIP_Bool* const nodes_isBorder = dpborder->nodes_isBorder;
   int* RESTRICT bordercharmap = dpborder->bordercharmap;
   int nbordernew = 0;

   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];

      if( nodes_isBorder[bordernode] )
         bordercharmap[i] = nbordernew++;
      else
         bordercharmap[i] = -1;
   }

   if( dpborder->nodes_outdeg[extnode] != 0 || (iteration == dpborder->nnodes - 1) )
   {
      nbordernew++;
   }

   /* now we set the delimiter */
   bordercharmap[StpVecGetSize(dpborder->bordernodes)] = nbordernew;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("char border map, old to new: \n");
   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      SCIPdebugMessage("%d->%d \n", i, bordercharmap[i]);
   }
   SCIPdebugMessage("delimiter: %d->%d \n", StpVecGetSize(dpborder->bordernodes),
         bordercharmap[StpVecGetSize(dpborder->bordernodes)]);

#endif
}

/** builds distances to extension node */
static
void borderBuildCharDists(
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   SCIP_Real* RESTRICT borderchardists = dpborder->borderchardists;
   const SCIP_Bool* const nodes_isBorder = dpborder->nodes_isBorder;
   const int* const bordernodes = dpborder->bordernodes;
   const CSR* const csr = graph->csr_storage;
   const int* const start_csr = csr->start;
   const int extnode = dpborder_getTopLevel(dpborder)->extnode;
   const int nbordernodes = StpVecGetSize(bordernodes);

   SCIPdebugMessage("setting up border distances for extnode=%d \n", extnode);

   for( int i = 0; i < nbordernodes; i++ )
      borderchardists[i] = FARAWAY;

#ifndef NDEBUG
   for( int i = nbordernodes; i < BPBORDER_MAXBORDERSIZE; i++ )
      borderchardists[i] = -BLOCKED;
#endif

   for( int e = start_csr[extnode]; e != start_csr[extnode + 1]; e++ )
   {
      const int head = csr->head[e];
      if( nodes_isBorder[head] )
      {
         int i;
         for( i = 0; i < BPBORDER_MAXBORDERSIZE; i++ )
         {
            const int bordernode = bordernodes[i];
            assert(bordernode != extnode);

            if( head == bordernode )
            {
               SCIPdebugMessage("setting edge distance for border node %d to %f \n", head, csr->cost[e]);
               assert(EQ(borderchardists[i], FARAWAY));
               borderchardists[i] = csr->cost[e];
               break;
            }
         }
         assert(i != BPBORDER_MAXBORDERSIZE);
      }
   }
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
   borderBuildCharDists(graph, dpborder);

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

#ifdef SCIP_DEBUG
   printf("iter=%d bordersize=%d\n", iteration, StpVecGetSize(dpborder->bordernodes));
   printf("global_npartitions=%d \n", dpborder->global_npartitions);
#endif

   borderBuildCharMap(iteration, dpborder);
   assert(!toplevel->bordernodesMapToOrg);

   /* remove outdated border nodes */
   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];
      assert(graph_knot_isInRange(graph, bordernode));

      if( nodes_isBorder[bordernode] )
      {
         dpborder->bordernodes[nbordernodes++] = dpborder->bordernodes[i];
         StpVecPushBack(scip, toplevel->bordernodesMapToOrg, bordernode);
      }
   }

   for( int i = StpVecGetSize(dpborder->bordernodes) - nbordernodes; i > 0; i-- )
      StpVecPopBack(dpborder->bordernodes);

   if( (nodes_outdegree[extnode] != 0) || (iteration == dpborder->nnodes - 1) )
   {
      dpborder->extborderchar = nbordernodes;
      nbordernodes++;
      StpVecPushBack(scip, dpborder->bordernodes, extnode);
      StpVecPushBack(scip, toplevel->bordernodesMapToOrg, extnode);
      nodes_isBorder[extnode] = TRUE;
   }
   else
   {
      dpborder->extborderchar = -1;
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
   const int delimiter_new = dpborder_getTopDelimiter(dpborder);

   partitionGetRangeGlobal(dpborder, globalposition, &part_start, &part_end);
   partition.partchars = &(global_partitions[part_start]);
   partition.partsize = (part_end - part_start);
   partition.delimiter = dpborder_getDelimiter(dpborder, StpVecGetSize(dpborder->borderlevels) - 2);

   candstarts = dpborder_partGetCandstarts(scip, &partition, dpborder);
   ncands = StpVecGetSize(candstarts);
   assert(ncands >= 0);
   powsize = (uint32_t) 1 << (uint32_t) ncands;
   assert(powsize == (uint32_t) pow(2.0, ncands));

   SCIP_CALL( SCIPallocBufferArray(scip, &subbuffer, BPBORDER_MAXBORDERSIZE) );

   /* make sure that partition storage is large enough */
   SCIP_CALL( partitionTryRealloc(scip, (powsize + 2) * BPBORDER_MAXBORDERSIZE * 2, dpborder) );
   global_partitions = dpborder->global_partitions;
   partition.partchars = &(global_partitions[part_start]);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("partition ncands=%d, isTerm=%d ...base partition: \n", ncands, toplevel->exnodeIsTerm);
   dpborder_partPrint(&partition);
#endif

   if( !toplevel->exnodeIsTerm )
   {
      /* try to add/update the partition that does not include the extension node */
      const int globalposition_new = dpborder_partGetIdxNewExclusive(scip, &partition, dpborder);
      if( globalposition_new != -1 && LT(part_cost, dpborder->global_partcosts[globalposition_new]) )
      {
         SCIPdebugMessage("...added partition at %d with cost=%f \n", globalposition_new, part_cost);
         dpborder->global_partcosts[globalposition_new] = part_cost;
         dpborder->global_predparts[globalposition_new] = globalposition;
         dpborder->global_partsUseExt[globalposition_new] = FALSE;
      }
   }

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
      cost_new = dpborder_partGetConnectionCost(dpborder, &partition, subbuffer, nsub);
      SCIPdebugMessage("connection cost: %f \n", cost_new);

      cost_new += part_cost;

      if( GT(dpborder->global_partcosts[globalposition_new], cost_new) )
      {
         assert(LT(cost_new, FARAWAY));
         SCIPdebugMessage("...added partition at %d with cost=%f \n", globalposition_new, cost_new);
         dpborder->global_partcosts[globalposition_new] = cost_new;
         dpborder->global_partsUseExt[globalposition_new] = TRUE;
         dpborder->global_predparts[globalposition_new] = globalposition;

         if( allTermsAreVisited )
         {
            if( 1 == dpborder_partglobalGetCard(globalposition_new, delimiter_new, dpborder) )
            {
               if( LT(cost_new, dpborder->global_obj) )
               {
                  dpborder->global_obj = cost_new;
                  dpborder->global_optposition = globalposition_new;
                  SCIPdebugMessage("updated global obj to %f \n", dpborder->global_obj);
               }
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

   if( iteration > 1 )
   {
      // todo extra method
      DPBHASHMAP* hashmap = &(dpborder->hashmap);
      const STP_Vectype(int) global_partstarts = dpborder->global_partstarts;
      const int levelstart = dpborder_getTopLevel(dpborder)->globalstartidx;
      const int levelend = dpborder->global_npartitions;

      assert(levelstart < levelend);

      for( int i = levelstart; i != levelend; i++ )
      {
         const int partstart = global_partstarts[i];
         const int partend = global_partstarts[i + 1];

         assert(partstart < partend);

#ifdef NDEBUG
         (void) hashmap_remove(hashmap, partstart, partend - partstart);
#else
         {
            int ret = hashmap_remove(hashmap, partstart, partend - partstart);
            assert(ret == 0);
         }
#endif
      }

      assert(hashmap_isEmpty(&dpborder->hashmap));
   }

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
   StpVecPushBack(scip, level->bordernodesMapToOrg, root);

   dpborder->global_partitions[0] = 0;
   StpVecPushBack(scip, dpborder->global_partstarts, 1);
   StpVecPushBack(scip, dpborder->global_partcosts, 0.0);
   StpVecPushBack(scip, dpborder->global_predparts, 0);
   StpVecPushBack(scip, dpborder->global_partsUseExt, TRUE);
   dpborder->global_npartitions = 1;

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

   dpborder->global_partcap = MAX(DPBORDER_GROWTH_FACTOR * maxnpartitions, BPBORDER_MAXBORDERSIZE);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpborder->global_partitions), dpborder->global_partcap) );

   hashmap_updateKeyarr(dpborder->global_partitions, &dpborder->hashmap);

   SCIP_CALL( addLevelFirst(scip, graph, dpborder) );

   return SCIP_OKAY;
}

// #define DPB_ORDERING_BFS
#ifdef DPB_ORDERING_BFS
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
#else
/** computes node ordering and updates if better */
static
SCIP_RETCODE computeOrderingFromNode(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   root,               /**< node to start from */
   DPBSEQUENCE*          dpbsequence         /**< sequence */
)
{
   STP_Vectype(int) frontier = NULL;
   SCIP_Bool* RESTRICT nodes_isVisited;
   SCIP_Bool* RESTRICT nodes_isFrontier;
   int* RESTRICT nodessquence = dpbsequence->nodessquence;
   const CSR* const csr = graph->csr_storage;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   int* RESTRICT nodes_outdegree;
   const int* const nodes_degree = graph->grad;
   const int nnodes = graph_get_nNodes(graph);
   int bordersize;
   int bordertermsize;

   assert(graph_knot_isInRange(graph, root));
   assert(Is_term(graph->term[root]));

   StpVecReserve(scip, frontier, nnodes);
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodes_isVisited, nnodes) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodes_isFrontier, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_outdegree, nnodes) );

   BMScopyMemoryArray(nodes_outdegree, nodes_degree, nnodes);

   dpbsequence->maxbordersize = 0;
   dpbsequence->maxnpartitions = 0;
   nodessquence[0] = root;
   nodes_isVisited[root] = TRUE;
   bordersize = 1;
   bordertermsize = 1;

   for( int e = start_csr[root]; e != start_csr[root + 1]; e++ )
   {
      const int head = head_csr[e];
      nodes_outdegree[head]--;

      StpVecPushBack(scip, frontier, head);
      nodes_isFrontier[head] = TRUE;
   }

   for( int s = 1; s < nnodes; s++ )
   {
      int priority_best = -nnodes * DPB_ORDERMULT_OUTDELTA;
      int node_best = -1;
      const int frontier_size = StpVecGetSize(frontier);

      assert(frontier_size > 0);

      for( int i = 0; i < frontier_size; i++ )
      {
         const int fnode = frontier[i];
         int priority_new = 0;
         int nprevs = 0;
         const int outdegree = nodes_outdegree[fnode];
         assert(!nodes_isVisited[fnode] && nodes_isFrontier[fnode]);

         /* count vertices that would be removed from border */
         for( int e = start_csr[fnode]; e != start_csr[fnode + 1]; e++ )
         {
            const int head = head_csr[e];
            assert(nodes_outdegree[head] >= 1);

            if( !nodes_isVisited[head] )
               continue;

            if( nodes_outdegree[head] == 1 )
               nprevs++;
         }
         assert(outdegree >= 0 && nodes_degree[fnode] >= outdegree);
         SCIPdebugMessage("(old) frontier node %d: outdegree=%d outdelta=%d nprevs=%d \n",
               fnode, outdegree, outdegree - nprevs, nprevs);

         priority_new += DPB_ORDERMULT_PREVS * nprevs;

         if( Is_term(graph->term[fnode]) )
            priority_new += DPB_ORDERMULT_TERM * 1;

         /* number of visited nodes in adjacency list of fnode */
         priority_new += DPB_ORDERMULT_OUTDEG * (nodes_degree[fnode] - outdegree);

         /* (negative) delta of number of edges going out of visited region */
         priority_new += DPB_ORDERMULT_OUTDELTA * (nprevs - outdegree);

         if( priority_new > priority_best )
         {
            node_best = fnode;
            priority_best = priority_new;
         }
      }
      assert(node_best >= 0);

      /* update frontier */
      for( int e = start_csr[node_best]; e != start_csr[node_best + 1]; e++ )
      {
         const int head = head_csr[e];
         nodes_outdegree[head]--;
         assert(nodes_outdegree[head] >= 0);

         if( nodes_isVisited[head] )
         {
            if( nodes_outdegree[head] == 0 )
            {
               if( Is_term(graph->term[head]) )
                  bordertermsize--;
               bordersize--;
            }
         }
         else if( !nodes_isFrontier[head] )
         {
            nodes_isFrontier[head] = TRUE;
            StpVecPushBack(scip, frontier, head);
         }
      }
      nodes_isFrontier[node_best] = FALSE;
      nodes_isVisited[node_best] = TRUE;
      nodessquence[s] = node_best;

      assert(StpVecGetSize(frontier) > 0);
      {
         int i;
         for( i = 0; i < StpVecGetSize(frontier); i++ )
         {
            if( frontier[i] == node_best )
               break;
         }
         assert(i < StpVecGetSize(frontier));
         SWAP_INTS(frontier[i], frontier[StpVecGetSize(frontier) - 1]);
         StpVecPopBack(frontier);
      }

      if( nodes_outdegree[node_best] != 0 )
      {
         if( Is_term(graph->term[node_best]) )
            bordertermsize++;
         bordersize++;
      }

      assert(bordersize >= bordertermsize);
      dpbsequence->maxnpartitions += bordersize * bordersize * (bordersize - bordertermsize);

      SCIPdebugMessage("ADDED node %d to border \n", node_best);
      SCIPdebugMessage("new size of border: %d (termsize=%d) \n", bordersize, bordertermsize);
      SCIPdebugMessage("new size of frontier: %d \n", StpVecGetSize(frontier));
      SCIPdebugMessage("new total number of partitions: %ld \n", dpbsequence->maxnpartitions);

      if( bordersize > dpbsequence->maxbordersize )
      {
         dpbsequence->maxbordersize = bordersize;
         if( bordersize >= BPBORDER_MAXBORDERSIZE )
         {
            SCIPdebugMessage("aborting early, border too large! \n");
            nodessquence[0] = -1;
            break;
         }
      }

      if( dpbsequence->maxnpartitions >= BPBORDER_MAXNPARTITIONS )
      {
         SCIPdebugMessage("aborting early, too many partitions! \n");
         nodessquence[0] = -1;
         break;
      }
   }

   StpVecFree(scip, frontier);
   SCIPfreeBufferArray(scip, &nodes_outdegree);
   SCIPfreeBufferArray(scip, &nodes_isFrontier);
   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}
#endif


/*
 * Interface methods
 */


// compute ordering


/** computes node ordering */
SCIP_RETCODE dpborder_coreComputeOrderingSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   SCIP_CALL( computeOrderingFromNode(scip, graph, graph->source, dpborder->dpbsequence) );

   return SCIP_OKAY;
}


/** updates given ordering */
SCIP_RETCODE dpborder_coreUpdateOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   DPBSEQUENCE* const sequence_base = dpborder->dpbsequence;
   DPBSEQUENCE* sequence_tmp;
   int* terms;
   const int nterms = graph->terms;
   const int root_org = sequence_base->nodessquence[0];

   assert(nterms > 1);
   assert(graph_knot_isInRange(graph, root_org));
   assert(Is_term(graph->term[root_org]));

   SCIP_CALL( dpborder_dpbsequenceInit(scip, graph, &sequence_tmp) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &terms, nterms) );
   SCIP_CALL( graph_getTermsRandom(scip, graph, terms) );

   for( int i = 0; i < MIN(nterms, DPB_ORDER_MAXNROOTS); i++ )
   {
      SCIP_Bool isBetter;
      const int root = terms[i];
      if( root == root_org )
         continue;

      SCIP_CALL( computeOrderingFromNode(scip, graph, root, sequence_tmp) );

      if( sequence_tmp->maxbordersize >= BPBORDER_MAXBORDERSIZE )
         continue;

      if( sequence_tmp->maxnpartitions >= BPBORDER_MAXNPARTITIONS )
         continue;

      isBetter = sequence_tmp->maxnpartitions < sequence_base->maxnpartitions
              && sequence_tmp->maxbordersize <= sequence_base->maxbordersize;

      // todo use parameter
      if( !isBetter )
         isBetter = (1.2 * (SCIP_Real) sequence_tmp->maxnpartitions < (SCIP_Real) sequence_base->maxnpartitions);

      if( isBetter )
      {
         SCIPdebugMessage("updating ordering: \n" );
         SCIPdebugMessage("---max n. partitions %ld->%ld \n", sequence_base->maxnpartitions, sequence_tmp->maxnpartitions);
         SCIPdebugMessage("---max bordersize %d->%d \n", sequence_base->maxbordersize, sequence_tmp->maxbordersize);

         dpborder_dpbsequenceCopy(sequence_tmp, sequence_base);
      }
   }

   /* NOTE: might happen that the partitions are actually 0 if only terminals are all borders */
   sequence_base->maxnpartitions = MAX(sequence_base->maxnpartitions, 2);

   SCIPdebugMessage("final ordering: \n" );
   SCIPdebugMessage("---max n. partitions %ld \n", sequence_base->maxnpartitions);
   SCIPdebugMessage("---max bordersize %d \n", sequence_base->maxbordersize);

   SCIPfreeMemory(scip, &terms);
   dpborder_dpbsequenceFree(scip, &sequence_tmp);

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
      assert(dpborder->global_optposition != -1);

      printf("solved with  partitions=%d, max. bordersize=%d \n", dpborder->global_npartitions, dpborder->dpbsequence->maxbordersize);
   }

   return SCIP_OKAY;
}
