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

/**@file   extreduce_data.c
 * @brief  plain data storages for extended reduction techniques for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements several (relatively) plain data storages needed for extended reduction techniques.
 * Basically data transfer objects with additional setup, free, and cleaning methods.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
// #define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"





/** helper function */
static inline
void postCleanSDs(
   MLDISTS*              sds                 /**< distance structure */
)
{
   const int nlevels = extreduce_mldistsNlevels(sds);

   assert(nlevels == 2 || nlevels == 1);

   extreduce_mldistsLevelRemoveTop(sds);

   if( nlevels == 2 )
      extreduce_mldistsLevelRemoveTop(sds);
}


/** helper function */
static inline
void postCleanMSTs(
    CSRDEPO*              msts                /**< CSR depository containing MSTs */
)
{
#ifndef NDEBUG
   const int nmsts = graph_csrdepo_getNcsrs(msts);

   assert(nmsts == 1 || nmsts == 2);
#endif

   graph_csrdepo_clean(msts);
}


/** cleans-up after trying to rule out a component */
void extreduce_extCompClean(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be cleaned for */
   SCIP_Bool             unhash,             /**< unhash component? */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   MLDISTS* const sds_vertical = reddata->sds_vertical;
   MLDISTS* const sds_horizontal = reddata->sds_horizontal;
   CSRDEPO* const msts_comp = reddata->msts_comp;
   CSRDEPO* const msts_levelbase = reddata->msts_levelbase;
   const int* const compedges = extcomp->compedges;
   const int ncompedges = extcomp->ncompedges;
   int* const tree_deg = extdata->tree_deg;
   int* const pseudoancestor_mark = extdata->reddata->pseudoancestor_mark;
   SCIP_Bool* const sdeq_edgesIsForbidden = extdata->sdeq_edgesIsForbidden;
   const STP_Vectype(int) sdeq_resetStack = extdata->sdeq_resetStack;
   const int sdeq_size = StpVecGetSize(sdeq_resetStack);

   assert(ncompedges >= 1);
   assert(sdeq_size >= 0);

   for( int i = 0; i < ncompedges; ++i )
   {
      const int edge = compedges[i];
      const int head = graph->head[edge];
      const int tail = graph->tail[edge];

      assert(graph_edge_isInRange(graph, edge));

      tree_deg[head] = 0;
      tree_deg[tail] = 0;

      if( unhash )
      {
         graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, edge, pseudoancestor_mark);
      }
   }

   if( unhash && extInitialCompIsGenStar(extdata) )
   {
      const int centeredge = extcomp->genstar_centeredge;
      assert(graph_edge_isInRange(graph, centeredge));

      graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, centeredge, pseudoancestor_mark);
   }

   for( int i = 0; i < sdeq_size; i++ )
   {
      const int edge = sdeq_resetStack[i];

      assert(graph_edge_isInRange(graph, 2 * edge));
      assert(sdeq_edgesIsForbidden[edge]);

      sdeq_edgesIsForbidden[edge] = FALSE;
   }
   extdata->sdeq_hasForbiddenEdges = FALSE;

   StpVecFree(scip, extdata->sdeq_resetStack);

   postCleanSDs(sds_vertical);
   postCleanSDs(sds_horizontal);

   if( extReddataHasBiasedSds(reddata) )
   {
      postCleanSDs(reddata->sdsbias_vertical);
      postCleanSDs(reddata->sdsbias_horizontal);
   }

   postCleanMSTs(msts_comp);
   postCleanMSTs(msts_levelbase);

   extreduce_extdataClean(extdata);
   extreduce_reddataClean(extdata->reddata);
   extreduce_pcdataClean(extdata->pcdata);

   assert(extreduce_extdataIsClean(graph, extdata));
   assert(extreduce_reddataIsClean(graph, extdata->reddata));
   assert(extreduce_pcdataIsClean(graph, extdata->pcdata));
}


/** initialize permanent extension data struct
 *  NOTE: Sets distdata and reddata entries to NULL, since non-owned */
SCIP_RETCODE extreduce_extPermaInit(
   SCIP*                 scip,               /**< SCIP */
   enum EXTRED_MODE      mode,               /**< mode */
   const GRAPH*          graph,              /**< graph data structure */
   STP_Bool*             edgedeleted,        /**< edge array to mark which directed edge can be removed */
   EXTPERMA**            extpermanent        /**< (uninitialized) extension data */
)
{
   EXTPERMA* extperm;
   STP_Vectype(int)* nodes_implied = NULL;
   SCIP_Bool* isterm = NULL;
   SCIP_Real* bottleneckDistNode = NULL;
   SCIP_Real* pcSdToNode = NULL;
   int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);
   const int msts_datasize = STP_EXT_MAXDFSDEPTH * STP_EXTTREE_MAXNLEAVES_GUARD * 2;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   const int msts_maxn = STP_EXT_MAXDFSDEPTH_GUARD;

#ifndef NDEBUG
   const SCIP_Bool sds_vertical_useids = TRUE;
#else
   const SCIP_Bool sds_vertical_useids = FALSE;
#endif

   assert(scip && extpermanent);
   assert(mode == extred_fast || mode == extred_full);

   SCIP_CALL( SCIPallocMemory(scip, extpermanent) );
   extperm = *extpermanent;

   SCIP_CALL( SCIPallocMemoryArray(scip, &nodes_implied, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bottleneckDistNode, nnodes) );

   if( pcmw )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &pcSdToNode, nnodes) );
   }

   reduce_impliedNodesGet(scip, graph, nodes_implied);
   SCIP_CALL( reduce_dcmstInit(scip, STP_EXTTREE_MAXNLEAVES_GUARD, &(extperm->dcmst)) );
   SCIP_CALL( graph_csrdepo_init(scip, msts_maxn, msts_datasize, &(extperm->msts_comp)) );
   SCIP_CALL( graph_csrdepo_init(scip, msts_maxn, msts_datasize, &(extperm->msts_levelbase)) );

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXTTREE_MAXNLEAVES_GUARD, 1, sds_vertical_useids, &(extperm->sds_vertical)) );

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXT_MAXGRAD, 1, TRUE, &(extperm->sds_horizontal)) );

   SCIP_CALL( extreduce_contractionInit(scip, STP_EXT_MAXDFSDEPTH_GUARD,
         STP_EXTTREE_MAXNLEAVES_GUARD, &(extperm->contration)) );

   extperm->randnumgen = NULL;
   extperm->sdsbias_horizontal = NULL;
   extperm->sdsbias_vertical = NULL;
   extperm->distdata_default = NULL;
   extperm->distdata_biased = NULL;
   extperm->redcostdata = NULL;
   extperm->solIsValid = FALSE;

   extperm->result = NULL;
   extperm->nodes_implications = nodes_implied;
   extperm->edgedeleted = edgedeleted;
   extperm->isterm = isterm;
   extperm->bottleneckDistNode = bottleneckDistNode;
   extperm->pcSdToNode = pcSdToNode;
   extperm->tree_deg = tree_deg;
   extperm->nnodes = nnodes;
   extperm->redcostEqualAllow = FALSE;
   extperm->useSdBias = FALSE;
   extperm->mode = mode;
   extperm->tree_maxdepth = extreduce_getMaxTreeDepth(graph, extperm);
   extperm->tree_maxnedges = STP_EXTTREE_MAXNEDGES;
   extperm->tree_maxnleaves = STP_EXTTREE_MAXNLEAVES;

   if( pcmw )
   {
      assert(pcSdToNode);

      for( int k = 0; k < nnodes; k++ )
         pcSdToNode[k] = -1.0;
   }

   graph_getIsTermArray(graph, isterm);

   for( int k = 0; k < nnodes; k++ )
   {
      bottleneckDistNode[k] = -1.0;

      if( graph->mark[k] )
         tree_deg[k] = 0;
      else
         tree_deg[k] = -1;
   }

   assert(extreduce_extPermaIsClean(graph, extperm));

   return SCIP_OKAY;
}


/** adds random number generator */
void extreduce_extPermaAddRandnumgen(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator to add (NON-OWNED!) */
   EXTPERMA*             extpermanent        /**< (initialized) extension data */
)
{
   assert(extpermanent && randnumgen);

   extpermanent->randnumgen = randnumgen;
}


/** adds biased ML distances containers */
SCIP_RETCODE extreduce_extPermaAddMLdistsbiased(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA*             extpermanent        /**< (initialized) extension data */
)
{
   const int msts_maxn = STP_EXT_MAXDFSDEPTH_GUARD;
#ifndef NDEBUG
   const SCIP_Bool sds_vertical_useids = TRUE;
#else
   const SCIP_Bool sds_vertical_useids = FALSE;
#endif

   assert(scip && extpermanent);
   assert(!extpermanent->sdsbias_vertical && !extpermanent->sdsbias_horizontal);

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXTTREE_MAXNLEAVES_GUARD, 1, sds_vertical_useids, &(extpermanent->sdsbias_vertical)) );

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXT_MAXGRAD, 1, TRUE, &(extpermanent->sdsbias_horizontal)) );

   return SCIP_OKAY;
}


/** initialize permanent extension data struct */
SCIP_Bool extreduce_extPermaIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperm             /**< extension data */
)
{
   const SCIP_Real* bottleneckDistNode = NULL;
   const SCIP_Real* pcSdToNode = NULL;
   const int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);

   assert(extperm);

   if( extperm->edgedeleted )
   {
      SCIPdebugMessage("use of edge-deleted array is deprecated! \n");
      return FALSE;
   }

   bottleneckDistNode = extperm->bottleneckDistNode;
   pcSdToNode = extperm->pcSdToNode;
   tree_deg = extperm->tree_deg;

   if( nnodes != extperm->nnodes )
   {
      SCIPdebugMessage("inconsistent number of nodes \n");
      return FALSE;
   }

   if( !graph_csrdepo_isEmpty(extperm->msts_comp) )
   {
      SCIPdebugMessage("msts_comp not empty! \n");
      return FALSE;
   }

   if( !graph_csrdepo_isEmpty(extperm->msts_levelbase) )
   {
      SCIPdebugMessage("msts_reduced not empty! \n");
      return FALSE;
   }

   if( !extreduce_mldistsIsEmpty(extperm->sds_vertical) )
   {
      SCIPdebugMessage("sds_vertical not empty! size=%d \n", extreduce_mldistsNlevels(extperm->sds_vertical));
      return FALSE;
   }

   if( !extreduce_mldistsIsEmpty(extperm->sds_horizontal) )
   {
      SCIPdebugMessage("sds_horizontal not empty! \n");
      return FALSE;
   }

   if( extperm->sdsbias_vertical && !extreduce_mldistsIsEmpty(extperm->sdsbias_vertical) )
   {
      SCIPdebugMessage("sdsbias_vertical not empty! size=%d \n", extreduce_mldistsNlevels(extperm->sds_vertical));
      return FALSE;
   }

   if( extperm->sdsbias_horizontal && !extreduce_mldistsIsEmpty(extperm->sdsbias_horizontal) )
   {
      SCIPdebugMessage("sdsbias_horizontal not empty! \n");
      return FALSE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !(tree_deg[i] == 0 || tree_deg[i] == -1) )
         return FALSE;

      if( bottleneckDistNode[i] > -0.5 )
         return FALSE;

      if( pcSdToNode && !EQ(pcSdToNode[i], -1.0) )
         return FALSE;
   }

   return TRUE;
}


/** frees  extension data */
void extreduce_extPermaFree(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA**            extpermanent        /**< extension data */
)
{
   EXTPERMA* extperm;

   assert(scip && extpermanent);

   extperm = *extpermanent;

   if( extperm->sdsbias_vertical )
      extreduce_mldistsFree(scip, &(extperm->sdsbias_vertical));

   if( extperm->sdsbias_horizontal )
      extreduce_mldistsFree(scip, &(extperm->sdsbias_horizontal));

   extreduce_contractionFree(scip, &(extperm->contration));
   extreduce_mldistsFree(scip, &(extperm->sds_horizontal));
   extreduce_mldistsFree(scip, &(extperm->sds_vertical));
   graph_csrdepo_free(scip, &(extperm->msts_comp));
   graph_csrdepo_free(scip, &(extperm->msts_levelbase));
   reduce_dcmstFree(scip, &(extperm->dcmst));

   SCIPfreeMemoryArrayNull(scip, &(extperm->pcSdToNode));
   SCIPfreeMemoryArray(scip, &(extperm->bottleneckDistNode));
   SCIPfreeMemoryArray(scip, &(extperm->tree_deg));
   SCIPfreeMemoryArray(scip, &(extperm->isterm));

   for( int i = extperm->nnodes - 1; i >= 0; i-- )
   {
      StpVecFree(scip, extperm->nodes_implications[i]);
   }
   SCIPfreeMemoryArray(scip, &(extperm->nodes_implications));

   SCIPfreeMemory(scip, extpermanent);
}


/** cleans extension data */
void extreduce_extdataClean(
   EXTDATA*              extdata             /**< extension data */
)
{
   assert(extdata);

   extdata->extstack_ncomponents = 0;
   extdata->tree_nDelUpArcs = 0;
   extdata->tree_nleaves = 0;
   extdata->tree_ninnerNodes = 0;
   extdata->tree_nedges = 0;
   extdata->tree_depth = 0;
   extdata->tree_root = -1;
   extdata->tree_starcenter = -1;
   extdata->tree_cost = 0.0;
   extdata->ncostupdatestalls = 0;
}


/** cleans reduction data */
void extreduce_reddataClean(
   REDDATA*              reddata             /**< reduction data */
)
{
   const int redcost_nlevels = reddata->redcost_nlevels;
   SCIP_Real* const redcost_treecosts = reddata->redcost_treecosts;

   assert(redcost_nlevels >= 1);
   assert(reddata->contration);

   for( int i = 0; i < redcost_nlevels; i++ )
   {
      redcost_treecosts[i] = 0.0;
   }
}


/** cleans PC data */
void extreduce_pcdataClean(
   PCDATA*               pcdata             /**< PC data */
)
{
   assert(pcdata);

   pcdata->tree_innerPrize = 0.0;
}


/** is the extension data clean? */
SCIP_Bool extreduce_extdataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   if( extdata->extstack_ncomponents != 0 )
   {
      printf("extdata->extstack_ncomponents %d \n", extdata->extstack_ncomponents);
      return FALSE;
   }

   if( extdata->tree_nDelUpArcs != 0 )
   {
      printf("extdata->tree_nDelUpArcs %d \n", extdata->tree_nDelUpArcs);
      return FALSE;
   }

   if( !EQ(extdata->tree_cost, 0.0) )
   {
      printf("extdata->tree_cost %f \n", extdata->tree_cost);
      return FALSE;
   }

   if( extdata->ncostupdatestalls != 0 )
   {
      printf("extdata->ncostupdatestalls %d \n", extdata->ncostupdatestalls);
      return FALSE;
   }

   if( extdata->tree_nleaves != 0 )
   {
      printf("extdata->tree_nleaves %d \n", extdata->tree_nleaves);
      return FALSE;
   }

   if( extdata->tree_ninnerNodes != 0 )
   {
      printf("extdata->tree_ninnerNodes %d \n", extdata->tree_ninnerNodes);
      return FALSE;
   }

   if( extdata->tree_root != -1 )
   {
      printf("extdata->tree_root %d \n", extdata->tree_root);
      return FALSE;
   }

   if( extdata->tree_nedges != 0 )
   {
      printf("extdata->tree_nedges %d \n", extdata->tree_nedges);
      return FALSE;
   }

   if( extdata->tree_starcenter != -1 )
   {
      printf("extdata->tree_starcenter %d \n", extdata->tree_starcenter);
      return FALSE;
   }

   if( extdata->tree_depth != 0 )
   {
      printf("extdata->tree_depth %d \n", extdata->tree_depth);
      return FALSE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !(extdata->tree_deg[i] == 0 || extdata->tree_deg[i] == -1) )
      {
         printf("extdata->tree_deg[i] %d \n", extdata->tree_deg[i]);
         return FALSE;
      }

      if( !EQ(extdata->tree_bottleneckDistNode[i], -1.0) )
      {
         printf("extdata->bottleneckDistNode[i] %f \n", extdata->tree_bottleneckDistNode[i]);
         return FALSE;
      }
   }

   for( int i = 0; i < nedges / 2; ++i )
   {
      if( extdata->sdeq_edgesIsForbidden[i] )
      {
         printf("extdata->sdeq_edgesIsForbidden[%d] not reset \n", i);
         return FALSE;
      }
   }

   return TRUE;
}


/** is the reduction data clean? */
SCIP_Bool extreduce_reddataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const REDDATA*        reddata             /**< reduction data */
)
{
   const int nancestors = graph_pseudoAncestorsGetHashArraySize(graph->pseudoancestors);
   const int redcost_nlevels = reddata->redcost_nlevels;
   const SCIP_Real* const redcost_treecosts = reddata->redcost_treecosts;

   assert(redcost_nlevels >= 1);

   for( int i = 0; i < nancestors; i++ )
   {
      if( reddata->pseudoancestor_mark[i] != 0 )
      {
         printf("%d pseudoancestor_mark %d \n", i, reddata->pseudoancestor_mark[i]);
         return FALSE;
      }
   }

   for( int i = 0; i < redcost_nlevels; i++ )
   {
      if( !EQ(redcost_treecosts[i], 0.0) )
      {
         printf("FAIL: tree redcost %f for level %d \n", redcost_treecosts[i], i);
         return FALSE;
      }
   }

   return TRUE;
}

/** is the reduction data clean? */
SCIP_Bool extreduce_pcdataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const PCDATA*         pcdata              /**< PC data */
)
{
   const int nnodes = graph_get_nNodes(graph);

   assert(pcdata);

   if( !graph_pc_isPc(graph) )
   {
      return TRUE;
   }

   if( pcdata->nPcSdCands != -1 )
   {
      return FALSE;
   }

   if( pcdata->pcSdStart != -1 )
   {
      return FALSE;
   }

   if( !EQ(pcdata->tree_innerPrize, 0.0) )
   {
      return FALSE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !EQ(pcdata->pcSdToNode[i], -1.0) )
      {
         printf("pcSdToNode[%d]=%f \n", i, pcdata->pcSdToNode[i]);
         return FALSE;
      }
   }


   return TRUE;
}
