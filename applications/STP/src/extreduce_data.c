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



/** initialize permanent extension data struct */
SCIP_RETCODE extreduce_extPermaInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   STP_Bool*             edgedeleted,        /**< edge array to mark which directed edge can be removed */
   EXTPERMA*             extperm             /**< (uninitialized) extension data */
)
{
   SCIP_Bool* isterm = NULL;
   SCIP_Real* bottleneckDistNode = NULL;
   SCIP_Real* pcSdToNode = NULL;
   int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);
   const int msts_datasize = STP_EXT_MAXDFSDEPTH * STP_EXTTREE_MAXNLEAVES_GUARD * 2;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   const int msts_maxn = STP_EXT_MAXDFSDEPTH + 1;

#ifndef NDEBUG
   const SCIP_Bool sds_vertical_useids = TRUE;
#else
   const SCIP_Bool sds_vertical_useids = FALSE;
#endif

   assert(scip && extperm);

   SCIP_CALL( SCIPallocMemoryArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bottleneckDistNode, nnodes) );

   if( pcmw )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &pcSdToNode, nnodes) );
   }

   SCIP_CALL( reduce_dcmstInit(scip, STP_EXTTREE_MAXNLEAVES_GUARD, &(extperm->dcmst)) );
   SCIP_CALL( graph_csrdepo_init(scip, msts_maxn, msts_datasize, &(extperm->msts_comp)) );
   SCIP_CALL( graph_csrdepo_init(scip, msts_maxn, msts_datasize, &(extperm->msts_levelbase)) );

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXTTREE_MAXNLEAVES_GUARD, 1, sds_vertical_useids, &(extperm->sds_vertical)) );

   SCIP_CALL( extreduce_mldistsInit(scip, msts_maxn, STP_EXT_MAXGRAD,
         STP_EXT_MAXGRAD, 1, TRUE, &(extperm->sds_horizontal)) );

   extperm->edgedeleted = edgedeleted;
   extperm->isterm = isterm;
   extperm->bottleneckDistNode = bottleneckDistNode;
   extperm->pcSdToNode = pcSdToNode;
   extperm->tree_deg = tree_deg;
   extperm->nnodes = nnodes;

   if( pcmw )
   {
      assert(pcSdToNode);

      for( int k = 0; k < nnodes; k++ )
         pcSdToNode[k] = -1.0;
   }

   graph_get_isTerm(graph, isterm);

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


/** frees members of extension data */
void extreduce_extPermaFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA*             extperm             /**< extension data */
)
{
   assert(scip && extperm);

   extreduce_mldistsFree(scip, &(extperm->sds_horizontal));
   extreduce_mldistsFree(scip, &(extperm->sds_vertical));
   graph_csrdepo_free(scip, &(extperm->msts_comp));
   graph_csrdepo_free(scip, &(extperm->msts_levelbase));
   reduce_dcmstFree(scip, &(extperm->dcmst));

   SCIPfreeMemoryArrayNull(scip, &(extperm->pcSdToNode));
   SCIPfreeMemoryArray(scip, &(extperm->bottleneckDistNode));
   SCIPfreeMemoryArray(scip, &(extperm->tree_deg));
   SCIPfreeMemoryArray(scip, &(extperm->isterm));

   extperm->nnodes = -1;
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
   extdata->tree_nedges = 0;
   extdata->tree_depth = 0;
   extdata->tree_root = -1;
   extdata->tree_redcost = 0.0;
   extdata->tree_cost = 0.0;
   extdata->ncostupdatestalls = 0;
}


/** cleans reduction data */
void extreduce_reddataClean(
   REDDATA*              reddata             /**< reduction data */
)
{
   assert(reddata);
}


/** is the extension data clean? */
SCIP_Bool extreduce_extdataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
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

   if( !EQ(extdata->tree_redcost, 0.0) )
   {
      printf("extdata->tree_redcost %f \n", extdata->tree_redcost);
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

   if( extdata->tree_depth != 0 )
   {
      printf("extdata->tree_depth %d \n", extdata->tree_depth);
      return FALSE;
   }

   for( int i = 0; i < graph->knots; i++ )
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

   return TRUE;
}


/** is the reduction data clean? */
SCIP_Bool extreduce_reddataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const REDDATA*        reddata             /**< reduction data */
)
{
   assert(graph && reddata);

   for( int i = 0; i < graph->knots; i++ )
   {
      if( reddata->pseudoancestor_mark[i] != 0 )
      {
         printf("pseudoancestor_mark %d \n", reddata->pseudoancestor_mark[i]);
         return FALSE;
      }
   }

   return TRUE;
}
