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

/**@file   extreduce_contract.c
 * @brief  extended-reduction tree contraction algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements rule-out algorithms for extended reduction techniques for Steiner problems.
 * Allows to compute and store special distance (SD) MSTs / SPGs between the leaves of extension tree,
 * after the contraction of certain parts of the tree.
 *
 * Similarly to the extmst methods, we keep two distance storages: One from all component vertices to the
 * lower leaves. One from the component root to the lower leafs and to the siblings.
 *
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
//#define STP_DEBUG_EXT

#include "extreduce.h"
#include "extreducedefs.h"



struct extension_tree_contraction
{
   CSR*                  mst_buffer;         /**< buffer that can keep at least leaves_maxn many nodes */
   SCIP_Real*            leafToCompRootDists;/**< stores distances from leaves to component root */
   SCIP_Real*            leafToCompUpDists;  /**< stores distances from leaves to all above components */
   SCIP_Real*            level_treecost;     /**< cost of tree per level */
   SCIP_Real*            leaves_mindists;    /**< minimum distances from leaves to contracted node (buffer, of size leaves_maxn) */
   int*                  leaves_start;       /**< start pointer per level */
   int                   leaves_maxn;        /**< maximum number of leaves */
   int                   level_maxn;         /**< maximum number of levels */
   int                   level_n;            /**< current number of levels */
#ifndef NDEBUG
   SCIP_Real*            leafToCompLeaves;   /**< stores leaves for leafToCompDists  */
#endif

};



/**@name Local methods
 *
 * @{
 */


/** gets position of distance value */
static inline
int compDistGetPosition(
   int                   leaf_id,            /**< id */
   int                   level,              /**< level */
   const CONTRACT*       contraction         /**< contraction data */
)
{
   const int maxnlevels = contraction->level_maxn;
   const int position = leaf_id * maxnlevels + level;

   assert(level >= 0);
   assert(0 <= leaf_id && leaf_id < contraction->leaves_maxn);

   assert(position >= 0);

   return position;
}


/** adds distance */
static inline
void compUpDistAddLeaf(
   int                   leaf_id,            /**< id */
   SCIP_Real             dist,               /**< distance */
   SCIP_Bool             override,           /**< override distance? (take in otherwise) */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const int level = contraction->level_n - 1;
   const int position = compDistGetPosition(leaf_id, level, contraction);

   assert(leaf_id < extdata->tree_nleaves);
   assert(extdata->tree_leaves[leaf_id] == contraction->leafToCompLeaves[leaf_id]);
   assert(GE(dist, 0.0));

   if( override )
   {
      contraction->leafToCompUpDists[position] = dist;
   }
   else if( LT(dist, contraction->leafToCompUpDists[position]) )
   {
      contraction->leafToCompUpDists[position] = dist;
   }

   assert(GE(contraction->leafToCompUpDists[position], 0.0));
}



/** updates distances from lower leafs to top component */
static inline
void compUpDistUpdateLeavesDists(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const int* const extstack_data = extdata->extstack_data;
   const int level = contraction->level_n;
   const int stackpos = extStackGetPosition(extdata);
   const int topedges_start = extStackGetTopOutEdgesStart(extdata, stackpos);
   const int topedges_end = extStackGetTopOutEdgesEnd(extdata, stackpos);
   const int leaves_end = contraction->leaves_start[level];

   assert(leaves_end >= 1);
   assert(extreduce_mldistsLevelNTopTargets(sds_vertical) == leaves_end);

   /* add the 'up' leaf to components distances */
   for( int i = topedges_start; i < topedges_end; i++ )
   {
      const SCIP_Bool override = (i == topedges_start);
      const int edge = extstack_data[i];
      const int topleaf = graph->head[edge];
      const SCIP_Real* const sdvertical_dists = extreduce_mldistsTopTargetDists(sds_vertical, topleaf);
#ifndef NDEBUG
      const int* const sdvertical_ids = extreduce_mldistsTopTargetIds(sds_vertical, topleaf);
#endif

      for( int j = 0; j < leaves_end; j++ )
      {
         assert(sdvertical_ids[j] == extdata->tree_leaves[j]);

         SCIPdebugMessage("top=%d leaf=%d, dist=%f \n", topleaf, extdata->tree_leaves[j], sdvertical_dists[j]);

         compUpDistAddLeaf(j, sdvertical_dists[j], override, extdata, contraction);
      }
   }
}


/** adds distance */
static inline
void compRootDistAddLeaf(
   int                   leaf_id,            /**< id */
   SCIP_Real             dist,               /**< distance */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const int level = contraction->level_n - 1;
   const int position = compDistGetPosition(leaf_id, level, contraction);

   assert(level >= 1);
   assert(leaf_id < extdata->tree_nleaves);
   assert(extdata->tree_leaves[leaf_id] == contraction->leafToCompLeaves[leaf_id]);
   assert(GE(dist, 0.0));

   contraction->leafToCompRootDists[position] = dist;
}


/** updates distances from lower leafs (w.r.t top component) to top component root */
static inline
void compRootDistsUpdateLeavesDists(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const MLDISTS* const sds_horizontal = extdata->reddata->sds_horizontal;
   const int level = contraction->level_n;
   const int leavesSame_start = contraction->leaves_start[level - 1];
   const int leavesSame_end = contraction->leaves_start[level];
   const int comproot = extStackGetTopRoot(graph, extdata);
   const SCIP_Bool hasSiblings = (leavesSame_start != leavesSame_end);
   const SCIP_Real* const sdvertical_dists = extreduce_mldistsTargetDists(sds_vertical, level - 1, comproot);
   const SCIP_Real* const sdhorizontal_dists = hasSiblings ? extreduce_mldistsTargetDists(sds_horizontal, level - 1, comproot) : NULL;
   const int* const sdhorizontal_ids = hasSiblings ? extreduce_mldistsTargetIds(sds_horizontal, level - 1, comproot) : NULL;
   const int sdhorizontal_ntargets = extreduce_mldistsLevelNTargets(sds_horizontal, level - 1);

#ifndef NDEBUG
   const int* const sdvertical_ids = extreduce_mldistsTargetIds(sds_vertical, level - 1, comproot);
#endif

   assert(level > 1);
   assert(leavesSame_end >= 1 && leavesSame_start >= 1);
   assert(extreduce_mldistsLevelNTargets(sds_vertical, level - 1) == leavesSame_start);
   assert(sdhorizontal_ntargets >= (leavesSame_end - leavesSame_start));

   SCIPdebugMessage("horizontal root component distances: \n");

   /* update leafs on same level as root */
   for( int i = leavesSame_start; i < leavesSame_end; i++ )
   {
      int j;
      assert(sdvertical_dists && sdhorizontal_dists);

      for( j = 0; j < sdhorizontal_ntargets; j++ )
      {
         if( sdhorizontal_ids[j] == extdata->tree_leaves[i] )
            break;
      }

      assert(j != sdhorizontal_ntargets);

      SCIPdebugMessage("comproot=%d leaf=%d, dist=%f \n", comproot, sdhorizontal_ids[j],
            sdhorizontal_dists[j]);

      compRootDistAddLeaf(i, sdhorizontal_dists[j], extdata, contraction);
   }

   SCIPdebugMessage("vertical root component distances: \n");

   /* update leafs below root */
   for( int i = 0; i < leavesSame_start; i++ )
   {
      assert(sdvertical_ids[i] == extdata->tree_leaves[i]);

      SCIPdebugMessage("comproot=%d leaf=%d, dist=%f \n", comproot, extdata->tree_leaves[i],
            sdvertical_dists[i]);

      compRootDistAddLeaf(i, sdvertical_dists[i], extdata, contraction);
   }
}


/** initializes distances from component to lower leaves */
static inline
void compUpDistInitMindists(
   int                   level,              /**< level from which to initialize */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const SCIP_Real* const leafToCompDists = contraction->leafToCompUpDists;
   SCIP_Real* RESTRICT leaves_mindists = contraction->leaves_mindists;
   const int leaves_end = contraction->leaves_start[level + 1];

   assert(0 <= level && level < contraction->level_n);
   assert(0 < leaves_end && leaves_end < extdata->tree_nleaves);

   for( int i = 0; i < leaves_end; i++ )
   {
      const int position = compDistGetPosition(i, level, contraction);
      leaves_mindists[i] = leafToCompDists[position];
   }

#ifndef NDEBUG
   for( int i = leaves_end; i < contraction->leaves_maxn; i++ )
      leaves_mindists[i] = -FARAWAY;
#endif
}


/** updates distances from component to lower leaves */
static inline
void compUpDistUpdateMindists(
   int                   level,              /**< level from which to initialize */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const SCIP_Real* const leafToCompDists = contraction->leafToCompUpDists;
   SCIP_Real* RESTRICT leaves_mindists = contraction->leaves_mindists;
   const int leaves_end = contraction->leaves_start[level + 1];

   assert(0 <= level && level < contraction->level_n);
   assert(0 < leaves_end && leaves_end < extdata->tree_nleaves);

   for( int i = 0; i < leaves_end; i++ )
   {
      const int position = compDistGetPosition(i, level, contraction);
      const SCIP_Real dist = leafToCompDists[position];
      assert(GE(dist, 0.0));

      if( LT(dist, leaves_mindists[i]) )
         leaves_mindists[i] = dist;
   }
}


/** updates leaf distances from component root to lower leaves */
static inline
void compRootDistUpdateMindists(
   int                   level,              /**< level from which to initialize */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const SCIP_Real* const leafToCompDists = contraction->leafToCompRootDists;
   SCIP_Real* RESTRICT leaves_mindists = contraction->leaves_mindists;
   const int leaves_end = contraction->leaves_start[level + 1];

   assert(1 <= level && level < contraction->level_n);
   assert(0 < leaves_end && leaves_end < extdata->tree_nleaves);

   for( int i = 0; i < leaves_end; i++ )
   {
      const int position = compDistGetPosition(i, level, contraction);
      const SCIP_Real dist = leafToCompDists[position];
      assert(GE(dist, 0.0));

      if( LT(dist, leaves_mindists[i]) )
         leaves_mindists[i] = dist;
   }
}


/** helper */
static inline
void addComponentUpdateTreeCosts(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   SCIP_Real treecost = extdata->tree_cost;
   const int* const extstack_data = extdata->extstack_data;
   const int stackpos = extStackGetPosition(extdata);
   const int topedges_start = extStackGetTopOutEdgesStart(extdata, stackpos);
   const int topedges_end = extStackGetTopOutEdgesEnd(extdata, stackpos);
   const int level = contraction->level_n;

   assert(level >= 1);

   for( int i = topedges_start; i < topedges_end; i++ )
   {
      const int edge = extstack_data[i];
      assert(graph_edge_isInRange(graph, edge));

      treecost -= graph->cost[edge];
   }

   if( extProbIsPc(graph, extdata) )
   {
      assert(extdata->pcdata);
      treecost -= extdata->pcdata->tree_innerPrize;
      assert(GE(treecost, 0.0));
   }

   SCIPdebugMessage("tree cost for level %d: %f \n", level - 1, treecost);
   contraction->level_treecost[level - 1] = treecost;
}


/** helper */
static inline
void addComponentUpdateLeavesStarts(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const int nleaves = extdata->tree_nleaves;
   const int stackpos = extStackGetPosition(extdata);
   const int topedges_start = extStackGetTopOutEdgesStart(extdata, stackpos);
   const int topedges_end = extStackGetTopOutEdgesEnd(extdata, stackpos);
   const int ntopedges = topedges_end - topedges_start;
   const int level = contraction->level_n;

   assert(ntopedges > 0);
   assert(ntopedges < nleaves);
   assert(level >= 1);

   contraction->leaves_start[level] = nleaves - ntopedges;

#ifndef NDEBUG
   for( int i = contraction->leaves_start[level - 1]; i < contraction->leaves_start[level]; i++ )
      contraction->leafToCompLeaves[i] = extdata->tree_leaves[i];
#endif

#ifdef SCIP_DEBUG
   {
      const int* const leaves = extdata->tree_leaves;

      SCIPdebugMessage("contraction lower leafs:  \n");
      for( int i = 0; i < contraction->leaves_start[level]; i++ )
         graph_knot_printInfo(graph, leaves[i]);
   }
#endif

}


/** helper */
static inline
void addComponentUpdateLeavesToCompDists(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   CONTRACT*             contraction         /**< contraction data */
)
{
   const int level = contraction->level_n;

   if( level > 1 )
   {
      compUpDistUpdateLeavesDists(graph, extdata, contraction);
      compRootDistsUpdateLeavesDists(graph, extdata, contraction);
   }
}


/** adds top component */
static
void addComponent(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
   )
{
   CONTRACT* const contraction = extdata->reddata->contration;

   SCIPdebugMessage("contraction: adding component for level %d \n", extdata->tree_depth);
   /* NOTE: will always be >= 1 */
   contraction->level_n = extdata->tree_depth;
   assert(contraction->level_n <= contraction->level_maxn);

   /* NOTE: calling order should not be changed */
   addComponentUpdateTreeCosts(graph, extdata, contraction);
   addComponentUpdateLeavesStarts(graph, extdata, contraction);
   addComponentUpdateLeavesToCompDists(graph, extdata, contraction);
}


/** can the current tree be ruled out? */
static
SCIP_Bool ruledOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
   )
{
   REDDATA* const reddata = extdata->reddata;
   CONTRACT* const contraction = reddata->contration;
   CSR* const mst_new = contraction->mst_buffer;
   DCMST* const dcmst = reddata->dcmst;
   CSRDEPO* const msts_levelbase = reddata->msts_levelbase;
   const int level = contraction->level_n;
   SCIP_Bool ruledOut = FALSE;
   const int mst_maxnedges = mst_new->nedges_max;
   const int mst_maxnnodes = mst_new->nnodes;

   for( int i = level - 1; i >= 1; i-- )
   {
      CSR mst_parent;
      SCIP_Real mstweight;
      const SCIP_Real tree_cost = contraction->level_treecost[i];
#ifdef SCIP_DEBUG
      const int leaves_end = contraction->leaves_start[i + 1];
#endif

      graph_csrdepo_getCSR(msts_levelbase, i + 1, &mst_parent);

      if( i == (level - 1) )
      {
         compUpDistInitMindists(i, extdata, contraction);
      }
      else
      {
         compUpDistUpdateMindists(i, extdata, contraction);
      }

      compRootDistUpdateMindists(i, extdata, contraction);

#ifdef SCIP_DEBUG
      printf("\n level=%d: \n", i);
      for( int j = 0; j < leaves_end; j++ )
      {
         printf("...distance from %d to contracted supernode: %f \n",
               extdata->tree_leaves[j], contraction->leaves_mindists[j]);
      }
#endif

      mst_new->nedges_max = mst_parent.nedges_max + 2;
      mst_new->nnodes = mst_parent.nnodes + 1;
      assert(mst_new->nedges_max <= mst_maxnedges && mst_new->nnodes <= mst_maxnnodes);

      SCIPdebugMessage("extending MST with n=%d, m=%d \n", mst_parent.nnodes ,mst_parent.nedges_max);

      reduce_dcmstAddNode(scip, &mst_parent, contraction->leaves_mindists, dcmst, mst_new);

      mstweight = reduce_dcmstGetWeight(scip, mst_new);

      SCIPdebugMessage("weigh of old %f \n", reduce_dcmstGetWeight(scip, &mst_parent));
      SCIPdebugMessage("weigh of new %f \n", mstweight);

      // todo also allow for equality if mst_new->nnodes > 3? */
      if( LT(mstweight, tree_cost) )
      {
         SCIPdebugMessage("ruled-out with %f < %f, contraction-level=%d \n", mstweight, tree_cost, i);
         ruledOut = TRUE;
         break;
      }

      if( mst_new->nnodes > 3 && EQ(mstweight, tree_cost) )
      {
         SCIPdebugMessage("ruled-out with equality %f <= %f, contraction-level=%d \n", mstweight, tree_cost, i);
         ruledOut = TRUE;
         break;
      }
   }

   mst_new->nedges_max = mst_maxnedges;
   mst_new->nnodes = mst_maxnnodes;

   return ruledOut;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */


/** initializes */
SCIP_RETCODE extreduce_contractionInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnlevels,         /**< maximum number of levels that can be handled */
   int                   maxnleaves,         /**< maximum number of leaves that can be handled */
   CONTRACT**            contraction         /**< to be initialized */
)
{
   CONTRACT* cont;

   assert(scip && contraction);
   assert(maxnlevels >= 1 && maxnleaves >= 1);

   SCIP_CALL( SCIPallocMemory(scip, contraction) );
   cont = *contraction;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->leaves_start), maxnlevels + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->leafToCompRootDists), maxnlevels * maxnleaves) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->leafToCompUpDists), maxnlevels * maxnleaves) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->level_treecost), maxnlevels) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->leaves_mindists), maxnleaves) );

   SCIP_CALL( graph_csr_alloc(scip, maxnleaves, 2 * (maxnleaves - 1), &(cont->mst_buffer)) );

   cont->leaves_start[0] = 0;

   cont->level_maxn = maxnlevels;
   cont->leaves_maxn = maxnleaves;

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cont->leafToCompLeaves), maxnleaves) );

   for( int i = 0; i < maxnlevels * maxnleaves; i++ )
   {
      cont->leafToCompRootDists[i] = -FARAWAY;
      cont->leafToCompUpDists[i] = -FARAWAY;
   }
#endif

   return SCIP_OKAY;
}


/** frees */
void extreduce_contractionFree(
   SCIP*                 scip,               /**< SCIP */
   CONTRACT**            contraction         /**< to be initialized */
)
{
   CONTRACT* cont;

   assert(scip && contraction);
   cont = *contraction;

#ifndef NDEBUG
   SCIPfreeMemoryArray(scip, &(cont->leafToCompLeaves));
#endif

   graph_csr_free(scip, &(cont->mst_buffer));

   SCIPfreeMemoryArray(scip, &(cont->leaves_mindists));
   SCIPfreeMemoryArray(scip, &(cont->level_treecost));
   SCIPfreeMemoryArray(scip, &(cont->leafToCompUpDists));
   SCIPfreeMemoryArray(scip, &(cont->leafToCompRootDists));
   SCIPfreeMemoryArray(scip, &(cont->leaves_start));

   SCIPfreeMemory(scip, contraction);
}



/** can current tree be peripherally ruled out by using contraction based arguments? */
SCIP_Bool extreduce_contractionRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   assert(scip && graph && extdata);
   assert(extdata->reddata && extdata->reddata->contration);

   addComponent(scip, graph, extdata);

   if( extdata->tree_nleaves <= 2 )
   {
      return FALSE;
   }

   if( ruledOut(scip, graph, extdata) )
   {
      SCIPdebugMessage("Rule-out periph (via contraction) \n");

      return TRUE;
   }

   return FALSE;
}

/**@} */
