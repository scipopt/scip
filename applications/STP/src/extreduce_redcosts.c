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

/**@file   extreduce_redcosts.c
 * @brief  reduced cost routines for extended reduction techniques for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements methods for using reduced costs in extended reduction techniques.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

 //#define SCIP_DEBUG
 // #define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"

#ifndef NDEBUG
/** recomputes simple reduced costs for current tree */
static
SCIP_Real getTreeRedcosts_dbg(
   const GRAPH*          graph,              /**< graph data structure */
   int                   redcostlevel,       /**< the reduced costs level */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   const REDCOST* const redcostdata = extdata->redcostdata;
   const PATH* const nodeTo3TermsPaths = redcosts_getNodeToTermsPaths(redcostdata, redcostlevel);
   const SCIP_Real* const rootToNodeDist = redcosts_getRootToNodeDist(redcostdata, redcostlevel);
   const SCIP_Real* const redcost = redcosts_getEdgeCosts(redcostdata, redcostlevel);
   const int* const tree_edges = extdata->tree_edges;
   const int* const tree_leaves = extdata->tree_leaves;
   const int tree_nedges = extdata->tree_nedges;
   const int nleaves = extdata->tree_nleaves;
   const int tree_root = extdata->tree_root;
   SCIP_Real tree_redcost;

   tree_redcost = rootToNodeDist[root];

   for( int i = 0; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      if( leaf == root )
         continue;

      tree_redcost += nodeTo3TermsPaths[leaf].dist;
   }

   for( int i = 0; i < tree_nedges; i++ )
   {
      const int edge = tree_edges[i];
      assert(edge >= 0 && edge < graph->edges);

      tree_redcost += redcost[edge];
      assert(LT(tree_redcost, FARAWAY));
   }

   for( int node = root; node != tree_root; node = extdata->tree_parentNode[node] )
   {
      int e;
      const int parent = extdata->tree_parentNode[node];
      assert(graph_knot_isInRange(graph, parent));

      for( e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == parent )
            break;
      }

      assert(e != EAT_LAST);
      tree_redcost += redcost[e];
      tree_redcost -= redcost[flipedge(e)];
   }

   return tree_redcost;
}
#endif


/** insertion sort; keyArr has sentinel value at position -1
 * : do something special: maybe sort index array
 * */
static inline
void sortDescendingIntRealReal(
   int* RESTRICT         keyArr,             /**< key array of size 'nentries' */
   SCIP_Real* RESTRICT   dataArr1,           /**< array of size 'nentries' */
   SCIP_Real* RESTRICT   dataArr2,           /**< array of size 'nentries' */
   int                   nentries            /**< number of entries */
)
{
   assert(keyArr && dataArr1 && dataArr2);
   assert(nentries >= 1);

#if 1
   for( int i = 1; i < nentries; i++ )
   {
      int j;
      const int currKey = keyArr[i];
      const SCIP_Real currData1 = dataArr1[i];
      const SCIP_Real currData2 = dataArr2[i];

      for( j = i - 1; currKey > keyArr[j]; j-- )
      {
         assert(j >= 0);
         keyArr[j + 1] = keyArr[j];
         dataArr1[j + 1] = dataArr1[j];
         dataArr2[j + 1] = dataArr2[j];
      }

      keyArr[j + 1] = currKey;
      dataArr1[j + 1] = currData1;
      dataArr2[j + 1] = currData2;
   }
#else
   for( int i = 0; i < nentries - 1; i++ )
   {
      int max_pos = i;

      for( int j = i + 1; j < nentries; j++ )
      {
         if( keyArr[j] > keyArr[max_pos] )
            max_pos = j;
      }

      if( max_pos != i )
      {
         int temp_int;
         SCIP_Real tmp_real;

         temp_int = keyArr[max_pos];
         keyArr[max_pos] = keyArr[i];
         keyArr[i] = temp_int;

         tmp_real = dataArr1[max_pos];
         dataArr1[max_pos] = dataArr1[i];
         dataArr1[i] = tmp_real;

         tmp_real = dataArr2[max_pos];
         dataArr2[max_pos] = dataArr2[i];
         dataArr2[i] = tmp_real;
      }
   }
#endif

#ifndef NDEBUG
   for( int i = 1; i < nentries; i++ )
      assert(keyArr[i - 1] >= keyArr[i] );
#endif
}


/** helper for rooted tree reduced cost computation */
static inline
SCIP_Real getMinDistCombination(
   const SCIP_Real*      firstTermDist,      /**< array of size 'nentries' */
   const SCIP_Real*      secondTermDist,     /**< array of size 'nentries' */
   int                   nentries            /**< number of entries to check */
)
{
   SCIP_Real min;

   assert(firstTermDist && secondTermDist);
   assert(nentries >= 1);

   if( nentries == 1 )
   {
      min = firstTermDist[0];
   }
   else
   {
      int i;
      SCIP_Real secondSum = 0.0;
      min = FARAWAY;

      for( i = 0; i < nentries; i++ )
      {
         assert(LE(firstTermDist[i], secondTermDist[i]));
         assert(LE(secondTermDist[i], FARAWAY));

         if( EQ(secondTermDist[i], FARAWAY) )
            break;

         secondSum += secondTermDist[i];
      }

      assert(LT(secondSum, FARAWAY));

      /* is there an index i with secondTermDist[i] == FARAWAY? */
      if( i < nentries )
      {
         assert(EQ(secondTermDist[i], FARAWAY));

         min = firstTermDist[i];
         for( int j = 0; j < nentries; j++ )
         {
            if( j == i )
               continue;

            min += secondTermDist[j];
         }
      }
      else
      {
         for( i = 0; i < nentries; i++ )
         {
            const SCIP_Real distCombination = secondSum + firstTermDist[i] - secondTermDist[i];

            assert(LT(secondTermDist[i], FARAWAY));

            if( distCombination < min )
               min = distCombination;
         }

         assert(LE(min, secondSum));
      }
   }

   return min;
}


/** gets reduced cost of current tree rooted at leave 'root' */
static inline
SCIP_Real extTreeGetDirectedRedcostProper(
   const GRAPH*          graph,              /**< graph data structure */
   int                   redcostlevel,       /**< the reduced costs level */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   int nearestTerms_x[STP_EXTTREE_MAXNLEAVES_GUARD + 1];
   int* nearestTerms;
   SCIP_Real firstTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   SCIP_Real secondTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const REDDATA* reddata = extdata->reddata;
   const REDCOST* const redcostdata = extdata->redcostdata;
   const PATH* const nodeTo3TermsPaths = redcosts_getNodeToTermsPaths(redcostdata, redcostlevel);
   const SCIP_Real* const rootToNodeDist = redcosts_getRootToNodeDist(redcostdata, redcostlevel);
   const int* const next3Terms = redcosts_getNodeToTermsBases(redcostdata, redcostlevel);
   const SCIP_Bool* const isterm = extdata->node_isterm;
   const int* tree_deg = extdata->tree_deg;
   const int nnodes = graph->knots;
   const SCIP_Real tree_redcost = reddata->redcost_treecosts[redcostlevel];
   const SCIP_Real swapcost = reddata->redcost_treenodeswaps[redcostlevel * nnodes + root];
   SCIP_Real redcost_directed = tree_redcost + rootToNodeDist[root] + swapcost;
   int leavescount = 0;
#ifndef NDEBUG
   SCIP_Real redcost_debug = redcost_directed;
#endif

   nearestTerms_x[0] = INT_MAX;
   nearestTerms = &(nearestTerms_x[1]);

#ifndef NDEBUG
   for( int i = 0; i < STP_EXTTREE_MAXNLEAVES_GUARD; i++ )
   {
      nearestTerms[i] = -1;
      firstTermDist[i] = -1.0;
      secondTermDist[i] = -1.0;
   }

   assert(LT(redcost_directed, FARAWAY));
   assert(GE(tree_redcost, 0.0));
   assert(graph_knot_isInRange(graph, root));
#endif



#ifdef STP_DEBUG_EXT
   SCIPdebugMessage("...reduced costs without leaves: %f (treecost=%f + rootdist=%f + swap=%f) \n",
         redcost_directed, tree_redcost, rootToNodeDist[root], swapcost);
#endif

   for( int j = 0; j < nleaves; j++ )
   {
      int i;
      int term;
      const int leaf = tree_leaves[j];

#ifdef STP_DEBUG_EXT
      SCIPdebugMessage("...checking leaf: %d \n", leaf);
#endif

      if( leaf == root || isterm[leaf] )
      {
         assert(leaf == root || EQ(0.0, nodeTo3TermsPaths[leaf].dist));
         continue;
      }

      /* find closest valid terminal for extension */
      for( i = 0; i < 3; i++ )
      {
         term = next3Terms[leaf + i * nnodes];

         if( term == UNKNOWN )
            break;

         assert(graph_pc_isPcMw(graph) || Is_term(graph->term[term]));
         assert(term >= 0 && term < graph->knots);
         assert(term != leaf);

         /* terminal not in current tree?*/
         if( tree_deg[term] == 0 )
            break;

         if( tree_deg[term] < 0 )
         {
            assert(graph_pc_isPcMw(graph) && tree_deg[term] == -1);
            assert(Is_pseudoTerm(graph->term[term]));
            break;
         }
      }

      /* no terminal reachable? */
      if( term == UNKNOWN )
      {
#ifdef STP_DEBUG_EXT
         SCIPdebugMessage("...no terminal reachable from leaf %d \n", leaf);
#endif
         assert(i < 3 && GE(nodeTo3TermsPaths[leaf + i * nnodes].dist, FARAWAY));
         return FARAWAY;
      }

      /* all terminals in current tree? */
      if( i == 3 )
         i = 2;

      nearestTerms[leavescount] = term;
      firstTermDist[leavescount] = nodeTo3TermsPaths[leaf + i * nnodes].dist;
      secondTermDist[leavescount] = (i < 2) ? nodeTo3TermsPaths[leaf + (i + 1) * nnodes].dist : firstTermDist[leavescount];

      assert(LE(nodeTo3TermsPaths[leaf].dist, firstTermDist[leavescount]));
      assert(LE(firstTermDist[leavescount], secondTermDist[leavescount]));

#ifdef STP_DEBUG_EXT
      SCIPdebugMessage("...closeterm=%d, dist1=%f dist2=%f distdef=%f \n", term,
            firstTermDist[leavescount], secondTermDist[leavescount], nodeTo3TermsPaths[leaf].dist);
#endif

      leavescount++;

#ifndef NDEBUG
      redcost_debug += nodeTo3TermsPaths[leaf].dist;
      assert(leavescount <= STP_EXTTREE_MAXNLEAVES_GUARD);
#endif
   }

   if( leavescount > 0 )
   {
      int first = 0;

      sortDescendingIntRealReal(nearestTerms, firstTermDist, secondTermDist, leavescount);

      for( int i = 1; i < leavescount; i++ )
      {
         assert(nearestTerms[i] >= 0 && firstTermDist[i] >= 0.0 && secondTermDist[i] >= 0.0);
         if( nearestTerms[i] != nearestTerms[i - 1] )
         {
            const int n = i - first;
            redcost_directed += getMinDistCombination(firstTermDist + first, secondTermDist + first, n);
            first = i;
         }
      }

      redcost_directed += getMinDistCombination(firstTermDist + first, secondTermDist + first, leavescount - first);
   }


#ifdef STP_DEBUG_EXT
   SCIPdebugMessage("...reduced costs with leaves: %f \n", redcost_directed);
#endif

 //  printf("redcost_directed=%f redcost_debug=%f  \n", redcost_directed, redcost_debug );

   assert(GE(redcost_directed, redcost_debug));

   return redcost_directed;
}


/** gets reduced cost of current tree rooted at leave 'root' */
static inline
SCIP_Real extTreeGetDirectedRedcost(
   const GRAPH*          graph,              /**< graph data structure */
   int                   redcostlevel,       /**< the reduced costs level */
   const EXTDATA*        extdata,            /**< extension data */
   const REDDATA*        reddata,            /**< reduction data */
   int                   root                /**< the root for the orientation */
)
{
   const SCIP_Real* const redcost_treenodeswaps = reddata->redcost_treenodeswaps;

   assert(extdata->tree_nleaves > 1 && extdata->tree_nleaves < STP_EXTTREE_MAXNLEAVES_GUARD);
   assert(extdata->tree_leaves[0] == extdata->tree_root);
   assert(graph_knot_isInRange(graph, root));

#ifdef STP_DEBUG_EXT
   SCIPdebugMessage("Check directed tree rooted at %d \n", root);
#endif

   /* are there any deleted arcs in the directed tree? */
   if( extdata->tree_nDelUpArcs > 0 && root == extdata->tree_root )
   {
      return FARAWAY;
   }

   /* is the rooting possible? */
   if( LT(redcost_treenodeswaps[redcostlevel * graph->knots + root], FARAWAY) )
   {
      return extTreeGetDirectedRedcostProper(graph, redcostlevel, extdata, root);
   }

#ifdef STP_DEBUG_EXT
   SCIPdebugMessage("Directed tree directly ruled-out \n");
#endif

   return FARAWAY;
}


/** checks for reduced-cost cutoff of current tree */
static inline
SCIP_Bool extTreeRedcostCutoff(
   const GRAPH*          graph,              /**< graph data structure */
   int                   redcostlevel,       /**< level */
   const EXTDATA*        extdata,            /**< extension data */
   const REDDATA*        reddata             /**< reduction data */
)
{
   const int* const tree_leaves = extdata->tree_leaves;
   const SCIP_Real cutoff = redcosts_getCutoff(extdata->redcostdata, redcostlevel);
   const int nleaves = extdata->tree_nleaves;
   const SCIP_Bool allowEquality = reddata->redcost_allowEquality;

#ifdef SCIP_DEBUG
   SCIP_Real tree_redcost = FARAWAY;
#endif

#ifdef STP_DEBUG_EXT
   SCIPdebugMessage("---Redcosts checking LEVEL %d--- \n", redcostlevel);
#endif

   /* take each leaf as root of the tree and check whether cutoff condition is violated */
   for( int i = 0; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      const SCIP_Real tree_redcost_new = extTreeGetDirectedRedcost(graph, redcostlevel, extdata, reddata, leaf);

     // printf("%f  >= %f \n", tree_redcost_new, getTreeRedcosts_dbg(graph, extdata, leaf));
      assert(GE_FEAS_EPS(tree_redcost_new, getTreeRedcosts_dbg(graph, redcostlevel, extdata, leaf), EPS_ZERO_HARD));

      if( allowEquality ? LT(tree_redcost_new, cutoff) : LE(tree_redcost_new, cutoff) )
      {
         SCIPdebugMessage("NO rule-out periph (red.cost<=%f from root=%d and cutoff=%f) \n", tree_redcost_new, leaf, cutoff);
         return FALSE;
      }

#ifdef SCIP_DEBUG
      tree_redcost = MIN(tree_redcost, tree_redcost_new);
#endif
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("Rule-out periph (with red.cost=%f and cutoff=%f, level=%d) \n", tree_redcost, cutoff, redcostlevel);
#endif

   return TRUE;
}


/** recomputes reduced cost tree information */
void extreduce_redcostTreeRecompute(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{

   REDDATA* const reddata = extdata->reddata;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   SCIP_Real* lvl_redcostsbuffer;
   const SCIP_Real* const lvl_redcosts = redcosts_getEdgeCosts(extdata->redcostdata, 0);
   const int* const tree_edges = extdata->tree_edges;
   SCIP_Real* const redcost_treecosts = extdata->reddata->redcost_treecosts;
   const int tree_nedges = extdata->tree_nedges;
   const int nlevels = reddata->redcost_nlevels;
   const int nedges = graph->edges;

   assert(!extreduce_treeIsFlawed(scip, graph, extdata));
   assert(nedges == redcosts_getNedges(extdata->redcostdata));

   SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &lvl_redcostsbuffer, nlevels) );

   for( int i = 0; i < nlevels; i++ )
      lvl_redcostsbuffer[i] = 0.0;

   for( int i = 0; i < tree_nedges; i++ )
   {
      const int edge = tree_edges[i];
      const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);
      assert(graph_edge_isInRange(graph, edge));

      if( !edgeIsDeleted )
      {
         /* NOTE: not clean, but better for better performance */
         for( int j = 0; j < nlevels; j++ )
         {
            lvl_redcostsbuffer[j] += lvl_redcosts[j * nedges + edge];
            assert(LT(lvl_redcostsbuffer[j], FARAWAY));
         }
      }
   }

   for( int i = 0; i < nlevels; i++ )
   {
      assert(SCIPisEQ(scip, redcost_treecosts[i], lvl_redcostsbuffer[i]));
      redcost_treecosts[i] = lvl_redcostsbuffer[i];
   }

   SCIPfreeBlockMemoryArray(scip, &lvl_redcostsbuffer, nlevels);
}


/** NOTE: call before adding edges for new expansion! */
void extreduce_redcostInitExpansion(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   const REDCOST* const redcostdata = extdata->redcostdata;
   SCIP_Bool* const noReverses = reddata->redcost_noReversedTree;
   const SCIP_Real* const treenodeswaps = reddata->redcost_treenodeswaps;
   const int nlevels = reddata->redcost_nlevels;
   const int comproot = extStackGetTopRoot(graph, extdata);
   const int tree_root = extdata->tree_root;
   const int nnodes = graph->knots;

   for( int i = 0; i < nlevels; i++ )
   {
      noReverses[i] = (redcosts_getRoot(redcostdata, i) == tree_root || GE(treenodeswaps[nnodes * i + comproot], FARAWAY));
   }
}


/** Updates reduced cost for added edge.
 *  NOTE: call extreduce_redcostInitExpansion before each new expansion! */
void extreduce_redcostAddEdge(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to be added */
   REDDATA*              reddata,            /**< reduction data */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Real* const redcosts = redcosts_getEdgeCosts(extdata->redcostdata, 0);
   SCIP_Real* const redcost_treenodeswaps = reddata->redcost_treenodeswaps;
   SCIP_Real* const redcost_treecosts = reddata->redcost_treecosts;
   SCIP_Bool* const redcost_noReverses = reddata->redcost_noReversedTree;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];
   const int nlevels = reddata->redcost_nlevels;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;

   assert(!edgedeleted && "deprecated!");
   assert(nedges == redcosts_getNedges(extdata->redcostdata));
   assert(nnodes == redcosts_getNnodes(extdata->redcostdata));

   for( int i = 0; i < nlevels; i++ )
   {
      const int offset_node = i * nnodes;
      const int offset_edge = i * nedges;

      if( redcost_noReverses[i] || (edgedeleted && edgedeleted[flipedge(edge)]) )
      {
         redcost_treenodeswaps[offset_node + head] = FARAWAY;

         if( head == extdata->tree_starcenter )
         {
            assert(extIsAtInitialStar(extdata) || extIsAtInitialGenStar(extdata));
            redcost_noReverses[i] = TRUE;
         }
      }
      else
      {
         redcost_treenodeswaps[offset_node + head] =
            redcost_treenodeswaps[offset_node + tail] + redcosts[offset_edge + flipedge(edge)];

         if( !edgeIsDeleted )
            redcost_treenodeswaps[offset_node + head] -= redcosts[offset_edge + edge];

         assert(LT(redcost_treenodeswaps[offset_node + head], FARAWAY));
      }

      if( !edgeIsDeleted )
      {
         redcost_treecosts[i] += redcosts[offset_edge + edge];
         assert(LT(redcost_treecosts[i], FARAWAY));
      }
   }

   if( edgeIsDeleted )
   {
      extdata->tree_nDelUpArcs++;
   }
}


/** update reduced cost for removed edge */
void extreduce_redcostRemoveEdge(
   int                   edge,               /**< edge to be added */
   const REDDATA*        reddata,            /**< reduction data */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Bool edgeIsDeleted = (reddata->edgedeleted && reddata->edgedeleted[edge]);

   if( !edgeIsDeleted )
   {
      const SCIP_Real* const redcosts = redcosts_getEdgeCosts(extdata->redcostdata, 0);
      const int nlevels = reddata->redcost_nlevels;
      const int nedges = redcosts_getNedges(extdata->redcostdata);
      SCIP_Real* const redcost_treecosts = reddata->redcost_treecosts;

      for( int i = 0; i < nlevels; i++ )
      {
         redcost_treecosts[i] -= redcosts[i * nedges + edge];
         assert(LT(redcost_treecosts[i], FARAWAY));
      }
   }
   else
   {
      extdata->tree_nDelUpArcs--;
      assert(extdata->tree_nDelUpArcs >= 0);
   }
}


/** Can current tree be peripherally ruled out by using reduced costs arguments? */
SCIP_Bool extreduce_redcostRuleOutPeriph(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const REDDATA* const reddata = extdata->reddata;
   const int nlevels = reddata->redcost_nlevels;

   assert(graph);
   assert(extdata->tree_nleaves > 1 && extdata->tree_leaves[0] == extdata->tree_root);
   assert(nlevels >= 1);

   for( int i = 0; i < nlevels; i++ )
   {
      if( extTreeRedcostCutoff(graph, i, extdata, reddata) )
      {
         return TRUE;
      }
   }

   return FALSE;
}
