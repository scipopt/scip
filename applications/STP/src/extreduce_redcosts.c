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

// #define SCIP_DEBUG
// #define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"



/** insertion sort; todo
 * : could be speed-up by use of sentinel value at position 0
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

   for( int i = 1; i < nentries; i++ )
   {
      int j;
      const int currKey = keyArr[i];
      const SCIP_Real currData1 = dataArr1[i];
      const SCIP_Real currData2 = dataArr2[i];

      for( j = i - 1; j >= 0 && currKey > keyArr[j]; j-- )
      {
         keyArr[j + 1] = keyArr[j];
         dataArr1[j + 1] = dataArr1[j];
         dataArr2[j + 1] = dataArr2[j];
      }

      keyArr[j + 1] = currKey;
      dataArr1[j + 1] = currData1;
      dataArr2[j + 1] = currData2;
   }

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


/** gets reduced cost of current tree rooted at leave 'root', called direct if tree cannot */
static inline
SCIP_Real extTreeGetDirectedRedcostProper(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   int nearestTerms[STP_EXTTREE_MAXNLEAVES_GUARD];
   SCIP_Real firstTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   SCIP_Real secondTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const SCIP_Real swapcost = extdata->tree_redcostSwap[root];
   const REDDATA* const reddata = extdata->reddata;
   const PATH* const nodeTo3TermsPaths = reddata->nodeTo3TermsPaths;
   const int* const next3Terms = reddata->nodeTo3TermsBases;
   const SCIP_Bool* const isterm = extdata->node_isterm;
   const int* tree_deg = extdata->tree_deg;

   SCIP_Real redcost_directed = extdata->tree_redcost + reddata->rootToNodeDist[root] + swapcost;
   const int nnodes = graph->knots;
   int leavescount = 0;

#ifndef NDEBUG
   SCIP_Real redcost_debug = redcost_directed;
   for( int i = 0; i < STP_EXTTREE_MAXNLEAVES_GUARD; i++ )
   {
      nearestTerms[i] = -1;
      firstTermDist[i] = -1.0;
      secondTermDist[i] = -1.0;
   }

   assert(LT(redcost_directed, FARAWAY));
#endif

   for( int j = 0; j < nleaves; j++ )
   {
      int i;
      int term;
      const int leaf = tree_leaves[j];

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

      //   printf("i %d \n", i);
 // printf("term=%d, first=%f second=%f def=%f \n", term, firstTermDist[leavescount], secondTermDist[leavescount], nodeTo3TermsPaths[leaf].dist);
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

 //  printf("redcost_directed=%f redcost_debug=%f  \n", redcost_directed, redcost_debug );

   assert(GE(redcost_directed, redcost_debug));

   return redcost_directed;
}


/** gets reduced cost of current tree rooted at leave 'root' */
static inline
SCIP_Real extTreeGetDirectedRedcost(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   const SCIP_Real* const tree_redcostSwap = extdata->tree_redcostSwap;

   assert(extdata->tree_nleaves > 1 && extdata->tree_nleaves < STP_EXTTREE_MAXNLEAVES_GUARD);
   assert(extdata->tree_leaves[0] == extdata->tree_root);
   assert(root >= 0 && root < graph->knots);

   /* are there any deleted arcs in the directed tree? */
   if( root == extdata->tree_root && extdata->tree_nDelUpArcs > 0 )
   {
      return FARAWAY;
   }

   /* is the rooting possible? */
   if( LT(tree_redcostSwap[root], FARAWAY) )
   {
      return extTreeGetDirectedRedcostProper(graph, extdata, root);
   }

   return FARAWAY;
}


/** gets reduced cost bound of current tree */
static inline
SCIP_Real extTreeGetRedcostBound(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Real tree_redcost;

   assert(graph && extdata);
   assert(nleaves > 1 && tree_leaves[0] == extdata->tree_root);

   tree_redcost = FARAWAY;

   // todo perhaps needs to be adapted for pseudo elimination tests...

   /* take each leaf as root of the tree */
   for( int i = 0; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      const SCIP_Real tree_redcost_new = extTreeGetDirectedRedcost(graph, extdata, leaf);
      int todo;
      // break early here!
      // move the entire redcost stuff to an extra file!

      tree_redcost = MIN(tree_redcost, tree_redcost_new);
   }

   return tree_redcost;
}


/** no reversed tree possible? */
SCIP_Bool extreduce_redcostReverseTreeRuledOut(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const REDDATA* const reddata = extdata->reddata;
   const int comproot = extStackGetTopRoot(graph, extdata);

   return (reddata->redCostRoot == extdata->tree_root || GE(extdata->tree_redcostSwap[comproot], FARAWAY));
}


/** update reduced cost for added edge */
void extreduce_redcostAddEdge(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to be added */
   SCIP_Bool             noReversedTree,     /**< don't consider reversed tree? */
   const REDDATA*        reddata,            /**< reduction data */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Real* const redcost = reddata->redCosts;
   SCIP_Real* const tree_redcostSwap = extdata->tree_redcostSwap;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);
   const int head = graph->head[edge];

   if( noReversedTree || (edgedeleted && edgedeleted[flipedge(edge)]) )
   {
      tree_redcostSwap[head] = FARAWAY;
   }
   else
   {
      const int tail = graph->tail[edge];

      tree_redcostSwap[head] = tree_redcostSwap[tail] + redcost[flipedge(edge)];

      if( !edgeIsDeleted )
         tree_redcostSwap[head] -= redcost[edge];

      assert(LT(tree_redcostSwap[head], FARAWAY));
   }

   if( !edgeIsDeleted )
   {
      extdata->tree_redcost += redcost[edge];
      assert(LT(extdata->tree_redcost, FARAWAY));
   }
   else
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
   const SCIP_Real* const redcost = reddata->redCosts;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);

   if( !edgeIsDeleted )
   {
      extdata->tree_redcost -= redcost[edge];
      assert(LT(extdata->tree_redcost, FARAWAY));
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
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcostBound(graph, extdata);
   const SCIP_Real cutoff = reddata->cutoff;

   if( reddata->equality ? GE(tree_redcost, cutoff) : GT(tree_redcost, cutoff) )
   {
      SCIPdebugMessage("Rule-out periph (with red.cost=%f) \n", tree_redcost);
      return TRUE;
   }

   return FALSE;
}
