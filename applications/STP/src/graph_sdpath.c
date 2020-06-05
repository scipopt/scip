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

/**@file   graph_sdpath.c
 * @brief  Special distance path and walk algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various special distance path (and walk) algorithms for
 * Steiner tree and related problems.
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"
#include "graphheaps.h"
#include "reduce.h"


/** internal data for clique-paths computations */
typedef struct special_distance_clique_paths
{
   SCIP_Real*            nodes_distBaseToMax;/**< distance from base to maximum profit */
   SCIP_Real*            nodes_distFromMax;  /**< distance from maximum profit to current node */
   SCIP_Real*            nodes_maxprofit;    /**< maximum profit on path */
   int*                  nodes_pred;         /**< per node: predecessor */
   int*                  nodes_base;         /**< per node: base */
   int*                  nodes_baseToClique; /**< per base: clique id
                                                  NOTE: undefined for non-base nodes! */
} CLIQUEPATHS;


/** initializes clique paths */
static
SCIP_RETCODE cliquePathsInitData(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   CLIQUEPATHS*          cliquepaths         /**< data */
)
{
   const int nnodes = graph_get_nNodes(g);

   assert(scip && cliquepaths);

   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_distBaseToMax), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_maxprofit), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_distFromMax), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_pred), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_base), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(cliquepaths->nodes_baseToClique), nnodes) );

   return SCIP_OKAY;
}


/** frees clique paths */
static
void cliquePathsFreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   CLIQUEPATHS*          cliquepaths         /**< data */
)
{
   assert(scip && cliquepaths);

   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_baseToClique));
   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_base));
   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_pred));
   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_distFromMax));
   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_maxprofit));
   SCIPfreeBufferArray(scip, &(cliquepaths->nodes_distBaseToMax));
}


/** gets distance for extension along edge e */
inline static
SCIP_Real sdwalkGetdistnewEdge(
   const int*            prevedges,          /**< previous edges per node */
   const int*            nprevedges,         /**< number of previous edges per node */
   const SCIP_Real*      cost,               /**< cost */
   const SCIP_Real*      dist,               /**< distance */
   int                   k,                  /**< previous node  */
   int                   e,                  /**< current outgoing edge  */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   const int nprevs = nprevedges[k];
   SCIP_Real dist_e;

   /* ancestor list not full? */
   if( nprevs != maxnprevs + 1 )
   {
      int i;
      const int e2 = e / 2;
      assert(nprevs <= maxnprevs);

      /* check whether m is contained in ancestor list */
      for( i = 0; i < nprevs; i++ )
      {
         const int prevedge = prevedges[maxnprevs * k + i];

         if( e2 == prevedge )
            break;
      }

      /* e2 in list? */
      if( i != nprevs )
      {
         assert(e2 == prevedges[maxnprevs * k + i]);
         dist_e = dist[k];
      }
      else
         dist_e = dist[k] + cost[e];
   }
   else
   {
      dist_e = dist[k] + cost[e];
   }

   return dist_e;
}


inline static
SCIP_Real sdwalkGetdistnewPrize(
   const int*            prevNPterms,        /**< previous np terminals per node */
   const int*            nprevNPterms,       /**< number of previous np terminals per node */
   const int*            termmark,           /**< terminal mark */
   const STP_Bool*       visited,            /**< visited */
   const SCIP_Real*      prize,              /**< prize */
   int                   k,                  /**< current node  */
   int                   m,                  /**< next node  */
   SCIP_Real             distnew,            /**< distance of m */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   SCIP_Real distnewP = distnew;

   assert(termmark[m] == 1 || termmark[m] == 2 );

   if( termmark[m] == 2 || !visited[m] )
   {
      distnewP = MAX(0.0, distnewP - prize[m]);
   }
   else
   {
      const int nprevs = nprevNPterms[k];

      /* ancestor list not full? */
      if( nprevs != maxnprevs + 1 )
      {
         int i;
         assert(nprevs <= maxnprevs);

         /* check whether m is contained in ancestor list */
         for( i = 0; i < nprevs; i++ )
         {
            const int prevterm = prevNPterms[maxnprevs * k + i];

            if( m == prevterm )
               break;
         }

         /* m not in list? */
         if( i == nprevs )
            distnewP = MAX(0.0, distnewP - prize[m]);
      }
   }

   return distnewP;
}



inline static
SCIP_Bool sdwalkHasConflict(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   const int*            prevterms,          /**< previous terminals */
   const int*            nprevterms,         /**< number of previous terminals */
   const SCIP_Bool       nodevisited         /**< visited node already? */
   )
{
   const int nprevs = nprevterms[prednode];
   SCIP_Bool conflict = FALSE;

   assert(Is_term(g->term[node]));

   if( !nodevisited )
      return FALSE;

   if( nprevs > maxnprevs )
   {
      assert(nprevs == maxnprevs + 1);
      return TRUE;
   }

   for( int i = 0; i < nprevs; i++ )
   {
      const int prevterm = prevterms[maxnprevs * prednode + i];
      assert(prevterm >= 0);

      if( prevterm == node )
      {
         conflict = TRUE;
         break;
      }
   }

   return conflict;
}

/** updates */
inline static
void sdwalkUpdate(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms          /**< number of previous terminals */
   )
{
   const int predsize = nprevterms[prednode];
   const SCIP_Bool isterm = Is_term(g->term[node]);

   assert(predsize <= maxnprevs + 1);

   if( predsize == maxnprevs + 1 || (isterm && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      for( int i = 0; i < predsize; i++ )
         prevterms[maxnprevs * node + i] = prevterms[maxnprevs * prednode + i];

      nprevterms[node] = predsize;

      if( isterm )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }
}

/** copies from of predecessor */
inline static
void sdwalkUpdateCopy(
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prev,               /**< previous data elements */
   int*                  nprev               /**< number of previous data elements */
   )
{
   const int predsize = nprev[prednode];

   assert(predsize <= maxnprevs);

   /* copy data from predecesseor */
   for( int i = 0; i < predsize; i++ )
      prev[maxnprevs * node + i] = prev[maxnprevs * prednode + i];

   nprev[node] = predsize;
}

/** update method for second version of SD walks */
static
void sdwalkUpdate2(
   const int*            termmark,           /**< terminal mark */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   edge,               /**< the edge to be updated */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Bool             clear,              /**< clear arrays */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges          /**< number of previous edges */
   )
{
   int predsize = nprevterms[prednode];

   /*** 1. proper terminals ***/

   /* not enough space? */
   if( predsize == maxnprevs + 1 || (termmark[node] == 2 && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      sdwalkUpdateCopy(node, prednode, maxnprevs, prevterms, nprevterms);

      if( termmark[node] == 2 )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }


   /*** 2. edges ***/

   if( clear )
   {
      nprevNPterms[node] = 0;
      nprevedges[node] = 0;
      return;
   }

   predsize = nprevedges[prednode];

   if( predsize >= maxnprevs )
   {
      assert(predsize == maxnprevs || predsize == maxnprevs + 1);

      nprevedges[node] = maxnprevs + 1;
      nprevNPterms[node] = maxnprevs + 1;
      return;
   }
   assert(predsize < maxnprevs);

   sdwalkUpdateCopy(node, prednode, maxnprevs, prevedges, nprevedges);

   prevedges[maxnprevs * node + predsize] = edge / 2;
   nprevedges[node]++;

   assert(nprevedges[node] <= maxnprevs);


   /*** 3. non-proper terminals ***/

   predsize = nprevNPterms[prednode];

   if( predsize == maxnprevs + 1 || (termmark[node] == 1 && predsize == maxnprevs) )
   {
      nprevNPterms[node] = maxnprevs + 1;
   }
   else
   {
      sdwalkUpdateCopy(node, prednode, maxnprevs, prevNPterms, nprevNPterms);

      if( termmark[node] == 1 )
      {
         assert(predsize < maxnprevs);
         prevNPterms[maxnprevs * node + predsize] = node;
         nprevNPterms[node]++;
      }

      assert(nprevNPterms[node] <= maxnprevs);
   }
}

/** resets temporary data */
inline static
void sdwalkReset(
   int                   nvisits,            /**< number of visited nodes */
   const int*            visitlist,          /**< stores all visited nodes */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   STP_Bool*             visited             /**< stores whether a node has been visited */
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      assert(node >= 0);

      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
   }
}


/** corrects heap entry */
inline static
void sdwalkCorrectHeap(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   SCIP_Real* RESTRICT pathdist,
   int node,
   SCIP_Real newcost
   )
{
   int    c;
   int    j;

   pathdist[node] = newcost;

   if (state[node] == UNKNOWN)
   {
      heap[++(*count)] = node;
      state[node]      = (*count);
   }

   /* Heap shift up
    */
   j = state[node];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      const int t    = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


/** gets position in Sd array
 *  todo: bad design, should be somewhere central...probably not too time expensive */
inline static
int sdCliqueStarGetSdPosition(
   int                   ncliquenodes,       /**< number of clique nodes */
   int                   newnode,            /**< new node */
   int                   prednode,           /**< predecessor node */
   const SCIP_Real*      nodes_dist,         /**< distance per node */
   const int*            nodes_base,         /**< base per node */
   const int*            baseToClique        /**< mapping to clique ID */
)
{
   const int base_new = nodes_base[newnode];
   const int base_pred = nodes_base[prednode];
   const int id_new = baseToClique[base_new];
   const int id_pred = baseToClique[base_pred];
   const int id_min = MIN(id_new, id_pred);
   const int id_max = MAX(id_new, id_pred);
   int pos = 0;

   assert(id_min >= 0);
   assert(id_max <= ncliquenodes);
   assert(id_min < id_max);

   for( int k = 1; k <= id_min; k++ )
   {
      pos += (ncliquenodes - k);
   }

   pos += (id_max - id_min - 1);
   assert(pos < ((ncliquenodes) * (ncliquenodes - 1)) / 2);

   return pos;
}


/** updates node data */
inline static
void sdCliqueStarUpdateNodeMaxDist(
   int                   newnode,            /**< new node */
   int                   prednode,           /**< predecessor node */
   const SCIP_Real*      nodes_dist,         /**< distance per node */
   SCIP_Real             newprofit,          /**< new profit */
   SCIP_Real             newedgecost,        /**< new edge cost */
   SCIP_Real* RESTRICT   nodes_distBaseToMax,/**< nodes to max */
   SCIP_Real* RESTRICT   nodes_distFromMax,  /**< maximum bias */
   SCIP_Real* RESTRICT   nodes_maxpathprofit /**< maximum profit */
   )
{
   assert(newnode >= 0 && prednode >= 0);
   assert(newnode != prednode);
   assert(GE(nodes_distBaseToMax[prednode], 0.0));
   assert(GE(nodes_maxpathprofit[prednode], 0.0));
   assert(GE(newprofit, 0.0));
   assert(GE(newedgecost, 0.0));

   if( newprofit > nodes_maxpathprofit[prednode] )
   {
      nodes_distBaseToMax[newnode] = nodes_dist[prednode];
      nodes_distFromMax[newnode] = newedgecost;
      nodes_maxpathprofit[newnode] = newprofit;
   }
   else
   {
      SCIP_Real bias = MIN(newedgecost, newprofit);
      bias = MIN(nodes_distFromMax[prednode], bias);

      nodes_distBaseToMax[newnode] = nodes_distBaseToMax[prednode];
      nodes_distFromMax[newnode] = nodes_distFromMax[prednode] + newedgecost - bias;
      nodes_maxpathprofit[newnode] = nodes_maxpathprofit[prednode];
   }
}


/** updates node data */
inline static
void sdCliqueStarUpdateNodeDist(
   int                   newnode,            /**< new node */
   int                   prednode,           /**< predecessor node */
   SCIP_Real             newdist,            /**< new distance */
   DHEAP* RESTRICT       dheap,              /**< heap */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance per node */
   int* RESTRICT         nodes_base,         /**< base per node */
   int* RESTRICT         nodes_pred          /**< predecessor per node */
   )
{
   assert(newnode >= 0 && prednode >= 0);
   assert(newnode != prednode);
   assert(GE(newdist, 0.0) && LT(newdist, FARAWAY));

   graph_heap_correct(newnode, newdist, dheap);

   nodes_pred[newnode] = prednode;
   nodes_dist[newnode] = newdist;
   nodes_base[newnode] = nodes_base[prednode];
}


/** gets biased special distance for two profits  */
static inline
SCIP_Real sdGet2ProfitsDist(
   SCIP_Real             distBaseToBase,
   SCIP_Real             distBaseToProfit1,
   SCIP_Real             distBaseToProfit2,
   SCIP_Real             profit1,
   SCIP_Real             profit2
   )
{
   const SCIP_Real distProfitToProfit = distBaseToBase - distBaseToProfit1 - distBaseToProfit2;
   const SCIP_Real sdAll = distBaseToBase - profit1 - profit2;
   const SCIP_Real sdProfit1 = distBaseToBase - distBaseToProfit2 - profit1;
   const SCIP_Real sdProfit2 = distBaseToBase - distBaseToProfit1 - profit2;
   const SCIP_Real sd_final = miscstp_maxReal((SCIP_Real[])
               { sdAll,
                 sdProfit1, sdProfit2,
                 distProfitToProfit,
                 distBaseToProfit1, distBaseToProfit2
               },
                 6);

   assert(GE(distBaseToProfit1, 0.0));
   assert(GE(distBaseToProfit2, 0.0));
   assert(GE(distProfitToProfit, 0.0));
   assert(LE(sd_final, distBaseToBase));

   return sd_final;
}


/** gets biased special distance for one profit  */
static inline
SCIP_Real sdGet1ProfitDist(
   SCIP_Real             distBaseToBase,
   SCIP_Real             distBaseToProfit,
   SCIP_Real             profit
   )
{
   const SCIP_Real sdAll = distBaseToBase - profit;
   const SCIP_Real distBaseToProfit2 = distBaseToBase - distBaseToProfit;
   SCIP_Real sd_final = MAX(distBaseToProfit, distBaseToProfit2);

   if( sdAll > sd_final )
      sd_final = sdAll;

   assert(GE(distBaseToProfit, 0.0));
   assert(GE(distBaseToProfit2, 0.0));
   assert(GE(sd_final, 0.0));
   assert(LE(sd_final, distBaseToBase));

   return sd_final;
}

/** gets node biased */
static inline
SCIP_Real sdCliqueStarGetNodeBias(
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   int                   node,
   int                   nextnode,
   int                   prevnode,
   SCIP_Real             edgecost,
   SCIP_Real             dist
   )
{
   SCIP_Real bias;
   if( node == prevnode )
   {
      bias = 0.0;
   }
   else
   {
      const SCIP_Real profit = reduce_sdprofitGetProfit(sdprofit, node, nextnode, prevnode);
      bias = MIN(edgecost, profit);

      if( dist < bias )
         bias = dist;
   }

   assert(GE(bias, 0.0));

   return bias;
}


#if 0

/** computes SD of path
 *  USE IN DEBUG MODE ONLY! */
static
void sdCliqueStarGetPathSd(
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const GRAPH*          g,                  /**< the graph */
   int                   startnode,          /**< start node */
   const int*            nodes_pred,         /**< predecessor per node */
   const int*            nodes_base,         /**< base per node */
   SCIP_Real*            endToMaxDist,
   SCIP_Real*            startToEndDist,
   SCIP_Real*            maxprofit
   )
{
   SCIP_Real dist = 0.0;
   const int endnode = nodes_base[startnode];
   int prevnode = startnode;
   *maxprofit = 0.0;
   *endToMaxDist = 0.0;

   printf("go from %d \n", startnode);

   for( int node = startnode; node != endnode; node = nodes_pred[node] )
   {
      int e;
      SCIP_Real profit;
      SCIP_Real offset;
      const int nextnode = nodes_pred[node];

      assert(nextnode != node);

      for( e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( g->head[e] == nextnode )
            break;
      }

      assert(e != EAT_LAST);

      profit = reduce_sdprofitGetProfit(sdprofit, node, prevnode, nextnode);

      printf("check node %d, profit=%f, edgecost=%f \n", node, profit, g->cost[e]);

      offset = MIN(profit, g->cost[e]);
      offset = MIN(offset, dist);

      dist += g->cost[e] - offset;

      if( node != startnode && profit - offset > *maxprofit )
      {
         *maxprofit = profit - offset;
         *endToMaxDist = dist;
      }

      prevnode = node;
   }

   assert(LE(*endToMaxDist, dist));
   *endToMaxDist = dist - *endToMaxDist;

   *startToEndDist = dist;
}


/** computes SD one combined path
  * USE IN DEBUG MODE ONLY! */
static
SCIP_Real sdCliqueStarGetRecomputedSd(
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const GRAPH*          g,                  /**< the graph */
   const int*            nodes_pred,         /**< predecessor per node */
   int                   edge,               /**< edge */
   const SCIP_Real*      nodes_dist,         /**< distance per node */
   const int*            nodes_base          /**< base per node */
   )
{
   SCIP_Real distBaseToBase;
   SCIP_Real distBaseToTail;
   SCIP_Real distBaseToHead;
   SCIP_Real maxprofit_tail;
   SCIP_Real maxprofit_head;
   SCIP_Real distBaseToMax_tail;
   SCIP_Real distBaseToMax_head;
   const int tail = g->tail[edge];
   const int head = g->head[edge];
   const SCIP_Real edgecost = g->cost[edge];
   SCIP_Real sd_final;

   printf("\n tail=%d, head=%d, edgecost=%f \n", tail, head, g->cost[edge]);

   sdCliqueStarGetPathSd(sdprofit, g, tail, nodes_pred, nodes_base, &distBaseToMax_tail, &distBaseToTail, &maxprofit_tail);
   sdCliqueStarGetPathSd(sdprofit, g, head, nodes_pred, nodes_base, &distBaseToMax_head, &distBaseToHead, &maxprofit_head);

   {
      SCIP_Real profit_tail = reduce_sdprofitGetProfit(sdprofit, tail, head, nodes_pred[tail]);
      SCIP_Real profit_head = reduce_sdprofitGetProfit(sdprofit, head, tail, nodes_pred[head]);
      const SCIP_Real bias_tail = sdCliqueStarGetNodeBias(sdprofit, tail, head, nodes_pred[tail], edgecost, distBaseToTail);
      const SCIP_Real bias_head = sdCliqueStarGetNodeBias(sdprofit, head, tail, nodes_pred[head], edgecost, distBaseToHead);

      if( bias_tail > bias_head )
      {
         profit_tail -= bias_tail;
         distBaseToBase = distBaseToTail + distBaseToHead + edgecost - bias_tail;
      }
      else
      {
         profit_head -= bias_head;
         distBaseToBase = distBaseToTail + distBaseToHead + edgecost - bias_head;
      }

      assert(GE(profit_tail, 0.0));
      assert(GE(profit_head, 0.0));

      if( GE(profit_head, maxprofit_head) || Is_term(g->term[head]) )
      {
         maxprofit_head = profit_head;
         distBaseToMax_head = distBaseToHead;
      }

      if( GE(profit_tail, maxprofit_tail) || Is_term(g->term[tail]) )
      {
         maxprofit_tail = profit_tail;
         distBaseToMax_tail = distBaseToTail;
      }
   }

   printf("maxprofit_tail=%f, maxprofit_head=%f \n", maxprofit_tail, maxprofit_head);

   if( GT(maxprofit_head, 0.0) && GT(maxprofit_tail, 0.0) )
   {
      sd_final = sdGet2ProfitsDist(distBaseToBase, distBaseToMax_tail, distBaseToMax_head, maxprofit_tail, maxprofit_head);
   }
   else if( GT(maxprofit_tail, 0.0) )
   {
      sd_final = sdGet1ProfitDist(distBaseToBase, distBaseToMax_tail, maxprofit_tail);
   }
   else if( GT(maxprofit_head, 0.0)  )
   {
      sd_final = sdGet1ProfitDist(distBaseToBase, distBaseToMax_head, maxprofit_head);
   }
   else
   {
      sd_final = distBaseToBase;
   }

   printf("sd_final=%f \n", sd_final);

   return sd_final;
}
#endif


/** helper */
inline static
void sdCliqueStarGetFinalProfitData(
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   int                   centernode,         /**< center node */
   int                   node,               /**< node */
   int                   neighbor,           /**< neighbor node */
   const SCIP_Real*      nodes_dist,         /**< distance */
   const int*            nodes_pred,         /**< predecessors */
   const SCIP_Real*      nodes_distBaseToMax,/**< distance */
   const SCIP_Real*      nodes_distFromMax,  /**< distance */
   const SCIP_Real*      nodes_maxpathprofit,/**< maximum profit */
   SCIP_Real*            distToMax,          /**< pointer */
   SCIP_Real*            distFromMax,        /**< pointer */
   SCIP_Real*            maxprofit_node,     /**< pointer */
   SCIP_Bool*            nodeHasMaxProfit    /**< pointer */
   )
{
   if( node == nodes_pred[node] || centernode == node )
   {
      *maxprofit_node = 0.0;
      assert(centernode == node || EQ(nodes_distBaseToMax[node], 0.0));
      assert(centernode == node || EQ(nodes_maxpathprofit[node], 0.0));
      assert(centernode == node || EQ(nodes_distFromMax[node], 0.0));
   }
   else
   {
      *maxprofit_node = reduce_sdprofitGetProfit(sdprofit, node, neighbor, nodes_pred[node]);
   }

   /* is the saved profit better? */
   if( *maxprofit_node < nodes_maxpathprofit[node] )
   {
      *distToMax = nodes_distBaseToMax[node];
      *maxprofit_node = nodes_maxpathprofit[node];
      *distFromMax = nodes_distFromMax[node];
      *nodeHasMaxProfit = FALSE;
   }
   else
   {
      *distToMax = nodes_dist[node];
      *distFromMax = 0.0;
      *nodeHasMaxProfit = TRUE;
   }
}


/** updates SD between nodes */
inline static
void sdCliqueStarUpdateSd(
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const int*            nodes_pred,         /**< predecessors */
   const SCIP_Real*      nodes_baseToMaxDist,/**< distance to max */
   const SCIP_Real*      nodes_distFromMax,  /**< distance from max */
   const SCIP_Real*      nodes_maxpathprofit,/**< maximum profit */
   int                   ncliquenodes,       /**< number of clique nodes */
   int                   centernode,         /**< center node */
   int                   newnode,            /**< new node */
   int                   prednode,           /**< predecessor node */
   SCIP_Real             newdist,            /**< the new distance to 'newnode' (along prednode) */
   SCIP_Real             edgecost,           /**< the edgecost from 'pred' to 'new' */
   const SCIP_Real*      nodes_dist,         /**< distance per node */
   const int*            nodes_base,         /**< base per node */
   const int*            baseToClique,       /**< mapping to clique ID */
   SCIP_Real* RESTRICT   sds                 /**< to be filled */
   )
{
   const int sdposition = sdCliqueStarGetSdPosition(ncliquenodes,
      newnode, prednode, nodes_dist, nodes_base, baseToClique);
   /* complete distance along new node: */
   SCIP_Real sdist = nodes_dist[newnode] + newdist;

   if( sdprofit != NULL )
   {
      SCIP_Bool nodeHasMaxProfit_pred;
      SCIP_Real distToMax_pred;
      SCIP_Real distFromMax_pred;
      SCIP_Real maxprofit_pred;
      SCIP_Bool nodeHasMaxProfit_new;
      SCIP_Real distToMax_new;
      SCIP_Real distFromMax_new;
      SCIP_Real maxprofit_new;
      SCIP_Real distBaseToBase;
      SCIP_Real sd_final;

      sdCliqueStarGetFinalProfitData(sdprofit, centernode, newnode, prednode, nodes_dist, nodes_pred, nodes_baseToMaxDist,
            nodes_distFromMax, nodes_maxpathprofit, &distToMax_new, &distFromMax_new, &maxprofit_new, &nodeHasMaxProfit_new);

      sdCliqueStarGetFinalProfitData(sdprofit, centernode, prednode, newnode, nodes_dist, nodes_pred, nodes_baseToMaxDist,
            nodes_distFromMax, nodes_maxpathprofit, &distToMax_pred, &distFromMax_pred, &maxprofit_pred, &nodeHasMaxProfit_pred);

      distBaseToBase = distToMax_pred + distToMax_new + distFromMax_pred + distFromMax_new + edgecost;

      if( !nodeHasMaxProfit_pred )
      {
         SCIP_Real bias = sdCliqueStarGetNodeBias(sdprofit, prednode, newnode, nodes_pred[prednode], edgecost, distFromMax_pred);
         distBaseToBase -= bias;
      }
      else if( !nodeHasMaxProfit_new )
      {
         SCIP_Real bias = sdCliqueStarGetNodeBias(sdprofit, newnode, prednode, nodes_pred[newnode], edgecost, distFromMax_new);
         distBaseToBase -= bias;
      }

      if( GT(maxprofit_new, 0.0) && GT(maxprofit_pred, 0.0) )
      {
         sd_final = sdGet2ProfitsDist(distBaseToBase, distToMax_pred, distToMax_new, maxprofit_pred, maxprofit_new);
      }
      else if( GT(maxprofit_pred, 0.0) )
      {
         sd_final = sdGet1ProfitDist(distBaseToBase, distToMax_pred, maxprofit_pred);
      }
      else if( GT(maxprofit_new, 0.0)  )
      {
         sd_final = sdGet1ProfitDist(distBaseToBase, distToMax_new, maxprofit_new);
      }
      else
      {
         sd_final = distBaseToBase;
      }

#if 0
         if( EQ(sd_final, 6.0) )
         {
            sdCliqueStarGetRecomputedSd(sdprofit, g, nodes_pred, edge, nodes_dist, nodes_base);
            printf("bases %d %d \n", nodes_base[prednode], nodes_base[newnode]);
            printf("my sd_final=%f \n", sd_final);
         }
#endif
      assert(LE(sd_final, sdist));
      sdist = sd_final;
   }

   assert(nodes_base[newnode] != nodes_base[prednode]);
   assert(GE(newdist, 0.0) && LT(newdist, FARAWAY));
   assert(GE(sdist, 0.0) && LT(sdist, FARAWAY));

   if( sdist < sds[sdposition] )
   {
      sds[sdposition] = sdist;
   }
}


/** initializes */
static
void sdCliqueStarInit(
   const GRAPH*          g,                  /**< graph data structure */
   SDCLIQUE*             cliquedata,         /**< data */
   CLIQUEPATHS*          cliquepaths         /**< paths data */
   )
{
   DIJK* RESTRICT dijkdata = cliquedata->dijkdata;
   SCIP_Real* RESTRICT nodes_dist = dijkdata->node_distance;
   int* RESTRICT visitlist = dijkdata->visitlist;
   DHEAP* RESTRICT dheap = dijkdata->dheap;
   const int* const cliquenodes = cliquedata->cliquenodes;
   const int ncliquenodes = cliquedata->ncliquenodes;
   SCIP_Real* RESTRICT nodes_distBaseToMax = cliquepaths->nodes_distBaseToMax;
   SCIP_Real* RESTRICT nodes_distFromMax = cliquepaths->nodes_distFromMax;
   SCIP_Real* RESTRICT nodes_maxprofit = cliquepaths->nodes_maxprofit;
   int* RESTRICT nodes_base = cliquepaths->nodes_base;
   int* RESTRICT nodes_pred = cliquepaths->nodes_pred;
   int* RESTRICT nodes_baseToClique = cliquepaths->nodes_baseToClique;
   STP_Bool* RESTRICT visited = dijkdata->node_visited;

   assert(nodes_dist && nodes_base && dheap && cliquenodes);
   assert(nodes_distBaseToMax && nodes_maxprofit && nodes_distFromMax);
   assert(graph_heap_isClean(dheap));
   assert(ncliquenodes >= 2);
   assert(dijkdata->edgelimit >= 0 && dijkdata->edgelimit < 20000);

#ifndef NDEBUG
   for( int i = 0; i < g->knots; i++ )
   {
      assert(UNKNOWN == dheap->position[i]);
      assert(EQ(FARAWAY, nodes_dist[i]));
      assert(!visited[i]);
   }
#endif

   /* add clique */
   for( int i = 0; i < ncliquenodes; i++ )
   {
      const int node = cliquenodes[i];

      assert(0 <= node && node < g->knots);

      visitlist[i] = node;
      visited[node] = TRUE;
      nodes_base[node] = node;
      nodes_pred[node] = node;
      nodes_dist[node] = 0.0;
      nodes_baseToClique[node] = i;
      nodes_distBaseToMax[node] = 0.0;
      nodes_maxprofit[node] = 0.0;
      nodes_distFromMax[node] = 0.0;

      graph_heap_correct(node, 0.0, dheap);
   }

   dijkdata->nvisits = ncliquenodes;
   assert(dheap->size == ncliquenodes);
}


/** returns distance limit */
static inline
SCIP_Real sdCliqueStarGetDistLimit(
   const SDCLIQUE*       cliquedata,         /**< data */
   const SCIP_Real*      sds                 /**< to be filled */
)
{
   const int ncliquenodes = cliquedata->ncliquenodes;
   const int nsds = (ncliquenodes * (ncliquenodes - 1)) / 2;
   SCIP_Real limit = 0.0;

   assert(nsds >= 1);

   for( int i = 0; i < nsds; i++ )
   {
      assert(GE(sds[i], 0.0) && LT(sds[i], FARAWAY));

      if( sds[i] > limit )
         limit = sds[i];
   }

   assert(GT(limit, 0.0));

   return limit;
}


/** computes SDs */
static
void sdCliqueStarComputeSds(
   const GRAPH*          g,                  /**< graph data structure */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   SDCLIQUE*             cliquedata,         /**< data */
   CLIQUEPATHS*          cliquepaths,        /**< paths data */
   SCIP_Real* RESTRICT   sds                 /**< to be filled */
)
{
   DIJK* RESTRICT dijkdata = cliquedata->dijkdata;
   SCIP_Real* RESTRICT nodes_dist = dijkdata->node_distance;
   DHEAP* RESTRICT dheap = dijkdata->dheap;
   STP_Bool* RESTRICT visited = dijkdata->node_visited;
   int* RESTRICT visitlist = dijkdata->visitlist;
   int* const state = dheap->position;
   const SCIP_Real* const gCost = g->cost;
   const int* const gOeat = g->oeat;
   const int* const gHead = g->head;
   const int* const nodes_baseToClique = cliquepaths->nodes_baseToClique;
   SCIP_Real* RESTRICT nodes_distBaseToMax = cliquepaths->nodes_distBaseToMax;
   SCIP_Real* RESTRICT nodes_distFromMax = cliquepaths->nodes_distFromMax;
   SCIP_Real* RESTRICT nodes_maxprofit = cliquepaths->nodes_maxprofit;
   int* RESTRICT nodes_base = cliquepaths->nodes_base;
   int* RESTRICT nodes_pred = cliquepaths->nodes_pred;
   const int ncliquenodes = cliquedata->ncliquenodes;
   int nvisits = dijkdata->nvisits;
   const SCIP_Real distlimit = sdCliqueStarGetDistLimit(cliquedata, sds);
   const SCIP_Bool useProfit = (sdprofit != NULL);
   const int centernode = cliquedata->centernode;
   int nchecks = 0;
   const int limit = dijkdata->edgelimit;

   assert(g->knots > 1);
   assert(dheap->size > 1);

   /* until the heap is empty */
   while( dheap->size > 0 && nchecks <= limit )
   {
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_base = nodes_base[k];
      const int k_predNode = nodes_pred[k];
      const SCIP_Real k_dist = nodes_dist[k];

      assert(0 <= k_base && k_base < g->knots);
      assert(CONNECT == state[k]);
      assert(CONNECT == state[k_base]);

      for( int e = g->outbeg[k]; e >= 0; e = gOeat[e] )
      {
         const int m = gHead[e];
         SCIP_Real profit = 0.0;
         SCIP_Real bias = 0.0;
         SCIP_Real newdist = k_dist + gCost[e];

         nchecks++;

         /* NOTE: need to make sure that we do not go over the center of the clique!
          * todo: Might be an issue if we pseudo-eliminate edges...probably need to block the edges as well */
         if( useProfit && k != centernode && k != k_predNode && m != k_predNode && m != centernode )
         {
            profit = reduce_sdprofitGetProfit(sdprofit, k, k_predNode, m);
            bias = MIN(gCost[e], profit);
            bias = MIN(k_dist, bias);
            newdist -= bias;

            assert(k != k_base || EQ(MIN(k_dist, bias), 0.0));
            assert(GE(newdist, 0.0));
         }


         if( GT(newdist, distlimit) )
            continue;

         /* first time visit of m? */
         if( UNKNOWN == state[m] )
         {
            assert(!visited[m]);

            sdCliqueStarUpdateNodeDist(m, k, newdist, dheap, nodes_dist, nodes_base, nodes_pred);
            sdCliqueStarUpdateNodeMaxDist(m, k, nodes_dist, profit, gCost[e], nodes_distBaseToMax,
                  nodes_distFromMax, nodes_maxprofit);

            visitlist[nvisits++] = m;
            visited[m] = TRUE;

            assert(k_base == nodes_base[m]);
            continue;
         }

         assert(visited[m]);

         /* can we update the special distances? */
         if( nodes_base[m] != k_base )
         {
            sdCliqueStarUpdateSd(sdprofit, nodes_pred, nodes_distBaseToMax, nodes_distFromMax, nodes_maxprofit,
                  ncliquenodes, centernode, m, k, newdist, gCost[e], nodes_dist, nodes_base, nodes_baseToClique, sds);
         }

         if( state[m] != CONNECT )
         {
            /* check whether the path (to m) including k is shorter than the so far best known */
            if( nodes_dist[m] > newdist )
            {
               sdCliqueStarUpdateNodeDist(m, k, newdist, dheap, nodes_dist, nodes_base, nodes_pred);
               sdCliqueStarUpdateNodeMaxDist(m, k, nodes_dist, profit, gCost[e], nodes_distBaseToMax,
                   nodes_distFromMax, nodes_maxprofit);
            }
         }
      }
   }

   dijkdata->nvisits = nvisits;
}


/*
 * Interface methods
 */

/** limited Dijkstra, stopping at terminals */
void graph_sdStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             with_zero_edges,    /**< telling name */
   int                   star_root,          /**< root of the start */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int*                  star_base,          /**< star base node, must be initially set to SDSTAR_BASE_UNSET */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool*            success             /**< will be set to TRUE iff at least one edge can be deleted */
   )
{
   int nchecks;
   int nstarhits;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int star_degree = range_csr[star_root].end - range_csr[star_root].start;
   SCIP_Real distlimit;
   /* NOTE: with zero edges case is already covered with state[k] = UNKNOWN if k == star_base[k] */
   const SCIP_Real eps = graph_pc_isPcMw(g) ? 0.0 : SCIPepsilon(scip);

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap && success);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->mark[star_root] && star_degree >= 1);
   assert(dheap->size == 0);
   assert(edgelimit >= 1);

   *nvisits = 0;
   *success = FALSE;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(dist[k] == FARAWAY);
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(state[k] == UNKNOWN);
   }
#endif

   distlimit = 0.0;
   dist[star_root] = 0.0;
   state[star_root] = CONNECT;
   visitlist[(*nvisits)++] = star_root;

   for( int e = range_csr[star_root].start, end = range_csr[star_root].end; e < end; e++ )
   {
      const int m = head_csr[e];

      assert(g->mark[m]);
      assert(!visited[m]);

      visitlist[(*nvisits)++] = m;
      visited[m] = TRUE;
      dist[m] = cost_csr[e];
      star_base[m] = m;

      /*  add epsilon to make sure that m is removed from the heap last in case of equality */
      graph_heap_correct(m, cost_csr[e] + eps, dheap);

      if( cost_csr[e] > distlimit )
         distlimit = cost_csr[e];
   }


   nchecks = 0;
   nstarhits = 0;

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != star_root);
      assert(state[k] == CONNECT);
      assert(LE(dist[k], distlimit));

      if( with_zero_edges && k == star_base[k] )
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         assert(g->mark[m] && star_base[k] >= 0);

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( GT(distnew, distlimit) )
               continue;

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               if( star_base[m] == m )
                  nstarhits++;

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }
            else if( EQ(distnew, dist[m]) && star_base[m] == m )
            {
               if( with_zero_edges && star_base[k] == star_base[m] )
                  continue;

               assert(visited[m]);
               nstarhits++;

               assert(star_base[m] != star_base[k]);

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }

            /* all star nodes hit already? */
            if( nstarhits == star_degree )
            {
               nchecks = edgelimit + 1;
               break;
            }
         }
         nchecks++;
      }
   }

  *success = (nstarhits > 0);
}


/** limited Dijkstra with node bias */
SCIP_RETCODE graph_sdStarBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   int                   star_root,          /**< root of the start */
   int*                  star_base,          /**< star base node, must be initially set to SDSTAR_BASE_UNSET */
   DIJK*                 dijkdata,           /**< Dijkstra data */
   SCIP_Bool*            success             /**< will be set to TRUE iff at least one edge can be deleted */
   )
{
   int nchecks;
   int nstarhits;
   const int nnodes = graph_get_nNodes(g);
   SCIP_Real* RESTRICT dist = dijkdata->node_distance;
   int* RESTRICT visitlist = dijkdata->visitlist;
   STP_Bool* RESTRICT visited = dijkdata->node_visited;
   int* node_predNode;
   DHEAP* dheap = dijkdata->dheap;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int star_degree = range_csr[star_root].end - range_csr[star_root].start;
   int nvisits = 0;
   SCIP_Real distlimit;
   const int edgelimit = dijkdata->edgelimit;
   /* NOTE: with zero edges case is already covered with state[k] = UNKNOWN if k == star_base[k] */
   const SCIP_Real eps = graph_pc_isPcMw(g) ? 0.0 : 2.0 * SCIPepsilon(scip);

   assert(dcsr && dist && visitlist && visited && dheap && success && sdprofit);
   assert(!g->extended);
   assert(g->mark[star_root] && star_degree >= 1);
   assert(dheap->size == 0);
   assert(edgelimit >= 1);

   nvisits = 0;
   *success = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &node_predNode, nnodes) );

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(dist[k] == FARAWAY);
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(state[k] == UNKNOWN);
      node_predNode[k] = UNKNOWN;
   }
#endif

   distlimit = 0.0;
   dist[star_root] = 0.0;
   state[star_root] = CONNECT;
   visitlist[(nvisits)++] = star_root;

   for( int e = range_csr[star_root].start, end = range_csr[star_root].end; e < end; e++ )
   {
      const int m = head_csr[e];

      assert(g->mark[m]);
      assert(!visited[m]);

      visitlist[(nvisits)++] = m;
      visited[m] = TRUE;
      dist[m] = cost_csr[e];
      star_base[m] = m;
      node_predNode[m] = star_root;

      /*  add epsilon to make sure that m is removed from the heap last in case of equality */
      graph_heap_correct(m, cost_csr[e] + eps, dheap);

      if( cost_csr[e] > distlimit )
         distlimit = cost_csr[e];
   }

   nchecks = 0;
   nstarhits = 0;

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;
      const int k_pred = node_predNode[k];

      assert(k != star_root);
      assert(k_pred >= 0 && k_pred < nnodes);
      assert(state[k] == CONNECT);
      assert(LE(dist[k], distlimit));

      if( k == star_base[k] )
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m] && star_base[k] >= 0);

         if( state[m] != CONNECT )
         {
            SCIP_Real distnew = dist[k] + cost_csr[e];
            if( m != k_pred )
            {
               SCIP_Real profitBias = reduce_sdprofitGetProfit(sdprofit, k, m, k_pred);
               profitBias = MIN(profitBias, cost_csr[e]);
               profitBias = MIN(profitBias, dist[k]);
               distnew  -= profitBias;
            }

            if( GT(distnew, distlimit) )
               continue;

            if( LT(distnew, dist[m]) )
            {
               if( !visited[m] )
               {
                  visitlist[(nvisits)++] = m;
                  visited[m] = TRUE;
               }

               if( star_base[m] == m )
                  nstarhits++;

               node_predNode[m] = k;
               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m && m != star_root);
            }
            else if( EQ(distnew, dist[m]) && star_base[m] == m )
            {
               if( star_base[k] == star_base[m] )
                  continue;

               assert(visited[m]);
               nstarhits++;

               assert(star_base[m] != star_base[k]);

               node_predNode[m] = k;
               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m && m != star_root);
            }

            /* all star nodes hit already? */
            if( nstarhits == star_degree )
            {
               nchecks = edgelimit + 1;
               break;
            }
         }
         nchecks++;
      }
   }

  dijkdata->nvisits = nvisits;
  *success = (nstarhits > 0);

  SCIPfreeBufferArray(scip, &node_predNode);

  return SCIP_OKAY;
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            sdwalkCorrectHeap(heap, state, &count, dist, m, newcost);
         }
         else
         {
            sdwalkCorrectHeap(heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               sdwalkCorrectHeap(heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT && m != start )
         {
            SCIP_Real distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksTriangle(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const int*            stateprev,          /**< state of previous run or NULL */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            prizeoffset,        /**< array for storing prize offset or NULL */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         if( stateprev && stateprev[m] == CONNECT )
            continue;

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);

            if( prizeoffset )
            {
               if( g->prize[m] > cost_csr[e] )
               {
                  prizeoffset[m] = cost_csr[e];
                  assert(SCIPisZero(scip, newcost));
               }
               else
                  prizeoffset[m] = g->prize[m];
            }
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT )
         {
            SCIP_Real distnew;

            assert(m != start);

            if( stateprev && stateprev[m] == CONNECT )
                continue;

            distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( distnew > distlimit )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( prizeoffset && termmark[m] != 0 )
               {
                  const SCIP_Real distnew0 = dist[k] + cost_csr[e];

                  if( g->prize[m] > distnew0 )
                  {
                     prizeoffset[m] = distnew0;
                     assert(SCIPisZero(scip, distnew));
                  }
                  else
                     prizeoffset[m] = g->prize[m];
               }

               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}

/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;
         sdwalkUpdate(g, m, start, maxnprevs, prevterms, nprevterms);

         if( Is_term(g->term[m]) )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            sdwalkCorrectHeap(heap, state, &count, dist, m, newcost);
         }
         else
         {
            sdwalkCorrectHeap(heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( Is_term(g->term[m]) )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               if( Is_term(g->term[m]) && sdwalkHasConflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalkUpdate(g, m, k, maxnprevs, prevterms, nprevterms);
               sdwalkCorrectHeap(heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt2(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges,         /**< number of previous edges */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         SCIP_Real distnew = cost[e];

         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
            distnew = MAX(distnew - g->prize[m], 0.0);

         sdwalkUpdate2(termmark, m, start, e, maxnprevs, SCIPisZero(scip, distnew),
               prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
         sdwalkCorrectHeap(heap, state, &count, dist, m, distnew);

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = sdwalkGetdistnewEdge(prevedges, nprevedges, cost, dist, k, e, maxnprevs);

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = sdwalkGetdistnewPrize(prevNPterms, nprevNPterms, termmark, visited, g->prize, k, m, distnew, maxnprevs);

            if( SCIPisLT(scip, distnew, dist[m]) )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               /* continue if m is proper terminals and is on the walk to k */
               if( termmark[m] == 2 && sdwalkHasConflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalkUpdate2(termmark, m, k, e, maxnprevs, SCIPisZero(scip, distnew),
                     prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
               sdwalkCorrectHeap(heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}


/** modified Dijkstra along walks for PcMw */
SCIP_Bool graph_sdWalksConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const SCIP_Real*      cost,               /**< edge costs */
   const STP_Bool*       endpoint,           /**< stores whether search should be ended at vertex */
   int                   start,              /**< start vertex */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool             resetarrays         /**< should arrays be reset? */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   int count;
   int nchecks;
   const SCIP_Real prize = g->prize[start];

   assert(heap != NULL);
   assert(state != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(Is_term(g->term[start]));
   assert(g->grad[start] > 0);
   assert(g->mark[start]);

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(state[k] == UNKNOWN);
      assert(visited[k] == FALSE);
      assert(dist[k] == FARAWAY);
   }
#endif

   *nvisits = 0;
   nchecks = 0;
   count = 1;
   heap[count] = start;
   state[start] = count;
   dist[start] = 0.0;
   visitlist[(*nvisits)++] = start;
   g->mark[start] = FALSE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(SCIPisLE(scip, dist[k], prize));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, prize) )
               continue;

            if( termmark[m] != 0 )
               distnew -= g->prize[m];

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( endpoint != NULL && endpoint[m] )
               {
                  g->mark[start] = TRUE;
                  if( resetarrays )
                     sdwalkReset(*nvisits, visitlist, dist, state, visited);

                  return TRUE;
               }

               sdwalkCorrectHeap(heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;

   if( resetarrays )
      sdwalkReset(*nvisits, visitlist, dist, state, visited);

   return FALSE;
}


/** computes (or rather updates) SDs between all */
SCIP_RETCODE graph_sdComputeCliqueStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   SDCLIQUE*             cliquedata          /**< data */
)
{
   CLIQUEPATHS cliquepaths = { NULL, NULL, NULL, NULL, NULL, NULL };
   SCIP_Real* sds = cliquedata->sds;

   assert(scip && g && cliquedata);
   assert(cliquedata->dijkdata && cliquedata->cliquenodes);
   assert(cliquedata->ncliquenodes >= 2);
   assert(g->knots >= cliquedata->ncliquenodes);

   SCIP_CALL( cliquePathsInitData(scip, g, &cliquepaths) );

   sdCliqueStarInit(g, cliquedata, &cliquepaths);
   sdCliqueStarComputeSds(g, sdprofit, cliquedata, &cliquepaths, sds);

   cliquePathsFreeData(scip, &cliquepaths);

   return SCIP_OKAY;
}
