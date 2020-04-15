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

/**@file   shortestpaths.c
 * @brief  Shortest path based algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various shortest path based algorithms.
 * Note: This file is supposed to replace graph_path.c in the long run, as it includes the faster implementations.
 *
 */

#include "shortestpath.h"
#include "portab.h"


#ifndef NDEBUG
/** all terminals reached? */
static
SCIP_Bool computeSteinerTree_allTermsAreReached(
   const GRAPH*          g,                  /**< graph data structure */
   const STP_Bool*       connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   const int nnodes = graph_get_nNodes(g);

   for( int k = 0; k < nnodes; ++k )
   {
      if( Is_term(g->term[k]) && !connected[k] )
      {
         SCIPdebugMessage("terminal %d not connected! \n", k);
         return FALSE;
      }
   }

   return TRUE;
}

/** all pseudo terminals reached? */
static
SCIP_Bool computeSteinerTree_allPseudoTermsAreReached(
   const GRAPH*          g,                  /**< graph data structure */
   const STP_Bool*       connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   const int nnodes = graph_get_nNodes(g);

   assert(graph_pc_isPcMw(g));

   for( int k = 0; k < nnodes; ++k )
   {
      if( Is_pseudoTerm(g->term[k]) && !connected[k] )
      {
         return FALSE;
      }
   }

   return TRUE;
}
#endif


/** initializes */
static inline
void computeSteinerTree_init(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   const int nnodes = graph_get_nNodes(g);
   const CSR* csr = spaths->csr;
   DHEAP* dheap = spaths->dheap;
   SCIP_Real* RESTRICT nodes_dist = spaths->nodes_dist;
   int* RESTRICT nodes_pred = spaths->nodes_pred;
   STP_Bool* RESTRICT connected = spaths->nodes_isConnected;

   assert(startnode >= 0 && startnode < g->knots);
   assert(nnodes >= 1);
   assert(csr->nnodes == nnodes);
   assert(csr->start[nnodes] <= g->edges);

   graph_heap_clean(TRUE, dheap);

#ifndef NDEBUG
   assert(dheap->size == 0);
   for( int k = 0; k < nnodes; k++ )
   {
      assert(dheap->position[k] == UNKNOWN);
      nodes_pred[k] = -1;
   }
#endif

   for( int k = 0; k < nnodes; k++ )
   {
      nodes_dist[k] = FARAWAY;
      connected[k] = FALSE;
   }

   nodes_dist[startnode] = 0.0;
   nodes_pred[startnode] = -1;
   graph_heap_correct(startnode, 0.0, dheap);
   connected[startnode] = TRUE;
}


/** connects node to current tree */
static inline
void computeSteinerTree_connectTerminal(
   const GRAPH*          g,                  /**< graph data structure */
   int                   k,                  /**< vertex to connect */
   const int*            nodes_pred,         /**< predecessor array (on vertices) */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   assert(k >= 0 && k < g->knots);
   assert(!connected[k]);

   connected[k] = TRUE;
   nodes_dist[k] = 0.0;

   SCIPdebugMessage("connect terminal %d \n", k);

   /* connect k to current solution */
   for( int node = nodes_pred[k]; !connected[node]; node = nodes_pred[node] )
   {
      assert(node >= 0 && node < g->knots);
      assert(!connected[node]);
      assert(!Is_term(g->term[node]));

      SCIPdebugMessage("connect path node %d \n", node);

      connected[node] = TRUE;
      nodes_dist[node] = 0.0;
      graph_heap_correct(node, 0.0, dheap);
   }
}


/** executes */
static inline
void computeSteinerTree_exec(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   const CSR* const csr = spaths->csr;
   DHEAP* const dheap = spaths->dheap;
   SCIP_Real* RESTRICT nodes_dist = spaths->nodes_dist;
   int* RESTRICT nodes_pred = spaths->nodes_pred;
   STP_Bool* RESTRICT connected = spaths->nodes_isConnected;
   int* const state = dheap->position;
   const SCIP_Real* const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   int termscount = 0;

   assert(dheap->size == 1);
   assert(connected[startnode]);

   if( Is_term(g->term[startnode]) )
      termscount++;

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = start_csr[k];
      const int k_end = start_csr[k + 1];
      register SCIP_Real k_dist;

      SCIPdebugMessage("take node %d from queue \n", k);

      assert(state[k] == CONNECT);
      /* NOTE: needs to be set for invariant of heap */
      state[k] = UNKNOWN;

      if( Is_term(g->term[k]) && k != startnode )
      {
         assert(!connected[k]);
         computeSteinerTree_connectTerminal(g, k, nodes_pred, nodes_dist, dheap, connected);

         assert(termscount < g->terms);

         /* have all terminals been reached? */
         if( ++termscount == g->terms )
         {
            break;
         }
      }

      k_dist = nodes_dist[k];

      for( int e = k_start; e != k_end; e++ )
      {
         const int m = head_csr[e];

         if( !connected[m] )
         {
            const SCIP_Real distnew = k_dist + cost_csr[e];

            /* closer to k than to current predecessor? */
            if( distnew < nodes_dist[m] )
            {
               nodes_pred[m] = k;
               nodes_dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }
}


/** connects node to current tree */
static inline
void computeSteinerTree_connectPseudoTerm(
   const GRAPH*          g,                  /**< graph data structure */
   int                   k,                  /**< vertex to connect */
   const int*            nodes_pred,         /**< predecessor array (on vertices) */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   SPATHSPC*             spaths_pc,          /**< PC/MW shortest paths data */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int* RESTRICT         nPseudoTerms        /**< pointer */
)
{
   assert(k >= 0 && k < g->knots);
   assert(!connected[k]);
   assert(Is_pseudoTerm(g->term[k]));

   connected[k] = TRUE;
   nodes_dist[k] = 0.0;
   shortestpath_pcConnectNode(g, connected, k, spaths_pc);
   (*nPseudoTerms)++;

   SCIPdebugMessage("connect pseudo terminal %d \n", k);

   /* connect k to current solution */
   for( int node = nodes_pred[k]; !connected[node]; node = nodes_pred[node] )
   {
      assert(node >= 0 && node < g->knots);
      assert(!connected[node]);
      assert(!Is_term(g->term[node]));

      SCIPdebugMessage("connect path node %d \n", node);

      if( Is_pseudoTerm(g->term[node]) )
      {
         shortestpath_pcConnectNode(g, connected, node, spaths_pc);
         (*nPseudoTerms)++;
      }

      connected[node] = TRUE;
      nodes_dist[node] = 0.0;
      graph_heap_correct(node, 0.0, dheap);
   }
}

/** connects node to current tree */
static inline
void computeSteinerTree_tryConnectNodePcMw(
   const GRAPH*          g,                  /**< graph data structure */
   int                   k,                  /**< vertex to connect */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Bool             costIsBiased,       /**< is cost biased? */
   const int*            nodes_pred,         /**< predecessor array (on vertices) */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SPATHSPC*             spaths_pc,          /**< PC/MW shortest paths data */
   SCIP_Bool*            isConnected,        /**< node connected to tree? */
   int* RESTRICT         nPseudoTerms        /**< pointer */
)
{
   SCIP_Bool connectK = FALSE;

   /* if k is positive vertex and close enough, connect k to current subtree */
   if( !connected[k] && Is_pseudoTerm(g->term[k]) )
   {
      connectK = (GE(prize[k], nodes_dist[k]));

      /* maybe if we count the prizes on the path, the extension becomes profitable? */
      if( !connectK )
      {
         SCIP_Real prizesum = 0.0;
         const SCIP_Bool isPc = graph_pc_isPc(g);

         for( int node = nodes_pred[k]; !connected[node]; node = nodes_pred[node] )
         {
            if( Is_pseudoTerm(g->term[node]) )
               prizesum += prize[node];
            else if( isPc && !costIsBiased && graph_pc_knotIsNonLeafTerm(g, node) )
               prizesum += prize[node];
         }

         assert(prizesum >= 0.0 && LT(prizesum, FARAWAY));
         connectK = GE(prize[k] + prizesum, nodes_dist[k]);
      }

      if( connectK )
      {
         computeSteinerTree_connectPseudoTerm(g, k, nodes_pred, nodes_dist, spaths_pc, dheap, connected, nPseudoTerms);
      }
   }

   *isConnected = connectK;
}


/** executes */
static inline
void computeSteinerTree_execPcMw(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Bool             costIsBiased,       /**< is cost biased? */
   SPATHSPC*             spaths_pc,          /**< PC/MW shortest paths data */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   const CSR* const csr = spaths->csr;
   DHEAP* const dheap = spaths->dheap;
   SCIP_Real* RESTRICT nodes_dist = spaths->nodes_dist;
   int* RESTRICT nodes_pred = spaths->nodes_pred;
   STP_Bool* RESTRICT connected = spaths->nodes_isConnected;
   int* const state = dheap->position;
   const SCIP_Real* const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   int termsposCount = 0;
   const int ntermspos = graph_pc_nNonFixedTerms(g);

   assert(g->extended);
   assert(dheap->size == 1);
   assert(connected[startnode]);
   assert(graph_pc_isPcMw(g) && !graph_pc_isRootedPcMw(g));

#ifndef NDEBUG
   {
      // todo deleteme
      int posCheck = 0;
      for( int k = 0; k < g->knots; k++ )
         if( Is_pseudoTerm(g->term[k]) )
            posCheck++;

      assert(posCheck == ntermspos);
   }
#endif

   shortestpath_pcStart(spaths_pc);

   if( Is_pseudoTerm(g->term[startnode]) )
   {
      termsposCount++;
      shortestpath_pcConnectNode(g, connected, startnode, spaths_pc);
   }

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = start_csr[k];
      const int k_end = start_csr[k + 1];
      register SCIP_Real k_dist;
      SCIP_Bool k_isConnected = FALSE;

      SCIPdebugMessage("take node %d from queue \n", k);

      assert(state[k] == CONNECT);
      state[k] = UNKNOWN;

      computeSteinerTree_tryConnectNodePcMw(g, k, prize, costIsBiased, nodes_pred, nodes_dist, dheap, connected,
            spaths_pc, &k_isConnected, &termsposCount);

      k_dist = nodes_dist[k];

      if( k_isConnected && termsposCount == ntermspos )
      {
         SCIPdebugMessage("all potential terminals reached \n");
         assert(computeSteinerTree_allPseudoTermsAreReached(g, connected));
         break;
      }
      else if( !k_isConnected && k_dist > spaths_pc->maxoutprize )
      {
         SCIPdebugMessage("further extension is not profitable \n");
         break;
      }

      for( int e = k_start; e != k_end; e++ )
      {
         const int m = head_csr[e];

         if( !connected[m] )
         {
            const SCIP_Real distnew = k_dist + cost_csr[e];

            /* closer to k than to current predecessor? */
            if( distnew < nodes_dist[m] )
            {
               nodes_pred[m] = k;
               nodes_dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   assert(termsposCount == ntermspos || !computeSteinerTree_allPseudoTermsAreReached(g, connected));
}


/*
 * Interface methods
 */


/** start computation */
void shortestpath_pcStart(
   SPATHSPC*             sppc                /**< PC/MW shortest path data */
)
{
   SCIP_Real* orderedprizes = sppc->orderedprizes;

   assert(orderedprizes);

   sppc->maxoutprize = orderedprizes[0];
   sppc->maxoutprize_idx = 0;
}

/** update maximum prize */
void shortestpath_pcConnectNode(
   const GRAPH*          g,                  /**< graph data structure */
   const STP_Bool*       connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int                   node_connected,     /**< node that is removed */
   SPATHSPC*             sppc                /**< PC data */
)
{
   const SCIP_Real* orderedprizes = sppc->orderedprizes;
   const int* orderedprizes_id = sppc->orderedprizes_id;
   int maxprizeidx = sppc->maxoutprize_idx;
   assert(maxprizeidx >= 0 && maxprizeidx <= g->terms);

   /* sentinel? */
   if( orderedprizes_id[maxprizeidx] < 0 )
   {
      assert(orderedprizes_id[maxprizeidx] == -1);
      assert(EQ(sppc->maxoutprize, 0.0));
      return;
   }

   assert(maxprizeidx < g->terms);

   /* is current node at the maximum? */
   if( node_connected == orderedprizes_id[maxprizeidx] )
   {
      while( orderedprizes_id[maxprizeidx] >= 0 && connected[orderedprizes_id[maxprizeidx]] )
      {
         maxprizeidx++;
         assert(maxprizeidx <= g->terms);
      }

      sppc->maxoutprize_idx = maxprizeidx;

      if( orderedprizes_id[maxprizeidx] < 0 )
         sppc->maxoutprize = 0.0;
      else
         sppc->maxoutprize = orderedprizes[maxprizeidx];
   }
}


/** shortest path based heuristic for computing a Steiner tree  */
void shortestpath_computeSteinerTree(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   assert(g && spaths);
   assert(spaths->csr && spaths->nodes_dist && spaths->nodes_pred && spaths->dheap && spaths->nodes_isConnected);
   assert(graph_typeIsSpgLike(g));

   computeSteinerTree_init(g, startnode, spaths);

   if( g->knots == 1 )
      return;

   computeSteinerTree_exec(g, startnode, spaths);

   assert(computeSteinerTree_allTermsAreReached(g, spaths->nodes_isConnected));
}


/** shortest path based heuristic for computing a Steiner tree in PC/MW case */
void shortestpath_computeSteinerTreePcMw(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Bool             costIsBiased,       /**< is cost biased? */
   SPATHSPC*             spaths_pc,          /**< PC/MW shortest paths data */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   assert(g && spaths);
   assert(spaths->csr && spaths->nodes_dist && spaths->nodes_pred && spaths->dheap && spaths->nodes_isConnected);
   assert(graph_pc_isPcMw(g) && !graph_pc_isRootedPcMw(g));
   assert(g->extended);
   assert(!graph_pc_knotIsDummyTerm(g, startnode));

   computeSteinerTree_init(g, startnode, spaths);

   if( g->knots == 1 )
      return;

   computeSteinerTree_execPcMw(g, startnode, prize, costIsBiased, spaths_pc, spaths);
}

