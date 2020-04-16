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


/** all fixed terminals reached? */
static
SCIP_Bool computeSteinerTree_allFixedTermsAreReached(
   const GRAPH*          g,                  /**< graph data structure */
   const STP_Bool*       connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   const int nnodes = graph_get_nNodes(g);

   assert(graph_pc_isPcMw(g));

   for( int k = 0; k < nnodes; ++k )
   {
      if( graph_pc_knotIsFixedTerm(g, k) && !connected[k] )
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
   DHEAP* dheap = spaths->dheap;
   SCIP_Real* RESTRICT nodes_dist = spaths->nodes_dist;
   int* RESTRICT nodes_pred = spaths->nodes_pred;
   STP_Bool* RESTRICT connected = spaths->nodes_isConnected;

   assert(startnode >= 0 && startnode < g->knots);
   assert(nnodes >= 1);
   assert(spaths->csr->nnodes == nnodes);
   assert(spaths->csr->start[nnodes] <= g->edges);

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


/** connects (also PC/MW potential) terminal to current tree */
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
      assert(!Is_term(g->term[node]) && !Is_pseudoTerm(g->term[node]));

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

   shortestpath_pcReset(spaths_pc);

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



/** executes */
static inline
void computeSteinerTree_execRpcMw(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
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
   const int nfixedterms = graph_pc_nFixedTerms(g);
   int fixedtermsCount = 0;

   assert(g->extended);
   assert(dheap->size == 1 && connected[startnode]);
   assert(graph_pc_isRootedPcMw(g));

   shortestpath_pcReset(spaths_pc);

   if( Is_pseudoTerm(g->term[startnode]) )
      shortestpath_pcConnectNode(g, connected, startnode, spaths_pc);

   if( Is_term(g->term[startnode]) )
   {
      assert(graph_pc_knotIsFixedTerm(g, startnode));
      fixedtermsCount++;
   }

   while( dheap->size > 0 )
   {
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = start_csr[k];
      const int k_end = start_csr[k + 1];
      register SCIP_Real k_dist;

      assert(state[k] == CONNECT);
      state[k] = UNKNOWN;

      /* if k is fixed terminal positive vertex and close enough, connect its path to current subtree */
      if( Is_anyTerm(g->term[k]) && (Is_term(g->term[k]) || GE(prize[k], nodes_dist[k])) && !connected[k] )
      {
         assert(!graph_pc_knotIsDummyTerm(g, k));
         assert(graph_pc_knotIsFixedTerm(g, k) || GE(prize[k], nodes_dist[k]));
         assert(graph_pc_knotIsFixedTerm(g, k) == Is_term(g->term[k]));

         if( Is_term(g->term[k]) )
            fixedtermsCount++;
         else if( Is_pseudoTerm(g->term[k]) )
            shortestpath_pcConnectNode(g, connected, k, spaths_pc);

         connected[k] = TRUE;
         nodes_dist[k] = 0.0;

         assert(nodes_pred[k] >= 0);

         for( int node = nodes_pred[k]; !connected[node]; node = nodes_pred[node] )
         {
            assert(node >= 0 && node < g->knots);
            SCIPdebugMessage("connect path node %d \n", node);

            if( Is_pseudoTerm(g->term[node]) )
               shortestpath_pcConnectNode(g, connected, node, spaths_pc);

            connected[node] = TRUE;
            nodes_dist[node] = 0.0;
            graph_heap_correct(node, 0.0, dheap);
         }
      }

      k_dist = nodes_dist[k];
      assert(!graph_pc_knotIsFixedTerm(g, k) || EQ(k_dist, 0.0));

      if( fixedtermsCount >= nfixedterms && k_dist > spaths_pc->maxoutprize )
      {
         SCIPdebugMessage("all fixed terminals reached \n");
         assert(fixedtermsCount == nfixedterms);
         break;
      }

      for( int e = k_start; e != k_end; e++ )
      {
         const int m = head_csr[e];

         if( !connected[m] )
         {
            const SCIP_Real distnew = k_dist + cost_csr[e];

            if( distnew < nodes_dist[m] )
            {
               nodes_pred[m] = k;
               nodes_dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   assert(computeSteinerTree_allFixedTermsAreReached(g, connected));
}


/** executes */
static inline
void computeSteinerTree_execPcMwFull(
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
   const int nterms = graph_pc_isRootedPcMw(g) ? g->terms : g->terms - 1;

   assert(g->extended && graph_pc_isPcMw(g));
   assert(dheap->size == 1);
   assert(connected[startnode]);

   if( Is_term(g->term[startnode]) || Is_pseudoTerm(g->term[startnode]) )
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
      state[k] = UNKNOWN;

      /* connect k to current subtree? */
      if( !connected[k] && (Is_term(g->term[k]) || Is_pseudoTerm(g->term[k])) )
      {
         computeSteinerTree_connectTerminal(g, k, nodes_pred, nodes_dist, dheap, connected);

         assert(termscount < nterms);

         /* have all terminals been reached? */
         if( ++termscount == nterms )
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

   assert(!graph_pc_isRootedPcMw(g) || computeSteinerTree_allFixedTermsAreReached(g, connected));
}


/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE shortestpath_pcInit(
    SCIP*                scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      costs,              /**< cost for edges (might be negative for MWCS or RMWCS) */
   const SCIP_Real*      prizes,             /**< prizes for all nodes */
   SPATHSPC**            sppc                /**< PC/MW shortest path data */
   )
{
   SCIP_Real* orderedprizes;
   int* orderedprizes_id;
   const int nnodes = graph_get_nNodes(graph);
   const int nterms = graph->terms;
   const SCIP_Bool usecosts = (costs && (graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP));
   int termcount;

   assert(graph_pc_isPcMw(graph));
   assert(prizes && graph && orderedprizes && orderedprizes_id);
   assert(graph->extended);
   assert(nterms >= 1);

   SCIP_CALL( SCIPallocMemoryArray(scip, &orderedprizes, nterms + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &orderedprizes_id, nterms + 1) );
   SCIP_CALL( SCIPallocMemory(scip, sppc) );

   termcount = 0;
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) )
      {
         orderedprizes[termcount] = prizes[k];

         /* consider incoming negative arcs */
         if( usecosts )
         {
            SCIP_Real mincost = 0.0;
            for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
               if( costs[e] < mincost )
                  mincost = costs[e];

            if( mincost < 0.0 )
               orderedprizes[termcount] -= mincost;
         }

         orderedprizes_id[termcount++] = k;
      }
   }

   for( int k = termcount; k < nterms; k++ )
   {
      orderedprizes[k] = 0.0;
      orderedprizes_id[k] = -1;
   }

   SCIPsortDownRealInt(orderedprizes, orderedprizes_id, nterms);

   /* set sentinel */
   orderedprizes[nterms] = 0.0;
   orderedprizes_id[nterms] = -1;

   assert(orderedprizes[0] >= orderedprizes[nterms - 1]);

   (*sppc)->orderedprizes = orderedprizes;
   (*sppc)->orderedprizes_id = orderedprizes_id;
   (*sppc)->maxoutprize_idx = -1;
   (*sppc)->maxoutprize = -FARAWAY;

   return SCIP_OKAY;
}


/** frees */
void shortestpath_pcFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SPATHSPC**            sppc                /**< PC/MW shortest path data */
   )
{
   SCIPfreeMemoryArray(scip, &((*sppc)->orderedprizes_id));
   SCIPfreeMemoryArray(scip, &((*sppc)->orderedprizes));

   SCIPfreeMemory(scip, sppc);
}


/** start computation */
void shortestpath_pcReset(
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


/** shortest path based heuristic for computing a Steiner tree in rooted PC/MW case */
void shortestpath_computeSteinerTreeRpcMw(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SPATHSPC*             spaths_pc,          /**< PC/MW shortest paths data */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   assert(g && spaths);
   assert(spaths->csr && spaths->nodes_dist && spaths->nodes_pred && spaths->dheap && spaths->nodes_isConnected);
   assert(graph_pc_isRootedPcMw(g) && g->extended);
   assert(!graph_pc_knotIsDummyTerm(g, startnode));

   computeSteinerTree_init(g, startnode, spaths);

   if( g->knots == 1 )
      return;

   computeSteinerTree_execRpcMw(g, startnode, prize, spaths_pc, spaths);
}


/** shortest path based heuristic for computing a Steiner tree in PC/MW case
 *  that contains all (potential and fixed) terminals */
void shortestpath_computeSteinerTreePcMwFull(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SPATHS*               spaths              /**< shortest paths data */
)
{
   assert(g && spaths);
   assert(spaths->csr && spaths->nodes_dist && spaths->nodes_pred && spaths->dheap && spaths->nodes_isConnected);
   assert(graph_pc_isPcMw(g) && g->extended);
   assert(!graph_pc_knotIsDummyTerm(g, startnode));

   computeSteinerTree_init(g, startnode, spaths);

   if( g->knots == 1 )
      return;

   computeSteinerTree_execPcMwFull(g, startnode, spaths);
}
