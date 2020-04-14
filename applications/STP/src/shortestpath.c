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

#ifndef NDEBUG
/** all terminals reached? */
static
SCIP_Bool computeSteinerTree_allTermsReached(
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
#endif


/** initializes */
static inline
void computeSteinerTree_init(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   int* RESTRICT         nodes_pred,         /**< predecessor edge array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT const state = dheap->position;

   assert(startnode >= 0 && startnode < g->knots);
   assert(nnodes >= 1);

   graph_heap_clean(TRUE, dheap);

#ifndef NDEBUG
   assert(dheap->size == 0);

   for( int k = 0; k < nnodes; k++ )
   {
      assert(state[k] == UNKNOWN);
   }
#endif

   for( int k = 0; k < nnodes; k++ )
   {
      nodes_dist[k] = FARAWAY;
      nodes_pred[k] = -1;
      connected[k] = FALSE;
   }

   nodes_dist[startnode] = 0.0;
   graph_heap_correct(startnode, 0.0, dheap);
}


/** connects node to current tree */
static inline
void computeSteinerTree_connectNode(
   const GRAPH*          g,                  /**< graph data structure */
   int                   k,                  /**< vertex to connect */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   int* RESTRICT         nodes_pred,         /**< predecessor edge array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   assert(k >= 0 && k < g->knots);
   assert(nodes_pred[k] >= 0 && !connected[k]);

   connected[k] = TRUE;
   nodes_dist[k] = 0.0;

   /* connect k to current solution */
   for( int node = g->tail[nodes_pred[k]]; !connected[node]; node = g->tail[nodes_pred[node]] )
   {
      assert(nodes_pred[node] != -1);
      assert(!Is_term(g->term[node]));

      connected[node] = TRUE;
      nodes_dist[node] = 0.0;
      graph_heap_correct(node, 0.0, dheap);
   }
}


/** initializes */
static inline
void computeSteinerTree_exec(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   int* RESTRICT         nodes_pred,         /**< predecessor edge array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   int termscount = 0;
   int *const state = dheap->position;
   const CSR *const csr = g->csr_storage;
   const SCIP_Real *const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;

   assert(dheap->size == 1);
   assert(state[startnode] == CONNECT);

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = start_csr[k];
      const int k_end = start_csr[k + 1];
      const SCIP_Real k_dist = nodes_dist[k];

      state[k] = UNKNOWN;

      if( Is_term(g->term[k]) && k != startnode )
      {
         computeSteinerTree_connectNode(g, k, nodes_dist, nodes_pred, dheap, connected);

         assert(termscount < g->terms);

         /* have all terminals been reached? */
         if( ++termscount == g->terms )
            break;
      }

      for( int e = k_start; e != k_end; e++ )
      {
         const int m = head_csr[e];

         /* not yet in the tree? */
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

/*
 * Interface methods
 */

/** shortest path based heuristic for computing a Steiner tree */
void shortestpath_computeSteinerTree(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   SCIP_Real* RESTRICT   nodes_dist,         /**< distance array (on vertices) */
   int* RESTRICT         nodes_pred,         /**< predecessor edge array (on vertices) */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool* RESTRICT    connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   assert(g && nodes_dist && nodes_pred && dheap && connected);
   assert(g->csr_storage);
   assert(graph_typeIsSpgLike(g));
   // todo assert that csr is valid for g!

   computeSteinerTree_init(g, startnode, nodes_dist, nodes_pred, dheap, connected);

   if( g->knots == 1 )
      return;

   computeSteinerTree_exec(g, startnode, nodes_dist, nodes_pred, dheap, connected);

   assert(computeSteinerTree_allTermsReached(g, connected));
}
