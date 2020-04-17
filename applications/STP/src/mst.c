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

/**@file   mst.c
 * @brief  minimum spanning tree based algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various minimum spanning tree based algorithms.
 * Note: This file is supposed to (partly) replace graph_path.c in the long run, as it includes the faster implementations.
 *
 */


#include "mst.h"
#include "portab.h"


/** initializes */
static inline
void computeOnMarked_init(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   MST*                  mst                 /**< MST data */
)
{
   const int nnodes = graph_get_nNodes(g);
   DHEAP* dheap = mst->dheap;
   SCIP_Real* RESTRICT nodes_dist = mst->nodes_dist;
   int* RESTRICT nodes_pred = mst->nodes_predEdge;

   assert(startnode >= 0 && startnode < g->knots);
   assert(g->mark[startnode]);
   assert(nnodes >= 1);
   assert(mst->csr->nnodes == nnodes);
   assert(mst->csr->start[nnodes] <= g->edges);

   graph_heap_clean(TRUE, dheap);

   for( int k = 0; k < nnodes; k++ )
   {
      nodes_pred[k] = UNKNOWN;
      nodes_dist[k] = FARAWAY;
   }

   nodes_dist[startnode] = 0.0;
   graph_heap_correct(startnode, 0.0, dheap);
}


/** executes */
static inline
void computeOnMarked_exec(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   MST*                  mst                 /**< MST data */
)
{
   const CSR* const csr = mst->csr;
   DHEAP* const dheap = mst->dheap;
   SCIP_Real* RESTRICT nodes_dist = mst->nodes_dist;
   int* RESTRICT nodes_pred = mst->nodes_predEdge;
   int* const state = dheap->position;
   const SCIP_Real* const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   const int* const nodeIsMarked = g->mark;

   assert(dheap->size == 1);

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = start_csr[k];
      const int k_end = start_csr[k + 1];

      SCIPdebugMessage("take node %d from queue \n", k);
      assert(state[k] == CONNECT);

      for( int e = k_start; e != k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT && nodeIsMarked[m] && cost_csr[e] < nodes_dist[m] )
         {
            nodes_pred[m] = e;
            nodes_dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }
      }
   }
}


/** update predecessors */
static inline
void computeOnMarked_computePredecessors(
   const GRAPH*          g,                  /**< graph data structure */
   MST*                  mst                 /**< MST data */
)
{
   int* RESTRICT nodes_pred = mst->nodes_predEdge;
   int* RESTRICT edgeid_csr = mst->csr->edge_id;
   const int nnodes = graph_get_nNodes(g);

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_pred[i] != UNKNOWN )
      {
         assert(g->mark[i]);
         nodes_pred[i] = edgeid_csr[nodes_pred[i]];
      }
   }
}


/*
 * Interface methods
 */


/** computes MST on marked vertices */
void mst_computeOnMarked(
   const GRAPH*          g,                  /**< graph data structure */
   int                   startnode,          /**< start vertex */
   MST*                  mst                 /**< MST data */
)
{
   assert(g && mst);
   assert(mst->csr && mst->dheap && mst->nodes_dist && mst->nodes_predEdge);
   assert(graph_typeIsSpgLike(g));

   computeOnMarked_init(g, startnode, mst);

   if( g->knots == 1 )
      return;

   computeOnMarked_exec(g, startnode, mst);
   computeOnMarked_computePredecessors(g, mst);
}
