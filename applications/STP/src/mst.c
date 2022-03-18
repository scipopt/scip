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
   const STP_Bool*       nodes_isMarked,     /**< should the node be considered? */
   int                   startnode,          /**< start vertex */
   MST*                  mst                 /**< MST data */
)
{
   const int nnodes = graph_get_nNodes(g);
   DHEAP* dheap = mst->dheap;
   SCIP_Real* RESTRICT nodes_dist = mst->nodes_dist;
   int* RESTRICT nodes_pred = mst->nodes_predEdge;

   assert(startnode >= 0 && startnode < g->knots);
   assert(nodes_isMarked[startnode]);
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
   const STP_Bool*       nodes_isMarked,     /**< should the node be considered? */
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

         if( state[m] != CONNECT && nodes_isMarked[m] && LT(cost_csr[e], nodes_dist[m]) )
         {
            nodes_pred[m] = e;
            nodes_dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }
      }
   }

   assert(nodes_pred[startnode] == UNKNOWN);
}



/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE mst_init(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   MST**                 minspantree         /**< MST data */
)
{
   MST* mst;
   CSR* csr;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   SCIP_CALL( SCIPallocMemory(scip, minspantree) );
   mst = *minspantree;

   SCIP_CALL( SCIPallocMemoryArray(scip, &mst->nodes_dist, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mst->nodes_predEdge, nnodes) );

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &mst->dheap) );
   SCIP_CALL( graph_csr_allocWithEdgeId(scip, nnodes, nedges, &csr) );
   graph_csr_build(g, g->cost, csr);
   mst->csr = csr;

   return SCIP_OKAY;
}


/** frees */
void mst_free(
   SCIP*                 scip,               /**< SCIP data structure */
   MST**                 minspantree         /**< MST data */
)
{
   MST* mst;
   CSR* csr;
   assert(scip && minspantree);

   mst = *minspantree;
   /* :/ */
   csr = (CSR*) mst->csr;

   graph_csr_free(scip, &csr);
   graph_heap_free(scip, TRUE, TRUE, &mst->dheap);
   SCIPfreeMemoryArray(scip, &mst->nodes_predEdge);
   SCIPfreeMemoryArray(scip, &mst->nodes_dist);

   SCIPfreeMemory(scip, minspantree);
}


/** gets objective value of MST */
SCIP_Real mst_getObj(
   const GRAPH*          g,                  /**< graph data structure */
   const MST*            mst                 /**< MST data */
)
{
   SCIP_Real obj = 0.0;
   const int nnodes = graph_get_nNodes(g);
   const int* const nodes_pred = mst->nodes_predEdge;
   const CSR* const csr = mst->csr;

   assert(nodes_pred && csr);

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_pred[i] != UNKNOWN )
      {
         assert(nodes_pred[i] >= 0);
         assert(GE(csr->cost[nodes_pred[i]], 0.0) && LE(csr->cost[nodes_pred[i]], FARAWAY));

         obj += csr->cost[nodes_pred[i]];
      }
   }

   return obj;
}


/** gets solution edges */
void mst_getSoledges(
   const GRAPH*          g,                  /**< graph data structure */
   const MST*            mst,                /**< MST data */
   int* RESTRICT         soledges            /**< to be filled (CONNECT/UNKNOWN) */
)
{
   const int* const nodes_pred = mst->nodes_predEdge;
   const int* const edgeid_csr = mst->csr->edge_id;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   for( int e = 0; e < nedges; e++ )
      soledges[e] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
   {
      const int edge_csr = nodes_pred[i];
      if( edge_csr != UNKNOWN )
      {
         assert(soledges[edgeid_csr[edge_csr]] == UNKNOWN);
         soledges[edgeid_csr[edge_csr]] = CONNECT;
      }
   }
}


/** computes MST on marked vertices */
void mst_computeOnMarked(
   const GRAPH*          g,                  /**< graph data structure */
   const STP_Bool*       nodes_isMarked,     /**< should the node be considered? */
   int                   startnode,          /**< start vertex */
   MST*                  mst                 /**< MST data */
)
{
   assert(g && mst);
   assert(mst->csr && mst->dheap && mst->nodes_dist && mst->nodes_predEdge);
   assert(g->stp_type == STP_MWCSP || g->stp_type == STP_PCSPG || startnode == g->source);

   computeOnMarked_init(g, nodes_isMarked, startnode, mst);

   if( g->knots == 1 )
      return;

   computeOnMarked_exec(g, nodes_isMarked, startnode, mst);
 //  computeOnMarked_computePredecessors(g, mst);
}
