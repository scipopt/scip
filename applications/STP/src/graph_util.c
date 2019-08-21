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

/**@file   graph_util.c
 * @brief  includes graph utility methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * Graph utility methods for Steiner problems, such as CSR data structure
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#define DHEAP_MAX_KEY  (1e20)
#define DHEAP_MIN_KEY  -(1e20)

#include "graph.h"
#include "portab.h"

/** clean the heap */
void
graph_heap_clean(
   SCIP_Bool             cleanposition,      /**< clean position array? */
   DHEAP*                heap                /**< the heap  */
   )
{
   int* const position = heap->position;
   int capacity = heap->capacity;

   assert(heap && position);

   heap->size = 0;

   if( cleanposition )
      for( int i = 0; i < capacity; i++ )
         position[i] = UNKNOWN;
}


/** creates new heap. If entries array is provided, it must be of size capacity + 2  */
SCIP_RETCODE graph_heap_create(
   SCIP*                 scip,               /**< SCIP */
   int                   capacity,           /**< heap capacity */
   int*                  position,           /**< heap position array or NULL */
   DENTRY*               entries,            /**< entries array or NULL */
   DHEAP**               heap                /**< the heap  */
   )
{
   int* position_heap;
   DENTRY* entries_heap;

   assert(scip && heap);
   assert(capacity >= 1);

   SCIP_CALL( SCIPallocMemory(scip, heap) );

   if( position )
      position_heap = position;
   else
      SCIP_CALL( SCIPallocMemoryArray(scip, &(position_heap), capacity) );

   if( entries )
      entries_heap = entries;
   else
      SCIP_CALL( SCIPallocMemoryArray(scip, &(entries_heap), capacity + 2) );

   (*heap)->capacity = capacity;
   (*heap)->position = position_heap;
   (*heap)->entries = entries_heap;

   /* sentinel */
   entries_heap[0].key = DHEAP_MIN_KEY;

   /* debug sentinel */
   entries_heap[capacity + 1].key = DHEAP_MAX_KEY;

   graph_heap_clean(TRUE, *heap);

   return SCIP_OKAY;
}

/** frees the heap */
void graph_heap_free(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Bool             freePositions,      /**< free positions array? */
   SCIP_Bool             freeEntries,        /**< free entries array? */
   DHEAP**               heap                /**< the heap  */
   )
{
   assert(scip && heap);

   if( freePositions )
      SCIPfreeMemoryArray(scip, &((*heap)->position));

   if( freeEntries )
      SCIPfreeMemoryArray(scip, &((*heap)->entries));

   SCIPfreeMemoryArray(scip, heap);
}

/** deletes heap minimum */
void graph_heap_deleteMin(
   int*                  node,               /**< pointer to value of minimum */
   SCIP_Real*            key,                /**< pointer to key of minimum */
   DHEAP*                heap                /**< the heap  */
   )
{
   *key = heap->entries[1].key;

   graph_heap_deleteMinGetNode(node, heap);
}

/** deletes heap minimum */
void graph_heap_deleteMinGetNode(
   int*                  node,               /**< pointer to node stored in minimum (set by method) */
   DHEAP*                heap                /**< the heap  */
   )
{
   assert(node);
   *node = graph_heap_deleteMinReturnNode(heap);
}


/** deletes heap minimum and returns corresponding node */
int graph_heap_deleteMinReturnNode(
   DHEAP*                heap                /**< the heap  */
   )
{
   int* const RESTRICT position = heap->position;
   DENTRY* const RESTRICT entries = heap->entries;
   SCIP_Real fill;
   int parent;
   int hole = 1;
   int child = 2;
   int node;
   const int lastentry = heap->size--;

   assert(heap && position && entries);
   assert(heap->size >= 0);

   node = entries[1].node;

   assert(position[node] == 1);

   position[node] = CONNECT;

   /* move down along min-path */
   while( child < lastentry )
   {
      const SCIP_Real key1 = entries[child].key;
      const SCIP_Real key2 = entries[child + 1].key;
      assert(hole >= 1);
      assert(key1 < DHEAP_MAX_KEY && key2 < DHEAP_MAX_KEY);

      /* second child with smaller key? */
      if( key1 > key2 )
      {
         entries[hole].key = key2;
         child++;
      }
      else
      {
         entries[hole].key = key1;
      }

      assert(entries[hole].node >= 0 && entries[hole].node < heap->capacity);

      entries[hole].node = entries[child].node;
      position[entries[hole].node] = hole;

      hole = child;
      child *= 2;
   }

   /* now hole is at last tree level, fill it with last heap entry and move it up */

   fill = entries[lastentry].key;
   parent = hole / 2;

   assert(fill < DHEAP_MAX_KEY && entries[parent].key < DHEAP_MAX_KEY);

   while( entries[parent].key > fill )
   {
      assert(hole >= 1);

      entries[hole] = entries[parent];

      assert(entries[hole].node >= 0 && entries[hole].node < heap->capacity);

      position[entries[hole].node] = hole;
      hole = parent;
      parent /= 2;

      assert(entries[parent].key < DHEAP_MAX_KEY);
   }

   /* finally, fill the hole */
   entries[hole].key = fill;
   entries[hole].node = entries[lastentry].node;

   assert(entries[hole].node >= 0 && entries[hole].node < heap->capacity);

   if( hole != lastentry )
      position[entries[hole].node] = hole;

#ifndef NDEBUG
   entries[lastentry].key = DHEAP_MAX_KEY;    /* set debug sentinel */
#endif

   return node;
}


/** corrects node position in heap according to new key (or newly inserts the node) */
void graph_heap_correct(
   int                   node,               /**< the node */
   SCIP_Real             newkey,             /**< the new key (needs to be smaller than current one) */
   DHEAP*                heap                /**< the heap  */
   )
{
   int* const RESTRICT position = heap->position;
   DENTRY* const RESTRICT entries = heap->entries;
   int hole;
   int parent;
   SCIP_Real parentkey;

   assert(heap && position && entries);
   assert(newkey < DHEAP_MAX_KEY && newkey > DHEAP_MIN_KEY);
   assert(heap->size <= heap->capacity);
   assert(node >= 0 && node <= heap->capacity);
   assert(position[node] != CONNECT);

   /* node not yet in heap? */
   if( position[node] == UNKNOWN )
   {
      assert(heap->size < heap->capacity);
      hole = ++(heap->size);
   }
   else
   {
      assert(position[node] >= 1);
      hole = position[node];

      assert(entries[hole].node == node);
      assert(GE(entries[hole].key, newkey));
   }

   parent = hole / 2;
   parentkey = entries[parent].key;

   assert(parentkey < DHEAP_MAX_KEY);

   /* move hole up */
   while( parentkey > newkey )
   {
      assert(hole >= 1);

      entries[hole].key = parentkey;
      entries[hole].node = entries[parent].node;
      position[entries[hole].node] = hole;
      hole = parent;
      parent /= 2;
      parentkey = entries[parent].key;
      assert(parentkey < DHEAP_MAX_KEY);
   }

   /* fill the hole */
   entries[hole].key = newkey;
   entries[hole].node = node;
   position[node] = hole;
}


/** initializes CSR storage */
SCIP_RETCODE graph_init_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   CSR* csr;
   int* start_csr;
   int* head_csr;
   SCIP_Real* cost_csr;
   const int nedges = g->edges;
   const int nnodes = g->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);

   assert(scip && g);
   assert(nnodes >= 1);
   assert(g->csr_storage == NULL);
   assert(!pcmw || !g->extended);

   SCIP_CALL( SCIPallocMemory(scip, &csr) );
   g->csr_storage = csr;

   csr->nedges = nedges;
   csr->nnodes = nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(start_csr), nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(head_csr), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cost_csr), nedges) );

   csr->start = start_csr;
   csr->head = head_csr;
   csr->cost = cost_csr;

   /* now fill the data in */

   start_csr[0] = 0;

   for( int k = 0; k < nnodes; k++ )
   {
      int pos = start_csr[k];

      if( !pcmw || g->mark[k] )
      {
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int ehead = g->head[e];

            if( pcmw && !g->mark[ehead] )
               continue;

            assert(g->cost[e] < FARAWAY && g->cost[flipedge(e)] < FARAWAY);

            head_csr[pos] = ehead;
            cost_csr[pos++] = g->cost[e];
         }
      }

      assert((pos == start_csr[k] + g->grad[k]) || pcmw);

      start_csr[k + 1] = pos;
   }

   assert(start_csr[nnodes] <= nedges);
   assert(graph_valid_csr(g, TRUE));

   return SCIP_OKAY;
}


/** frees dynamic CSR storage */
void graph_free_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   CSR* csr = g->csr_storage;

   assert(scip && g);
   assert(csr != NULL && csr->nnodes >= 1);

   SCIPfreeMemoryArray(scip, &(csr->cost));
   SCIPfreeMemoryArray(scip, &(csr->head));
   SCIPfreeMemoryArray(scip, &(csr->start));

   SCIPfreeMemoryArray(scip, &(g->csr_storage));
   assert(g->csr_storage == NULL);
}

/** is CSR storage of graph valid? */
SCIP_Bool graph_valid_csr(
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool             verbose             /**< be verbose? */
)
{
   const CSR* csr = g->csr_storage;
   const int* start = csr->start;
   const int nedges = g->edges;
   const int nnodes = g->knots;

   assert(g && csr && start && csr->head && csr->cost);

   if( nnodes != g->knots || nedges != g->edges )
   {
      if( verbose )
         printf("CSR: wrong node/edge cound \n");
      return FALSE;
   }

   if( start[0] != 0 )
      return FALSE;

   for( int i = 0; i < nnodes; i++ )
   {
      if( start[i] > start[i + 1] )
      {
         if( verbose )
            printf("CSR: ranges corrupted \n");

         return FALSE;
      }
   }

   return TRUE;
}


/** initializes dynamic CSR storage */
SCIP_RETCODE graph_init_dcsr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   int* edgeid_csr;
   int* id2csr_csr;
   SCIP_Real* cost_csr;
   const int nedges = g->edges;
   const int nnodes = g->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);

   assert(scip && g);
   assert(nnodes >= 1);
   assert(g->dcsr_storage == NULL);
   assert(!pcmw || !g->extended);

   SCIP_CALL( SCIPallocMemory(scip, &dcsr) );
   g->dcsr_storage = dcsr;

   dcsr->nedges = nedges;
   dcsr->nnodes = nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(range_csr), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(head_csr), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(edgeid_csr), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(id2csr_csr), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cost_csr), nedges) );

   dcsr->range = range_csr;
   dcsr->head = head_csr;
   dcsr->edgeid = edgeid_csr;
   dcsr->id2csredge = id2csr_csr;
   dcsr->cost = cost_csr;
   dcsr->cost2 = NULL;
   dcsr->cost3 = NULL;

   /* now fill the data in */

   for( int e = 0; e < nedges; e++ )
      id2csr_csr[e] = -1;

   range_csr[0].start = 0;

   for( int k = 0; k < nnodes; k++ )
   {
      int pos = range_csr[k].start;

      if( !pcmw || g->mark[k] )
      {
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int ehead = g->head[e];

            if( pcmw && !g->mark[ehead] )
               continue;

            assert(g->cost[e] < FARAWAY && g->cost[flipedge(e)] < FARAWAY);

            id2csr_csr[e] = pos;
            head_csr[pos] = ehead;
            edgeid_csr[pos] = e;
            cost_csr[pos++] = g->cost[e];
         }
      }

      assert(pos == range_csr[k].start + g->grad[k] || pcmw);

      range_csr[k].end = pos;

      if( k != nnodes - 1 )
         range_csr[k + 1].start = range_csr[k].end;
   }

   assert(range_csr[nnodes - 1].end <= nedges);
   assert(graph_valid_dcsr(g, TRUE));

   return SCIP_OKAY;
}

/** frees dynamic CSR storage */
void graph_free_dcsr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   DCSR* dcsr = g->dcsr_storage;

   assert(scip && g);
   assert(dcsr != NULL && dcsr->nnodes >= 1);

   SCIPfreeMemoryArray(scip, &(dcsr->cost));
   SCIPfreeMemoryArray(scip, &(dcsr->id2csredge));
   SCIPfreeMemoryArray(scip, &(dcsr->edgeid));
   SCIPfreeMemoryArray(scip, &(dcsr->head));
   SCIPfreeMemoryArray(scip, &(dcsr->range));

   SCIPfreeMemoryArray(scip, &(g->dcsr_storage));

   assert(g->dcsr_storage == NULL);
}

/** updates dynamic CSR storage */
void graph_update_dcsr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   assert(scip && g);

   assert(0 && "implement"); // check whether enough edges and nodes, otherwise reallocate
}

/** deletes CSR indexed edge */
void graph_dcsr_deleteEdge(
   DCSR*                 dcsr,               /**< DCSR container */
   int                   tail,               /**< tail of edge */
   int                   e_csr               /**< CSR indexed edge */
)
{
   RANGE* const range = dcsr->range;
   int* const head = dcsr->head;
   int* const edgeid = dcsr->edgeid;
   int* const id2csredge = dcsr->id2csredge;
   SCIP_Real* const cost = dcsr->cost;
   SCIP_Real* const cost2 = dcsr->cost2;
   SCIP_Real* const cost3 = dcsr->cost3;
   int last;

   assert(dcsr);
   assert(tail >= 0 && tail < dcsr->nnodes);
   assert(e_csr >= 0 && e_csr < dcsr->nedges);
   assert(range[tail].start <= e_csr && e_csr < range[tail].end);

   last = --(range[tail].end);

#ifndef NDEBUG
   id2csredge[edgeid[e_csr]] = -1;
#endif

   /* e_csr not already deleted? */
   if( e_csr != last )
   {
      head[e_csr] = head[last];
      edgeid[e_csr] = edgeid[last];
      id2csredge[edgeid[last]] = e_csr;
      cost[e_csr] = cost[last];

      if( cost2 )
         cost2[e_csr] = cost2[last];

      if( cost3 )
         cost3[e_csr] = cost3[last];
   }

#ifndef NDEBUG
   head[last] = -1;
   edgeid[last] = -1;
   cost[last] = -FARAWAY;
   if( cost2 )
      cost2[last] = -FARAWAY;
   if( cost3 )
      cost3[last] = -FARAWAY;
#endif

}

/** deletes CSR indexed edge and anti-parallel one */
void graph_dcsr_deleteEdgeBi(
   SCIP*                 scip,               /**< SCIP data structure */
   DCSR*                 dcsr,               /**< DCSR container */
   int                   e_csr               /**< CSR indexed edge */
)
{
   int* const head = dcsr->head;
   int* const edgeid = dcsr->edgeid;
   int* const id2csredge = dcsr->id2csredge;
   const int erev_csr = id2csredge[flipedge(edgeid[e_csr])];
   const int i1 = head[erev_csr];
   const int i2 = head[e_csr];

   assert(scip && dcsr);
   assert(e_csr >= 0 && edgeid[e_csr] >= 0);
   assert(erev_csr >= 0);
   assert(SCIPisEQ(scip, dcsr->cost[e_csr], dcsr->cost[erev_csr]));

   SCIPdebugMessage("delete %d %d \n", i1, i2);

   graph_dcsr_deleteEdge(dcsr, i2, erev_csr);
   graph_dcsr_deleteEdge(dcsr, i1, e_csr);
}

/** is DCSR storage of graph valid? */
SCIP_Bool graph_valid_dcsr(
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool             verbose             /**< be verbose? */
)
{
   const DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const range = dcsr->range;
   const int* const head = dcsr->head;
   const int* const edgeid = dcsr->edgeid;
   const int* const id2csredge = dcsr->id2csredge;
   const int nnodes = dcsr->nnodes;
   const int nedges = dcsr->nedges;

   assert(g && dcsr && range && head && edgeid && id2csredge);

   if( nnodes != g->knots || nedges != g->edges )
   {
      if( verbose )
         printf("DCSR: wrong node/edge cound \n");
      return FALSE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      const int start = range[i].start;
      const int end = range[i].end;

      if( start > end )
      {
         if( verbose )
            printf("DCSR: ranges corrupted \n");

         return FALSE;
      }

      for( int e = start; e < end; e++ )
      {
         const int ordedge = edgeid[e];

         assert(ordedge >= 0 && ordedge < nedges);

         if( id2csredge[ordedge] != e )
         {
            if( verbose )
               printf("DCSR: id assignment corrupted \n");

            return FALSE;
         }

         if( head[e] != g->head[ordedge] || i != g->tail[ordedge] )
         {
            if( verbose )
               printf("DCSR: edge assignment corrupted \n");

            printf(" %d == %d %d == %d \n",  head[e], g->head[ordedge], i, g->tail[ordedge]);

            return FALSE;
         }
      }
   }

   return TRUE;
}


/** initializes (allocates and fills) limited Dijkstra structure members */
SCIP_RETCODE graph_dijkLimited_init(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   DIJK*                 dijkdata            /**< data for limited Dijkstra */
)
{
   const int nnodes = g->knots;
   SCIP_Real* distance;
   STP_Bool* visited;

   assert(scip && g && dijkdata);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(dijkdata->distance), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dijkdata->visitlist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dijkdata->visited), nnodes) );

   graph_heap_create(scip, nnodes, NULL, NULL, &(dijkdata->dheap));

   dijkdata->nvisits = -1;

   distance = dijkdata->distance;
   visited = dijkdata->visited;

   for( int k = 0; k < nnodes; k++ )
   {
      visited[k] = FALSE;
      distance[k] = FARAWAY;
   }

   return SCIP_OKAY;
}

/** cleans limited Dijkstra structure members */
void graph_dijkLimited_clean(
   const GRAPH*          g,                  /**< the graph */
   DIJK*                 dijkdata            /**< data for limited Dijkstra */
)
{
   const int nnodes = g->knots;
   STP_Bool* const visited = dijkdata->visited;
   SCIP_Real* const distance = dijkdata->distance;
   dijkdata->nvisits = -1;

   for( int k = 0; k < nnodes; k++ )
   {
      visited[k] = FALSE;
      distance[k] = FARAWAY;
   }

   graph_heap_clean(TRUE, dijkdata->dheap);
}


/** reset data of limited Dijkstra structure */
void graph_dijkLimited_reset(
   const GRAPH*          g,                  /**< the graph */
   DIJK*                 dijkdata            /**< data for limited Dijkstra */
)
{
   STP_Bool* const visited = dijkdata->visited;
   int* const visitlist = dijkdata->visitlist;
   int* const state = dijkdata->dheap->position;
   SCIP_Real* const distance = dijkdata->distance;
   const int nvisits = dijkdata->nvisits;

   assert(dijkdata && g);

   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      assert(node >= 0 && node < g->knots);

      visited[node] = FALSE;
      distance[node] = FARAWAY;
      state[node] = UNKNOWN;
   }

   graph_heap_clean(FALSE, dijkdata->dheap);

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(visited[k] == FALSE);
      assert(state[k] == UNKNOWN);
      assert(distance[k] == FARAWAY);
   }
#endif
}


/** frees limited Dijkstra structure member */
void graph_dijkLimited_freeMembers(
   SCIP*                 scip,               /**< SCIP */
   DIJK*                 dijkdata            /**< data for limited Dijkstra */
)
{
   SCIPfreeMemoryArray(scip, &(dijkdata->distance));
   SCIPfreeMemoryArray(scip, &(dijkdata->visitlist));
   SCIPfreeMemoryArray(scip, &(dijkdata->visited));

   graph_heap_free(scip, TRUE, TRUE, &(dijkdata->dheap));
}
