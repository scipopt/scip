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

/**@file   graph_tpath.c
 * @brief  Graph algorithms for Steiner problems related to shortest paths to terminals.
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file encompasses various algorithms that are used to compute shortest paths from
 * non-terminals to terminals.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "graph.h"
#include "graphdefs.h"
#include "stpvector.h"
#include "reduce.h"

#define STP_TPATHS_NTERMBASES 4
#define STP_TPATHS_RESERVESIZE 16

/** reset single path struct */
#define pathResetSingle(entry) \
   do                          \
   {                           \
      entry.dist = FARAWAY;    \
      entry.edge = UNKNOWN;    \
   } while( 0 )


/**  resets node */
#define tpathsRepairResetNode(scip, node, resetnodes, stack, nodes_isvisited)        \
   do                                                                                 \
   {                                                                 \
      assert(scip && resetnodes && stack && nodes_isvisited); \
      assert(!nodes_isvisited[node]);                              \
      StpVecPushBack(scip, resetnodes, node);                        \
      StpVecPushBack(scip, stack, node);                             \
      nodes_isvisited[node] = TRUE;                                  \
   } while( 0 )


/** Steiner nodes to terminal paths
 * NOTE: all arrays are of size STP_TPATHS_NTERMBASES * nnodes */
struct nodes_to_terminal_paths
{
   PATH*                 termpaths;         /**< path data (leading to first, second, ... terminal) */
   int*                  termbases;         /**< terminals to each non terminal */
   int*                  termbases_org;     /**< original terminals to each non terminal; only needed for repair */
   int*                  state;             /**< array to mark the state of each node during calculation */
   int                   nnodes;            /**< number of nodes of underlying graph */
};


/** data needed to repair node to terminal paths for edge-deletion */
typedef struct tpaths_repair
{
   STP_Vectype(int)*     resetnodes;         /**< nodes to be reseted; for different levels (1st to 4th) */
   STP_Vectype(int)      stack;              /**< temporary vector */
   TPATHS*               tpaths;             /**< the terminal paths */
   SCIP_Bool*            nodes_isvisited;    /**< visited nodes during repair process */
   int                   edge;               /**< edge about to be eliminated */
   int                   nHeapElems;         /**< number of heap elements */
   SCIP_Bool             withEdgeReplacement;/**< with edge replacement? */
} TREPAIR;




/*
 * Local methods
 */



/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static
int pathheapGetNearest(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   const PATH* path)
{
   int   k;
   int   t;
   int   c;
   int   j;

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   if ((*count) > 2)
      if (LT(path[heap[3]].dist, path[heap[2]].dist))
         c++;

   while((c <= (*count)) && GT(path[heap[j]].dist, path[heap[c]].dist))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
         if (LT(path[heap[c + 1]].dist, path[heap[c]].dist))
            c++;
   }
   return(k);
}



/** set new distance and predeccesor */
inline static
void pathheapCorrectDist(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   PATH* RESTRICT path,
   int    l,
   int    k,
   int    edge,
   SCIP_Real newdist
   )
{
   int    t;
   int    c;
   int    j;

   path[l].dist = newdist;
   path[l].edge = edge;

   /* new node? */
   if( state[l] == UNKNOWN )
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up */
   j = state[l];
   c = j / 2;
   while( (j > 1) && path[heap[c]].dist > path[heap[j]].dist )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


/** terminal to terminal distance */
inline
static void utdist(
   SCIP*            scip,
   const GRAPH*  g,
   PATH* path,
   SCIP_Real ecost,
   int* vbase,
   int k,
   int l,
   int k2,
   int shift,
   int nnodes
   )
{
   SCIP_Real dist;
   int vbk;
   int vbk2;

   if( Is_term(g->term[k]) )
      vbk = k;
   else
      vbk = vbase[k];

   if( l == 0 )
   {
      assert(shift == 0);

      dist = ecost;
      if( !Is_term(g->term[k]) )
         dist += path[k].dist;

      if( !Is_term(g->term[k2]) )
      {
         dist += path[k2].dist;
         vbk2 = vbase[k2];
      }
      else
      {
         vbk2 = k2;
      }

      if( SCIPisLT(scip, dist, path[vbk].dist) )
      {
         path[vbk].dist = dist;
         vbase[vbk] = vbk2;
         return;
      }
   }
   else
   {
      int max;
      int pos;
      int r;
      int s;
      int t;

      pos = vbk + shift;
      max = MIN((l + 1), 3);

      for( r = 0; r <= max; r++ )
      {
         if( Is_term(g->term[k2]) )
         {
            if( r == 0 )
               t = k2;
            else
               break;
         }
         else
         {
            t = vbase[k2 + (r * nnodes)];
         }
         for( s = 0; s < l; s++ )
            if( vbase[vbk + s * nnodes] == t )
               break;
         if( s < l || vbk == t )
            continue;

         dist = ecost;
         if( !Is_term(g->term[k]) )
            dist += path[k].dist;
         if( !Is_term(g->term[k2]) )
            dist += path[k2 + (r * nnodes)].dist;

         if( SCIPisLT(scip, dist, path[pos].dist) )
         {
            path[pos].dist = dist;
            vbase[pos] = t;
            return;
         }
      }
   }
   return;
}


/** prints terminal path from given node */
static
void tpathsPrintPath(
   const GRAPH*          g,                  /**< graph data structure */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   const TPATHS*         tpaths,             /**< storage for terminal paths */
   int                   node,               /**< node to start from */
   int                   level               /**< between 1 and 4 */
)
{
   const int* const termbases = tpaths->termbases;
   const int nnodes = graph_get_nNodes(g);
   const int offset = (level - 1) * nnodes;
   const int base = termbases[node + offset];
   const PATH* const termpaths = tpaths->termpaths;
   int node_pred;
   int node_curr;
   int node_next;
   int edge_curr;
   assert(1 <= level && level <= STP_TPATHS_NTERMBASES);
   assert(graph_knot_isInRange(g, base));
   assert(base != UNKNOWN);

   node_pred = -1;
   node_curr = node;
   edge_curr = termpaths[node_curr + offset].edge;
   node_next = (edge_curr >= 0) ? g->tail[edge_curr] : -1;

   while( node_curr != base )
   {
      assert(graph_knot_isInRange(g, node_curr));

      if( !graph_edge_isInRange(g, edge_curr) )
      {
         printf("found deleted edge, break \n");
         break;
      }

      if( sdprofit && node_curr != node )
      {
         const SCIP_Real profit = reduce_sdprofitGetProfit(sdprofit, node_curr, node_pred, node_next);
         printf("node %d profit: %f \n", node_curr, profit);
      }

      graph_edge_printInfo(g, edge_curr);

      if( graph_edge_isDeleted(g, edge_curr) )
      {
         printf("...edge killed! \n");
      }

      node_pred = node_curr;
      node_curr = node_next;
      edge_curr = termpaths[node_curr + offset].edge;
      node_next = (edge_curr >= 0) ? g->tail[edge_curr] : -1;
   }
}


/** allocates TPATHS data */
static
SCIP_RETCODE tpathsAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TPATHS**              tpaths              /**< the terminal paths */
)
{

   TPATHS* tp;
   const int nnodes = graph_get_nNodes(g);
   assert(nnodes >= 1);

   SCIP_CALL( SCIPallocMemory(scip, tpaths) );

   tp = *tpaths;
   tp->nnodes = nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->termpaths), nnodes * STP_TPATHS_NTERMBASES) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->termbases), nnodes * STP_TPATHS_NTERMBASES) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tp->state), nnodes * STP_TPATHS_NTERMBASES) );
   tp->termbases_org = NULL;

   return SCIP_OKAY;
}


/** resets TPATHS members */
static
void tpathsResetMembers(
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level between 2 and 4 */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;
   const int nnodes = graph_get_nNodes(g);
   const int offset = (level - 1) * nnodes;

   assert(2 <= level && level <= 4);
   assert(offset == nnodes || offset == 2 * nnodes || offset == 3 * nnodes);

   for( int i = 0; i < nnodes; i++ )
   {
      /* copy of node i */
      const int k = i + offset;

      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      pathResetSingle(path[k]);
   }

   for( int i = 0; i < offset; i++ )
   {
      state[i] = CONNECT;
   }
}


/** gets (up to) four close terminals to given node i */
static inline
void tpathsGetKCloseTerms(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   maxnumber,          /**< number of close terms */
   int                   node,               /**< node */
   SCIP_Real             maxdist,            /**< maximum valid distance  */
   SCIP_Bool             distIsStrict,       /**< is 'maxdist' strict? */
   int* RESTRICT         closeterms,         /**< four close terminals */
   int* RESTRICT         firstedges,         /**< corresponding first edge (can be NULL) */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< four close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const PATH* const termpaths = tpaths->termpaths;
   const int* const termbases = tpaths->termbases;
   const int nnodes = tpaths->nnodes;
   int pos = node;
   int nnterms = 0;
   const SCIP_Bool withFirstEdges = (firstedges != NULL);

   assert(closeterms && closeterms_dist && ncloseterms && g);
   assert(maxnumber <= STP_TPATHS_NTERMBASES);
   assert(0 <= node && node < nnodes);
   assert(nnodes == g->knots);
   assert(GE(maxdist, 0.0));

   if( Is_term(g->term[node]))
   {
      *ncloseterms = 1;
      closeterms[0] = node;
      closeterms_dist[0] = 0.0;

      if( withFirstEdges )
         firstedges[0] = UNKNOWN;

      return;
   }

   for( int k = 0; k < maxnumber; k++ )
   {
      const SCIP_Bool addTerm =
         distIsStrict ?
            LT(termpaths[pos].dist, maxdist)
         :  LE(termpaths[pos].dist, maxdist);

      if( addTerm )
      {
         closeterms[nnterms] = termbases[pos];

         if( withFirstEdges )
         {
            firstedges[nnterms] = termpaths[pos].edge;
         }
         closeterms_dist[nnterms++] = termpaths[pos].dist;
      }
      else
      {
         assert(distIsStrict || !LE(termpaths[pos].dist, maxdist));
      }

      pos += nnodes;
   }


#ifndef NDEBUG
   for( int k = 0; k < nnterms; k++ )
   {
      const int term = closeterms[k];
      assert(graph_knot_isInRange(g, term));
      assert(Is_term(g->term[term]) || g->grad[term] == 0);
   }
#endif

   *ncloseterms = nnterms;
}


/** gets new distance for extension */
inline static
SCIP_Real tpathsGetDistNew(
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const int*            vbase,              /**< bases */
   int                   node,               /**< node to extend from */
   int                   nextnode,           /**< node to extend to */
   int                   edge,               /**< the edge */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   const PATH*           path                /**< the actual terminal paths */
)
{
   SCIP_Real distnew;
   const SCIP_Real ecost = ((g->source == vbase[node]) ? cost[edge] : costrev[edge]);

#ifndef NDEBUG
   assert(node % g->knots == g->tail[edge]);
   assert(nextnode % g->knots == g->head[edge]);
#endif

   if( sdprofit && path[node].edge >= 0 )
   {
      const int nnodes = graph_get_nNodes(g);
      const int node_pred = g->tail[path[node].edge];

      assert(!Is_term(g->term[node % g->knots]));

      if( node_pred == nextnode % nnodes )
      {
         distnew = path[node].dist + ecost;
      }
      else
      {
         distnew = reduce_sdprofitGetBiasedDist(sdprofit, node % nnodes,
               ecost, path[node].dist, nextnode % nnodes, node_pred);
      }
   }
   else
   {
      distnew = path[node].dist + ecost;
   }

   assert(graph_pc_isMw(g) || GE(distnew, 0.0));
   assert(LE(distnew, path[node].dist + ecost));

   return distnew;
}


/** allocates TPATHS data */
static inline
void tpathsBuild(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   assert(STP_TPATHS_NTERMBASES == 4);

   graph_tpathsSetAll4(g, g->cost, g->cost, NULL, tpaths);
}

/** allocates TPATHS data */
static inline
void tpathsBuildBiased(
   const SDPROFIT*       sdprofit,           /**< SD profit */
   GRAPH*                g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   assert(STP_TPATHS_NTERMBASES == 4);

   graph_tpathsSetAll4(g, g->cost, g->cost, sdprofit, tpaths);
}


/** sub-method to find closest terminal to each non terminal (Voronoi diagram) */
static
void tpathsScan1st(
   SCIP*                 scip,               /**< SCIP or NULL */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   int                   heapsize,           /**< size of heap */
   TPATHS*               tpaths,             /**< storage for terminal paths */
   TREPAIR*              repair              /**< data for repairing or NULL */
)
{
   int nHeapElems = heapsize;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;
   const SCIP_Bool withProfit = (sdprofit != NULL);
   const SCIP_Bool extend = (repair != NULL);
   SCIP_Bool* nodes_isvisited = extend ? repair->nodes_isvisited : NULL;

   assert(nHeapElems >= 0);
   assert(repair == NULL || repair->withEdgeReplacement);

   while( nHeapElems > 0 )
   {
      const int k = pathheapGetNearest(heap, state, &nHeapElems, path);
      const SCIP_Real k_dist = path[k].dist;



      /* mark vertex k as scanned */
      state[k] = CONNECT;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( state[m] != CONNECT && g->mark[m]  )
         {
            SCIP_Real distnew;

            if( withProfit && vbase[k] != k )
            {
               int k_pred;
               assert(path[k].edge >= 0);
               k_pred = g->tail[path[k].edge];
               distnew = reduce_sdprofitGetBiasedDist(sdprofit, k, cost[e], k_dist, m, k_pred);
            }
            else
            {
               assert(!withProfit || Is_term(g->term[k]));
               distnew = k_dist + cost[e];
            }

            assert(GE(distnew, 0.0) && LE(distnew, k_dist + cost[e]));

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( GT(path[m].dist, distnew) )
            {
               if( extend && !nodes_isvisited[m] )
               {
                  StpVecPushBack(scip, repair->resetnodes[0], m);
                  nodes_isvisited[m] = TRUE;
               }

               pathheapCorrectDist(heap, state, &nHeapElems, path, m, k, e, distnew);
               vbase[m] = vbase[k];
            }
         }
      }
   }
}


/** sub-method to find 2nd closest terminal to each non terminal */
static
void tpathsScan2nd(
   SCIP*                 scip,               /**< SCIP or NULL */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reverse edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   int                   heapsize,           /**< size of heap */
   TPATHS*               tpaths,             /**< storage for terminal paths */
   TREPAIR*              repair              /**< data for repairing or NULL */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nheapElems = heapsize;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;
   const SCIP_Bool extend = (repair != NULL);
   SCIP_Bool* nodes_isvisited = extend ? repair->nodes_isvisited : NULL;

   assert(nheapElems >= 0);
   assert(repair == NULL || repair->withEdgeReplacement);

   while( nheapElems > 0 )
   {
      const int k2 = pathheapGetNearest(heap, state, &nheapElems, path);
      const int k = k2 - nnodes;

      assert(graph_knot_isInRange(g, k));

      /* mark vertex k as removed from heap */
      state[k2] = UNKNOWN;

      /* iterate over all outgoing edges of vertex (ancestor of) k */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( !Is_term(g->term[head]) && vbase[head] != vbase[k2] && g->mark[head] )
         {
            const int head2 = head + nnodes;
            const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, k2, head2, e, sdprofit, path);

            /* check whether the path (to j) including k is shorter than the so far best known */
            if( GT(path[head2].dist, distnew) )
            {
               if( extend && !nodes_isvisited[head] )
               {
                  StpVecPushBack(scip, repair->resetnodes[1], head);
                  nodes_isvisited[head] = TRUE;
               }

               pathheapCorrectDist(heap, state, &nheapElems, path, head2, k2, e, distnew);
               vbase[head2] = vbase[k2];
            }
         }
      }
   }
}


/** sub-method to find 3rd closest terminal to each non terminal */
static
void tpathsScan3rd(
   SCIP*                 scip,               /**< SCIP or NULL */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reverse edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   int                   heapsize,           /**< size of heap */
   TPATHS*               tpaths,             /**< storage for terminal paths */
   TREPAIR*              repair              /**< data for repairing or NULL */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int dnnodes = 2 * nnodes;
   int nheapElems = heapsize;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;
   const SCIP_Bool extend = (repair != NULL);
   SCIP_Bool* nodes_isvisited = extend ? repair->nodes_isvisited : NULL;

   assert(nheapElems >= 0);
   assert(repair == NULL || repair->withEdgeReplacement);

   while( nheapElems > 0 )
   {
      /* get the next (i.e. a nearest) vertex of the heap */
      const int k3 = pathheapGetNearest(heap, state, &nheapElems, path);
      const int k = k3 - dnnodes;
      assert(0 <= k && k < nnodes);

      /* mark vertex k as removed from heap */
      state[k3] = UNKNOWN;

      /* iterate over all outgoing edges of vertex k */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( !Is_term(g->term[head]) && vbase[head] != vbase[k3] && vbase[head + nnodes] != vbase[k3] && g->mark[head] )
         {
            const int head3 = head + dnnodes;
            const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, k3, head3, e, sdprofit, path);

            /* check whether the path (to j) including k is shorter than the so far best known */
            if( GT(path[head3].dist, distnew) )
            {
               if( extend && !nodes_isvisited[head] )
               {
                  StpVecPushBack(scip, repair->resetnodes[2], head);
                  nodes_isvisited[head] = TRUE;
               }

               pathheapCorrectDist(heap, state, &nheapElems, path, head3, k3, e, distnew);
               vbase[head3] = vbase[k3];
            }
         }
      }
   }
}


/** sub-method to find 4th closest terminal to each non terminal */
static
void tpathsScan4th(
   SCIP*                 scip,               /**< SCIP or NULL */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reverse edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   int                   heapsize,           /**< size of heap */
   TPATHS*               tpaths,             /**< storage for terminal paths */
   TREPAIR*              repair              /**< data for repairing or NULL */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int dnnodes = 2 * nnodes;
   const int tnnodes = 3 * nnodes;
   int nHeapElems = heapsize;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;
   const SCIP_Bool extend = (repair != NULL);
   SCIP_Bool* nodes_isvisited = extend ? repair->nodes_isvisited : NULL;

   assert(nHeapElems >= 0);
   assert(repair == NULL || repair->withEdgeReplacement);

   /* until the heap is empty */
   while( nHeapElems > 0 )
   {
      /* get the next (i.e. a nearest) vertex of the heap */
      const int k4 = pathheapGetNearest(heap, state, &nHeapElems, path);
      const int k = k4 - tnnodes;
      assert(0 <= k && k < nnodes);

      /* mark vertex k as removed from heap */
      state[k4] = UNKNOWN;

      /* iterate over all outgoing edges of vertex k */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( !Is_term(g->term[head])
            && vbase[head] != vbase[k4] && vbase[head + nnodes] != vbase[k4] && vbase[head + dnnodes] != vbase[k4] && g->mark[head] )
         {
            const int head4 = head + tnnodes;
            const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, k4, head4, e, sdprofit, path);

            /* check whether the path4 (to j) including k is shorter than the so far best known */
            if( GT(path[head4].dist, distnew) )
            {
               if( extend && !nodes_isvisited[head] )
               {
                  StpVecPushBack(scip, repair->resetnodes[3], head);
                  nodes_isvisited[head] = TRUE;
               }

               pathheapCorrectDist(heap, state, &nHeapElems, path, head4, k4, e, distnew);
               vbase[head4] = vbase[k4];
            }
         }
      }
   }
}


/** initializes */
static inline
SCIP_RETCODE tpathsRepairInitLevel(
   SCIP*                 scip,               /**< SCIP */
   int                   level,              /**< 0 - 3 */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
   )
{
   const int nnodes = graph_get_nNodes(g);

#ifndef NDEBUG
   const int* const state = repair->tpaths->state;
   for( int i = 0; i < nnodes; ++i )
   {
      assert(state[i] == UNKNOWN);
   }
#endif
   assert(0 <= level && level <= 3);
   assert(!repair->nodes_isvisited);
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(repair->nodes_isvisited), nnodes) );

   repair->nHeapElems = -1;

   return SCIP_OKAY;
}


/** cleans and frees */
static inline
void tpathsRepairExitLevel(
   SCIP*                 scip,               /**< SCIP */
   int                   level,              /**< 0 - 3 */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
   )
{
   const STP_Vectype(int) resetnodes_level = repair->resetnodes[level];
   int* RESTRICT state = repair->tpaths->state;
   SCIP_Bool* RESTRICT nodes_isvisited = repair->nodes_isvisited;
   const int size = StpVecGetSize(resetnodes_level);
   const int shift = g->knots * level;

   assert(0 <= level && level <= 3);

   /* clear visited list */
   for( int i = 0; i < size; i++ )
   {
      const int node = resetnodes_level[i];
      assert(graph_knot_isInRange(g, node));

      nodes_isvisited[node] = FALSE;
      state[node + shift] = UNKNOWN;
   }

#ifndef NDBEUG
   for( int i = 0; i < g->knots; ++i )
   {
      assert(!nodes_isvisited[i]);
   }

   for( int i = 0; i < g->knots * STP_TPATHS_NTERMBASES; ++i )
   {
      assert(UNKNOWN == state[i]);
   }
#endif

   SCIPfreeCleanBufferArray(scip, &(repair->nodes_isvisited));
}


/** traverses graph to find nodes that need to be reseted */
static
SCIP_RETCODE tpathsRepairTraverse1st(
   SCIP*                 scip,               /**< SCIP */
   int                   start,              /**< node to start from */
   int                   pred,               /**< predecessor node */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const TPATHS* tpaths = repair->tpaths;
   const PATH* termpaths = tpaths->termpaths;
   SCIP_Bool* RESTRICT nodes_isvisited = repair->nodes_isvisited;

   if( termpaths[start].edge >= 0 && g->tail[termpaths[start].edge] == pred )
   {
      STP_Vectype(int) resetnodes1st = repair->resetnodes[0];
      STP_Vectype(int) stack = repair->stack;

      assert(!Is_term(g->term[start]));
      assert(resetnodes1st);

      tpathsRepairResetNode(scip, start, resetnodes1st, stack, nodes_isvisited);

      SCIPdebugMessage("starting DFS from %d ... \n", start);

      /* DFS loop */
      while( !StpVecIsEmpty(stack) )
      {
         const int node = stack[StpVecGetSize(stack) - 1];
         StpVecPopBack(stack);

         assert(graph_knot_isInRange(g, node));

         for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];

            if( nodes_isvisited[head] )
               continue;

            if( termpaths[head].edge >= 0 && g->tail[termpaths[head].edge] == node )
            {
               assert(tpaths->termbases[node] == tpaths->termbases[head]);
               assert(head != pred);

               SCIPdebugMessage("add %d to reset-nodes \n", head);

               tpathsRepairResetNode(scip, head, resetnodes1st, stack, nodes_isvisited);
            }
         }
      }

      repair->stack = stack;
      repair->resetnodes[0] = resetnodes1st;
   }

   return SCIP_OKAY;
}


/** DFS on level */
static
void tpathsRepairTraverseLevelWithStack(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level from 1-3 */
   TREPAIR*              repair              /**< data for repairing */
)
{
   TPATHS* tpaths = repair->tpaths;
   const PATH* termpaths = tpaths->termpaths;
   const int* vbase_org = tpaths->termbases_org;
   SCIP_Bool* RESTRICT nodes_isvisited = repair->nodes_isvisited;
   const int nnodes = graph_get_nNodes(g);
   const int shift = level * nnodes;
   STP_Vectype(int) stack = repair->stack;
   STP_Vectype(int) resetnodes_level = repair->resetnodes[level];

   assert(scip && nodes_isvisited && resetnodes_level);
   assert(1 <= level && level <= 3);

   SCIPdebugMessage("starting DFS for level %d ... \n", level);

   /* DFS loop */
   while( !StpVecIsEmpty(stack) )
   {
      const int node = stack[StpVecGetSize(stack) - 1];
      const int nodebaseorg_top = vbase_org[node + shift];
      StpVecPopBack(stack);

      assert(graph_knot_isInRange(g, node));

      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( !nodes_isvisited[head] )
         {
            const int head_shifted = head + shift;

            if( termpaths[head_shifted].edge >= 0 && g->tail[termpaths[head_shifted].edge] == node )
            {
               if( nodebaseorg_top != vbase_org[head_shifted] )
                  continue;

               SCIPdebugMessage("level%d: add %d to reset-nodes \n", level, head);

               tpathsRepairResetNode(scip, head, resetnodes_level, stack, nodes_isvisited);
            }
         }
      }
   }

   repair->stack = stack;
   repair->resetnodes[level] = resetnodes_level;
}


/** adds outdated nodes whose predecessor points along edge to be deleted */
static
void tpathsRepairTraverseStackAddEdge(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level (1-3, level 0 is handled separately) */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const TPATHS* tpaths = repair->tpaths;
   const PATH* termpaths = tpaths->termpaths;
   const int nnodes = graph_get_nNodes(g);
   const int shift = level * nnodes;
   const int edge = repair->edge;
   const int tail = g->tail[edge];
   const int head = g->head[edge];

   if( termpaths[tail + shift].edge >= 0 && g->tail[termpaths[tail + shift].edge] == head )
   {
      tpathsRepairResetNode(scip, tail, repair->resetnodes[level], repair->stack, repair->nodes_isvisited);
   }

   if( termpaths[head + shift].edge >= 0 && g->tail[termpaths[head + shift].edge] == tail )
   {
      tpathsRepairResetNode(scip, head, repair->resetnodes[level], repair->stack, repair->nodes_isvisited);
   }
}


/** adds outdated nodes whose predecessor point to lower level */
static
void tpathsRepairTraverseStackAddBelow(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level (1-3, level 0 is handled separately) */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const TPATHS* tpaths = repair->tpaths;
   const PATH* termpaths = tpaths->termpaths;
   const int* vbase = tpaths->termbases;
   const int* vbase_org = tpaths->termbases_org;
   SCIP_Bool* RESTRICT nodes_isvisited = repair->nodes_isvisited;
   const int nnodes = graph_get_nNodes(g);
   const int shift = level * nnodes;

   assert(repair->resetnodes[level] && vbase && nodes_isvisited);
   assert(1 <= level && level <= 3);

   for( int d = 0; d < level; d++ )
   {
      const STP_Vectype(int) resetnodes_d = repair->resetnodes[d];
      const int size = StpVecGetSize(resetnodes_d);

      assert(size >= 0);

      for( int i = 0; i < size; i++ )
      {
         const int node = resetnodes_d[i];
         const int nodebase_d = vbase[node + d * nnodes];
         const int nodebaseorg_d = vbase_org[node + d * nnodes];
         assert(graph_knot_isInRange(g, node));

         /* duplicate? */
         if( !nodes_isvisited[node] && nodebase_d == vbase[node + shift] )
         {
            tpathsRepairResetNode(scip, node, repair->resetnodes[level], repair->stack, nodes_isvisited);
         }

         for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];

            if( !nodes_isvisited[head] )
            {
               const int head_shifted = head + shift;

               if( termpaths[head_shifted].edge >= 0 && g->tail[termpaths[head_shifted].edge] == node )
               {
                  if( nodebaseorg_d != vbase[head_shifted] )
                     continue;

                  assert(!Is_term(g->term[head]));

                  SCIPdebugMessage("level%d: add %d to reset-nodes \n", level, head);

                  tpathsRepairResetNode(scip, head, repair->resetnodes[level], repair->stack, nodes_isvisited);
               }
            }
         }
      }
   }
}


/** Traverses graph to find nodes on given level that need to be reseted,
 *  by updating from previous levels. Also considers duplicates. */
static
void tpathsRepairTraverseLevel(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level (1-3, level 0 is handled separately) */
   TREPAIR*              repair              /**< data for repairing */
)
{
   assert(1 <= level && level <= 3);

   tpathsRepairTraverseStackAddEdge(scip, g, level, repair);
   tpathsRepairTraverseStackAddBelow(scip, g, level, repair);

   tpathsRepairTraverseLevelWithStack(scip, g, level, repair);
}


/** updates reset nodes from non-reset nodes */
static
void tpathsRepairUpdate1st(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const STP_Vectype(int) resetnodes1st = repair->resetnodes[0];
   const SCIP_Bool* nodes_isvisited = repair->nodes_isvisited;
   TPATHS* tpaths = repair->tpaths;
   PATH* RESTRICT termpaths = tpaths->termpaths;
   int* RESTRICT termbases = tpaths->termbases;
   int* RESTRICT heap = g->path_heap;
   int* RESTRICT state = tpaths->state;
   int nHeapElems = 0;
   const int edge = repair->edge;
   const int edge_rev = flipedge(edge);
   const int size = StpVecGetSize(resetnodes1st);

   for( int i = 0; i < size; i++ )
   {
      const int node = resetnodes1st[i];
      assert(graph_knot_isInRange(g, node));
      assert(!Is_term(g->term[node]));

      pathResetSingle(termpaths[node]);
      termbases[node] = UNKNOWN;

      for( int e = g->inpbeg[node]; e != EAT_LAST; e = g->ieat[e] )
      {
         SCIP_Real dist;
         const int tail = g->tail[e];

         if( nodes_isvisited[tail] || e == edge || e == edge_rev )
            continue;

         dist = g->cost[e] + termpaths[tail].dist;

         if( dist < termpaths[node].dist )
         {
            termpaths[node].dist = dist;
            termpaths[node].edge = e;
            termbases[node] = termbases[tail];

            assert(graph_knot_isInRange(g, termbases[node]));
         }
      }

      if( termbases[node] != UNKNOWN )
      {
         SCIPdebugMessage("update node %d: base=%d, pred=%d dist=%f \n", node,
               termbases[node], g->tail[termpaths[node].edge], termpaths[node].dist);

         graph_pathHeapAdd(termpaths, node, heap, state, &nHeapElems);
      }
   }

   repair->nHeapElems = nHeapElems;
}




/** updates reset nodes from non-reset nodes of given level */
static
void tpathsRepairUpdateLevel(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   int                   level,              /**< level (1-3, 0 is handled separately) */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const STP_Vectype(int) resetnodes_top = repair->resetnodes[level];
   const SCIP_Bool* nodes_isvisited = repair->nodes_isvisited;
   TPATHS* tpaths = repair->tpaths;
   PATH* RESTRICT termpaths = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT heap = g->path_heap;
   int* RESTRICT state = tpaths->state;
   int nHeapElems = 0;
   const int nnodes = graph_get_nNodes(g);
   const int shift = level * nnodes;
   const int edge = repair->edge;
   const int edge_rev = flipedge(edge);
   const int size = StpVecGetSize(resetnodes_top);

   assert(1 <= level && level <= 3);

   for( int i = 0; i < size; i++ )
   {
      const int node = resetnodes_top[i];
      const int node_shift = node + shift;
      assert(graph_knot_isInRange(g, node));
      assert(!Is_term(g->term[node]));
      assert(nodes_isvisited[node]);

      pathResetSingle(termpaths[node_shift]);
      vbase[node_shift] = UNKNOWN;

      for( int e = g->inpbeg[node]; e != EAT_LAST; e = g->ieat[e] )
      {
         const int tail = g->tail[e];

         if( e == edge || e == edge_rev )
            continue;

         /* loop over all levels of tail */
         for( int d = 0; d <= level; d++ )
         {
            SCIP_Real dist;
            const int tail_d = tail + d * nnodes;
            const int tailbase_d = vbase[tail_d];
            int h;

            if( d == level && nodes_isvisited[tail] )
               break;

            for( h = 0; h < level; h++ )
               if( vbase[node + h * nnodes] == tailbase_d )
                  break;

            /* has tail_d a base from node? */
            if( h != level )
               continue;

            dist = g->cost[e] + termpaths[tail_d].dist;

            if( dist < termpaths[node_shift].dist )
            {
               termpaths[node_shift].dist = dist;
               termpaths[node_shift].edge = e;
               vbase[node_shift] = tailbase_d;

               assert(graph_knot_isInRange(g, vbase[node_shift]));
            }
         }
      }

      if( vbase[node_shift] != UNKNOWN )
      {
         SCIPdebugMessage("update node %d: base=%d, pred=%d dist=%f \n", node,
               vbase[node_shift], g->tail[termpaths[node_shift].edge], termpaths[node_shift].dist);

         graph_pathHeapAdd(termpaths, node_shift, heap, state, &nHeapElems);
      }
   }

   repair->nHeapElems = nHeapElems;
}


/** repairs TPATHS structure for imminent edge deletion (1st level) */
static
SCIP_RETCODE tpathsRepair1st(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const int edge = repair->edge;
   const int tail = g->tail[edge];
   const int head = g->head[edge];

   SCIP_CALL( tpathsRepairInitLevel(scip, 0, g, repair) );

   /* find nodes that need to be reseted */
   SCIP_CALL( tpathsRepairTraverse1st(scip, head, tail, g, repair) );
   SCIP_CALL( tpathsRepairTraverse1st(scip, tail, head, g, repair) );

   /* update newly found nodes from remaining nodes */
   tpathsRepairUpdate1st(scip, g, repair);

   /* complete the repair process for this level */
   tpathsScan1st(scip, g, g->cost, NULL, repair->nHeapElems, repair->tpaths,
       repair->withEdgeReplacement ? repair : NULL);

   tpathsRepairExitLevel(scip, 0, g, repair);

   return SCIP_OKAY;
}


/** repairs TPATHS structure for imminent edge deletion (2nd level) */
static
SCIP_RETCODE tpathsRepair2nd(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const int level = 1; /* NOTE: level is 0-indexed */

   SCIP_CALL( tpathsRepairInitLevel(scip, level, g, repair) );

   tpathsRepairTraverseLevel(scip, g, level, repair);
   tpathsRepairUpdateLevel(scip, g, level, repair);
   tpathsScan2nd(scip, g, g->cost,  g->cost, NULL, repair->nHeapElems, repair->tpaths,
         repair->withEdgeReplacement ? repair : NULL);

   tpathsRepairExitLevel(scip, level, g, repair);

   return SCIP_OKAY;
}

/** repairs TPATHS structure for imminent edge deletion (3rd level) */
static
SCIP_RETCODE tpathsRepair3rd(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const int level = 2; /* NOTE: level is 0-indexed */

   SCIP_CALL( tpathsRepairInitLevel(scip, level, g, repair) );

   tpathsRepairTraverseLevel(scip, g, level, repair);
   tpathsRepairUpdateLevel(scip, g, level, repair);
   tpathsScan3rd(scip, g, g->cost,  g->cost, NULL, repair->nHeapElems, repair->tpaths,
         repair->withEdgeReplacement ? repair : NULL);

   tpathsRepairExitLevel(scip, level, g, repair);

   return SCIP_OKAY;
}


/** repairs TPATHS structure for imminent edge deletion (4th level) */
static
SCIP_RETCODE tpathsRepair4th(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   const int level = 3; /* NOTE: level is 0-indexed */

   SCIP_CALL( tpathsRepairInitLevel(scip, level, g, repair) );

   tpathsRepairTraverseLevel(scip, g, level, repair);
   tpathsRepairUpdateLevel(scip, g, level, repair);
   tpathsScan4th(scip, g, g->cost,  g->cost, NULL, repair->nHeapElems, repair->tpaths,
         repair->withEdgeReplacement ? repair : NULL);

   tpathsRepairExitLevel(scip, level, g, repair);

   return SCIP_OKAY;
}


/** initializes */
static
void tpathsRepairInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{

#ifndef NDEBUG
   TPATHS* tpaths = repair->tpaths;
   for( int i = 0; i < STP_TPATHS_NTERMBASES * g->knots; i++ )
   {
      assert(tpaths->termbases_org[i] == tpaths->termbases[i]);
   }
#endif

   StpVecReserve(scip, repair->stack, STP_TPATHS_RESERVESIZE);

   for( int i = 0; i < STP_TPATHS_NTERMBASES; i++ )
   {
      assert(!repair->resetnodes[i]);
      StpVecReserve(scip, repair->resetnodes[i], STP_TPATHS_RESERVESIZE);
      assert(repair->resetnodes[i]);
   }
}


/** frees and resets */
static
void tpathsRepairExit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   TREPAIR*              repair              /**< data for repairing */
)
{
   TPATHS* tpaths = repair->tpaths;
   int* RESTRICT vbase_org = tpaths->termbases_org;
   const int* const vbase = tpaths->termbases;
   const int nnodes = graph_get_nNodes(g);

   for( int i = STP_TPATHS_NTERMBASES - 1; i >= 0; i-- )
   {
      STP_Vectype(int) resetnodes_i = repair->resetnodes[i];
      const int size = StpVecGetSize(resetnodes_i);
      const int shift = nnodes * i;

      assert(resetnodes_i);

      for( int j = 0; j < size; j++ )
      {
         const int node = resetnodes_i[j];
         const int node_shift = node + shift;
         assert(graph_knot_isInRange(g, node));

         vbase_org[node_shift] = vbase[node_shift];
      }

      StpVecFree(scip, repair->resetnodes[i]);
   }

   StpVecFree(scip, repair->stack);

#ifndef NDEBUG
   for( int i = 0; i < STP_TPATHS_NTERMBASES * nnodes; i++ )
   {
      assert(vbase_org[i] == vbase[i]);

   }
#endif
}


/*
 * Interface methods
 */


/** builds a Voronoi region w.r.t. shortest paths, for all terminals
 *  NOTE: serves as a wrapper */
void graph_add1stTermPaths(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path,               /**< path data structure (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path, .termbases = vbase, .state = state, .nnodes = nnodes };

   assert(path && vbase && state);

   graph_tpathsAdd1st(g, cost, NULL, &tpaths);
}

/** 2nd next terminal to all non terminal nodes
 *  NOTE: legacy wrapper */
void graph_add2ndTermPaths(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path2,              /**< path data structure (leading to first and second nearest terminal) */
   int*                  vbase2,             /**< first and second nearest terminal to each non terminal */
   int*                  state2              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path2, .termbases = vbase2, .state = state2, .nnodes = nnodes };

   assert(path2 && vbase2 && state2);

   graph_tpathsAdd2nd(g, cost, costrev, NULL, &tpaths);
}

/** 3rd next terminal to all non terminal nodes
 *  NOTE: legacy wrapper */
void graph_add3rdTermPaths(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path3,              /**< path data structure (leading to first, second and third nearest terminal) */
   int*                  vbase3,             /**< first, second and third nearest terminal to each non terminal */
   int*                  state3              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path3, .termbases = vbase3, .state = state3, .nnodes = nnodes };

   assert(path3 && vbase3 && state3);

   graph_tpathsAdd3rd(g, cost, costrev, NULL, &tpaths);
}


/* 4th next terminal to all non terminal nodes
 * NOTE: legacy wrapper */
void graph_add4thTermPaths(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path4,              /**< path data structure (leading to first, second, third and fourth nearest terminal) */
   int*                  vbase4,             /**< first, second, third, and fourth nearest terminal to each non terminal */
   int*                  state4              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path4, .termbases = vbase4, .state = state4, .nnodes = nnodes };

   assert(path4 && vbase4 && state4);

   graph_tpathsAdd4th(g, cost, costrev, NULL, &tpaths);
}



/** gets non-terminal shortest paths to 2 closest terminal for each non-terminal
 *  NOTE: legacy wrapper */
void graph_get2nextTermPaths(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path2,              /**< path data structure (leading to first, second terminal) */
   int*                  vbase2,             /**< first, second and third nearest terminal to each non terminal */
   int*                  state2              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path2, .termbases = vbase2, .state = state2, .nnodes = nnodes };

   assert(path2 && vbase2 && state2);

   graph_tpathsSetAll2(g, cost, costrev, NULL, &tpaths);
}


/** gets non-terminal shortest paths to 3 closest terminal for each non-terminal
 *  NOTE: legacy wrapper */
void graph_get3nextTermPaths(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path3,              /**< path data structure (leading to first, second and third nearest terminal) */
   int*                  vbase3,             /**< first, second and third nearest terminal to each non terminal */
   int*                  state3              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path3, .termbases = vbase3, .state = state3, .nnodes = nnodes };

   assert(path3 && vbase3 && state3);

   graph_tpathsSetAll3(g, cost, costrev, NULL, &tpaths);
}


/** gets non-terminal shortest paths to 4 closest terminal for each non-terminal
 *  NOTE: legacy wrapper */
void graph_get4nextTermPaths(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path4,              /**< path data struture (leading to first, second, third and fouth nearest terminal) */
   int*                  vbase4,             /**< first, second and third nearest terminal to each non terminal */
   int*                  state4              /**< array to mark the state of each node during calculation */
   )
{
   const int nnodes = graph_get_nNodes(g);
   TPATHS tpaths = { .termpaths = path4, .termbases = vbase4, .state = state4, .nnodes = nnodes };

   assert(path4 && vbase4 && state4);

   graph_tpathsSetAll4(g, cost, costrev, NULL, &tpaths);
}


/** get 4 close terminals to each terminal */
SCIP_RETCODE graph_get4nextTTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path,               /**< path data structure (leading to first, second, third and fourth nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int* boundedges;
   int k;
   int e;
   int l;
   int k2;
   int bedge;
   int shift;
   int nnodes;
   int nboundedges;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &boundedges, g->edges) );

   shift = 0;
   nnodes = g->knots;

   if( !graph_pc_isPcMw(g) )
      graph_mark(g);

   nboundedges = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         k2 = g->head[e];
         if( !g->mark[k2] || k2 < k )
            continue;

         /* is e a boundary edge? */
         if( vbase[k] != vbase[k2] )
         {
            boundedges[nboundedges++] = e;
         }
      }
      if( Is_term(g->term[k]) )
      {
         path[k].dist = FARAWAY;
         vbase[k] = UNKNOWN;
      }

   }

   for( l = 0; l < 4; l++ )
   {
      for( e = 0; e < nboundedges; e++ )
      {
         bedge = boundedges[e];
         k = g->tail[bedge];
         k2 = g->head[bedge];
         utdist(scip, g, path, cost[bedge], vbase, k, l, k2, shift, nnodes);
         utdist(scip, g, path, cost[bedge], vbase, k2, l, k, shift, nnodes);
      }
      shift += nnodes;
   }

   SCIPfreeBufferArray(scip, &boundedges);

   return SCIP_OKAY;
}



/** initializes TPATHS structure */
SCIP_RETCODE graph_tpathsInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   TPATHS**              tpaths              /**< the terminal paths */
)
{
   assert(scip && g);
   assert(STP_TPATHS_NTERMBASES == 4);

   SCIP_CALL( tpathsAlloc(scip, g, tpaths) );
   tpathsBuild(scip, g, *tpaths);

   return SCIP_OKAY;
}


/** initializes biased TPATHS structure */
SCIP_RETCODE graph_tpathsInitBiased(
   SCIP*                 scip,               /**< SCIP */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   TPATHS**              tpaths              /**< the terminal paths */
)
{
   assert(scip && g && sdprofit);
   assert(STP_TPATHS_NTERMBASES == 4);

   SCIP_CALL( tpathsAlloc(scip, g, tpaths) );
   tpathsBuildBiased(sdprofit, g, *tpaths);

   return SCIP_OKAY;
}


/** sets up TPATHS for subsequent edge deletions */
SCIP_RETCODE graph_tpathsRepairSetUp(
   const GRAPH*          g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT state = tpaths->state;

   if( !tpaths->termbases_org )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tpaths->termbases_org), nnodes * STP_TPATHS_NTERMBASES) );
   }

   BMScopyMemoryArray(tpaths->termbases_org, tpaths->termbases, nnodes * STP_TPATHS_NTERMBASES);

   for( int i = 0; i < nnodes * STP_TPATHS_NTERMBASES; i++ )
   {
      state[i] = UNKNOWN;
   }

   return SCIP_OKAY;
}


/** repairs TPATHS structure for imminent edge deletion */
SCIP_RETCODE graph_tpathsRepair(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge about to be eliminated */
   SCIP_Bool             withEdgeReplacement,/**< with edge replacement? */
   const GRAPH*          g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   STP_Vectype(int) resetnodes[STP_TPATHS_NTERMBASES] = { NULL, NULL, NULL, NULL };
   TREPAIR repair = { .resetnodes = resetnodes, .stack = NULL, .tpaths = tpaths, .nodes_isvisited = NULL,
         .edge = edge, .nHeapElems = 0, .withEdgeReplacement = withEdgeReplacement };

   assert(scip && g && tpaths);
   assert(STP_TPATHS_NTERMBASES == 4);
   assert(graph_edge_isInRange(g, edge));
   assert(graph_isMarked(g));
   assert(tpaths->termbases_org);

   tpathsRepairInit(scip, g, &repair);

   SCIP_CALL( tpathsRepair1st(scip, g, &repair) );
   SCIP_CALL( tpathsRepair2nd(scip, g, &repair) );
   SCIP_CALL( tpathsRepair3rd(scip, g, &repair) );
   SCIP_CALL( tpathsRepair4th(scip, g, &repair) );

   tpathsRepairExit(scip, g, &repair);

   return SCIP_OKAY;
}


/** recomputes biased TPATHS structure */
SCIP_RETCODE graph_tpathsRecomputeBiased(
   const SDPROFIT*       sdprofit,           /**< SD profit */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const
                                                 :( terrible design */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   assert(tpaths && g && sdprofit);
   assert(STP_TPATHS_NTERMBASES == 4);

   tpathsBuildBiased(sdprofit, g, tpaths);

   return SCIP_OKAY;
}


/** prints terminal paths from given node */
void graph_tpathsPrintForNode(
   const GRAPH*          g,                  /**< graph data structure */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   const TPATHS*         tpaths,             /**< storage for terminal paths */
   int                   node                /**< node to start from */
)
{
   const int* termbases;
   const int nnodes = graph_get_nNodes(g);
   int pos;

   assert(0 && "currently does not work, need to use vbase to find correct ancestor");
   assert(tpaths);
   assert(graph_knot_isInRange(g, node));

   printf("printing terminal paths for ");
   graph_knot_printInfo(g, node);

   if( Is_term(g->term[node]))
   {
      printf("single node path only: \n");
      graph_knot_printInfo(g, node);

      return;
   }

   termbases = tpaths->termbases;
   pos = node;

   for( int k = 0; k < STP_TPATHS_NTERMBASES; k++ )
   {
      const int base = termbases[pos];

      if( base != UNKNOWN )
      {
         assert(graph_knot_isInRange(g, base));

         printf("Path to terminal (distance=%f) ", tpaths->termpaths[pos].dist);
         graph_knot_printInfo(g, base);

         tpathsPrintPath(g, sdprofit, tpaths, node, k + 1);
      }
      pos += nnodes;
   }

}


/** gets nodes with positive profit on path */
void graph_tpathsGetProfitNodes(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const TPATHS*         tpaths,             /**< storage for terminal paths */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes */
   int                   start,              /**< start vertex */
   int                   base,               /**< base vertex (terminal) */
   STP_Vectype(int)      profitnodes         /**< NOTE: needs to be of capacity at least g->knots! */
)
{
   const int nnodes = graph_get_nNodes(g);
   int level;

   assert(tpaths && profitnodes && sdprofit);
   assert(graph_knot_isInRange(g, start));
   assert(graph_knot_isInRange(g, base));
   assert(Is_term(g->term[base]));
   assert(StpVecGetcapacity(profitnodes) >= g->knots);

   StpVecClear(profitnodes);

   if( base != start )
   {
      const int* const termbases = tpaths->termbases;
      const PATH* const termpaths = tpaths->termpaths;
      int node_pred;
      int node_curr;
      int node_next;
      int edge_curr;

      for( level = 0; termbases[start + nnodes * level] != base; level++ )
      {
         assert(level <= STP_TPATHS_NTERMBASES);
      }
      assert(termbases[start + nnodes * level] == base);
      assert(0 <= level && level <= STP_TPATHS_NTERMBASES);

      node_pred = -1;
      node_curr = start;
      edge_curr = termpaths[node_curr + level * nnodes].edge;
      assert(edge_curr >= 0);
      node_next = g->tail[edge_curr];

   //   printf("%d->%d \n", start, base);

      while( node_curr != base )
      {
         assert(graph_knot_isInRange(g, node_curr));
         assert(graph_edge_isInRange(g, edge_curr));
         assert(termbases[node_curr + level * nnodes] == base);

       //  printf("%d \n", node_curr);

         if( node_curr != start )
         {
            const SCIP_Real profit = reduce_sdprofitGetProfit(sdprofit, node_curr, node_pred, node_next);
            if( GT(profit, 0.0) )
            {
               StpVecPushBack(scip, profitnodes, node_curr);
            }
         }

         node_pred = node_curr;
         node_curr = node_next;

         while( termbases[node_curr + level * nnodes] != base )
         {
            level--;
            assert(level >= 0);
         }

         edge_curr = termpaths[node_curr + level * nnodes].edge;
         node_next = (edge_curr >= 0) ? g->tail[edge_curr] : -1;
      }
   }

#ifndef NDEBUG
   for( int k = 0; k < StpVecGetSize(profitnodes); k++ )
   {
      const int node = profitnodes[k];
      assert(graph_knot_isInRange(g, node));
      assert(!Is_term(g->term[node]));
   }
#endif
}



/** Computes next terminal to all non-terminal nodes.
 *  Essentially a Voronoi diagram */
void graph_tpathsAdd1st(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   int nHeapElems = 0;
   int nbases = 0;
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;

   assert(path != NULL);
   assert(cost != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(nnodes == tpaths->nnodes);

   /* initialize */
   for( int i = 0; i < nnodes; i++ )
   {
      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
         {
            heap[++nHeapElems] = i;
         }

         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = nHeapElems;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i] = UNKNOWN;
      }
   }

   if( nbases == 0 || nnodes <= 1 )
      return;

   tpathsScan1st(NULL, g, cost, sdprofit, nHeapElems, tpaths, NULL);
}


/** computes 2nd next terminal to all non-terminal nodes */
void graph_tpathsAdd2nd(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   int nheapElems;
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;

   assert(nnodes == tpaths->nnodes);
   assert(path != NULL);
   assert(cost != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(costrev != NULL);

   nheapElems = 0;
   tpathsResetMembers(g, 2, tpaths);

   /* scan original nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( !Is_term(g->term[head]) && vbase[k] != vbase[head] && g->mark[head] )
         {
            const int head2 = head + nnodes;
            const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, k, head2, e, sdprofit, path);

            if( GT(path[head2].dist, distnew) )
            {
               pathheapCorrectDist(heap, state, &nheapElems, path, head2, k, e, distnew);
               vbase[head2] = vbase[k];
            }
         }
      }
   }

   if( nnodes <= 1 )
      return;

   tpathsScan2nd(NULL, g, cost, costrev, sdprofit, nheapElems, tpaths, NULL);
}


/** computes 3rd next terminal to all non-terminal nodes */
void graph_tpathsAdd3rd(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   int nheapElems = 0;
   const int nnodes = graph_get_nNodes(g);
   const int dnnodes = 2 * nnodes;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;

   assert(nnodes == tpaths->nnodes);
   assert(path != NULL);
   assert(cost != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(costrev != NULL);

   tpathsResetMembers(g, 3, tpaths);

   /* scan original nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         const int head3 = head + dnnodes;

         if( !Is_term(g->term[head]) && g->mark[head] )
         {
            int kn = k;

            for( int level = 0; level < 2; level++ )
            {
               if( vbase[kn] != vbase[head] && vbase[kn] != vbase[head + nnodes] )
               {
                  const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, kn, head3, e, sdprofit, path);

                  if( GT(path[head3].dist, distnew) )
                  {
                     pathheapCorrectDist(heap, state, &nheapElems, path, head3, kn, e, distnew);
                     vbase[head3] = vbase[kn];
                  }
               }

               kn += nnodes;
            }
         }
      }
   }

   if( nnodes <= 1 )
      return;

   tpathsScan3rd(NULL, g, cost, costrev, sdprofit, nheapElems, tpaths, NULL);
}



/** computes 4th next terminal to all non-terminal nodes */
void graph_tpathsAdd4th(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   int nHeapElems = 0;
   const int nnodes = graph_get_nNodes(g);
   const int dnnodes = 2 * nnodes;
   const int tnnodes = 3 * nnodes;
   int* RESTRICT heap = g->path_heap;
   PATH* RESTRICT path = tpaths->termpaths;
   int* RESTRICT vbase = tpaths->termbases;
   int* RESTRICT state = tpaths->state;

   assert(nnodes == tpaths->nnodes);
   assert(path != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(costrev != NULL);

   tpathsResetMembers(g, 4, tpaths);

   /* scan original nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         const int head4 = head + tnnodes;

         if( !Is_term(g->term[head]) && g->mark[head] )
         {
            int kn = k;

            for( int level = 0; level < 3; level++ )
            {
               if( vbase[kn] != vbase[head] && vbase[kn] != vbase[head + nnodes] && vbase[kn] != vbase[head + dnnodes] )
               {
                  const SCIP_Real distnew = tpathsGetDistNew(g, cost, costrev, vbase, kn, head4, e, sdprofit, path);

                  if( GT(path[head4].dist, distnew))
                  {
                     pathheapCorrectDist(heap, state, &nHeapElems, path, head4, kn, e, distnew);
                     vbase[head4] = vbase[kn];
                  }
               }
               kn += nnodes;
            }
         }
      }
   }

   if( nnodes <= 1 )
      return;

   tpathsScan4th(NULL, g, cost, costrev, sdprofit, nHeapElems, tpaths, NULL);
}


/** computes 2 next terminal to all non-terminal nodes */
void graph_tpathsSetAll2(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(tpaths      != NULL);

   if( !graph_pc_isPcMw(g) )
      graph_mark(g);

   graph_tpathsAdd1st(g, cost, sdprofit, tpaths);
   graph_tpathsAdd2nd(g, cost, costrev, sdprofit, tpaths);

#ifndef NDEBUG
   if( !sdprofit )
   {
      const PATH* RESTRICT path2 = tpaths->termpaths;
      const int nnodes = graph_get_nNodes(g);

      for( int level = 0; level < 1; level++ )
      {
         for( int k = 0; k < nnodes; ++k )
         {
            assert(LE_FEAS_EPS(path2[level * nnodes + k].dist, path2[(level + 1) * nnodes + k].dist, EPSILON));
         }
      }
   }
#endif

   return;
}

/** computes 3 next terminal to all non-terminal nodes */
void graph_tpathsSetAll3(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(tpaths      != NULL);

   if( !graph_pc_isPcMw(g) )
      graph_mark(g);

   graph_tpathsAdd1st(g, cost, sdprofit, tpaths);
   graph_tpathsAdd2nd(g, cost, costrev, sdprofit, tpaths);
   graph_tpathsAdd3rd(g, cost, costrev, sdprofit, tpaths);

#ifndef NDEBUG
   if( !sdprofit )
   {
      const PATH* RESTRICT path3 = tpaths->termpaths;
      const int nnodes = graph_get_nNodes(g);

      for( int level = 0; level < 2; level++ )
      {
         for( int k = 0; k < nnodes; ++k )
         {
            assert(LE_FEAS_EPS(path3[level * nnodes + k].dist, path3[(level + 1) * nnodes + k].dist, EPSILON));
         }
      }
   }
#endif

   return;
}


/** computes 4 next terminal to all non-terminal nodes */
void graph_tpathsSetAll4(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const SDPROFIT*       sdprofit,           /**< SD bias for nodes or NULL */
   TPATHS*               tpaths              /**< storage for terminal paths */
)
{
   assert(g         != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);

   graph_mark(g);

   graph_tpathsAdd1st(g, cost, sdprofit, tpaths);
   graph_tpathsAdd2nd(g, cost, costrev, sdprofit, tpaths);
   graph_tpathsAdd3rd(g, cost, costrev, sdprofit, tpaths);
   graph_tpathsAdd4th(g, cost, costrev, sdprofit, tpaths);

#ifndef NDEBUG
   if( !sdprofit )
   {
      PATH* RESTRICT path4 = tpaths->termpaths;
      const int nnodes = graph_get_nNodes(g);

      for( int level = 0; level < 3; level++ )
      {
         for( int k = 0; k < nnodes; ++k )
         {
            assert(LE(path4[level * nnodes + k].dist, path4[(level + 1) * nnodes + k].dist));
         }
      }
   }
#endif

   return;
}


/** frees TPATHS structure */
void graph_tpathsFree(
   SCIP*                 scip,               /**< SCIP */
   TPATHS**              tpaths              /**< the terminal paths */
)
{
   TPATHS* tp;
   assert(scip && tpaths);

   tp = *tpaths;
   assert(tp);

   SCIPfreeMemoryArrayNull(scip, &(tp->termbases_org));
   SCIPfreeMemoryArray(scip, &(tp->state));
   SCIPfreeMemoryArray(scip, &(tp->termbases));
   SCIPfreeMemoryArray(scip, &(tp->termpaths));

   SCIPfreeMemory(scip, tpaths);
}


/** gets (up to) four close terminals to given node i;
 *  with strict upper bound on allowed distances */
void graph_tpathsGetClosestTerm(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   node,               /**< node */
   int* RESTRICT         closeterm,          /**< terminal */
   int* RESTRICT         firstedge,          /**< corresponding first edge (can be NULL) */
   SCIP_Real* RESTRICT   closeterm_dist      /**< terminal distance */
)
{
   assert(g && tpaths && closeterm && closeterm_dist);
   assert(graph_knot_isInRange(g, node));

   if( Is_term(g->term[node]) )
   {
      *closeterm = node;
      *closeterm_dist = 0.0;
      if( firstedge )
         *firstedge = UNKNOWN;
   }
   else
   {
      const PATH* const termpaths = tpaths->termpaths;
      const int* const termbases = tpaths->termbases;

      *closeterm = termbases[node];
      *closeterm_dist = termpaths[node].dist;
      if( firstedge )
         *firstedge = termpaths[node].edge;
   }
}


/** gets (up to) three close terminals to given node i;
 *  with strict upper bound on allowed distances */
void graph_tpathsGet3CloseTerms(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   node,               /**< node */
   SCIP_Real             maxdist_strict,     /**< maximum valid distance (strict) */
   int* RESTRICT         closeterms,         /**< three close terminals */
   int* RESTRICT         firstedges,         /**< corresponding first edge (can be NULL) */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< three close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const SCIP_Bool isStrict = TRUE;
   tpathsGetKCloseTerms(g, tpaths, 3, node, maxdist_strict, isStrict, closeterms, firstedges, closeterms_dist, ncloseterms);
}


/** gets (up to) four close terminals to given node i;
 *  with strict upper bound on allowed distances */
void graph_tpathsGet4CloseTerms(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   node,               /**< node */
   SCIP_Real             maxdist_strict,     /**< maximum valid distance (strict) */
   int* RESTRICT         closeterms,         /**< four close terminals */
   int* RESTRICT         firstedges,         /**< corresponding first edge (can be NULL) */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< four close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const SCIP_Bool isStrict = TRUE;
   tpathsGetKCloseTerms(g, tpaths, 4, node, maxdist_strict, isStrict, closeterms, firstedges, closeterms_dist, ncloseterms);
}


/** gets (up to) four close terminals to given node i
  *  with non-strict upper bound on allowed distances */
void graph_tpathsGet4CloseTermsLE(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   node,               /**< node */
   SCIP_Real             maxdist_nonstrict,  /**< maximum valid distance (not strict) */
   int* RESTRICT         closeterms,         /**< four close terminals */
   int* RESTRICT         firstedges,         /**< corresponding first edge (can be NULL) */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< four close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const SCIP_Bool isStrict = FALSE;
   tpathsGetKCloseTerms(g, tpaths, 4, node, maxdist_nonstrict, isStrict, closeterms, firstedges, closeterms_dist, ncloseterms);
}
