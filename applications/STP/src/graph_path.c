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

/**@file   graph_path.c
 * @brief  Shortest path based graph algorithms for Steiner problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file encompasses various (heap-based) shortest path based algorithms including
 * Dijkstra's algorithm.
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
#include "shortestpath.h"

#define STP_TPATHS_NTERMBASES 4



/** Steiner nodes to terminal paths
 * NOTE: all arrays are of size STP_TPATHS_NTERMBASES * nnodes */
struct nodes_to_terminal_paths
{
   PATH*                 termpaths;         /**< path data (leading to first, second, ... terminal) */
   int*                  termbases;         /**< terminals to each non terminal */
   int*                  state;             /**< array to mark the state of each node during calculation */
   int                   nnodes;            /**< number of nodes of underlying graph */
};



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


/*---------------------------------------------------------------------------*/
/*--- Name     : CORRECT heap                                             ---*/
/*--- Function : Setzt ein neues Element auf den Heap, bzw. korrigiert    ---*/
/*---            die Position eines vorhandenen Elementes                 ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzten Kanten,            ---*/
/*---            Neuer Knoten, Vorgaengerknoten, Kante von der man aus    ---*/
/*---            den neuen Knoten erreicht, Kosten der Kante,             ---*/
/*---            sowie Betriebsmodus                                      ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
inline static
void pathheapCorrect(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   PATH* RESTRICT path,
   int    l,
   int    k,
   int    edge,
   SCIP_Real cost,
   int    mode)
{
   int    t;
   int    c;
   int    j;

   path[l].dist = (mode == MST_MODE) ? cost : (path[k].dist + cost);
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

/** resets node */
inline static
void heapReset(
   PATH* RESTRICT path,
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,
   int    node
   )
{
   int    t;
   int    c;
   int    j;

   path[node].dist = 0.0;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* heap shift up */
   j = state[node];
   c = j / 2;

   while( (j > 1) && GT(path[heap[c]].dist, path[heap[j]].dist) )
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


/** connect given node to tree */
static inline
void stPcmwConnectNode(
   int                   k,                  /**< the vertex */
   const GRAPH*          g,                  /**< graph data structure */
   SPATHSPC*             spaths_pc,          /**< shortest paths data */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  count,              /**< for the heap */
   int*                  nterms              /**< terminal count */
)
{
   assert(k >= 0);

   connected[k] = TRUE;
   pathdist[k] = 0.0;
   shortestpath_pcConnectNode(g, connected, k, spaths_pc);
   (*nterms)++;

   assert(pathedge[k] != -1);

   /* connect k to current subtree */
   for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
   {
      connected[node] = TRUE;
      resetX(pathdist, g->path_heap, g->path_state, count, node, 0.0);

      if( Is_pseudoTerm(g->term[node]) )
      {
         shortestpath_pcConnectNode(g, connected, node, spaths_pc);
         (*nterms)++;
      }

      assert(pathedge[node] != -1);
   }
}


/** initialize */
static inline
void stPcmwInit(
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  npseudoterms        /**< number of pseudo terminals */
)
{
   const int nnodes = graph_get_nNodes(g);
   int ntermspos = 0;
   int* const RESTRICT state = g->path_state;

   for( int k = 0; k < nnodes; k++ )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k] = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( Is_pseudoTerm(g->term[k]) )
         ntermspos++;
   }

   if( npseudoterms )
      *npseudoterms = ntermspos;
}


/** initialize */
static inline
void stRpcmwInit(
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  nrealterms          /**< number of real terminals */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nrterms = 0;
   int* const RESTRICT state = g->path_state;

   /* unmark dummy terminals */
   graph_pc_markOrgGraph(g);
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   for( int k = 0; k < nnodes; k++ )
   {
      state[k] = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(g->mark[k]);
         nrterms++;
         assert(!Is_pseudoTerm(g->term[k]));
      }
   }

   if( nrealterms )
      *nrealterms = nrterms;
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
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   for( int i = 0; i < offset; i++ )
   {
      state[i] = CONNECT;
   }
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

   if( path[node].edge >= 0 )
   {
      assert(!Is_term(g->term[node % g->knots]));
   }
   else
   {
      assert(Is_term(g->term[node % g->knots]));
   }
#endif

   if( sdprofit && path[node].edge >= 0 )
   {
      const int nnodes = graph_get_nNodes(g);
      const int node_pred = g->tail[path[node].edge];

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
   assert(GE(distnew, 0.0) && LE(distnew, path[node].dist + ecost));

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
   SCIP*                 scip,               /**< SCIP */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   GRAPH*                g,                  /**< graph */
   TPATHS*               tpaths              /**< the terminal paths */
)
{
   assert(STP_TPATHS_NTERMBASES == 4);

   graph_tpathsSetAll4(g, g->cost, g->cost, sdprofit, tpaths);
}


/*
 * Interface methods
 */

/** adds element 'node' to heap */
void graph_pathHeapAdd(
   const PATH* path,
   int node,               /* the node to be added */
   int* RESTRICT heap,     /* heap array */
   int* RESTRICT state,
   int* count             /* pointer to store the number of elements on the heap */
   )
{
   int c;
   int j;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* Heap shift up */
   j = state[node];
   c = j / 2;

   while((j > 1) && GT(path[heap[c]].dist, path[heap[j]].dist))
   {
      const int t = heap[c];
      heap[c] = heap[j];
      heap[j] = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j = c;
      c = j / 2;
   }
}


/*---------------------------------------------------------------------------*/
/*--- Name     : INIT shortest PATH algorithm                             ---*/
/*--- Function : Initialisiert den benoetigten Speicher fuer die          ---*/
/*---            kuerzeste Wege Berechnung                                ---*/
/*--- Parameter: Graph in dem berechnet werden soll.                      ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
SCIP_RETCODE graph_path_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   assert(g != NULL);
   assert(g->path_heap  == NULL);
   assert(g->path_state == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->path_heap), g->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->path_state), g->knots) );

   return SCIP_OKAY;
}

/*---------------------------------------------------------------------------*/
/*--- Name     : EXIT shortest PATH algorithm                             ---*/
/*--- Function : Gibt den bei der initialisierung angeforderten Speicher  ---*/
/*---            wieder frei                                              ---*/
/*--- Parameter: Keine                                                    ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
void graph_path_exit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   assert(g != NULL);
   assert(g->path_heap  != NULL);
   assert(g->path_state != NULL);

   SCIPfreeMemoryArray(scip, &(g->path_state));
   SCIPfreeMemoryArray(scip, &(g->path_heap));
}

/*---------------------------------------------------------------------------*/
/*--- Name     : FIND shortest PATH / minimum spanning tree               ---*/
/*--- Function : Dijkstras Algorithmus zur bestimmung eines kuerzesten    ---*/
/*---            Weges in einem gerichteten Graphen.                      ---*/
/*--- Parameter: Graph, Nummer des Startknotens und Betriebsmodus         ---*/
/*---            (Kuerzeste Wege oder minimal spannender Baum)            ---*/
/*--- Returns  : Liefert einen Vektor mit einem PATH Element je Knotem.   ---*/
/*----           Setzt das .dist Feld zu jedem Knoten auf die Entfernung  ---*/
/*---            zum Startknoten, sowie das .prev Feld auf den Knoten von ---*/
/*---            dem man zum Knoten gekommen ist. Im MST Modus steht im   ---*/
/*---            .dist Feld die Entfernung zum naechsten Knoten.          ---*/
/*---------------------------------------------------------------------------*/
void graph_path_exec(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          p,                  /**< graph data structure */
   const int             mode,               /**< shortest path (FSP_MODE) or minimum spanning tree (MST_MODE)? */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path                /**< shortest paths data structure */
   )
{
   int* heap;
   int* state;
   int k;
   const int nnodes = p->knots;
   int count;

   assert(scip      != NULL);
   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) || (mode == MST_MODE));
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(p->mark[start]);

   /* no nodes?, return*/
   if( nnodes == 0 )
      return;

   heap = p->path_heap;
   state = p->path_state;

   /* initialize */
   for( int i = 0; i < nnodes; i++ )
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY + 1;
      path[i].edge = UNKNOWN;
   }
   /* add first node to heap */
   k            = start;
   path[k].dist = 0.0;

   if( nnodes > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         /* get nearest labeled node */
         k = pathheapGetNearest(heap, state, &count, path);

         /* mark as scanned */
         state[k] = CONNECT;

         for( int i = p->outbeg[k]; i >= 0; i = p->oeat[i] )
         {
            const int m = p->head[i];

            assert(i != EAT_LAST);

            /* node not scanned and valid? */
            if( state[m] )
            {
               /* closer than previously and valid? */
               if( path[m].dist > ((mode == MST_MODE) ? cost[i] : (path[k].dist + cost[i])) && p->mark[m] )
                  pathheapCorrect(heap, state, &count, path, m, k, i, cost[i], mode);
            }
         }
      }
   }
}

/** limited Dijkstra, stopping at terminals */
void graph_sdPaths(
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  memlbl,             /**< array to save labelled nodes */
   int*                  nlbl,               /**< number of labelled nodes */
   int                   tail,               /**< tail of the edge */
   int                   head,               /**< head of the edge */
   int                   limit               /**< maximum number of edges to consider during execution */
   )
{
   int count;
   int nchecks;
   const int limit1 = limit / 2;

   assert(g      != NULL);
   assert(heap   != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(nlbl   != NULL);
   assert(memlbl != NULL);
   assert(limit1 >= 0);

   *nlbl = 0;

   if( g->grad[tail] == 0 || g->grad[head] == 0 )
      return;

   assert(g->mark[head] && g->mark[tail]);

   count = 0;
   nchecks = 0;
   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;

   /* NOTE: for MW we do not consider the edge between tail and head */
   if( !graph_pc_isMw(g) )
      g->mark[head] = FALSE;

   for( int e = g->outbeg[tail]; e >= 0; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && (GE(distlimit, cost[e])) )
      {
         assert(GT(path[m].dist, path[tail].dist + cost[e]));

         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         pathheapCorrect(heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }

   g->mark[head] = TRUE;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      const int k = pathheapGetNearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* distance limit reached? */
      if( GT(path[k].dist, distlimit) )
         break;

      /* stop at terminals */
      if( Is_term(g->term[k]) || k == head )
         continue;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( state[m] && g->mark[m] && GE(distlimit, cost[e]) && (path[m].dist > path[k].dist + cost[e]) )
         {
            /* m labelled for the first time? */
            if( state[m] == UNKNOWN )
               memlbl[(*nlbl)++] = m;
            pathheapCorrect(heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
         if( nchecks++ > limit )
            break;
      }
   }
}


/** limited Dijkstra for PCSTP, taking terminals into account */
void graph_path_PcMwSd(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int*                  pathmaxnode,        /**< maximum weight node on path */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  stateblock,         /**< array to indicate whether a node has been scanned during previous SP calculation */
   int*                  memlbl,             /**< array to save labelled nodes */
   int*                  nlbl,               /**< number of labelled nodes */
   int                   tail,               /**< tail of the edge */
   int                   head,               /**< head of the edge */
   int                   limit               /**< maximum number of edges to consider during execution */
   )
{
   const int limit1 = limit / 2;
   int count;
   int nchecks;
   const SCIP_Bool block = (stateblock != NULL);

   assert(limit > 0);
   assert(g      != NULL);
   assert(heap   != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(nlbl   != NULL);
   assert(memlbl != NULL);
   assert(g->prize != NULL);
   assert(pathmaxnode != NULL);

   *nlbl = 0;

   if( g->grad[tail] == 0 || g->grad[head] == 0 )
      return;

   assert(g->mark[head] && g->mark[tail]);

   nchecks = 0;
   count = 0;
   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;

   if( g->stp_type != STP_MWCSP )
      g->mark[head] = FALSE;

   for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(SCIPisGT(scip, path[m].dist, path[tail].dist + cost[e]));

         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         pathheapCorrect(heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }

   g->mark[head] = TRUE;

   /* main loop */
   while( count > 0 )
   {
      const int k = pathheapGetNearest(heap, state, &count, path);
      SCIP_Real maxweight = pathmaxnode[k] >= 0 ? g->prize[pathmaxnode[k]] : 0.0;

      assert(k != tail);
      assert(maxweight >= 0);
      assert(SCIPisLE(scip, path[k].dist - maxweight, distlimit));

      /* scanned */
      state[k] = CONNECT;

      /* stop at other end */
      if( k == head )
         continue;

      if( Is_term(g->term[k]) && g->prize[k] > maxweight && distlimit >= path[k].dist )
      {
         pathmaxnode[k] = k;
         maxweight = g->prize[k];
      }

      /* stop at node scanned in first run */
      if( block && stateblock[k] == CONNECT )
         continue;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( state[m] && g->mark[m] && path[m].dist > (path[k].dist + cost[e])
               && distlimit >= (path[k].dist + cost[e] - maxweight) )
         {
            if( state[m] == UNKNOWN ) /* m labeled for the first time? */
               memlbl[(*nlbl)++] = m;

            pathmaxnode[m] = pathmaxnode[k];
            pathheapCorrect(heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
         if( nchecks++ > limit )
            break;
      }
   }
}



/** Dijkstra's algorithm starting from node 'start' */
void graph_path_execX(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge            /**< predecessor edge array (on vertices) */
   )
{
   int   k;
   int   m;
   int   i;
   int   count;
   int   nnodes;
   int* heap;
   int* state;

   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cost   != NULL);

   nnodes = g->knots;

   if( nnodes == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;

   for( i = nnodes - 1; i >= 0; --i )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k = start;
   pathdist[k] = 0.0;

   if( nnodes > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         k = nearestX(heap, state, &count, pathdist);

         state[k] = CONNECT;

         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( state[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[i])) )
               correctX(heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}

/** Dijkstra on incoming edges until root is reached */
void graph_path_invroot(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge            /**< predecessor edge array (on vertices) */
   )
{
   SCIP_Real rootdist;
   int   k;
   int   m;
   int   i;
   int   count;
   int   nnodes;
   int* heap;
   int* state;

   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cost   != NULL);

   nnodes = g->knots;

   if( nnodes == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;
   rootdist = FARAWAY;

   for( i = nnodes - 1; i >= 0; --i )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k = start;
   pathdist[k] = 0.0;

   if( nnodes > 1 )
   {
      int root = g->source;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         k = nearestX(heap, state, &count, pathdist);

         state[k] = CONNECT;

         if( k == root )
            rootdist = pathdist[k];
         else if( SCIPisGT(scip, pathdist[k], rootdist) )
            break;

         for( i = g->inpbeg[k]; i != EAT_LAST; i = g->ieat[i] )
         {
            m = g->tail[i];

            if( state[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[i])) )
               correctX(heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}


/** extension heuristic */
void graph_path_st_pcmw_extendOut(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Real*            dist,               /**< distances array */
   int*                  pred,               /**< predecessor node */
   STP_Bool*             connected_out,      /**< array for internal stuff */
   DHEAP*                dheap,              /**< Dijkstra heap */
   SCIP_Bool*            success             /**< extension successful? */
   )
{
   int* const state = dheap->position;
   const CSR* const csr = g->csr_storage;
   const int* const start_csr = csr->start;
   const int* const head_csr = csr->head;
   const SCIP_Real* const cost_csr = csr->cost;
   SCIP_Real outerprofit;
   SCIP_Real outermaxprize;
   const int nnodes = g->knots;

   assert(csr && g && dist && dheap && pred && connected && connected_out);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(!connected[start]);

   *success = FALSE;
   outermaxprize = 0.0;

   for( int k = 0; k < nnodes; k++ )
   {
      state[k] = UNKNOWN;
      dist[k] = FARAWAY;
      connected_out[k] = FALSE;
#ifndef NDEBUG
      pred[k] = -1;
#endif

      if( !connected[k] && Is_term(g->term[k]) && g->prize[k] > outermaxprize && k != start )
         outermaxprize = g->prize[k];
   }

   graph_heap_clean(FALSE, dheap);

   dist[start] = 0.0;
   graph_heap_correct(start, 0.0, dheap);
   connected_out[start] = TRUE;
   outerprofit = g->prize[start];

   for( int rounds = 0; rounds < 2 && *success == FALSE; rounds++  )
   {
      if( rounds == 1 )
      {
         /* no improvement in last round? */
         if( !SCIPisGT(scip, outerprofit, g->prize[start]) )
            break;

         if( dheap->size > 0 )
            graph_heap_clean(TRUE, dheap);

         assert(dheap->size == 0);

         /* insert outer tree vertices into heap */
         for( int k = 0; k < nnodes; k++ )
         {
            if( connected_out[k] )
            {
               dist[k] = 0.0;
               graph_heap_correct(k, 0.0, dheap);
            }
            else
               dist[k] = FARAWAY;
         }
      }

      while( dheap->size > 0 )
      {
         /* get nearest labelled node */
         const int k = graph_heap_deleteMinReturnNode(dheap);
         const int k_start = start_csr[k];
         const int k_end = start_csr[k + 1];
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect k to current subtree */
         if( (connected[k] && SCIPisGT(scip, outerprofit, dist[k]))
               || (!connected[k] && !connected_out[k] && Is_term(g->term[k]) && SCIPisGE(scip, g->prize[k], dist[k])) )
         {
            assert(k != start);
            assert(pred[k] != -1);
            assert(!connected_out[k] || !connected[k]);

            outerprofit += g->prize[k] - dist[k];
            connected_out[k] = TRUE;
            dist[k] = 0.0;

            assert(SCIPisGE(scip, outerprofit, g->prize[start]) || connected[k]);

            /* connect k to current subtree */
            for( int node = pred[k]; !connected_out[node]; node = pred[node] )
            {
               connected_out[node] = TRUE;
               dist[node] = 0.0;
               graph_heap_correct(node, 0.0, dheap);

               if( Is_term(g->term[node]) )
                  outerprofit += g->prize[node];

               assert(state[node]);
               assert(pred[node] >= 0);
            }

            if( connected[k] )
            {
               *success = TRUE;
               break;
            }
         }
         else if( outerprofit + outermaxprize < dist[k] )
         {
            assert(*success == FALSE);
            break;
         }

         /* correct incident nodes */
         for( int e = k_start; e < k_end; e++ )
         {
            const int m = head_csr[e];
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( distnew < dist[m] )
            {
               dist[m] = distnew;
               pred[m] = k;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   if( *success )
   {
      for( int k = 0; k < nnodes; k++ )
         if( connected_out[k] )
            connected[k] = TRUE;
   }
}

/** Find a directed tree rooted in node 'start' and spanning all terminals */
void graph_path_st(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edgecosts */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int   count;
   int   nnodes;
   int*  heap;
   int*  state;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   count = 0;
   nnodes = g->knots;
   heap = g->path_heap;
   state = g->path_state;

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int nterms = 0;

      if( Is_term(g->term[start]) )
         nterms++;

      /* add start vertex to heap */
      count = 1;
      heap[count] = start;
      state[start] = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* k is terminal an not connected yet? */
         if( Is_term(g->term[k]) && k != start )
         {
            assert(pathedge[k] >= 0 && !connected[k]);

            connected[k] = TRUE;
            pathdist[k] = 0.0;

            /* connect k to current solution */
            for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
            {
               assert(pathedge[node] != -1);
               assert(!Is_term(g->term[node]));

               connected[node] = TRUE;
               resetX(pathdist, heap, state, &count, node, 0.0);
            }

            /* have all terminals been reached? */
            if( ++nterms == g->terms )
               break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) && g->mark[m] )
               correctX(heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}



/** !!!LEGACY CODE!!!
 *  Find a tree rooted in node 'start' and connecting
 *  positive vertices as long as this is profitable.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw(
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            orderedprizes,      /**< legacy code */
   int*                  orderedprizes_id,   /**< legacy code */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Bool             costIsBiased,       /**< is cost biased? */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   const int nnodes = g->knots;
   int* const RESTRICT heap = g->path_heap;
   int* const RESTRICT state = g->path_state;
   int ntermspos = -1;
   SPATHSPC spaths_pc = { .orderedprizes = orderedprizes, .orderedprizes_id = orderedprizes_id,
                          .maxoutprize = -FARAWAY, .maxoutprize_idx = -1};

   assert(g && cost && prize && pathdist && pathedge && connected);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->extended);
   assert(graph_pc_isPcMw(g) && !graph_pc_isRootedPcMw(g));

   /* initialize */
   stPcmwInit(g, pathdist, pathedge, connected, &ntermspos);

   assert(g->mark[start]);
   assert(ntermspos >= 0);

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int count = 1;
      int nterms = 0;
      const SCIP_Bool isPc = graph_pc_isPc(g);

      shortestpath_pcReset(&spaths_pc);

      if( Is_pseudoTerm(g->term[start]) )
      {
         nterms++;
         shortestpath_pcConnectNode(g, connected, start, &spaths_pc);
      }

      /* add start vertex to heap */
      heap[count] = start;
      state[start] = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         SCIP_Bool connectK = FALSE;
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect k to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) )
         {
            connectK = (prize[k] >= pathdist[k]);

            assert(k != start);

            if( !connectK )
            {
               SCIP_Real prizesum = 0.0;

               for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
               {
                  if( Is_pseudoTerm(g->term[node]) )
                  {
                     prizesum += prize[node];
                  }
                  else if( isPc && !costIsBiased && graph_pc_knotIsNonLeafTerm(g, node) )
                  {
                     prizesum += prize[node];
                  }
               }

               assert(prizesum >= 0.0 && LT(prizesum, FARAWAY));

               connectK = (prize[k] + prizesum >= pathdist[k]);
            }

            if( connectK )
            {
               stPcmwConnectNode(k, g, &spaths_pc, pathdist, pathedge, connected, &count, &nterms);

               assert(nterms <= ntermspos);

               /* have all biased terminals been connected? */
               if( nterms == ntermspos )
               {
                  SCIPdebugMessage("all terms reached \n");
                  break;
               }
            }
         }

         if( !connectK && pathdist[k] > spaths_pc.maxoutprize )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m] && e != EAT_LAST);

            /* is m not connected, allowed and closer? */
            if( g->mark[m] && !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) )
               correctX(heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      } /* while( count > 0 ) */
   } /* nnodes > 1*/
}

/** Reduce given solution
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw_reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            tmpnodeweight,      /**< node weight array */
   int*                  result,             /**< incoming arc array */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   assert(tmpnodeweight != NULL);
   assert(result   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      if( result[e] == CONNECT )
      {
         const int head = g->head[e];

         if( Is_term(g->term[head]) )
            continue;

         graph_path_st_pcmw_reduce(scip, g, cost, tmpnodeweight, result, head, connected);

         assert(connected[head]);

         if( SCIPisGE(scip, cost[e], tmpnodeweight[head]) )
         {
            connected[head] = FALSE;
            result[e] = UNKNOWN;
         }
         else
         {
            tmpnodeweight[start] += tmpnodeweight[head] - cost[e];
         }
      }
   }

   SCIPfreeBufferArray(scip, &tmpnodeweight);
}



/** Find a tree rooted in node 'start' and connecting
 *  all positive vertices.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw_full(
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   const int nnodes = g->knots;
   const int nterms = graph_pc_isRootedPcMw(g) ? g->terms : g->terms - 1;

   assert(start >= 0 && start < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->extended);

   if( graph_pc_isRootedPcMw(g) )
      stRpcmwInit(g, pathdist, pathedge, connected, NULL);
   else
      stPcmwInit(g, pathdist, pathedge, connected, NULL);

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int heapsize = 1;
      int termscount = 0;

      /* add start vertex to heap */
      heap[heapsize] = start;
      state[start] = heapsize;

      if( Is_term(g->term[start]) || Is_pseudoTerm(g->term[start]) )
         termscount++;

      /* repeat until heap is empty */
      while( heapsize > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &heapsize, pathdist);
         state[k] = UNKNOWN;

         /* if k is unconnected proper terminal, connect its path to current subtree */
         if( !connected[k] && (Is_term(g->term[k]) || Is_pseudoTerm(g->term[k])) )
         {
            connected[k] = TRUE;
            pathdist[k] = 0.0;

            assert(k != start);
            assert(pathedge[k] != -1);

            for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
            {
               connected[node] = TRUE;
               resetX(pathdist, heap, state, &heapsize, node, 0.0);

               assert(!Is_term(g->term[node]) && !Is_pseudoTerm(g->term[node]));
               assert(pathedge[node] != -1);
            }

            /* have all terminals been reached? */
            if( ++termscount == nterms )
            {
               break;
            }
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && g->mark[m] && GT(pathdist[m], (pathdist[k] + cost[e])) )
               correctX(heap, state, &heapsize, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }

#ifndef NDEBUG
   if( graph_pc_isRootedPcMw(g) )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( graph_pc_knotIsFixedTerm(g, k) )
            assert(connected[k]);
      }
   }
#endif
}


/** greedy extension of a given tree for PC or MW problems */
void graph_path_st_pcmw_extend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Bool             breakearly,         /**< finish computation early if no profitable extension possible? */
   PATH*                 path,               /**< shortest paths data structure */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            extensions          /**< extensions performed? */
   )
{
   SCIP_Real maxprize;
   int   count;
   int   nnodes;
   int   nstnodes;
   int   outerterms;
   int*  heap;
   int*  state;

   assert(path   != NULL);
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);
   assert(g->extended);

   maxprize = 0.0;
   count = 0;
   nstnodes = 0;
   nnodes = g->knots;
   heap = g->path_heap;
   state = g->path_state;
   *extensions = FALSE;
   outerterms = 0;

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      if( connected[k] && g->mark[k] )
      {
         /* add node to heap */
         nstnodes++;
         if( nnodes > 1 )
            heap[++count] = k;

         state[k]     = count;
         path[k].dist = 0.0;
         assert(path[k].edge != UNKNOWN || k == g->source);
      }
      else
      {
         state[k]     = UNKNOWN;
         path[k].dist = FARAWAY;

         if( Is_pseudoTerm(g->term[k]) && g->mark[k] )
         {
            outerterms++;
            if( g->prize[k] > maxprize )
            {
               maxprize = g->prize[k];
            }
         }
      }

      if( !connected[k] )
         path[k].edge = UNKNOWN;
   }

   /* nothing to extend? */
   if( nstnodes == 0 )
      return;

   if( nnodes > 1 )
   {
      int nterms = 0;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = pathheapGetNearest(heap, state, &count, path);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough (or fixnode), connect its path to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) && SCIPisGE(scip, g->prize[k], path[k].dist) )
         {
            int node;

            nterms++;
            *extensions = TRUE;
            connected[k] = TRUE;
            path[k].dist = 0.0;

            assert(path[k].edge >= 0);
            node = g->tail[path[k].edge];

            while( !connected[node] )
            {
               assert(path[node].edge != UNKNOWN);
               connected[node] = TRUE;
               heapReset(path, heap, state, &count, node);
               assert(state[node]);

               if( Is_pseudoTerm(g->term[node]) )
                  nterms++;

               node = g->tail[path[node].edge];
            }

            assert(path[node].dist == 0.0);
            assert(nterms <= outerterms);

            /* have all terminals been reached? */
            if( nterms == outerterms )
               break;
         }
         else if( breakearly && SCIPisGT(scip, path[k].dist, maxprize) )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m]);
            assert(e != EAT_LAST);

            /* is m not connected, allowed and closer (as close)? */

            if( !connected[m] && path[m].dist > (path[k].dist + cost[e]) && g->mark[m] )
               pathheapCorrect(heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
      }
   }
}


/** greedy extension of a given tree for PC or MW problems; path[i].edge needs to be initialized */
void graph_path_st_pcmw_extendBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   PATH*                 path,               /**< shortest paths data structure with .edge initialized */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            extensions          /**< extensions performed? */
   )
{
   SCIP_Real maxprize;
   int count;
   int nstnodes;
   int outermscount;
   const int nnodes = g->knots;
   int* const heap = g->path_heap;
   int* const state = g->path_state;

   assert(path   != NULL);
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);
   assert(g->extended);

   maxprize = 0.0;
   count = 0;
   nstnodes = 0;
   outermscount = 0;

   *extensions = FALSE;

   /* unmark dummy terminals */
   graph_pc_markOrgGraph(g);
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      path[k].dist = FARAWAY;

      if( !g->mark[k] )
         continue;

      if( connected[k] )
      {
         /* add node to heap */
         nstnodes++;
         if( nnodes > 1 )
            heap[++count] = k;

         state[k]     = count;
         path[k].dist = 0.0;
         assert(path[k].edge != UNKNOWN || k == g->source);
      }
      else if( Is_pseudoTerm(g->term[k]) )
      {
         assert(g->mark[k]);
         outermscount++;

         if( prize[k] > maxprize )
            maxprize = prize[k];

      }
   }

   /* with at least two nodes and at least one in the solution? */
   if( nnodes > 1 && nstnodes > 0 )
   {
      int nterms = 0;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = pathheapGetNearest(heap, state, &count, path);
         state[k] = UNKNOWN;

         assert(g->mark[k]);

         /* if k is positive vertex and close enough (or fixnode), connect its path to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) && SCIPisGE(scip, prize[k], path[k].dist) )
         {
            int node;

            nterms++;
            *extensions = TRUE;
            connected[k] = TRUE;
            path[k].dist = 0.0;

            assert(path[k].edge >= 0);
            node = g->tail[path[k].edge];

            while( !connected[node] )
            {
               assert(g->mark[node]);
               assert(path[node].edge >= 0);
               connected[node] = TRUE;
               heapReset(path, heap, state, &count, node);
               assert(state[node]);

               if( Is_pseudoTerm(g->term[node]) )
                  nterms++;

               node = g->tail[path[node].edge];
            }

            assert(path[k].dist == 0.0);

            assert(nterms <= outermscount);

            /* have all terminals been reached? */
            if( nterms == outermscount )
               break;
         }
         else if( path[k].dist > maxprize )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];
            assert(state[m]);

            if( connected[m] )
               continue;

            /* is m allowed and closer? */
            if( path[m].dist > (path[k].dist + cost[e]) && g->mark[m] )
               pathheapCorrect(heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
      }
   }
}


/** !!!LEGACY CODE!!!
 *  Shortest path heuristic for the RMWCSP and RPCSPG
 *  Find a directed tree rooted in node 'start' and connecting all terminals as well as all
 *  positive vertices (as long as this is profitable).
 *  */
void graph_path_st_rpcmw(
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            orderedprizes,      /**< legacy code */
   int*                  orderedprizes_id,   /**< legacy code */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   const int nnodes = g->knots;
   int nrterms = -1;
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   SPATHSPC spaths_pc = { .orderedprizes = orderedprizes, .orderedprizes_id = orderedprizes_id,
                          .maxoutprize = -FARAWAY, .maxoutprize_idx = -1};

   assert(g && cost && prize && pathdist && pathedge && connected);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->extended);
   assert(graph_pc_isRootedPcMw(g));

   stRpcmwInit(g, pathdist, pathedge, connected, &nrterms);

   assert(nrterms >= 1);
   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int count;
      const int nterms = g->terms;
      int termscount = 0;
      int rtermscount = 0;

      shortestpath_pcReset(&spaths_pc);

      /* add start vertex to heap */
      count = 1;
      heap[count] = start;
      state[start] = count;

      if( Is_anyTerm(g->term[start]) )
      {
         shortestpath_pcConnectNode(g, connected, start, &spaths_pc);

         termscount++;
      }

      if( Is_term(g->term[start]) )
      {
         assert(graph_pc_knotIsFixedTerm(g, start));
         rtermscount++;
      }

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is fixed terminal positive vertex and close enough, connect its path to current subtree */
         if( Is_anyTerm(g->term[k]) && (Is_term(g->term[k]) || GE(prize[k], pathdist[k]))
            && !connected[k] )
         {
            int node;

            assert(k != start);
            assert(pathedge[k] != -1);
            assert(!graph_pc_knotIsDummyTerm(g, k));
            assert(graph_pc_knotIsFixedTerm(g, k) || GE(prize[k], pathdist[k]));

            if( !graph_pc_knotIsNonLeafTerm(g, k) )
               termscount++;

            if( Is_term(g->term[k]) )
            {
               assert(graph_pc_knotIsFixedTerm(g, k));
               rtermscount++;
            }
            else if( Is_pseudoTerm(g->term[k]) )
            {
               shortestpath_pcConnectNode(g, connected, k, &spaths_pc);
            }

            connected[k] = TRUE;
            pathdist[k] = 0.0;

            node = k;

            while( !connected[node = g->tail[pathedge[node]]] )
            {
               assert(pathedge[node] != -1);
               assert(!Is_term(g->term[node]));
               assert(!graph_pc_knotIsFixedTerm(g, node));
               assert(g->mark[node]);

               connected[node] = TRUE;
               resetX(pathdist, heap, state, &count, node, 0.0);

               if( Is_pseudoTerm(g->term[node]) )
               {
                  termscount++;
                  shortestpath_pcConnectNode(g, connected, node, &spaths_pc);
               }
            }

            assert(termscount <= nterms);

            /* have all terminals been reached? */
            if( termscount == nterms )
            {
               SCIPdebugMessage("all terminals reached \n");
               break;
            }
         }
         else if( rtermscount >= nrterms && pathdist[k] > spaths_pc.maxoutprize )
         {
            SCIPdebugMessage("all fixed terminals reached \n");

            assert(rtermscount == nrterms);
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int head = g->head[e];

            assert(state[head]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[head] && g->mark[head] && pathdist[head] > (pathdist[k] + cost[e]) )
               correctX(heap, state, &count, pathdist, pathedge, head, k, e, cost[e]);
         }
      }
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(connected[k]);
      }
   }
#endif
}


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

/** gets non-terminal shortest paths to 4 closest terminal for each non-terminal
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
   SCIP_CALL( tpathsBuild(scip, g, *tpaths) );

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
   SCIP_CALL( tpathsBuildBiased(scip, sdprofit, g, *tpaths) );

   return SCIP_OKAY;
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
   const SCIP_Bool withProfit = (sdprofit != NULL);

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
               pathheapCorrectDist(heap, state, &nHeapElems, path, m, k, e, distnew);
               vbase[m] = vbase[k];
            }
         }
      }
   }
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

   while( nheapElems > 0 )
   {
      const int k2 = pathheapGetNearest(heap, state, &nheapElems, path);
      const int k = k2 - nnodes;
      assert(0 <= k && k < nnodes);

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
               pathheapCorrectDist(heap, state, &nheapElems, path, head2, k2, e, distnew);
               vbase[head2] = vbase[k2];
            }
         }
      }
   }


   return;
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
               pathheapCorrectDist(heap, state, &nheapElems, path, head3, k3, e, distnew);
               vbase[head3] = vbase[k3];
            }
         }
      }
   }

   return;
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
               pathheapCorrectDist(heap, state, &nHeapElems, path, head4, k4, e, distnew);
               vbase[head4] = vbase[k4];
            }
         }
      }
   }

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
            assert(LE(path3[level * nnodes + k].dist, path3[(level + 1) * nnodes + k].dist));
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

   if( !graph_pc_isPcMw(g) )
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

   SCIPfreeMemoryArray(scip, &(tp->state));
   SCIPfreeMemoryArray(scip, &(tp->termbases));
   SCIPfreeMemoryArray(scip, &(tp->termpaths));

   SCIPfreeMemory(scip, tpaths);
}


/** gets (up to) four close terminals to given node i */
void graph_tpathsGet4CloseTerms(
   const GRAPH*          g,                  /**< graph */
   const TPATHS*         tpaths,             /**< the terminal paths */
   int                   node,               /**< node */
   SCIP_Real             maxdist_strict,     /**< maximum valid distance (strict) */
   int*                  closeterms,         /**< four close terminals */
   SCIP_Real*            closeterms_dist,    /**< four close terminal distance */
   int*                  ncloseterms         /**< number of close terminals found */
)
{
   const PATH* const termpaths = tpaths->termpaths;
   const int* const termbases = tpaths->termbases;
   const int nnodes = tpaths->nnodes;
   int pos = node;
   int nnterms = 0;

   assert(closeterms && closeterms_dist && ncloseterms && g);
   assert(4 <= STP_TPATHS_NTERMBASES);
   assert(0 <= node && node < nnodes);
   assert(nnodes == g->knots);
   assert(GE(maxdist_strict, 0.0));

   if( Is_term(g->term[node]))
   {
      *ncloseterms = 1;
      closeterms[0] = node;
      closeterms_dist[0] = 0.0;
      return;
   }

   for( int k = 0; k < 4; k++ )
   {
      if( LT(termpaths[pos].dist, maxdist_strict) )
      {
         closeterms[nnterms] = termbases[pos];
         closeterms_dist[nnterms++] = termpaths[pos].dist;
      }
      else
      {
         break;
      }

      pos += nnodes;
   }

   *ncloseterms = nnterms;
}
