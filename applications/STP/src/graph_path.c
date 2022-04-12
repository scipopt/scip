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

/**@file   graph_path.c
 * @brief  Shortest path based graph algorithms for Steiner problems
 * @author Thorsten Koch
 * @author Henriette Franz
 * @author Daniel Rehfeldt
 *
 * This file encompasses various (heap-based) shortest path based algorithms including
 * Dijkstra's algorithm.
 *
 */
//#define SCIP_DEBUG
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"
#include "graphheaps.h"
#include "shortestpath.h"
#include "solstp.h"



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

/** resets node */
inline static
void pathheapReset(
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


/** For DCSTP can terminal be connected to current tree?  */
static inline
SCIP_Bool stDcTermIsConnectable(
   const GRAPH*          g,                  /**< graph data structure */
   int                   term,               /**< terminal to be checked */
   const int*            deg_free,           /**< free degree of each node */
   const int*            pathedge,           /**< predecessor edge array (on vertices) */
   const STP_Bool*       connected           /**< array to mark whether a vertex is part of computed Steiner tree */
)
{
   int node;

   assert(!connected[term]);
   assert(g->stp_type == STP_DCSTP);
   assert(deg_free[term] >= 0);

   if( deg_free[term] == 0 )
   {
      SCIPdebugMessage("failed to add %d-path because free degree 0 end \n", term);
      return FALSE;
   }

   /* NOTE: intermediary path vertices need to be of degree at least 2 */
   for( node = g->tail[pathedge[term]]; !connected[node]; node = g->tail[pathedge[node]] )
   {
      assert(pathedge[node] != -1);
      assert(!Is_term(g->term[node]));
      assert(deg_free[node] >= 0);

      if( deg_free[node] <= 1 )
      {
         SCIPdebugMessage("failed to add %d-path due to node %d (free degree: %d) \n", term, node, deg_free[node]);
         return FALSE;
      }
   }

   assert(connected[node]);

   if( deg_free[node] == 0 )
   {
      SCIPdebugMessage("failed to add %d-path because free degree 0 start \n", node);
      return FALSE;
   }

   return TRUE;
}



/** For DCSTP: extends (implicitly) given tree */
static
void sdDcExtendTree(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edgecosts */
   int*                  heapsize,           /**< heap size */
   int*                  pathdeg_free,       /**< currently available degree of path */
   int*                  nodedeg_free,       /**< currently available degree of node */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  result,             /**< solution array */
   int*                  soldegfree,         /**< per node: solution degree */
   int*                  nsolterms           /**< number of solution terminals */
   )
{
   int* RESTRICT heap = g->path_heap;
   int* RESTRICT state = g->path_state;

   assert(*heapsize > 0);
   assert(*soldegfree > 0);

   /* repeat until heap is empty */
   while( *heapsize > 0 )
   {
      const int k = nearestX(heap, state, heapsize, pathdist);
      state[k] = UNKNOWN;

      /* k is terminal and not connected yet? */
      if( Is_term(g->term[k]) && !connected[k] )
      {
         /* no degree problem? */
         if( stDcTermIsConnectable(g, k, nodedeg_free, pathedge, connected) )
         {
            int node;
            assert(pathedge[k] >= 0);
            assert(nodedeg_free[k] >= 1);

            nodedeg_free[k] -= 1;
            result[pathedge[k]] = CONNECT;
            connected[k] = TRUE;
            pathdist[k] = 0.0;
            pathdeg_free[k] = 0;
            SCIPdebugMessage("connecting terminal %d \n", k);

            /* connect k to current solution */
            for( node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
            {
               assert(pathedge[node] != -1);
               assert(!Is_term(g->term[node]));
               assert(nodedeg_free[node] >= 2);

               SCIPdebugMessage("...adding path node %d \n", node);

               pathdeg_free[node] = 0;
               nodedeg_free[node] -= 2;
               connected[node] = TRUE;
               result[pathedge[node]] = CONNECT;
               resetX(pathdist, heap, state, heapsize, node, 0.0);
            }

            assert(nodedeg_free[node] >= 1);
            nodedeg_free[node] -= 1;

            /* have all terminals been reached? */
            if( ++(*nsolterms) == g->terms )
               break;
         }
      }

      assert((*nsolterms) < g->terms);

      if( !connected[k]  )
      {
         /* NOTE: if k is not in tree yet, adding it as a non-leaf reduces its free degree by two */
         if( nodedeg_free[k] <= 1 )
            continue;

         /* NOTE: if path to terminal could not be connect now, it can also not be used as sub-path later */
         if( Is_term(g->term[k]) )
            continue;
      }
      else if( nodedeg_free[k] == 0 )
      {
         continue;
      }

      /* update adjacent vertices */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( connected[m] )
            continue;

         assert(state[m]);

         if( GT(pathdist[m], (pathdist[k] + cost[e]))
            || (EQ(pathdist[m], (pathdist[k] + cost[e])) && (pathdeg_free[k] + nodedeg_free[m] - 2 > pathdeg_free[m]) )
         )
         {
            assert(nodedeg_free[m] > 0);
            correctX(heap, state, heapsize, pathdist, pathedge, m, k, e, cost[e]);

            pathdeg_free[m] = pathdeg_free[k] + nodedeg_free[m] - 2;
         }
      }
   }
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
   int* count              /* pointer to store the number of elements on the heap */
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


/** existing? */
SCIP_Bool graph_path_exists(
   const GRAPH*          g                   /**< graph data structure */
   )
{
   assert(g);
   assert((g->path_heap != NULL) == (g->path_state != NULL));

   return (g->path_heap != NULL) ;
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
   assert(g);
#ifndef WITH_UG
   assert(g->path_heap && g->path_state);
#endif

   SCIPfreeMemoryArrayNull(scip, &(g->path_state));
   SCIPfreeMemoryArrayNull(scip, &(g->path_heap));
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
   const GRAPH*          g,                  /**< graph data structure */
   int                   mode,               /**< shortest path (FSP_MODE) or minimum spanning tree (MST_MODE)? */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path                /**< shortest paths data structure */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT heap = g->path_heap;
   int* RESTRICT state = g->path_state;

   assert(g != NULL);
   assert(start >= 0);
   assert(start < g->knots);
   assert((mode == FSP_MODE) || (mode == MST_MODE));
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);
   assert(path != NULL);
   assert(cost != NULL);
   assert(g->mark[start]);

   /* initialize */
   for( int i = 0; i < nnodes; i++ )
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY + 1.0;
      path[i].edge = UNKNOWN;
   }

   path[start].dist = 0.0;

   if( nnodes > 1 )
   {
      /* add first node to heap */
      int nheapElems = 1;
      heap[nheapElems] = start;
      state[start] = nheapElems;

      while( nheapElems > 0 )
      {
         /* get nearest labeled node */
         const int k = pathheapGetNearest(heap, state, &nheapElems, path);

         /* mark as scanned */
         state[k] = CONNECT;

         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];
            assert(e != EAT_LAST);

            /* node not scanned and valid? */
            if( state[m] )
            {
               /* closer than previously and valid? */
               if( path[m].dist > ((mode == MST_MODE) ? cost[e] : (path[k].dist + cost[e])) && g->mark[m] )
               {
                  pathheapCorrect(heap, state, &nheapElems, path, m, k, e, cost[e], mode);
               }
            }
         }
      }
   }
}


/** limited Dijkstra on incoming edges */
void graph_pathInLimitedExec(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      edges_cost,         /**< edge cost */
   const SCIP_Bool*      nodes_abort,        /**< nodes to abort at */
   int                   startnode,          /**< start */
   DIJK*                 dijkdata,           /**< Dijkstra data */
   SCIP_Real*            abortdistance       /**< distance at which abort happened */
)
{
   SCIP_Real* RESTRICT nodes_dist = dijkdata->node_distance;
   int* RESTRICT visitlist = dijkdata->visitlist;
   STP_Bool* RESTRICT visited = dijkdata->node_visited;
   DHEAP* dheap = dijkdata->dheap;
   int* const state = dheap->position;
   int nvisits = 0;

   assert(nodes_dist && visitlist && visited && dheap && edges_cost && nodes_abort && abortdistance);
   assert(dheap->size == 0);
   assert(graph_knot_isInRange(g, startnode));

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(nodes_dist[k] == FARAWAY);
      assert(state[k] == UNKNOWN);
   }
#endif

   *abortdistance = -FARAWAY;
   nodes_dist[startnode] = 0.0;
   visitlist[(nvisits)++] = startnode;
   visited[startnode] = TRUE;
   graph_heap_correct(startnode, 0.0, dheap);

   while( dheap->size > 0 )
   {
      const int k = graph_heap_deleteMinReturnNode(dheap);
      assert(state[k] == CONNECT);

      if( nodes_abort[k] )
      {
         *abortdistance = nodes_dist[k];
         break;
      }

      for( int e = g->inpbeg[k]; e >= 0; e = g->ieat[e] )
      {
         const int m = g->tail[e];

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = nodes_dist[k] + edges_cost[e];

            if( !visited[m] )
            {
               assert(nvisits < g->knots);
               visitlist[(nvisits)++] = m;
               visited[m] = TRUE;
            }

            if( LT(distnew, nodes_dist[m]) )
            {
               nodes_dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   dijkdata->nvisits = nvisits;
}


/** limited Dijkstra, stopping at terminals */
void graph_sdPaths(
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   const SCIP_Real*      cost,               /**< edge costs */
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
   const SCIP_Bool isDirected = !graph_typeIsUndirected(g);

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
         if( isDirected && GE(cost[e], FARAWAY) )
            continue;

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

      if( !isDirected )
      {
         if( Is_term(g->term[k]) || k == head )
            continue;
      }

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( state[m] && g->mark[m] && GE(distlimit, cost[e]) && (GT(path[m].dist, path[k].dist + cost[e])) )
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




/** For DCSTP: Find a directed tree rooted in node 'start' and spanning all terminals, while respecting degree constraints */
SCIP_RETCODE graph_path_st_dc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edgecosts */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  result,             /**< solution */
   STP_Bool*             solFound            /**< pointer to store whether solution was found */
   )
{
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT heap = g->path_heap;
   int* RESTRICT state = g->path_state;
   int* nodedeg_free;
   int* pathdeg_free;

   assert(cost && pathdist && pathedge && connected && result && solFound);
   assert(heap && state);
   assert(graph_knot_isInRange(g, start));
   assert(g->stp_type == STP_DCSTP);
   assert(g->maxdeg);

   *solFound = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodedeg_free, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdeg_free, nnodes) );

   BMScopyMemoryArray(nodedeg_free, g->maxdeg, nnodes);

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   for( int k = 0; k < nnodes; k++ )
   {
      pathdeg_free[k] = -1;
      state[k] = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int heapsize;
      int nsolterms = 0;
      int nsolterms_old;
      int soldegfree;

      if( Is_term(g->term[start]) )
         nsolterms++;

#ifdef SCIP_DEBUG
      printf("TM start=%d \n", start);
      graph_knot_printInfo(g, start);
#endif

      soldegfree = nodedeg_free[start];

      /* add start vertex to heap */
      heapsize = 1;
      heap[heapsize] = start;
      state[start] = heapsize;
      pathdeg_free[start] = 0;

      /* repeat as long as the tree can be extended */
      do
      {
         SCIPdebugMessage("extension loop with nsolterms=%d  \n", nsolterms);


         nsolterms_old = nsolterms;
         sdDcExtendTree(g, cost, &heapsize, pathdeg_free, nodedeg_free, pathdist, pathedge, connected, result, &soldegfree, &nsolterms);

         assert(nsolterms <= g->terms);
         assert(nsolterms >= nsolterms_old);

         /* will we continue? */
         if( nsolterms_old != nsolterms && nsolterms != g->terms && soldegfree > 0 )
         {
            heapsize = 0;

            for( int k = 0; k < nnodes; k++ )
            {
               if( connected[k] )
               {
                  heap[++heapsize] = k;
                  state[k] = heapsize;
                  assert(pathdeg_free[k] == 0);
               }
               else
               {
                  pathdeg_free[k] = 0;
                  state[k] = UNKNOWN;
                  pathdist[k] = FARAWAY;
                  pathedge[k] = -1;
               }
            }

            assert(heapsize > 0);
         }

         /* finished? */
         if( nsolterms == g->terms )
            break;

      } while( nsolterms_old != nsolterms  && soldegfree > 0 );


      *solFound = (g->terms == nsolterms);

      SCIPdebugMessage("g->terms == nterms? %d==%d \n", g->terms, nsolterms);
   }

   SCIPfreeBufferArray(scip, &pathdeg_free);
   SCIPfreeBufferArray(scip, &nodedeg_free);

   return SCIP_OKAY;
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
               pathheapReset(path, heap, state, &count, node);
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
               pathheapReset(path, heap, state, &count, node);
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



/**
 * todo refactor, was copied from Henriette's branch
 */
SCIP_RETCODE graph_path_st_brmwcs(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            solfound            /**< could a solution be found? */
   )
{

   /* node utilities */
   const SCIP_Real* nodeweight = prize;
   /* node costs */
   const SCIP_Real* nodebudget = g->costbudget;
   /* given budget */
   const SCIP_Real budget = g->budget;
   /* number of nodes */
   const int nnodes = g->knots;

   /* number of terminals which are not part of the solution yet */
   int numberTerminals;
   /* current available budget */
   SCIP_Real newBudget;
   /* positive sum of all negative utilities */
   double smallestWeight;
   /* alpha (in [0,1]) is used by the key of the heap */
   double alpha;

   /* upper bound for node costs (pi) and the corresponding node utilities (my) */
   double* pi;
   double* my;

   int lookNeighbour;

   int* alreadyContained;

  // printf("\n\n start = %d \n", start+1);

   *solfound = TRUE;

   graph_pc_markOrgGraph(g);

   assert( pathdist != NULL );
   assert( pathedge != NULL );
   assert( g != NULL );
   assert( start >= 0 );
   assert( start < g->knots );
   assert( prize != NULL );
   assert( connected != NULL );
   assert( g->mark[start] );
   assert( Is_term( g->term[start] ) );

   SCIP_CALL( SCIPallocBufferArray( scip, &alreadyContained, nnodes ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &my, nnodes ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &pi, nnodes ) );

   numberTerminals = g->terms;
   newBudget = budget - nodebudget[start];
   smallestWeight = 1.0;

   for( int k = 0; k < nnodes; k++ )
   {
      /* node costs have to be positive */
      if( (nodebudget[k] < 0.0) && g->mark[k] )
      {
         printf("all costs have to be non-negativ! \n");
         SCIPfreeBufferArray(scip, &pi);
         SCIPfreeBufferArray(scip, &my);
         SCIPfreeBufferArray(scip, &alreadyContained);
         return SCIP_ERROR;
      }

      if( (nodeweight[k] < 0) && g->mark[k] )
         smallestWeight = smallestWeight + ((-1.0) * nodeweight[k]);
   }

   /* add start vertex to the solution */
   connected[start] = TRUE;

   /* at the beginning set alpha equal one */
   alpha = 1.0;

  /* first find a solution containing all terminals
   * second expand the solution by adding additional nodes until no further budget is available */
   while( ( (alpha >= 0.0) && (numberTerminals > 0) ) || ( (alpha >= 0.0) && (newBudget >= 0.0) ) )
   {

      /* construct a heap */
      int* heap = g->path_heap;
      int* state = g->path_state;
      int count = 0;

      /* calculate number of terminals again if alpha + 0.1 does not work for the key
       * and no feasible solution is found yet */
      if( numberTerminals > 0 )
            numberTerminals = g->terms;

      /* initialize heap and define node start as start solution
       * if no feasible solution is found yet */
      for( int k = 0; k < nnodes; k++ )
      {
         if( (k != start) && (numberTerminals > 0) )
            connected[k] = FALSE;

         alreadyContained[k] = FALSE;

         if( !g->mark[k] && (numberTerminals > 0) )
            numberTerminals--;

         state[k] = UNKNOWN;
         pathedge[k] = -1;

         if( connected[k] && (numberTerminals <= 0) )
         {
            pathdist[k] = 0.0;
         }else
         {
            pathdist[k] = FARAWAY;
         }
      }

      /* reduce number of searched terminals if start is a terminal and
       * no feasible solution is found yet */
      if( g->mark[start] && Is_term( g->term[start] ) && (numberTerminals > 0) )
            numberTerminals--;

      /* initialize pi and my and add nodes contained by the solution to the heap;
       * if a feasible solution is already found, then add all nodes to the heap */
      for( int k = 0; k < nnodes; k++ )
      {
         if( g->mark[k] )
         {
            if( connected[k] )
            {
               pi[k] = 0.0;
               my[k] = 0.0;

               resetX(pathdist, heap, state, &count, k, 0.0 );
            }else
            {
               if( Is_term( g->term[k] ) && (numberTerminals > 0) ){
                  my[k] = 0.0;
               }else
               {
                  my[k] = nodeweight[k];
               }
               pi[k] = DBL_MAX;
               if( numberTerminals <= 0 ){
                  resetX(pathdist, heap, state, &count, k, pathdist[k] );
               }
            }
         }
      }

      /* return if the heap is empty */
      if( count == 0 )
      {
         SCIPfreeBufferArray( scip, &pi );
         SCIPfreeBufferArray( scip, &my );
         SCIPfreeBufferArray( scip, &alreadyContained );
         return SCIP_OKAY;
      }

           /* TRUE if we want to consider the neighbor nodes of node k */
      lookNeighbour = TRUE;

      /* repeat until the heap is empty */
      while( count > 0 )
      {
         /* get first node  */
         int k = nearestX( heap, state, &count, pathdist );
         alreadyContained[k] = TRUE;

         /* TRUE if we want to consider the neighbor nodes of node k */
         lookNeighbour = TRUE;

         /*if( numberTerminals <= 0 )
            state[k] = UNKNOWN;*/

         /* if k is a terminal node that is  not contained in the current solution,
          * add k and its path to the solution */
         if( ( numberTerminals > 0 ) && !connected[k] && Is_term( g->term[k] ) && g->mark[k] )
         { //numberTerminals>0 braucht man nicht eigentlich

            lookNeighbour = FALSE;

            /* if the costs of the path are greater than the current budget,
            * then no feasible solution exists for alpha, decrease alpha by 0.1 */
            if( ( newBudget - pi[k] ) < 0.0 )
            {
               alpha = alpha - 0.1;
               count = 0;

               /* if alpha smaller than zero, no feasible solution could be found */
               if( alpha < 0.0 )
               {
          //        printf("No feasible solution found! \n" );

                  *solfound = FALSE;

                  SCIPfreeBufferArray( scip, &pi );
                  SCIPfreeBufferArray( scip, &my );
                  SCIPfreeBufferArray( scip, &alreadyContained );

                  return SCIP_OKAY;
               }
            /* otherwise, add k and its path to the solution and update the current budget */
            }else
            {
               int v = k;

               newBudget = newBudget - pi[k];

               while( !connected[v] )
               {
                  connected[v] = TRUE;
                  pi[v] = 0.0;
                  my[v] = 0.0;
                  if( pathedge[v] != -1 )
                     v = g->tail[pathedge[v]];
               }


               /* update the number of terminals */
               numberTerminals = numberTerminals - 1;

               /* neu */
               heap = g->path_heap;
               state = g->path_state;
               count = 0;

               /* initialize heap and define node start as start solution
                * if no feasible solution is found yet */
               for( int n = 0; n < nnodes; n++ )
               {
                  state[n] = UNKNOWN;
                  pathedge[n] = -1;
                  alreadyContained[n] = FALSE;

                  if( connected[n] )
                  {
                     pathdist[n] = 0.0;
                     resetX(pathdist, heap, state, &count, n, 0.0);
                  }
                  else
                  {
                     pathdist[n] = FARAWAY;
                     if( Is_term(g->term[n]) )
                     {
                        my[n] = 0.0;
                     }
                     else
                     {
                        my[n] = nodeweight[n];
                     }
                     pi[n] = DBL_MAX;
                  }
               } /* bis hier ist der heap jetzt leer? */

               /* break the while-loop if all terminals are contained in the solution */
               if( numberTerminals == 0 )
                  count = 0;
            }
         }
         /* iterate over incident edges */
         for( int e = g->outbeg[k]; lookNeighbour && e != EAT_LAST; e = g->oeat[e]  )
         {
            /* get neighbor head of node k */
              const int head = g->head[e];
            assert( state[head] );

            /* calculate new key value for head if neighbor head is not contained in the solution
             * and no feasible solution is found yet or neighbor head is still part of the heap
             */
            //if( g->mark[head] && !connected[head] && ( ( numberTerminals > 0 ) || ( state[head] != UNKNOWN ) ) )
//if( g->mark[head] && !connected[head] && ( ( numberTerminals > 0 ) || !alreadyContained[head] ) )
            if( g->mark[head] && !connected[head] && ( ( numberTerminals > 0 && pathedge[head]==-1 ) || !alreadyContained[head] ) )
            {
               int j = k;
               int isConnected = TRUE;

               double possibleUtility;
               if( Is_term( g->term[head] ) && (numberTerminals > 0) )
               {
                  if( my[k] >= 0.0 )
                  {
                     possibleUtility = ( 1.0 - alpha ) * ( pi[k] + nodebudget[head] ) + ( alpha * ( pi[k] + nodebudget[head] ) /( my[k] + 1.0 ) );
                  }else
                  {
                     possibleUtility = ( 1.0 - alpha ) * ( pi[k] + nodebudget[head] ) + ( alpha * ( pi[k] + nodebudget[head] ) / ( ( my[k] + smallestWeight ) / smallestWeight ) );
                  }
               }else
               {
                  if( ( my[k] + nodeweight[head] ) >= 0.0 )
                  {
                     possibleUtility = ( 1.0 - alpha ) * ( pi[k] + nodebudget[head] )
                               + ( alpha * ( pi[k] + nodebudget[head] ) / ( my[k] + nodeweight[head] + 1.0 ) );
                  }else
                  {
                     possibleUtility = ( 1.0 - alpha ) * ( pi[k] + nodebudget[head] )
                               + ( alpha * ( pi[k] + nodebudget[head] ) / ( ( my[k] + nodeweight[head] + smallestWeight ) / smallestWeight ) );
                  }
               }

               while( !connected[j] && isConnected && ( pathedge[j] != -1 ) )
               {
                  if( j == head )
                     isConnected = FALSE;

                  if( pathedge[j] != -1 )
                     j = g->tail[pathedge[j]];
               }

               /* if this the new key value is smaller than the current key value,
                * then update the pu, my, and the heap */
               if( ( pathdist[head] > possibleUtility ) && isConnected )
               {
                  pi[head] = pi[k] + nodebudget[head];



                  if( Is_term(g->term[head] ) )
                  {
                     my[head] = my[k];
                  }else
                  {
                     my[head] = my[k] + nodeweight[head];
                  }
                  correctX( heap, state, &count, pathdist, pathedge, head, k, e, (-1.0) * pathdist[k] + possibleUtility );
               }
            }
         }
      }
      /* if a feasible solution is already found, then expand the solution by a suitable node and its path
       * If such a node cannot be found, then decrease alpha by 0.1 */
      if( numberTerminals <= 0 && lookNeighbour ) //lookNeighbor neu
      {
         int foundVertex = FALSE;

         /* add all nodes to the heap which are not contained in the solution */
         for( int n = 0; n < nnodes; n++ )
         {
            if( !connected[n] && g->mark[n] )
               resetX(pathdist, heap, state, &count, n, pathdist[n]);
         }
         /* find a suitable node */

         while( count > 0 && !foundVertex )
         {
           /* get first node n from the heap */
            int n = nearestX( heap, state, &count, pathdist );

            /* if the cost of node n are smaller than the current budget,
             * then add the node n and its path to the solution and update the
             * available budget */
            if( g->mark[n] && ( newBudget - pi[n] >= 0.0 ) && ( my[n] > 0.0 ) )
            {
               int v = n;

//printf("budget=%f\n", newBudget);
 //printf("pi[%d]=%f\n",n, pi[n]);
               newBudget = newBudget - pi[n];
               foundVertex = TRUE;

               while( !connected[v] )
               {
//printf("nodebudget[%d]=%f \n", v, nodebudget[v]);
                  connected[v] = TRUE;
                  if( pathedge[v] != -1 )
                     v = g->tail[ pathedge[v] ];
               }
            }
         }
         /* decrease alpha by 0.1 if such a node cannot be found */
         if( !foundVertex )
            alpha = alpha - 0.1;
      }
   }
   SCIPfreeBufferArray( scip, &pi );
   SCIPfreeBufferArray( scip, &my );
   SCIPfreeBufferArray( scip, &alreadyContained );

   return SCIP_OKAY;
}
