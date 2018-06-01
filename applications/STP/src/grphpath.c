/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphpath.c
 * @brief  Shortest path based graph algorithms for Steiner problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file encompasses various (heap-based) shortest path based algorithms including
 * Dijkstra's algorithm and Voronoi diagram algorithms
 *
 * The underlying heap routines can be found in Jon Bentley, Programming Pearls, Addison-Wesley 1989
 *
 * The heap array is initialized with n elements (nodes), but only at most n-1 nodes can be included
 * in the array, since element 0 is not used for storing.
 */
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "grph.h"

/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearest(
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
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
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearestX(
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   const SCIP_Real* pathdist)
{
   int   k;
   int   t;
   int   c;
   int   j;
   int   dcount;

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   dcount = *count;

   if (dcount > 2)
      if (LT(pathdist[heap[3]], pathdist[heap[2]]))
         c++;

   while((c <= dcount) && GT(pathdist[heap[j]], pathdist[heap[c]]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= dcount)
         if (LT(pathdist[heap[c + 1]], pathdist[heap[c]]))
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
inline static void correct(
   SCIP* scip,
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   PATH*  path,
   int    l,
   int    k,
   int    e,
   SCIP_Real cost,
   int    mode)
{
   int    t;
   int    c;
   int    j;

   path[l].dist = (mode == MST_MODE) ? cost : (path[k].dist + cost);
   path[l].edge = e;

   /* new node? */
   if( state[l] == UNKNOWN )
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up */
   j = state[l];
   c = j / 2;
   while( (j > 1) && SCIPisGT(scip, path[heap[c]].dist, path[heap[j]].dist) )
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
inline static void correctX(
   SCIP* scip,
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   SCIP_Real*  pathdist,
   int*   pathedge,
   int    l,
   int    k,
   int    e,
   SCIP_Real cost
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[l] = (pathdist[k] + cost);
   pathedge[l] = e;

   /* Ist der Knoten noch ganz frisch ?
    */
   if (state[l] == UNKNOWN)
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up
    */
   j = state[l];
   c = j / 2;

   while( (j > 1) && SCIPisGT(scip, pathdist[heap[c]], pathdist[heap[j]]) )
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

void heap_add(
   int* heap,     /* heaparray */
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   int    node,    /* the node to be added */
   PATH*  path
   )
{
   int    t;
   int    c;
   int    j;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* Heap shift up */
   j = state[node];
   c = j / 2;

   while((j > 1) && GT(path[heap[c]].dist, path[heap[j]].dist))
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
inline static void resetX(
   SCIP* scip,
   SCIP_Real*  pathdist,
   int* heap,
   int* state,
   int* count,
   int    node
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[node] = 0.0;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* heap shift up */
   j = state[node];
   c = j / 2;

   while( (j > 1) && SCIPisGT(scip, pathdist[heap[c]], pathdist[heap[j]]) )
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

inline static void reset(
   SCIP* scip,
   PATH* path,
   int* heap,
   int* state,
   int* count,
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

   while( (j > 1) && SCIPisGT(scip, path[heap[c]].dist, path[heap[j]].dist) )
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

inline static void utdist(
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
         k = nearest(heap, state, &count, path);

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
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], mode);
            }
         }
      }
   }
}

/** limited Dijkstra, stopping at terminals */
void graph_sdPaths(
   SCIP*                 scip,               /**< SCIP data structure */
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

   if( g->stp_type != STP_MWCSP )
      g->mark[head] = FALSE;

   for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && (distlimit >= cost[e]) )
      {
         assert(SCIPisGT(scip, path[m].dist, path[tail].dist + cost[e]));

         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         correct(scip, heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }
   g->mark[head] = TRUE;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      const int k = nearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* distance limit reached? */
      if( SCIPisGT(scip, path[k].dist, distlimit) )
         break;

      /* stop at terminals */
      if( Is_term(g->term[k]) || k == head )
         continue;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];
         if( state[m] && g->mark[m] && (distlimit >= cost[e]) && (path[m].dist > path[k].dist + cost[e]) )
         {
            /* m labelled for the first time? */
            if( state[m] == UNKNOWN )
               memlbl[(*nlbl)++] = m;
            correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
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
         correct(scip, heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }

   g->mark[head] = TRUE;

   /* main loop */
   while( count > 0 )
   {
      const int k = nearest(heap, state, &count, path);
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
            correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
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
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
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
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}

/** Find a directed tree rooted in node 'start' and spanning all terminals */
void graph_path_st(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edgecosts */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int   k;
   int   m;
   int   e;
   int   count;
   int   nnodes;
   int*  heap;
   int*  state;

   assert(randnumgen != NULL);
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
   for( k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   /* add start vertex to heap */
   k            = start;
   pathdist[k] = 0.0;
   connected[k] = TRUE;

   if( nnodes > 1 )
   {
      int node;
      int nterms = 0;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is terminal, connect its path to current subtree */
         if( Is_term(g->term[k]) )
         {
            assert(k == start || !connected[k]);
            connected[k] = TRUE;
            pathdist[k] = 0.0;
            node = k;

            if( k != start )
            {
               assert(pathedge[k] != - 1);

               while( !connected[node = g->tail[pathedge[node]]] )
               {
                  assert(pathedge[node] != - 1);
                  connected[node] = TRUE;
                  resetX(scip, pathdist, heap, state, &count, node);
               }
            }
            /* have all terminals been reached? */
            if( ++nterms == g->terms )
               break;
         }

         /* update adjacent vertices */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) && g->mark[m] )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}


/** For rooted prize-collecting problem find a tree rooted in node 'start' and connecting
 *  positive vertices as long as this is profitable.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_rpc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   SCIP_Real maxprize;
   int   k;
   int   m;
   int   e;
   int   root;
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
   assert(g->stp_type == STP_RPCSPG);

   root = g->source;
   heap = g->path_heap;
   count = 0;
   state = g->path_state;
   nnodes = g->knots;
   maxprize = 0.0;

   /* initialize */
   for( k = nnodes - 1; k >= 0; --k )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( Is_pterm(g->term[k]) && SCIPisGT(scip, g->prize[k], maxprize) && k != start )
      {
         maxprize = g->prize[k];
      }
   }

   g->mark[root] = TRUE;

   /* add start vertex to heap */
   k            = start;
   pathdist[k] = 0.0;
   connected[k] = TRUE;

   if( nnodes > 1 )
   {
      int node;
      int nterms = 0;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect its path to current subtree */
         if( ((Is_pterm(g->term[k]) && SCIPisGE(scip, g->prize[k], pathdist[k])) || k == root) && !connected[k] )
         {
            connected[k] = TRUE;
            pathdist[k] = 0.0;
            node = k;

            if( k != start )
            {
               assert(pathedge[k] != - 1);

               while( !connected[node = g->tail[pathedge[node]]] )
               {
                  assert(pathedge[node] != - 1);
                  connected[node] = TRUE;
                  resetX(scip, pathdist, heap, state, &count, node);
               }
            }
            /* have all terminals been reached? */
            if( ++nterms == g->terms - 1 )
               break;
         }
         else if( SCIPisGT(scip, pathdist[k], maxprize) && connected[root] )
         {
            break;
         }

         /* update adjacent vertices */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) && g->mark[m] )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}


/** Find a tree rooted in node 'start' and connecting
 *  positive vertices as long as this is profitable.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   SCIP_Real maxprize;
   int k;
   int count;
   const int nnodes = g->knots;
   int* const heap = g->path_heap;
   int* const state = g->path_state;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   maxprize = 0.0;
   count = 0;

   /* initialize */
   for( k = nnodes - 1; k >= 0; --k )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( Is_pterm(g->term[k]) && (g->prize[k] > maxprize) && k != start )
         maxprize = g->prize[k];
   }

   /* add start vertex to heap */
   k            = start;
   pathdist[k] = 0.0;
   connected[k] = TRUE;

   if( nnodes > 1 )
   {
      int nterms = 0;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect its path to current subtree */
         if( !connected[k] && Is_pterm(g->term[k]) && (g->prize[k] >= pathdist[k]) )
         {
            int node = k;

            connected[k] = TRUE;
            pathdist[k] = 0.0;

            if( k != start )
            {
               assert(pathedge[k] != - 1);

               while( !connected[node = g->tail[pathedge[node]]] )
               {
                  assert(pathedge[node] != - 1);
                  connected[node] = TRUE;
                  resetX(scip, pathdist, heap, state, &count, node);
                  assert(state[node]);
               }
            }
            /* have all terminals been reached? */
            if( ++nterms == g->terms - 1 )
               break;
         }
         else if( SCIPisGT(scip, pathdist[k], maxprize) )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m] && e != EAT_LAST);
            /* is m not connected, allowed and closer (as close)? */

            if( !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) && g->mark[m] )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}

#if 1

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
#endif


/** Find a tree rooted in node 'start' and connecting
 *  all positive vertices.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw_full(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   int k;
   const int nnodes = g->knots;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   /* initialize */
   for( k = nnodes - 1; k >= 0; --k )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   /* add start vertex to heap */
   k            = start;
   pathdist[k] = 0.0;
   connected[k] = TRUE;

   if( nnodes > 1 )
   {
      int count = 1;
      int nterms = 0;

      heap[count] = k;
      state[k]    = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is positive vertex, connect its path to current subtree */
         if( !connected[k] && Is_pterm(g->term[k]) )
         {
            connected[k] = TRUE;
            pathdist[k] = 0.0;

            if( k != start )
            {
               int node = k;
               assert(pathedge[k] != - 1);

               while( !connected[node = g->tail[pathedge[node]]] )
               {
                  assert(pathedge[node] != - 1);
                  connected[node] = TRUE;
                  resetX(scip, pathdist, heap, state, &count, node);
                  assert(state[node]);
               }
            }
            /* have all terminals been reached? */
            if( ++nterms == g->terms - 1 )
               break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            int m = g->head[e];

            assert(state[m]);
            /* is m not connected, allowed and closer (as close)? */

            if( !connected[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[e])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}


/** greedy extension of a given tree for PC or MW problems */
void graph_path_st_pcmw_extend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path,               /**< shortest paths data structure */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            extensions          /**< extensions performed? */
   )
{
   SCIP_Real maxprize;
   int   k;
   int   count;
   int   nnodes;
   int   nstnodes;
   int*  heap;
   int*  state;

   assert(path   != NULL);
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);

   maxprize = 0.0;
   count = 0;
   nstnodes = 0;
   nnodes = g->knots;
   heap = g->path_heap;
   state = g->path_state;
   *extensions = FALSE;

   /* initialize */
   for( k = 0; k < nnodes; k++ )
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

         if( Is_pterm(g->term[k]) && SCIPisGT(scip, g->prize[k], maxprize) && g->mark[k] )
         {
            maxprize = g->prize[k];
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
      int node;
      int nterms = 0;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearest(heap, state, &count, path);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough (or fixnode), connect its path to current subtree */
         if( !connected[k] && Is_pterm(g->term[k]) &&
                SCIPisGE(scip, g->prize[k], path[k].dist) )
         {
            connected[k] = TRUE;
            *extensions = TRUE;
            path[k].dist = 0.0;
            node = k;

            assert(path[k].edge != UNKNOWN);

            while( !connected[node = g->tail[path[node].edge]] )
            {
               assert(path[node].edge != UNKNOWN);
               connected[node] = TRUE;
               reset(scip, path, heap, state, &count, node);
               assert(state[node]);
            }

            /* have all terminals been reached? */
            if( ++nterms == g->terms - 1 )
               break;
         }
         else if( SCIPisGT(scip, path[k].dist, maxprize) )
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
               correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
      }
   }
}


/** Shortest path heuristic for the RMWCSP
 * Find a directed tree rooted in node 'start' and connecting all terminals as well as all
 *  positive vertices (as long as this is profitable).
 *  */
void graph_path_st_rmw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   SCIP_Real maxprize;
   int   k;
   int   e;
   int   root;
   int   count;
   int   nrterms;
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

   root = g->source;
   count = 0;
   nrterms = 0;
   state = g->path_state;
   heap = g->path_heap;
   nnodes = g->knots;
   maxprize = 0.0;

   for( k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   for( e = g->outbeg[root]; e != EAT_LAST; e = g->oeat[e] )
   {
      if( SCIPisGT(scip, g->cost[e], 0.0) && Is_term(g->term[g->head[e]]) )
      {
         k = g->head[e];
         g->mark[k] = FALSE;
         assert(g->grad[k] == 2);
      }
   }

   for( k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( Is_term(g->term[k]) && g->mark[k] )
         nrterms++;

      if( Is_pterm(g->term[k]) )
         if( SCIPisGT(scip, g->prize[k], maxprize) && k != start )
            maxprize = g->prize[k];
   }

   /* add start vertex to heap */
   k            = start;
   pathdist[k] = 0.0;
   connected[k] = TRUE;

   if( nnodes > 1 )
   {
      int node;
      int nterms = g->terms;
      int termscount = 0;
      int rtermscount = 0;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;
         /* if k is positive vertex and close enough, connect its path to current subtree */
         if( Is_gterm(g->term[k]) && (Is_term(g->term[k]) || SCIPisGE(scip, g->prize[k], pathdist[k])) && !connected[k] )
         {
            if( Is_term(g->term[k]) )
               rtermscount++;
            connected[k] = TRUE;
            pathdist[k] = 0.0;
            node = k;

            if( k != start )
            {
               assert(pathedge[k] != - 1);
               while( !connected[node = g->tail[pathedge[node]]] )
               {
                  assert(pathedge[node] != - 1);

                  connected[node] = TRUE;
                  resetX(scip, pathdist, heap, state, &count, node);
               }
            }

            /* have all terminals been reached? */
            if( ++termscount == nterms )
            {
               break;
            }
         }
         else if( SCIPisGT(scip, pathdist[k], maxprize) && rtermscount >= nrterms )
         {
            break;
         }

         /* update adjacent vertices */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            int m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[e])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}


/** extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE graph_voronoiExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edgecosts */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real**           distarr,            /**< array to store distance from each node to its base */
   int**                 basearr,            /**< array to store the bases */
   int**                 edgearr,            /**< array to store the ancestor edge */
   STP_Bool*             termsmark,          /**< array to mark terminal */
   int*                  reachednodes,       /**< array to mark reached nodes */
   int*                  nreachednodes,      /**< pointer to number of reached nodes */
   int*                  nodenterms,         /**< array to store number of terimals to each node */
   int                   nneighbterms,       /**< number of neighbouring terminals */
   int                   base,               /**< voronoi base */
   int                   countex             /**< count of heap elements */
   )
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count = countex;

   assert(g      != NULL);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   if( g->knots == 0 || nneighbterms <= 0 )
      return SCIP_OKAY;

   heap = g->path_heap;
   state = g->path_state;

   if( g->knots > 1 )
   {
      while( count > 0 && nneighbterms > 0 )
      {
         k = nearest(heap, state, &count, path);
         state[k] = CONNECT;
	 reachednodes[(*nreachednodes)++] = k;
         distarr[k][nodenterms[k]] = path[k].dist;
         edgearr[k][nodenterms[k]] = path[k].edge;
         basearr[k][nodenterms[k]] = base;

         nodenterms[k]++;

         if( termsmark[k] == TRUE )
         {
            termsmark[k] = FALSE;
            if( --nneighbterms == 0 )
            {
               while( count > 0 )
                  reachednodes[(*nreachednodes)++] = heap[count--];

               return SCIP_OKAY;
            }
         }
         else if( !Is_term(g->term[k]) )
         {
            /* iterate over all outgoing edges of vertex k */
            for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
            {
               m = g->head[i];

               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (g->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
            }
         }
      }
      assert(nneighbterms == 0);
   }
   return SCIP_OKAY;
}


/** build a voronoi region, w.r.t. shortest paths, for a given set of bases */
void graph_voronoi(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   STP_Bool*             base,               /**< array to indicate whether a vertex is a Voronoi base */
   int*                  vbase,              /**< voronoi base to each vertex */
   PATH*                 path                /**< path data struture (leading to respective Voronoi base) */
   )
{
   int k;
   int m;
   int i;
   int* heap;
   int* state;
   int count = 0;
   int root;
   int nbases = 0;

   assert(g != NULL);
   assert(scip != NULL);
   assert(path != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);

   if( g->knots == 0 )
      return;

   root = g->source;
   heap = g->path_heap;
   state = g->path_state;

   /* initialize */
   for( i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( base[i] )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY + 1;
         path[i].edge = UNKNOWN;
         state[i] = UNKNOWN;
      }
   }
   assert(nbases > 0);

   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + ((vbase[k] == root)? cost[i] : costrev[i])) )
            {
               assert(!base[m]);
               correct(scip, heap, state, &count, path, m, k, i, ((vbase[k] == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
}

/** 2th next terminal to all non terminal nodes */
void graph_get2next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to first and second nearest terminal) */
   int*                  vbase,              /**< first and second nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int e;
   int root;
   int count;
   int nnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   count = 0;
   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;

      /* copy of node i */
      k = i + nnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + nnodes;
         if( !Is_term(g->term[j]) && SCIPisGT(scip, path[k].dist, path[i].dist + ((root == vbase[i])? cost[e] : costrev[e])) &&
            vbase[i] != vbase[j] && g->mark[j] )
         {
            correct(scip, heap, state, &count, path, k, i, e, ((root == vbase[i])? cost[e] : costrev[e]), FSP_MODE);
            vbase[k] = vbase[i];
         }
      }
   }

   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - nnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - nnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + nnodes;

            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) )
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}

/* 3th next terminal to all non terminal nodes */
void graph_get3next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*       costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int l;
   int v;
   int e;
   int root;
   int count = 0;
   int nnodes;
   int dnnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   nnodes = g->knots;
   dnnodes = 2 * nnodes;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* copy of node i */
      k = i + dnnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;
      state[i + nnodes] = CONNECT;
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + dnnodes;
         if( !Is_term(g->term[j]) && g->mark[j] )
         {
            v = i;
            for( l = 0; l < 2; l++ )
            {
               if( SCIPisGT(scip, path[k].dist, path[v].dist + ((root == vbase[v])? cost[e] : costrev[e])) &&
                  vbase[v] != vbase[j] && vbase[v] != vbase[j + nnodes] )
               {
                  correct(scip, heap, state, &count, path, k, v, e, ((root == vbase[v])? cost[e] : costrev[e]), FSP_MODE);
                  vbase[k] = vbase[v];
               }
               v += nnodes;
            }
         }
      }
   }
   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - dnnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - dnnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + dnnodes;
            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && vbase[j + nnodes] != vbase[k]
               && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) /*TODO(state[jc])??*/
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}


/* 4th next terminal to all non terminal nodes */
void graph_get4next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data struture (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int l;
   int v;
   int e;
   int root;
   int count = 0;
   int nnodes;
   int dnnodes;
   int tnnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   nnodes = g->knots;
   dnnodes = 2 * nnodes;
   tnnodes = 3 * nnodes;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* copy of node i */
      k = i + tnnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;
      state[i + nnodes] = CONNECT;
      state[i + dnnodes] = CONNECT;
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + tnnodes;
         if( !Is_term(g->term[j]) && g->mark[j] )
         {
            v = i;
            for( l = 0; l < 3; l++ )
            {
               if( SCIPisGT(scip, path[k].dist, path[v].dist + ((root == vbase[v])? cost[e] : costrev[e])) &&
                  vbase[v] != vbase[j] && vbase[v] != vbase[j + nnodes] && vbase[v] != vbase[j + dnnodes] )
               {
                  correct(scip, heap, state, &count, path, k, v, e, ((root == vbase[v])? cost[e] : costrev[e]), FSP_MODE);
                  vbase[k] = vbase[v];
               }
               v += nnodes;
            }
         }
      }
   }
   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - tnnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - tnnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + tnnodes;
            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && vbase[j + nnodes] != vbase[k] && vbase[j + dnnodes] != vbase[k]
               && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) /*TODO(state[jc])??*/
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}

/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void graph_get3nextTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   PATH*                 path3,              /**< path data struture (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   assert(g      != NULL);
   assert(path3   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   if( g->stp_type != STP_PCSPG && g->stp_type != STP_RPCSPG )
      for( k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   /* build voronoi diagram */
   graph_voronoiTerms(scip, g, cost, path3, vbase, heap, state);

   /* get 2nd nearest terms */
   graph_get2next(scip, g, cost, costrev, path3, vbase, heap, state);

   /* get 3th nearest terms */
   graph_get3next(scip, g, cost, costrev, path3, vbase, heap, state);

   return;
}

/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void graph_get4nextTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data struture (leading to first, second, third and fouth nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;

   assert(g         != NULL);
   assert(path      != NULL);
   assert(cost      != NULL);
   assert(heap      != NULL);
   assert(state     != NULL);
   assert(costrev   != NULL);

   if( g->stp_type != STP_PCSPG && g->stp_type != STP_RPCSPG )
      for( k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   /* build voronoi diagram */
   graph_voronoiTerms(scip, g, cost, path, vbase, heap, state);

   /* get 2nd nearest terms */
   graph_get2next(scip, g, cost, costrev, path, vbase, heap, state);

   /* get 3th nearest terms */
   graph_get3next(scip, g, cost, costrev, path, vbase, heap, state);

   /* get 4th nearest terms */
   graph_get4next(scip, g, cost, costrev, path, vbase, heap, state);

   return;
}


/** get 4 close terminals to each terminal */
SCIP_RETCODE graph_get4nextTTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
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

   if( g->stp_type != STP_PCSPG && g->stp_type != STP_RPCSPG )
      for( k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   nboundedges = 0;
   for( k = 0; k < nnodes; k ++ )
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

/** build a Voronoi region in presolving, w.r.t. shortest paths, for all terminals */
void graph_voronoiTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path,               /**< path data structure (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int m;
   int i;
   int count = 0;
   int nbases = 0;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   /* initialize */
   for( i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }
   if( nbases == 0 )
      return;
   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) && g->mark[m] )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}


/** build a Voronoi region, w.r.t. shortest paths, for all positive vertices */
void graph_voronoiMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int m;
   int i;
   int count = 0;
   int nbases = 0;
   int nnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(costrev   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = -g->prize[i];
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }
   if( nbases == 0 )
      return;
   if( nnodes > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + costrev[i]) && g->mark[m] && !Is_term(g->term[m]) )
            {
               assert(vbase[m] != m);
               correct(scip, heap, state, &count, path, m, k, i, costrev[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}


/** build a voronoi region, w.r.t. shortest paths, for all terminal and the distance */
SCIP_RETCODE graph_voronoiWithDist(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            distance,           /**< array storing path from a terminal over shortest
                                                incident edge to nearest terminal */
   int*                  minedgepred,        /**< shortest edge predecessor array */
   int*                  vbase,              /**< array containing Voronoi base to each node */
   int*                  minarc,             /**< array to mark whether an edge is one a path corresponding to 'distance' */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array indicating state of each vertex during calculation of Voronoi regions */
   int*                  distnode,           /**< array to store terminal corresponding to distance stored in distance array */
   PATH*                 path                /**< array containing Voronoi paths data */
   )
{
   SCIP_Real new;
   int e;
   int k;
   int m;
   int i;
   int pred;
   int count;
   int nbases;
   int nnodes;
   int prededge;

   assert(g        != NULL);
   assert(path     != NULL);
   assert(cost     != NULL);
   assert(heap     != NULL);
   assert(state    != NULL);
   assert(distance != NULL);

   count = 0;
   nbases = 0;
   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      distance[i] = FARAWAY;
      if( distnode != NULL )
         distnode[i] = UNKNOWN;

      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         state[i] = count;
         path[i].dist = 0.0;
      }
      else
      {
         vbase[i] = UNKNOWN;
         state[i] = UNKNOWN;
         path[i].dist = FARAWAY;
      }
      path[i].edge = UNKNOWN;
   }

   /* initialize predecessor array */

   for( e = 0; e < g->edges; e++ )
      minedgepred[e] = FALSE;

   for( k = 0; k < nbases; k++ )
      if( minarc[k] != UNKNOWN )
         minedgepred[minarc[k]] = TRUE;

   /* if graph contains more than one vertices, start main loop */
   if( nnodes > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         prededge = path[k].edge;

         if( prededge != UNKNOWN )
         {
            pred = g->tail[prededge];

            assert(vbase[k] != UNKNOWN);
            assert(vbase[pred] != UNKNOWN);
            assert(vbase[pred] == vbase[k]);
            assert(g->head[prededge] == k);

            if( !Is_term(g->term[pred]) && g->mark[pred] )
            {
               assert(path[pred].edge != UNKNOWN);
               minedgepred[prededge] = minedgepred[path[pred].edge];
            }

         }

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( state[m] == CONNECT && vbase[m] != vbase[k] && g->mark[m] )
            {
               if( minedgepred[i] || (prededge != UNKNOWN && minedgepred[prededge] ) )
               {
                  new = path[k].dist + g->cost[i] + path[m].dist;
                  if( SCIPisLT(scip, new, distance[vbase[k]]) )
                  {
                     if( distnode != NULL )
                        distnode[vbase[k]] = vbase[m];

                     distance[vbase[k]] = new;
                  }
               }
               if( minedgepred[Edge_anti(i)] || (path[m].edge != UNKNOWN && minedgepred[path[m].edge] ) )
               {
                  new = path[m].dist + g->cost[i] + path[k].dist;
                  if( SCIPisLT(scip, new, distance[vbase[m]]) )
                  {
                     if( distnode != NULL )
                        distnode[vbase[m]] = vbase[k];

                     distance[vbase[m]] = new;
                  }
               }
            }

            /* Check whether the path (to m) including k is shorter than the so far best known.
               In case of equality, also update if k is sucessor of minedge. */
            if( state[m] && g->mark[m] &&
               (SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) ||
                  (prededge != UNKNOWN && minedgepred[prededge] && SCIPisEQ(scip, path[m].dist, path[k].dist + cost[i]))) )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** build voronoi regions, w.r.t. shortest paths, for all terminals and compute the radii */
SCIP_RETCODE graph_voronoiWithRadius(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   GRAPH*                adjgraph,           /**< graph data structure */
   PATH*                 path,               /**< array containing Voronoi paths data */
   SCIP_Real*            rad,                /**< array storing shortest way from a terminal out of its Voronoi region */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark state of each node during calculation */
   )
{
   int* nodesid;
   int i;
   int k;
   int m;
   int vbm;
   int vbk;
   int root;
   int count = 0;
   int nnodes;
   int nterms = 0;
   STP_Bool pc;
   STP_Bool mw;

   assert(graph != NULL);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(rad != NULL);
   assert(vbase != NULL);

   nnodes = graph->knots;

   if( graph->terms == 0 )
      return SCIP_OKAY;

   root = graph->source;
   mw = (graph->stp_type == STP_MWCSP);
   pc = ((graph->stp_type == STP_PCSPG) || (graph->stp_type == STP_RPCSPG));

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      rad[i] = FARAWAY;

      /* set the base of vertex i */
      if( Is_term(graph->term[i]) && graph->mark[i] )
      {
         if( nnodes > 1 )
            heap[++count] = i;

         if( !mw )
         {
            adjgraph->mark[adjgraph->knots] = TRUE;
            graph_knot_add(adjgraph, -1);
         }

         nodesid[i] = nterms++;
         vbase[i] = i;

         if( mw )
            assert(SCIPisGE(scip, graph->prize[i], 0.0));
         if( mw )
            path[i].dist = -graph->prize[i];
         else
            path[i].dist = 0.0;

         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         nodesid[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }

   if( nnodes > 1 )
   {
      SCIP_Real c1;
      SCIP_Real c2;
      SCIP_Real ecost;
      int ne;

      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = graph->outbeg[k]; i != EAT_LAST; i = graph->oeat[i] )
         {
            m = graph->head[i];
            vbm = vbase[m];
            vbk = vbase[k];

            if( state[m] == CONNECT && vbm != vbk && graph->mark[m] )
            {
               assert(graph->mark[vbm]);
               assert(graph->mark[vbk]);
               if( mw )
               {
                  if( SCIPisGT(scip, rad[vbk], path[k].dist) )
                     rad[vbk] = path[k].dist;
                  if( SCIPisGT(scip, rad[vbm], path[m].dist) )
                     rad[vbm] = path[m].dist;
#if 0
                  if( SCIPisGT(scip, path[m].dist + graph->prize[vbm], path[k].dist + graph->prize[vbk]) )
                  ecost = graph->prize[vbm] - path[k].dist;
                  else
                  ecost = graph->prize[vbk] - path[m].dist;
                  if( SCIPisLT(scip, ecost, 0.0) )
                  ecost = 0.0;
                  /* find edge in adjgraph */
                  for( ne = adjgraph->outbeg[nodesid[vbk]]; ne != EAT_LAST; ne = adjgraph->oeat[ne] )
                  if( adjgraph->head[ne] == nodesid[vbm] )
                  break;

                  /* edge exists? */
                  if( ne != EAT_LAST )
                  {
                     assert(ne >= 0);
                     assert(adjgraph->head[ne] == nodesid[vbm]);
                     assert(adjgraph->tail[ne] == nodesid[vbk]);
                     if( SCIPisLT(scip, adjgraph->cost[ne], ecost) )
                     {
                        adjgraph->cost[ne] = ecost;
                        adjgraph->cost[Edge_anti(ne)] = ecost;
                     }
                  }
                  else
                  {
                     graph_edge_add(scip, adjgraph, nodesid[vbm], nodesid[vbk], ecost, ecost);
                  }
#endif
               }
               else
               {
                  if( graph->stp_type == STP_DHCSTP )
                  {
                     if( m == root )
                        c1 = path[m].dist + costrev[i];
                     else
                        c1 = path[m].dist + cost[i];
                     if( k == root )
                        c2 = path[k].dist + cost[i];
                     else
                        c2 = path[k].dist + costrev[i];

                     if( SCIPisGT(scip, c1, c2) )
                        ecost = c2;
                     else
                        ecost = c1;
                  }
                  else
                  {
                     if( SCIPisGT(scip, path[m].dist, path[k].dist) )
                        ecost = path[k].dist + cost[i];
                     else
                        ecost = path[m].dist + cost[i];
                  }

                  if( pc )
                  {
                     if( SCIPisGT(scip, ecost, graph->prize[vbm])
                           && root != vbm )
                        ecost = graph->prize[vbm];
                     if( SCIPisGT(scip, ecost, graph->prize[vbk])
                           && root != vbk )
                        ecost = graph->prize[vbk];
                  }

                  /* find edge in adjgraph */
                  for( ne = adjgraph->outbeg[nodesid[vbk]]; ne != EAT_LAST; ne = adjgraph->oeat[ne] )
                     if( adjgraph->head[ne] == nodesid[vbm] )
                        break;

                  /* edge exists? */
                  if( ne != EAT_LAST )
                  {
                     assert(ne >= 0);
                     assert(adjgraph->head[ne] == nodesid[vbm]);
                     assert(adjgraph->tail[ne] == nodesid[vbk]);
                     if( SCIPisGT(scip, adjgraph->cost[ne], ecost) )
                     {
                        adjgraph->cost[ne] = ecost;
                        adjgraph->cost[Edge_anti(ne)] = ecost;
                     }
                  }
                  else
                  {
                     graph_edge_add(scip, adjgraph, nodesid[vbm], nodesid[vbk],
                           ecost, ecost);
                  }

                  if( SCIPisGT(scip, rad[vbk],
                        path[k].dist + ((vbk == root) ? cost[i] : costrev[i])) )
                     rad[vbk] = path[k].dist
                           + ((vbk == root) ? cost[i] : costrev[i]);

                  if( SCIPisGT(scip, rad[vbm],
                        path[m].dist + ((vbm == root) ? costrev[i] : cost[i])) )
                     rad[vbm] = path[m].dist
                           + ((vbm == root) ? costrev[i] : cost[i]);
               }
            }

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && graph->mark[m] && !Is_term(graph->term[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + ((vbk == root)? cost[i] : costrev[i])) )
            {
               correct(scip, heap, state, &count, path, m, k, i, ((vbk == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbk;
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &nodesid);

   return SCIP_OKAY;
}

/** Build partition of graph for MWCSP into regions around the positive vertices.
 * Store the weight of a minimum weight center-boundary path for each region
 * in the radius array (has to be reverted to obtain the final r() value).
 */
void graph_voronoiWithRadiusMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< array containing graph decomposition data */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            radius,             /**< array storing shortest way from a positive vertex out of its region */
   int*                  vbase,              /**< array containing base of each node */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark state of each node during calculation */
   )
{
   SCIP_Real* nodeweight;
   int i;
   int k;
   int m;
   int count;
   int nbases;
   int nnodes;
   int nterms;

   assert(g != NULL);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(radius != NULL);
   assert(vbase != NULL);

   count = 0;
   nbases = 0;
   nnodes = g->knots;
   nterms = g->terms - 1;
   nodeweight = g->prize;

   assert(nodeweight != NULL);

   if( nnodes == 0 || nterms <= 0 )
      return;

   assert(!g->mark[g->source]);

   /* initialize data */
   for( i = 0; i < nnodes; i++ )
   {
      radius[i] = FARAWAY;

      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }

   if( nbases == 0 )
      return;

   if( nnodes > 1 )
   {
      int basem;
      int basek;

      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;
         basek = vbase[k];

         assert(g->mark[basek]);

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( !(g->mark[m]) )
               continue;

            basem = vbase[m];

            /* are both m and k finally in different regions? */
            if( basem != basek && (state[m] == CONNECT || vbase[m] == m ||
               SCIPisGE(scip, path[k].dist, nodeweight[basek])) )
            {
               if( SCIPisGT(scip, radius[basek], path[k].dist) )
                  radius[basek] = path[k].dist;
               if( (state[m] == CONNECT || vbase[m] == m) && SCIPisGT(scip, radius[basem], path[m].dist) )
               {
                  assert(g->mark[basem]);
                  radius[basem] = path[m].dist;
               }
            }

            /* is the distance of vertex k smaller than the weight of its base node and smaller than the temp. radius? */
            if( SCIPisLT(scip, path[k].dist, nodeweight[basek]) && SCIPisLT(scip, path[k].dist, radius[basek]) )
            {
               /* check whether the path (to m) including k is shorter than the so far best known */
               if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) )
               {
                  assert(vbase[m] != m);
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
                  vbase[m] = basek;
               }
            }
         }
      }
   }
   return;
}

/** repair a Voronoi diagram for a given set of base nodes */
void graph_voronoiRepair(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   int*                  count,              /**< pointer to number of heap elements */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   PATH*                 path,               /**< Voronoi paths data struture */
   int*                  newedge,            /**< the new edge */
   int                   crucnode,           /**< the current crucial node */
   UF*                   uf                  /**< union find data structure */
   )
{
   int k;
   int m;
   int i;
   int e;
   int* heap;
   int* state;
   int node1;
   int node2;

   *newedge = UNKNOWN;
   e = UNKNOWN;
   assert(g != NULL);
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);
   assert(path != NULL);
   assert(cost != NULL);

   if( g->knots == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;

   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( *count > 0 )
      {
         /* get the next (i.e. a nearest) vertex from the heap */
         k = nearest(heap, state, count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;
         assert(g->mark[k]);
         assert(g->mark[vbase[k]]);
         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && (SCIPisGT(scip, path[m].dist, path[k].dist + cost[i])) )/*&& g->mark[m] ) */
            {
               assert(g->mark[m]);
               correct(scip, heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }

            /* check whether there is a better new boundary edge adjacent to vertex k */
            else
            {
               node1 = SCIPStpunionfindFind(uf, vbase[m]);
               node2 = SCIPStpunionfindFind(uf, vbase[k]);
               if( state[m] == CONNECT && ((node1 == crucnode) != (node2 == crucnode)) && g->mark[m] && g->mark[vbase[m]]
                  && g->mark[node1] && g->mark[node2] && ((e == UNKNOWN) || SCIPisGT(scip,
                        (path[g->tail[e]].dist + cost[e] + path[g->head[e]].dist), path[k].dist + cost[i] + path[m].dist)) )
                  e = i;
            }
         }
      }
   }
   *newedge = e;
}


/** repair the Voronoi diagram for a given set nodes */
void graph_voronoiRepairMult(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   int*                  count,              /**< pointer to number of heap elements */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   int*                  boundedges,         /**< boundary edges */
   int*                  nboundedges,        /**< number of boundary edges */
   STP_Bool*             nodesmark,          /**< array to mark temporarily discarded nodes */
   UF*                   uf,                 /**< union find data structure */
   PATH*                 path                /**< Voronoi paths data structure */
   )
{
   int k;
   int m;
   int i;
   int* heap;
   int* state;
   int node1;
   int node2;

   assert(g != NULL);
   assert(g->mark != NULL);
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);
   assert(path != NULL);
   assert(cost != NULL);

   if( g->knots == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;

   assert(heap != NULL);
   assert(state != NULL);

   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( (*count) > 0 )
      {
         /* get the next (i.e. a nearest) vertex from the heap */
         k = nearest(heap, state, count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];
#if 1
            if( !(g->mark[m]) )
               continue;
#endif
            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && (SCIPisGT(scip,path[m].dist, path[k].dist + cost[i])) )
            {
               correct(scip, heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
            /* check whether there is a new boundary edge adjacent to vertex k */
            else if( (state[m] == CONNECT) && ((node1 = SCIPStpunionfindFind(uf, vbase[m])) != (node2 = SCIPStpunionfindFind(uf, vbase[k])))
               && g->mark[node1] && g->mark[node2] && (nodesmark[node1] || nodesmark[node2]) )
            {
               boundedges[(*nboundedges)++] = i;
            }
         }
      }
   }
}
