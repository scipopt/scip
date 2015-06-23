/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphpath.c
 * @brief  Shortest path based graph algorithms for Steiner problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * Das Nachfolgende ist eine Implementierung von Dijkstras Algorithmus
 * mit einem Heap zur Verwaltung der aktiven Knoten.
 *
 * Heaproutinen siehe:
 *
 *       Jon Bentley, Programming Pearls, Addison-Wesley 1989
 *
 * Die Heaptabelle wird mit n (Knoten) Elementen initialisiert, aber erst ab
 * Element 1 benutzt, da aber Knoten 0 als Erster in die Tabelle aufgenommen
 * und dann sofort wieder entfernt wird, koennen maximal n-1 Knoten
 * in der Tabelle sein.
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

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   if ((*count) > 2)
      if (LT(pathdist[heap[3]], pathdist[heap[2]]))
         c++;

   while((c <= (*count)) && GT(pathdist[heap[j]], pathdist[heap[c]]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
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
   assert(g->path_heap  != NULL);
   assert(g->path_state != NULL);

   SCIPfreeMemoryArray(scip, &(g->path_heap));
   SCIPfreeMemoryArray(scip, &(g->path_state));
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
   int                   mode,               /**< shortest path (FSP_MODE) or minimum spanning tree (MST_MODE)? */
   int                   start,              /**< start vertex */
   SCIP_Real*            cost,               /**< edge costs */
   PATH*                 path                /**< shortest paths data structure */
   )
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count;

   assert(scip      != NULL);
   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) ||  (mode == MST_MODE));
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   /* no nodes?, return*/
   if (p->knots == 0)
      return;

   heap = p->path_heap;
   state = p->path_state;

   /* Erstmal alles auf null, unbekannt und weit weg
    */
   for(i = 0; i < p->knots; i++)
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY + 1;
      path[i].edge = -1;
   }
   /* Startknoten in den Heap
    */
   k            = start;
   path[k].dist = 0.0;

   /* Wenn nur ein Knoten drin ist funktioniert der Heap nicht sehr gut,
    * weil dann fuer genau 0 Elemente Platz ist.
    */
   if( p->knots > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* Wenn nichts mehr auf dem Heap ist, sind wir fertig
       * und jetzt erstmal Hula Loop
       */
      while( count > 0 )
      {
         /* Na, wer ist der Naechste ?
          */
         k = nearest(heap, state, &count, path);

         /* Wieder einen erledigt
          */
         state[k] = CONNECT;

         /* Verbunden Knoten berichtigen ...
          *
          * Wenn ein Knoten noch nicht erledigt ist
          * werden wir dann auf diesem Wege besser ?
          */
         for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
         {
            m = p->head[i];
            /* node not scanned and valid? */
            if( (state[m]) && (p->mark[m]) )
            {
	       /* closer than previously? */
               if( SCIPisGT(scip, path[m].dist, (mode == MST_MODE) ? cost[i] : (path[k].dist + cost[i])) )
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], mode);
            }
         }
      }
   }
}

/** limited Dijkstra, stoping at terminals */
void sdpaths(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real*            cost,               /**< edge costs */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  memlbl,             /**< array to save labelled nodes */
   int*                  nlbl,               /**< number of labelled nodes */
   int                   tail,               /**< tail of the edge */
   int                   head,               /**< head of the edge */
   int                   limit               /**< maximum number of edges to consider during execution */
   )
{
   int   k;
   int   m;
   int   e;
   int   limit1;
   int count;
   int   nchecks;

   assert(g      != NULL);
   assert(heap   != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(nlbl   != NULL);
   assert(memlbl != NULL);

   limit1 = limit - limit / 3;
   *nlbl = 0;
   nchecks = 0;

   if( g->grad[tail] == 0 || g->grad[head] == 0 )
      return;

   assert(g->mark[head] && g->mark[tail]);
   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;
   count = 0;
   g->mark[head] = FALSE;
   for( e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
   {
      m = g->head[e];
      if( g->mark[m] )
      {
	 assert(SCIPisGT(scip, path[m].dist, path[tail].dist + cost[e]));
         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         correct(scip, heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);
      }
      if( nchecks++ > limit1 )
         break;
   }
   g->mark[head] = TRUE;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      k = nearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* stop at terminals */
      if( Is_term(g->term[k]) || k == head )
         continue;

      /* correct incident nodes */
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         m = g->head[e];
         if( state[m] && g->mark[m] && SCIPisGT(scip, path[m].dist, path[k].dist + cost[e]) )
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


/** Dijkstra's algorithm starting from node 'start' */
void graph_path_execX(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start vertex */
   SCIP_Real*            cost,               /**< edgecosts */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge            /**< predecessor edge array (on vertices) */
   )
{
   int   k;
   int   m;
   int   i;
   int   count;
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

   if (g->knots == 0)
      return;

   heap = g->path_heap;
   state = g->path_state;

   for( i = 0; i < g->knots; i++ )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k            = start;
   pathdist[k] = 0.0;

   if( g->knots > 1 )
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
#if 0
	    if( Is_term(g->term[k]) )
	       continue;
#endif
            m = g->head[i];
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
   char*                 connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int   k;
   int   m;
   int   e;
   int count;
   int* heap;
   int* state;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   heap = g->path_heap;
   state = g->path_state;
   count = 0;

   /* initialize */
   for( k = 0; k < g->knots; k++ )
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

   if( g->knots > 1 )
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

	    /* is m not connected, allowed and closer? */
            if( !connected[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[e])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}
#if 0
/* computes the shortest path from each terminal to every other vertex */
void calculate_distances(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH**                path,
   double*               cost,
   int                   mode
   )
{
   int i;

   assert(mode == FSP_MODE || mode == BSP_MODE);

   SCIPdebug(fputc('C', stdout));
   SCIPdebug(fflush(stdout));

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]) && (g->grad[i] > 0))
      {
         if (path[i] == NULL)
            path[i] = malloc((size_t)g->knots * sizeof(PATH));

         assert(path[i] != NULL);

         graph_path_exec(scip, g, mode, i, cost, path[i]);
      }
      else
      {
         if (path[i] != NULL)
         {
            free(path[i]);

            path[i] = NULL;
         }
      }
   }
}
#endif
/** extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE voronoi_extend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edgecosts */
   PATH*                 path,               /**< shortest paths data structure */
   VLIST**               adjterms,           /**< structure for adjacent terminal data */
   char*                 termsmark,          /**< array to mark terminal */
   int*                  reachednodes,       /**< array to mark reached nodes */
   int*                  nreachednodes,      /**< pointer to number of reached nodes */
   int*                  nodenterms,         /**< array to store number of terimals to each node */
   int                   nneighbterms,       /**< number of neighbouring terminals */
   int                   base,               /**< voronoi base */
   int                   countex             /**< number of heap elements */
   )
{
   int* heap;
   int* state;
   int   k;
   int   m;
   int   i;
   int count = countex;
   VLIST* curr;

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
      /*
       */
      while( count > 0 && nneighbterms > 0 )
      {

         k = nearest(heap, state, &count, path);

         state[k] = CONNECT;
	 reachednodes[(*nreachednodes)++] = k;
         SCIP_CALL( SCIPallocMemory(scip, &curr) );
         curr->dist = path[k].dist;
         curr->edge = path[k].edge;
         curr->base = base;
         curr->next = adjterms[k];
         adjterms[k] = curr;
         nodenterms[k]++;

         if( termsmark[k] == TRUE )
         {
            termsmark[k] = FALSE;
            if( --nneighbterms == 0 )
            {
               while( count > 0 )
               {
                  reachednodes[(*nreachednodes)++] = heap[count--];
               }
               return SCIP_OKAY;
            }

         }
         else
         {
            /* iterate over all outgoing edges of vertex k */
            for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
            {
               m = g->head[i];


               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (g->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
               {
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               }
            }
         }
      }
      assert(nneighbterms == 0);
   }
   return SCIP_OKAY;
}


/** extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE voronoi_extend2(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edgecosts */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real**           distarr,            /**< array to store distance from each node to its base */
   int**                 basearr,            /**< array to store the bases */
   int**                 edgearr,            /**< array to store the ancestor edge */
   char*                 termsmark,          /**< array to mark terminal */
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
void voronoi(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   char*                 base,               /**< array to indicate whether a vertex is a Voronoi base */
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

   assert(scip != NULL);
   assert(g      != NULL);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   if( g->knots == 0 )
      return;

   root = g->source[0];
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
         state[i]     = UNKNOWN;
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

/* 2th next terminal to all non terminal nodes */
void get2next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data struture (leading to first and second nearest terminal) */
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
   int count = 0;
   int nnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source[0];
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
void get3next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
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

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source[0];
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

/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void getnext3terms(
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

   if( g->stp_type != STP_PRIZE_COLLECTING && g->stp_type != STP_ROOTED_PRIZE_COLLECTING )
   for( k = 0; k < g->knots; k++ )
      g->mark[k] = (g->grad[k] > 0);

   /* build voronoi diagram */
   voronoi_terms(scip, g, cost, path3, vbase, heap, state);

   /* get 2nd nearest terms */
   get2next(scip, g, cost, costrev, path3, vbase, heap, state);

   /* get 3th nearest terms */
   get3next(scip, g, cost, costrev, path3, vbase, heap, state);

   return;
}


/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void voronoi_terms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   PATH*                 path,               /**< path data struture (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex*/
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


/** build a voronoi region, w.r.t. shortest paths, for all terminal and the distance */
SCIP_RETCODE voronoi_dist(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            distance,           /**< array storing path from a terminal over shortest
                                                incident edge to nearest terminal */
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
   int count = 0;
   int nbases = 0;
   int nnodes;
   char* minedgepred;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(distance   != NULL);

   nnodes = g->knots;
   SCIP_CALL( SCIPallocBufferArray(scip, &minedgepred, g->edges) );

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

   for( e = 0; e < g->edges; e++ )
      minedgepred[e] = FALSE;

   for( k = 0; k < nbases; k++ )
   {
      if( minarc[k] != UNKNOWN )
         minedgepred[minarc[k]] = TRUE;
   }
   if( nnodes > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;
         if( path[k].edge != UNKNOWN )
         {
            assert(g->head[path[k].edge] == k);
            pred = g->tail[path[k].edge];
            assert(vbase[k] != UNKNOWN);
            assert(vbase[pred] != UNKNOWN);
            assert(vbase[pred] == vbase[k]);
            if( !Is_term(g->term[pred]) && g->mark[pred] )
            {
               assert(path[pred].edge != UNKNOWN);
               minedgepred[path[k].edge] = minedgepred[path[pred].edge];
            }
         }

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( state[m] == CONNECT && vbase[m] != vbase[k] && g->mark[m] )
            {
               if( minedgepred[i] || (path[k].edge != UNKNOWN && minedgepred[path[k].edge] ) )
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

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) && g->mark[m] )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &minedgepred);

   return SCIP_OKAY;
}

/** build voronoi regions, w.r.t. shortest paths, for all terminals and compute the radii */
SCIP_RETCODE voronoi_radius(
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

   assert(graph != NULL);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(rad != NULL);
   assert(vbase != NULL);

   nnodes = graph->knots;
   if( nnodes == 0 || graph->terms == 0 )
      return SCIP_OKAY;
   root = graph->source[0];
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

         graph_knot_add(adjgraph, -1);
         nodesid[i] = nterms++;
         vbase[i] = i;
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
      SCIP_Real ecost;
      SCIP_Real c1;
      SCIP_Real c2;
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
               if( graph->stp_type == STP_HOP_CONS )
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

               if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
               {
                  if( SCIPisGT(scip, ecost, graph->prize[vbm]) && root != vbm )
                     ecost = graph->prize[vbm];
                  if( SCIPisGT(scip, ecost, graph->prize[vbk]) && root != vbk )
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
                     adjgraph->cost[ne]            = ecost;
                     adjgraph->cost[Edge_anti(ne)] = ecost;
                  }
               }
               else
               {
                  graph_edge_add(scip, adjgraph, nodesid[vbm], nodesid[vbk], ecost, ecost);
               }

               if( SCIPisGT(scip, rad[vbk], path[k].dist + ((vbk == root)? cost[i] : costrev[i])) )
                  rad[vbk] = path[k].dist + ((vbk == root)? cost[i] : costrev[i]);

               if( SCIPisGT(scip, rad[vbm], path[m].dist + ((vbm == root)? costrev[i] : cost[i])) )
                  rad[vbm] = path[m].dist + ((vbm == root)? costrev[i] : cost[i]);
            }

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && graph->mark[m] && SCIPisGT(scip, path[m].dist, path[k].dist + ((vbk == root)? cost[i] : costrev[i])) )
            {
               assert(!Is_term(graph->term[m]));
               correct(scip, heap, state, &count, path, m, k, i, ((vbk == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbk;
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &nodesid);
   return SCIP_OKAY;
}

/** repair the voronoi diagram for SL reduction */
void voronoi_slrepair(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   PATH*                 path,               /**< array to store */
   int*                  vbase,              /**< array containing Voronoi base of each node*/
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array indicating state of each vertex during calculation of Voronoi regions */
   int                   start,              /**< node to start repair from */
   int                   adjnode             /**< adjacent node */
   )
{
   int k;
   int m;
   int i;
   int count = 0;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   /* initialize */
   for( i = 0; i < g->knots; i++ )
      state[i] = UNKNOWN;

   if( SCIPisGT(scip, path[start].dist, path[adjnode].dist) )
   {
      if( vbase[adjnode] == adjnode )
      {
         assert(Is_term(g->term[start]));
         vbase[start] = start;
      }
      else
      {
         vbase[start] = vbase[adjnode];
      }
      path[start].dist = path[adjnode].dist;
      path[start].edge = path[adjnode].edge;
   }
   k = start;

   if( g->knots > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

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
            if( (state[m] != CONNECT) && SCIPisGE(scip, path[m].dist, path[k].dist + cost[i]) && g->mark[m] )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}

/** repair the voronoi diagram for a given set nodes */
void voronoi_repair(
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
               node1 = SCIPunionfindFind(uf, vbase[m]);
               node2 = SCIPunionfindFind(uf, vbase[k]);
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


/** repair the voronoi diagram for a given set nodes */
void voronoi_repair_mult(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   int*                  count,              /**< pointer to number of heap elements */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   int*                  boundedges,         /**< boundary edges */
   int*                  nboundedges,        /**< number of boundary edges */
   char*                 nodesmark,          /**< array to mark temporarily discarded nodes */
   UF*                   uf,                 /**< union find data structure */
   PATH*                 path                /**< Voronoi paths data struture */
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
            /*   printf("scrut edge %d->%d vbases:  %d %d \n ", k,m, vbase[k], vbase[m]);
                 printf("          uf : %d %d \n ", SCIPunionfindFind(uf, vbase[k]), SCIPunionfindFind(uf, vbase[m]) ); */
            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && (SCIPisGT(scip,path[m].dist, path[k].dist + cost[i])) && g->mark )
            {
               correct(scip, heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
            /* check whether there is a new boundary edge adjacent to vertex k */
            else if( (state[m] == CONNECT) && ((node1 = SCIPunionfindFind(uf, vbase[m])) != (node2 = SCIPunionfindFind(uf, vbase[k])))
               && g->mark[node1] && g->mark[node2] && (nodesmark[node1] || nodesmark[node2]) )
            {
               boundedges[(*nboundedges)++] = i;
               /*printf("adding new boundaryedge [%d] %d_%d \n", (*nboundedges) - 1, k,m);*/
            }
         }
      }
   }
}

/* ARGSUSED */
void graph_path_length(
   const GRAPH*          g,                  /**< graph data structure */
   const PATH*           p                   /**< path data structure */
   )
{
#ifndef NDEBUG
   int    len  = 0;
   SCIP_Real dist = 0.0;
   int    e;
   int    i;

   for(i = 0; i < g->knots; i++)
   {
      if ((e = p[i].edge) < 0)
         continue;

      len++;
      dist += g->cost[e];
   }
   (void)printf("Graph path length: %d  Distance=%g\n", len, dist);

#endif /* NDEBUG */
}
