/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: grphpath.c                                                    */
/*   Name....: Find Shortest Path / Minimum Spanning Tree                    */
/*   Author..: Thorsten Koch, Daniel Rehfeldt                                */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,GRPHPATH_C) -esym(766,stdlib.h)                             */

/* Das Nachfolgende ist eine Implementierung von Dijkstras Algorithmus
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
#define GRPHPATH_C
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "scip/scip.h"

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
   path[l].hops = path[k].hops + 1;

   /* new node? */
   if( state[l] == UNKNOWN )
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up */
   j = state[l];
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

   while((j > 1) && GT(pathdist[heap[c]], pathdist[heap[j]]))
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

   //assert(state[node] == CONNECT);
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
void graph_path_init(
   GRAPH* p)
{
   assert(p->path_heap  == NULL);
   assert(p->path_state == NULL);

   p->path_heap  = malloc((size_t)(p->knots + 1) * sizeof(int));
   p->path_state = malloc((size_t)p->knots * sizeof(int));

   assert(p->path_heap  != NULL);
   assert(p->path_state != NULL);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : EXIT shortest PATH algorithm                             ---*/
/*--- Function : Gibt den bei der initialisierung angeforderten Speicher  ---*/
/*---            wieder frei                                              ---*/
/*--- Parameter: Keine                                                    ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
void graph_path_exit(
   GRAPH* p)
{
   assert(p->path_heap  != NULL);
   assert(p->path_state != NULL);

   free(p->path_heap);
   free(p->path_state);

   p->path_heap  = NULL;
   p->path_state = NULL;
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
   const GRAPH*  p,
   int           mode,
   int           start,
   SCIP_Real*    cost,
   PATH*         path)
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   double pathdist;
   double pathhops;
   int count;

   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) || (mode == BSP_MODE) || (mode == MST_MODE));
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   /* Kein Baum ohne Knoten
    */
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
      path[i].hops = 0;
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
         for(i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i])
         {

            m = p->head[i];
            /* 1. Ist der Knoten noch nicht festgelegt ?
             *    Ist der wohlmoeglich tabu ?
             */
            if ((state[m]) && (p->mark[m]))
            {
               /* 2. Ist es ueberhaupt eine Verbesserung diesen Weg zu nehmen ?
                *    Wenn ja, dann muss die Entferung kuerzer sein.
                */
               pathdist = ((mode == MST_MODE) ? cost[i] : (mode == FSP_MODE) ? (path[k].dist + cost[i]) :
                  (path[k].dist + cost[Edge_anti(i)]));
               pathhops = path[k].hops + 1;
               if( GT(path[m].dist, pathdist) ||
                  (mode != MST_MODE && EQ(path[m].dist, pathdist) && GT(path[m].hops, pathhops)) )
	       {
                  correct(heap, state, &count, path, m, k, i, (mode == BSP_MODE) ? cost[Edge_anti(i)] : cost[i], mode);
	       }
            }
         }
      }
   }
}

/* limited dijkstra, stoping at terminals */
void sdtail(
   SCIP* scip,
   const GRAPH*  g,
   PATH*         path,
   SCIP_Real*    cost,
   int*          heap,
   int*          state,
   int*          memlbl,
   int*          nlbl,
   int           tail,
   int           head,
   int           limit
   )
{
   int   k;
   int   m;
   int   i;
   int count;

   int   nnodes;
   int   nchecks;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   *nlbl = 0;
   nchecks = 0;
   nnodes = g->knots;

   assert(nnodes > 1);

   /* TODO: outsource */
   for( i = 0; i < nnodes; i++ )
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY;
      path[i].edge = -1;
   }

   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;
   count       = 0;
   g->mark[head] = FALSE;
   for( i = g->outbeg[tail]; i != EAT_LAST; i = g->oeat[i] )
   {
      m = g->head[i];
      if( nchecks++ > limit / 2 )
         break;
      if( g->mark[m] && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) )
      {
         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
      }
   }
   g->mark[head] = TRUE;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      k = nearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* stop at terminals */
      if( Is_term(g->term[k]) )
         continue;

      /* correct incident nodes */
      for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
      {
         m = g->head[i];
         if( nchecks++ > limit )
            break;
         if( state[m] && g->mark[m] && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) )
         {
            /* m labelled for the first time? */
            if( state[m] == UNKNOWN )
	       memlbl[(*nlbl)++] = m;
            correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
         }
      }
   }
}


/* limited dijkstra, stoping at terminals */
SCIP_Real sdhead(
   SCIP* scip,
   const GRAPH*  g,
   PATH*         path,
   PATH*         pathtail,
   SCIP_Real*    cost,
   int*          heap,
   int*          state,
   int*          statetail,
   int*          memlbl,
   int*          nlbl,
   int           tail,
   int           head,
   int           limit
   )
{
   SCIP_Real    dist;
   SCIP_Real    sdist;
   int   k;
   int   m;
   int   i;
   int count;
   int   nnodes;
   int   nchecks;

   assert(scip      != NULL);
   assert(g      != NULL);
   assert(path   != NULL);
   assert(pathtail   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(statetail   != NULL);
   assert(memlbl   != NULL);
   assert(cost   != NULL);

   (*nlbl) = 0;
   nchecks = 0;
   nnodes = g->knots;
   if( statetail[head] != UNKNOWN )
   {
      sdist = pathtail[head].dist;
      assert(SCIPisGT(scip, FARAWAY, sdist));
   }
   else
      sdist = FARAWAY;

   /* TODO: outsource */
   for( i = 0; i < nnodes; i++ )
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY;
      path[i].edge = UNKNOWN;
   }

   k            = head;
   path[k].dist = 0.0;
   memlbl[(*nlbl)++] = head;

   count       = 1;
   heap[count] = k;
   state[k]    = count;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      k = nearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* stop at terminals */
      if( Is_term(g->term[k]) )
         continue;

      /* correct incident nodes */
      for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
      {
         m = g->head[i];
         if( nchecks++ > limit )
            break;
         if( state[m] && g->mark[m] && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) )
         {
            /* m labelled for the first time? */
            if( state[m] == UNKNOWN )
	       memlbl[(*nlbl)++] = m;
            correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);

            if( statetail[m] != UNKNOWN )
            {
	       assert(SCIPisGT(scip, FARAWAY, pathtail[m].dist));
	       assert(SCIPisGT(scip, FARAWAY, path[m].dist));
               if( Is_term(g->term[m]) )
               {
		  dist = 0.0;
		  if( SCIPisLT(scip, dist, path[m].dist) )
		     dist = path[m].dist;
                  if( SCIPisLT(scip, dist, pathtail[m].dist) )
		     dist = pathtail[m].dist;
		  if( SCIPisGT(scip, sdist, dist) )
		     sdist = dist;
		  /* PC*/
               }
               else
               {
                  if( SCIPisGT(scip, sdist, path[m].dist + pathtail[m].dist) )
		     sdist = dist;
               }
            }
         }
      }
   }
   return sdist;
}


/* Dijkstra's algorithm starting from node 'start' */
void graph_path_execX(
   SCIP* scip,
   const GRAPH*  p,
   int           start,
   SCIP_Real*    cost,
   SCIP_Real*   pathdist,
   int*         pathedge
   )
{
   int   k;
   int   m;
   int   i;
   int   count;
   int* heap;
   int* state;

   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cost   != NULL);

   if (p->knots == 0)
      return;

   heap = p->path_heap;
   state = p->path_state;

   for( i = 0; i < p->knots; i++ )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k            = start;
   pathdist[k] = 0.0;

   if( p->knots > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         k = nearestX(heap, state, &count, pathdist);

         state[k] = CONNECT;

         for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
         {
            m = p->head[i];
            if( state[m] && p->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[i])) )
               correctX(heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}

/* Find a directed tree rooted in node 'start' and spanning all terminals */
void graph_path_st(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real* cost,
   SCIP_Real*   pathdist,
   int*   pathedge,
   int           start,
   char*         connected
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
               correctX(heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}

/* computes the shortest path from each terminal to every other vertex */
void calculate_distances(
   const GRAPH* g,
   PATH** path,
   double* cost,
   int mode)
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

         graph_path_exec(g, mode, i, cost, path[i]);
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




/* extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE voronoi_extend(
   SCIP* scip,
   const GRAPH*  p,
   SCIP_Real*    cost,
   PATH*         path,
   VLIST**       adjterms,
   char*         termsmark,
   int*          reachednodes,
   int*          nreachednodes,
   int*          nodenterms,
   int           nneighbterms,
   int           base,
   int           countex)
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count = countex;
   VLIST* curr;

   assert(p      != NULL);
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);


   if( p->knots == 0 || nneighbterms <= 0 )
      return SCIP_OKAY;

   heap = p->path_heap;
   state = p->path_state;


   if( p->knots > 1 )
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
            for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
            {
               m = p->head[i];


               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (p->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
               {
                  correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               }
            }
         }
      }
      assert( nneighbterms == 0);
   }
   return SCIP_OKAY;
}




/* extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE voronoi_extend2(
   SCIP* scip,
   const GRAPH*  p,
   SCIP_Real*    cost,
   PATH*         path,
   SCIP_Real**   distarr,
   int**         basearr,
   int**         edgearr,
   char*         termsmark,
   int*          reachednodes,
   int*          nreachednodes,
   int*          nodenterms,
   int           nneighbterms,
   int           base,
   int           countex
   )
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count = countex;

   assert(p      != NULL);
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   if( p->knots == 0 || nneighbterms <= 0 )
      return SCIP_OKAY;

   heap = p->path_heap;
   state = p->path_state;

   if( p->knots > 1 )
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
         else if( !Is_term(p->term[k]) )
         {
#if 0
	    if( Is_term(p->term[k]) )
               printf("terminal reached! \n");
#endif
            /* iterate over all outgoing edges of vertex k */
            for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
            {
               m = p->head[i];

               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (p->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
               {
                  correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               }
            }
         }
      }
      assert( nneighbterms == 0);
   }
   return SCIP_OKAY;
}


#if 0
/* extend a voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE voronoi_extend3(
   SCIP* scip,
   const GRAPH*  p,
   const SCIP_Real* cost,
   PATH*         path,
   GNODE***      nodedist,
   int**         basearr,
   int**         edgearr,
   char*         termsmark,
   int*          reachednodes,
   int*          nreachednodes,
   int*          nodenterms,
   int           nneighbterms,
   int           base,
   int           countex)
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count = countex;


   assert(p      != NULL);
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);


   if( p->knots == 0 || nneighbterms <= 0 )
      return SCIP_OKAY;

   heap = p->path_heap;
   state = p->path_state;


   if( p->knots > 1 )
   {

      while( count > 0 && nneighbterms > 0 )
      {

         k = nearest(heap, state, &count, path);

         state[k] = CONNECT;
	 reachednodes[(*nreachednodes)++] = k;
         SCIP_CALL( SCIPallocBuffer(scip, &nodedist[k][nodenterms[k]]) );
         nodedist[k][nodenterms[k]]->dist = path[k].dist;
         nodedist[k][nodenterms[k]]->number = k;
         edgearr[k][nodenterms[k]] = path[k].edge;
         basearr[k][nodenterms[k]] = base;

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
            for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
            {
               m = p->head[i];

               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (p->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
               {
                  correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               }
            }
         }
      }
      assert( nneighbterms == 0);
   }
   return SCIP_OKAY;
}
#endif
/*** build a voronoi region, w.r.t. shortest paths, for a given set of bases ***/
void voronoi(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   char*         base,
   int*          vbase,
   PATH*         path
   )
{
   int e;
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

   for( e = 0; e < g->edges; e++)
   {
      assert(SCIPisGE(scip, g->cost[e], 0));
   }

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
               correct(heap, state, &count, path, m, k, i, ((vbase[k] == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
}

/* 2th next terminal to all non terminal nodes */
void get2next(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   PATH*         path,
   int*          vbase,
   int*          heap,
   int*          state
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
	    correct(heap, state, &count, path, k, i, e, ((root == vbase[i])? cost[e] : costrev[e]), FSP_MODE);
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
            if( vbase[j] != vbase[k] && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) //(state[jc])??
            {
               correct(heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}


/* 3th next terminal to all non terminal nodes */
void get3next(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   PATH*         path,
   int*          vbase,
   int*          heap,
   int*          state
   )
{
   int k;
   int j;
   int i;
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
	 if( !Is_term(g->term[j]) && SCIPisGT(scip, path[k].dist, path[i].dist + ((root == vbase[i])? cost[e] : costrev[e])) &&
	    vbase[i] != vbase[j] && vbase[i] != vbase[j + nnodes] && g->mark[j] )
	 {
	    correct(heap, state, &count, path, k, i, e, ((root == vbase[i])? cost[e] : costrev[e]), FSP_MODE);
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
               && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) //(state[jc])??
            {
               correct(heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}

/*** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals***/
void getnext3terms(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   PATH*         path3,
   int*          vbase,
   int*          heap,
   int*          state
   )
{
   assert(g      != NULL);
   assert(path3   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   /* build voronoi diagram */
   voronoi_terms(scip, g, cost, path3, vbase, heap, state);

   /* get 2nd nearest terms */
   get2next(scip, g, cost, costrev, path3, vbase, heap, state);

   /* get 3th nearest terms */
   get3next(scip, g, cost, costrev, path3, vbase, heap, state);

   return;
}


/*** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals***/
void voronoi_terms(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   PATH*         path,
   int*          vbase,
   int*          heap,
   int*          state
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
      if( Is_term(g->term[i]) && g->mark[i] ) // g->grad[i] > 0
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
               correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}


/*** build a voronoi region, w.r.t. shortest paths, for all terminal and the distance ***/
SCIP_RETCODE voronoi_dist(
   SCIP*         scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   double*       distance,
   int*          vbase,
   int*          minarc,
   int*          heap,
   int*          state,
   int*          distnode,
   PATH*         path
   )
{
   int e;
   int k;
   int m;
   int i;
   int pred;
   int new;
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

   for( e = 0; e < g->edges; e++)
   {
      minedgepred[e] = FALSE;
      assert(g->cost[i] == g->cost[Edge_anti(i)]);
      assert(GE(g->cost[e], 0.0));
   }

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
         // printf("edge: %d \n", path[k].edge);
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

         //printf("get node %d \n ", k);
         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            //    printf("m: %d\n", m);
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
               correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &minedgepred);

   return SCIP_OKAY;
}

/*** build voronoi regions, w.r.t. shortest paths, for all terminals and compute the radii ***/
SCIP_RETCODE voronoi_radius(
   SCIP* scip,
   const GRAPH*  graph,
   GRAPH*        adjgraph,
   PATH*         path,
   SCIP_Real*    rad,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   int*          vbase,
   int*          heap,
   int*          state
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
                  graph_edge_add(adjgraph, nodesid[vbm], nodesid[vbk], ecost, ecost);
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
               correct(heap, state, &count, path, m, k, i, ((vbk == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbk;
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &nodesid);
   return SCIP_OKAY;
}

#if 0
/*** build the voronoi regions for a directed graph. This calculates the inward and outward voronoi regions ***/
void voronoi_inout(
   const GRAPH*   g
   )
{
   int   e;
   int   i;
   int   j;
   int   k;
   int*  heap;
   int*  state;
   int*  inpred;
   int*  outpred;
   int*  terms;
   int* minArc1;
   int* minArc2;
   PATH* inpath;
   double* distance;
   double* radius;
   PATH* outpath;

   distance = malloc((size_t)g->knots * sizeof(double));
   radius = malloc((size_t)g->knots * sizeof(double));
   inpath = malloc((size_t)g->knots * sizeof(PATH));
   outpath = malloc((size_t)g->knots * sizeof(PATH));
   inpred = malloc((size_t)g->edges * sizeof(int));
   outpred = malloc((size_t)g->edges * sizeof(int));
   terms = malloc((size_t)g->terms * sizeof(int));
   minArc1 = malloc((size_t)g->knots * sizeof(int));
   minArc2 = malloc((size_t)g->knots * sizeof(int));
   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   voronoi_term(g, g->cost, distance, radius, inpath, heap, state, inpred, 1);
   voronoi_term(g, g->cost, distance, radius, outpath, heap, state, outpred, 0);

   k = 0;
   for( i = 0; i < g->knots; i++ )
   {
      printf("Node %d - In region: %d, Out region: %d\n", i, g->in_vregion[i], g->out_vregion[i]);
      printf("Node %d - In path dist: %g, Out path dist: %g\n", i, inpath[i].dist, outpath[i].dist);
      printf("Node %d - In pred: %d, Out pred: %d\n", i, inpred[i], outpred[i]);
      if( Is_term(g->term[i]) )
         terms[k++] = i;

      minArc1[i] = -1;
      minArc2[i] = -1;
   }

   for( e = 0; e < g->edges; e++ )
   {
      i = g->head[e];
      j = g->tail[e];
      if( g->in_vregion[i] != g->in_vregion[j] )
      {
         if( minArc1[g->in_vregion[i]] < 0 )
            minArc1[g->in_vregion[i]] = e;
         else if( g->cost[e] < g->cost[minArc1[g->in_vregion[i]]] )
         {
            minArc2[g->in_vregion[i]] = minArc1[g->in_vregion[i]];
            minArc1[g->in_vregion[i]] = e;
         }
      }
   }

   for( i = 0; i < g->terms; i++ )
      printf("Shortest Distance [%d]: %g %g %g\n", i, g->cost[minArc1[terms[i]]], g->cost[minArc2[terms[i]]],
         outpath[g->tail[minArc1[terms[i]]]].dist + g->cost[minArc1[terms[i]]] + inpath[g->head[minArc1[terms[i]]]].dist);

   free(state);
   free(heap);
   free(minArc2);
   free(minArc1);
   free(terms);
   free(outpred);
   free(inpred);
   free(outpath);
   free(inpath);
}
#endif

/*** build a voronoi region, w.r.t. shortest paths, for all terminals ***/
void voronoi_term(
   const GRAPH*   g,
   double*        cost,
   double*        distance,
   double*        radius,
   PATH*          path,
   int*           vregion,
   int*           heap,
   int*           state,
   int*           predecessor,
   int            inward
   )
{
   int k;
   int m;
   int i;
   int curr_edge;
   int nv;
   int nvedge;
   int count = 0;
   int nbases = 0;

   assert(g          != NULL);
   assert(path       != NULL);
   assert(cost       != NULL);
   assert(distance   != NULL);
   assert(radius     != NULL);
   assert(vregion    != NULL);
   assert(heap       != NULL);
   assert(state      != NULL);
   assert(predecessor != NULL);

   if( g->knots == 0 )
      return;

   /* initialize */
   for( i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( g->term[i] >= 0 )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vregion[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
         distance[i] = FARAWAY;
         predecessor[i] = UNKNOWN;
         radius[i] = FARAWAY;
      }
      else
      {
         vregion[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
         distance[i] = FARAWAY;
         predecessor[i] = UNKNOWN;
         radius[i] = FARAWAY;
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

         /* iterate over all ingoing edges of vertex k.
          * working in reverse from the terminal nodes.
          * For the undirected case this does not make any difference,
          * for the directed case, this finds the voronoi region for paths
          * entering a terminal.
          */
         if( inward > 0 )
            nvedge = g->inpbeg[k];
         else
            nvedge = g->outbeg[k];

         for( i = g->inpbeg[k]; i != EAT_LAST; i = g->ieat[i] )
         {

            m = g->tail[i];
            if( inward > 0 )
               curr_edge = i;
            else
               curr_edge = Edge_anti(i);

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) )
            {
               if( GT(path[m].dist, path[k].dist + cost[curr_edge]) )
               {
                  assert(g->term[m] < 0);
                  correct(heap, state, &count, path, m, k, curr_edge, cost[curr_edge], FSP_MODE);
                  predecessor[m] = predecessor[k];
                  vregion[m] = vregion[k];
               }


               if( g->term[k] >= 0 )
               {
                  if( LE(cost[curr_edge], cost[nvedge]) )
                     nvedge = curr_edge;
               }
            }
            else if( state[m] == CONNECT )
            {
               // updating the shortest distance between two terminals. This is used for the nearest vertex test.
               if( vregion[m] != vregion[k] )
               {
                  // this only gives an upper bound on the shortest distance from k to the nearest terminal through the
                  // nearest vertex. The conditions should be checked.
                  if( predecessor[k] == vregion[k] && predecessor[m] != vregion[k] )
                  {
                     assert(predecessor[k] != vregion[m]);
                     assert(predecessor[m] != vregion[k]);

                     assert(g->grad[vregion[m]] > 0);

                     //printf("predecessor[path[k].edge]: %d, predecessor[path[m].edge]: %d, VR[k]: %d, VR[m]: %d, "
                     //"k: %d, path[k].dist: %f, m: %d, path[m].dist: %f, cost[curr_edge]: %f, distance: %f, "
                     //"(%d, %d)\n", predecessor[path[k].edge], predecessor[path[m].edge], vregion[k], vregion[m],
                     //k, path[k].dist, m, path[m].dist, cost[curr_edge], distance[predecessor[path[k].edge]],
                     //g->head[curr_edge], g->tail[curr_edge]);

                     if( predecessor[k] != UNKNOWN )
                        distance[predecessor[k]] = MIN(distance[predecessor[k]],
                           path[k].dist + cost[curr_edge] + path[m].dist);

#if 0
                     if( predecessor[path[m].edge] != UNKNOWN )
                        distance[predecessor[path[m].edge]] = MIN(distance[predecessor[path[m].edge]],
                           path[k].dist + cost[curr_edge] + path[m].dist);
#endif
                  }

                  radius[vregion[k]] = MIN(radius[vregion[k]], path[k].dist + cost[curr_edge]);
                  radius[vregion[m]] = MIN(radius[vregion[m]], path[m].dist + cost[curr_edge]);
               }
            }
         }

         if( g->term[k] >= 0 && g->grad[k] > 0 && nvedge != EAT_LAST )
         {
            if( inward > 0 )
            {
               assert(k == g->head[nvedge]);
               nv = g->tail[nvedge];
            }
            else
            {
               assert(k == g->tail[nvedge]);
               nv = g->head[nvedge];
            }

            predecessor[nv] = k;
         }
      }
   }
}




/*** build a voronoi region for hop constrained problem, w.r.t. shortest paths, for all terminals ***/
/* The voronoi region surronding the root is an outward region and the voronoi region surronding each of the terminals
 * is an inward region. */
void voronoi_hop(
   const GRAPH*   g,
   double*        cost,
   double*        distance,
   double*        radius,
   PATH*          path,
   int*           vregion,
   int*           heap,
   int*           state,
   int*           predecessor,
   int*           radiushops
   )
{
   int k;
   int m;
   int i;
   int curr_edge;
   int nv;
   int nvedge;
   int count = 0;
   int nbases = 0;
   int source;
   double edgecost;

   assert(g          != NULL);
   assert(path       != NULL);
   assert(cost       != NULL);
   assert(distance   != NULL);
   assert(radius     != NULL);
   assert(vregion    != NULL);
   assert(heap       != NULL);
   assert(state      != NULL);
   assert(predecessor != NULL);
   assert(radiushops != NULL);

   if( g->knots == 0 )
      return;

   /* initialize */
   for( i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( g->term[i] >= 0 )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vregion[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         path[i].hops = 0;
         state[i] = count;
         distance[i] = FARAWAY;
         predecessor[i] = UNKNOWN;
         radius[i] = FARAWAY;
         radiushops[i] = 0;
      }
      else
      {
         vregion[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         path[i].hops = -1;
         state[i]     = UNKNOWN;
         distance[i] = FARAWAY;
         predecessor[i] = UNKNOWN;
         radius[i] = FARAWAY;
         radiushops[i] = 0;
      }

   }
   assert(nbases > 0);

   source = g->source[0];

   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all ingoing edges of vertex k.
          * We traverse the graph in two different directions, out from the source and in to the terminals.
          * This is required for the hop constained problems because each of the terminals are leaves of the steiner
          * tree.
          */
         if( vregion[k] == source )
            nvedge = g->outbeg[k];
         else
            nvedge = g->inpbeg[k];

         for( i = g->inpbeg[k]; i != EAT_LAST; i = g->ieat[i] )
         {
            m = g->tail[i];
            if( vregion[k] == source )
               curr_edge = Edge_anti(i);
            else
               curr_edge = i;

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) )
            {
               if( GT(path[m].dist, path[k].dist + cost[curr_edge]) )
               {
                  assert(g->term[m] < 0);
                  correct(heap, state, &count, path, m, k, curr_edge, cost[curr_edge], FSP_MODE);
                  predecessor[m] = predecessor[k];
                  vregion[m] = vregion[k];
               }


               if( g->term[k] >= 0 )
               {
                  if( LE(cost[curr_edge], cost[nvedge]) )
                     nvedge = curr_edge;
               }
            }
            else if( state[m] == CONNECT )
            {
               // updating the shortest distance between two terminals. This is used for the nearest vertex test.
               if( vregion[m] != vregion[k] )
               {
                  // this only gives an upper bound on the shortest distance from k to the nearest terminal through the
                  // nearest vertex. The conditions should be checked.
                  if( predecessor[k] == vregion[k] && predecessor[m] != vregion[k] )
                  {
                     assert(predecessor[k] != vregion[m]);
                     assert(predecessor[m] != vregion[k]);

                     assert(g->grad[vregion[m]] > 0);

                     //printf("predecessor[path[k].edge]: %d, predecessor[path[m].edge]: %d, VR[k]: %d, VR[m]: %d, "
                     //"k: %d, path[k].dist: %f, m: %d, path[m].dist: %f, cost[curr_edge]: %f, distance: %f, "
                     //"(%d, %d)\n", predecessor[path[k].edge], predecessor[path[m].edge], vregion[k], vregion[m],
                     //k, path[k].dist, m, path[m].dist, cost[curr_edge], distance[predecessor[path[k].edge]],
                     //g->head[curr_edge], g->tail[curr_edge]);

                     if( predecessor[k] != UNKNOWN )
                        distance[predecessor[k]] = MIN(distance[predecessor[k]],
                           path[k].dist + cost[curr_edge] + path[m].dist);

#if 0
                     if( predecessor[path[m].edge] != UNKNOWN )
                        distance[predecessor[path[m].edge]] = MIN(distance[predecessor[path[m].edge]],
                           path[k].dist + cost[curr_edge] + path[m].dist);
#endif
                  }

                  /* For the directed case, which is the situation for the hop constrained problems, the radius is a
                   * little difficult to calculate. In the situation where the edge (k, m) crosses a vregion boundary
                   * and m is a terminal, the cost of the edge will be FARAWAY. This is because all terminals are leaves
                   * in this problem type. So the radius must take into account this different cost and add it to the
                   * distance from the vregion[k]. */
                  if( Is_term(g->term[m]) && !LT(cost[curr_edge], FARAWAY) )
                     edgecost = cost[Edge_anti(curr_edge)];
                  else
                     edgecost = cost[curr_edge];

                  if( GT(radius[vregion[k]], path[k].dist + edgecost) ||
                     (EQ(radius[vregion[k]], path[k].dist + edgecost) &&
                        GT(radiushops[vregion[k]], path[k].hops + 1)) )
                  {
                     radius[vregion[k]] = path[k].dist + edgecost;
                     radiushops[vregion[k]] = path[k].hops + 1;
                  }

                  if( GT(radius[vregion[m]], path[m].dist + edgecost) ||
                     (EQ(radius[vregion[m]], path[m].dist + edgecost) &&
                        GT(radiushops[vregion[m]], path[m].hops + 1)) )
                  {
                     radius[vregion[m]] = path[m].dist + edgecost;
                     radiushops[vregion[m]] = path[m].hops + 1;
                  }
               }
            }
         }

         if( g->term[k] >= 0 && g->grad[k] > 0 && nvedge != EAT_LAST )
         {
            if( vregion[k] == source )
            {
               assert(k == g->tail[nvedge]);
               nv = g->head[nvedge];
            }
            else
            {
               assert(k == g->head[nvedge]);
               nv = g->tail[nvedge];
            }

            predecessor[nv] = k;
         }
      }
   }
}

/*** repair the voronoi diagram for SL reduction ***/
void voronoi_slrepair(
   SCIP* scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   PATH*         path,
   int*          vbase,
   int*          heap,
   int*          state,
   int           start,
   int           adjnode
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
               correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}

/*** repair the voronoi diagram for a given set nodes ***/
void voronoi_repair(
   SCIP*         scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   int*          count,
   int*          vbase,
   PATH*         path,
   int*          newedge,
   int           crucnode,
   UF*           uf
   )
{
   int k;
   int m;
   int i;
   int* heap;
   int* state;
   int node1;
   int node2;

   *newedge = UNKNOWN;
   assert(g != NULL);
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);
   assert(path != NULL);
   assert(cost != NULL);

   if (g->knots == 0)
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
               correct(heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }

            /* check whether there is a better new boundary edge adjacent to vertex k */
            else if( (state[m] == CONNECT) && (((node1 = SCIPunionfindFind(uf, vbase[m])) == crucnode) ^ ((node2 = SCIPunionfindFind(uf, vbase[k])) == crucnode))
               && g->mark[m] && g->mark[vbase[m]]
               &&(g->mark[node1]) && (g->mark[node2]) && (SCIPisGT(scip, (*newedge == UNKNOWN)? FARAWAY :
                     (path[g->tail[*newedge]].dist + cost[*newedge] + path[g->head[*newedge]].dist), path[k].dist + cost[i] + path[m].dist) ) )
	    {
	       if ( !g->mark[m])
	       {
		  printf("notnode: %d \n", m);
		  assert(0);
	       }
               *newedge = i;
	    }
         }
      }
   }
}


/*** repair the voronoi diagram for a given set nodes ***/
void voronoi_repair_mult(
   SCIP*         scip,
   const GRAPH*  g,
   SCIP_Real*    cost,
   int*          count,
   int*          vbase,
   int*          boundedges,
   int*          nboundedges,
   char*         nodesmark,
   UF*           uf,
   PATH*         path
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
               correct(heap, state, count, path, m, k, i, cost[i], FSP_MODE);
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
   const GRAPH* g,
   const PATH*  p)
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
