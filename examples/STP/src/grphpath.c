/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: grphpath.c                                                    */
/*   Name....: Find Shortest Path / Minimum Spanning Tree                    */
/*   Author..: Thorsten Koch                                                 */
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

#include <stdio.h>
#include <assert.h>
#include "portab.h"
#include <stdlib.h>

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
   double cost,
   int    mode)
{
   int    t;
   int    c;
   int    j;

   path[l].dist = (mode == MST_MODE) ? cost : (path[k].dist + cost);
   path[l].edge = e;

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
inline static void correct2(
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   PATH*  path,
   int    l,
   int    k,
   int    e,
   double cost,
   int    mode)
{
   int    t;
   int    c;
   int    j;

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

   p->path_heap  = malloc((size_t)p->knots * sizeof(int));
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
   const double* cost,
   PATH*         path)
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count;

   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) || (mode == MST_MODE));
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
   count = 0;

   /* Erstmal alles auf null, unbekannt und weit weg
    */
   for(i = 0; i < p->knots; i++)
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY;
      path[i].edge = -1;
   }
   /* Startknoten in den Heap
    */
   k            = start;
   path[k].dist = 0.0;

   /* Wenn nur ein Knoten drin ist funktioniert der Heap nicht sehr gut,
    * weil dann fuer genau 0 Elemente Platz ist.
    */
   if (p->knots > 1)
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* Wenn nichts mehr auf dem Heap ist, sind wir fertig
       * und jetzt erstmal Hula Loop
       */
      while(count > 0)
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
            if ((state[m]) && (p->mark[m])

               /* 2. Ist es ueberhaupt eine Verbesserung diesen Weg zu nehmen ?
                *    Wenn ja, dann muss die Entferung kuerzer sein.
                */
               && (GT(path[m].dist, (mode == MST_MODE) ? cost[i] : (path[k].dist + cost[i]))))
               correct(heap, state, &count, path, m, k, i, cost[i], mode);
         }
      }
   }
}


/* Find a tree conntaining all terminals */
void graph_path_exec2(
   const GRAPH*  p,
   int           mode,
   int           start,
   const double* cost,
   PATH*         path,
   char*         connected,
   int* cluster,
   int* csize
   )
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count;

   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) || (mode == MST_MODE));
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);
   /* Kein Baum ohne Knoten
    */
   if (p->knots == 0)
      return;

   heap = p->path_heap;
   state = p->path_state;
   count = 0;

   /* Erstmal alles auf null, unbekannt und weit weg
    */
   for(i = 0; i < p->knots; i++)
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY;
      path[i].edge = -1;
      connected[i] = FALSE;
   }
   /* Startknoten in den Heap
    */
   k            = start;
   path[k].dist = 0.0;
   connected[k] = TRUE;

   /* Wenn nur ein Knoten drin ist funktioniert der Heap nicht sehr gut,
    * weil dann fuer genau 0 Elemente Platz ist.
    */
   if (p->knots > 1)
   {
      int termsn = 0;
      int node;
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      /* Wenn nichts mehr auf dem Heap ist, sind wir fertig
       * und jetzt erstmal Hula Loop
       */
      while(count > 0)
      {
         /* Na, wer ist der Naechste ?
          */
         k = nearest(heap, state, &count, path);

         /* Wieder einen erledigt
          */
         state[k] = CONNECT;

         if( p->term[k] != -1 )
         {
            /* add path to the terminal */


            ++termsn;
            connected[k] = TRUE;

            if (cluster[*csize-1] != k)

               cluster[(*csize)++] = k;
            path[k].dist = 0.0;


	    node = k;
            while( path[k].edge != -1 && connected[node = p->tail[path[node].edge]] == FALSE )
            {
               connected[node] = TRUE;
	       if (cluster[*csize-1] != k)
                  cluster[(*csize)++] = node;
               path[node].dist = 0.0;
               state[node] = UNKNOWN;
               correct2(heap, state, &count, path, node, 0, 0, 0, mode);
            }
            /* have all terminals been reached? */
            if( termsn == p->terms )
            {
               break;
            }
         }

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
            if( (state[m]) /*  && (p->mark[m]) */

               /* 2. Ist es ueberhaupt eine Verbesserung diesen Weg zu nehmen ?
                *    Wenn ja, dann muss die Entferung kuerzer sein.
                */
               && (GT(path[m].dist, (path[k].dist + cost[i]))) )
               correct(heap, state, &count, path, m, k, i, cost[i], mode);
         }
      }
   }
}


/*** build a voronoi region, w.r.t. shortest paths, for a given set of bases ***/
void voronoi(
   const GRAPH*  p,
   const double* cost,
   char*         base,
   int*          vbase,
   PATH*         path
   )
{
   int k;
   int m;
   int i;
   int* heap;
   int* state;
   int count = 0;
   int nbases = 0;

   assert(p      != NULL);
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   if( p->knots == 0 )
      return;

   heap = p->path_heap;
   state = p->path_state;

   /* init */
   for( i = 0; i < p->knots; i++ )
   {

      /* set the base of vertex i */
      if( base[i] )
      {
         nbases++;
         if( p->knots > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0;
         path[i].edge = UNKNOWN;
         state[i] = CONNECT; /* = CONNECT? TODO */
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }

   }

   /* add the starting vertex to the heap */
   assert(nbases > 0);

   if( p->knots > 1 && nbases < p->knots )
   {

      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex from the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = p->outbeg[k]; i != EAT_LAST; i = p->oeat[i] )
         {
            m = p->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && (GT(path[m].dist, path[k].dist + cost[i])) ) /* TODO: && mark[m] ?? */
            {
               correct(heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
}
/*** repair the voronoi diagram for a given set nodes ***/
void voronoi_repair(
   SCIP*         scip,
   const GRAPH*  g,
   const double* cost,
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

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && (SCIPisGT(scip,path[m].dist, path[k].dist + cost[i])) ) /* TODO: && mark[m] ?? */
            {
               correct(heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
            /* check whether there is a better new boundary edge adjacent to vertex k */
            else if( (state[m] == CONNECT) && ((UF_find(uf, vbase[m]) == crucnode) ^ (UF_find(uf, vbase[k]) == crucnode))  && (SCIPisGT(scip, (*newedge == -1)? FARAWAY :
                     (path[g->tail[*newedge]].dist + cost[*newedge] + path[g->head[*newedge]].dist), path[k].dist + cost[i] + path[m].dist) ) )
	    {
               *newedge = i;
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
   double dist = 0.0;
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
