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
inline static void correct2(
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


/* Find a tree conntaining all terminals */
void graph_path_exec2(
   const GRAPH*  p,
   int           mode,
   int           start,
   const SCIP_Real* cost,
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
   if( p->knots == 0 )
      return;

   heap = p->path_heap;
   state = p->path_state;
   count = 0;

   /* Erstmal alles auf null, unbekannt und weit weg
    */
   for( i = 0; i < p->knots; i++ )
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
   if ( p->knots > 1 )
   {
      int termsn = 0;
      int node;
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
         path[i].dist = FARAWAY;
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
         //printf("get node %d \n ", k);
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
