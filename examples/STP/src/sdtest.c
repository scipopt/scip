/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: sdtest.c                                                      */
/*   Name....: Special Distance Test                                         */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"

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

typedef struct sd_path
{
   double dist;
   double tran;
} SDPTH;


static int compare(
   const SDPTH* path,
   int          a,
   int          b)
{
   if (NE(path[a].dist, path[b].dist))
      return(LT(path[a].dist, path[b].dist) ? -1 : 1);

   if (EQ(path[a].tran, path[b].tran))
      return(0);

   return(LT(path[a].tran, path[b].tran) ? -1 : 1);
}

static int issmaller(
   const SDPTH* path,
   int          a,
   int          b)
{
   return (LT(path[a].dist, path[b].dist) || (!GT(path[a].dist, path[b].dist) && LT(path[a].tran, path[b].tran)));
}
static int islarger(
   const SDPTH* path,
   int          a,
   int          b)
{
   return (GT(path[a].dist, path[b].dist) || (!LT(path[a].dist, path[b].dist) && GT(path[a].tran, path[b].tran)));
}
/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearest(
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   const SDPTH* path)
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
      if (islarger(path, heap[2], heap[3]))
         c++;

   while((c <= (*count)) && islarger(path, heap[j], heap[c]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
         if (issmaller(path, heap[c + 1], heap[c]))
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
/*---            den neuen Knoten erreicht, Kosten der Kante.             ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
inline static void correct(
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   SDPTH* path,
   int    l)
{
   int   t;
   int   c;
   int   j;

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

   while((j > 1) && (islarger(path, heap[c], heap[j])))
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
static void compute_sd(
   const GRAPH*  p,
   int           start,
   const double* cost,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   SDPTH*        path)
{
   int    k;
   int    m;
   int    i;
   int    done = 0;
   double tran;
   double dist;

   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(count != NULL);
   assert(*count >= 0);

   /* Kein Baum ohne Knoten
    */
   if (p->knots == 0)
      return;

   (*count) = 0;

   /* Erstmal alles auf null, unbekannt und weit weg
    */
   for(i = 0; i < p->knots; i++)
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY;
      path[i].tran = FARAWAY;
   }
   /* Startknoten in den Heap
    */
   k            = start;
   path[k].dist = 0.0;
   path[k].tran = 0.0;

   /* Wenn nur ein Knoten drin ist funktioniert der Heap nicht sehr gut,
    * weil dann fuer genau 0 Elemente Platz ist.
    */
   if (p->knots > 1)
   {
      (*count)       = 1;
      heap[(*count)] = k;
      state[k]    = (*count);

      /* Wenn nichts mehr auf dem Heap ist, sind wir fertig
       * und jetzt erstmal Hula Loop
       */
      while((*count) > 0)
      {
         /* Na, wer ist der Naechste ?
          */
         k = nearest(heap, state, count, path);

         /* Wieder einen erledigt
          */
         state[k] = CONNECT;

         if (p->mark[k] == 2)
            if (++done >= p->grad[start])
               break;

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
               /* The special distance is the length of the longest path between two terminals (elementary path)
                * contained in the path between knots i and j.
                * - tran measures the distance between two terminals.
                * - dist stores the current longest elementary path.
                */
               tran = Is_term(p->term[m]) ? 0.0 : path[k].tran + cost[i];
               dist = Max(path[k].dist, path[k].tran + cost[i]);

               if (LT(dist, path[m].dist)
                || (EQ(dist, path[m].dist) && LT(tran, path[m].tran)))
               {
                  path[m].dist = dist;
                  path[m].tran = tran;

                  correct(heap, state, count, path, m);
               }
            }
         }
      }
   }
}

/* C. W. Duin and A. Volganant
 *
 * "An Edge Elimination Test for the Steiner Problem in Graphs"
 *
 * Operations Research Letters 8, (1989), 79-83
 *
 * Special Distance Test
 */
int sd_reduction(
   GRAPH* g)
{
   SDPTH*  sd;
   double* cost;
   int*    heap;
   int*    state;
   int     count = 0;
   int     i;
   int     e;
   int     j;
   int     elimins = 0;

   SCIPdebugMessage("SD-Reduktion: ");
   fflush(stdout);

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   assert(heap  != NULL);
   assert(state != NULL);

   sd = malloc((size_t)g->knots * sizeof(SDPTH));

   assert(sd != NULL);

   cost  = malloc((size_t)g->edges * sizeof(double));

   assert(cost != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = (g->grad[i] > 0);

   for(i = 0; i < g->edges; i++)
      cost[i] = g->cost[i] * 1000.0 + (double)(rand() % 512);

   for(i = 0; i < g->knots; i++)
   {
      if (!(i % 100))
      {
         SCIPdebug(fputc('.', stdout));
         SCIPdebug(fflush(stdout));
      }
      if (g->grad[i] == 0)
         continue;

      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ( (g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source[0] )
         continue;

      for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
      {
         assert(g->mark[g->head[e]] == 1);

         g->mark[g->head[e]] = 2;
      }

      compute_sd(g, i, cost, heap, state, &count, sd);

      for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
      {
         assert(g->mark[g->head[e]] == 2);
         /* assert(sd[g->head[e]].dist < FARAWAY); */

         g->mark[g->head[e]] = 1;
      }

      for(e = g->outbeg[i]; e != EAT_LAST; e = j)
      {
         assert(g->tail[e] == i);

         j = g->oeat[e];

         if (LT(cost[e], FARAWAY) && LT(sd[g->head[e]].dist, cost[e]))
         {
            graph_edge_del(g, e);

            elimins++;
         }
      }
   }
   free(heap);
   free(state);

   heap  = NULL;
   state = NULL;

   free(sd);
   free(cost);

   assert(graph_valid(g));

   SCIPdebugMessage("%d Edges deleted\n", elimins * 2);

   return(elimins);
}

static void redirect_edge(
   GRAPH* g,
   int    eki,
   int    k,
   int    j,
   double cost)
{
   int e;

   graph_edge_del(g, eki);

   for(e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e])
      if ((g->tail[e] == k) && (g->head[e] == j))
         break;

   /* Gibt es die Kante schon ?
    */
   if (e != EAT_LAST)
   {
      /* Ja, dann nur Kosten korrigieren.
       */
      if (GT(g->cost[e], cost))
      {
         g->cost[e]            = cost;
         g->cost[Edge_anti(e)] = cost;
      }
   }
   else
   {
      assert(g->oeat[eki] == EAT_FREE);

      e = eki;

      g->grad[k]++;
      g->grad[j]++;

      g->cost[e]   = cost;
      g->head[e]   = j;
      g->tail[e]   = k;
      g->ieat[e]   = g->inpbeg[j];
      g->oeat[e]   = g->outbeg[k];
      g->inpbeg[j] = e;
      g->outbeg[k] = e;

      e = Edge_anti(eki);

      g->cost[e]   = cost;
      g->head[e]   = k;
      g->tail[e]   = j;
      g->ieat[e]   = g->inpbeg[k];
      g->oeat[e]   = g->outbeg[j];
      g->inpbeg[k] = e;
      g->outbeg[j] = e;
   }
}

/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Bottleneck Degree 3 Test
 */
int bd3_reduction(
   GRAPH* g)
{
   SDPTH* path1;
   SDPTH* path2;
   int*    heap;
   int*    state;
   int     count = 0;
   int    i;
   int    k1;
   int    k2;
   int    k3;
   int    e1;
   int    e2;
   int    e3;
   int    elimins = 0;
   double c1;
   double c2;
   double c3;
   double c123;

   SCIPdebugMessage("BD3-Reduction: ");
   fflush(stdout);

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   assert(heap  != NULL);
   assert(state != NULL);

   path1 = malloc((size_t)g->knots * sizeof(SDPTH));
   path2 = malloc((size_t)g->knots * sizeof(SDPTH));

   assert(path1 != NULL);
   assert(path2 != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = (g->grad[i] > 0);

   for(i = 0; i < g->knots; i++)
   {
      if (!(i % 100))
      {
         SCIPdebug(fputc('.', stdout));
         SCIPdebug(fflush(stdout));
      }
      if (g->grad[i] != 3)
         continue;

      if (Is_term(g->term[i]))
         continue;

      e1 = g->outbeg[i];

      assert(e1 != EAT_LAST);

      k1 = g->head[e1];
      c1 = g->cost[e1];
      e2 = g->oeat[e1];

      assert(e2 != EAT_LAST);

      k2 = g->head[e2];
      c2 = g->cost[e2];
      e3 = g->oeat[e2];

      assert(e3 != EAT_LAST);

      k3 = g->head[e3];
      c3 = g->cost[e3];

      assert(g->oeat[e3] == EAT_LAST);

      compute_sd(g, k1, g->cost, heap, state, &count, path1);
      compute_sd(g, k2, g->cost, heap, state, &count, path2);
#if 0
      graph_path_exec(g, FSP_MODE, k1, g->cost, path1);
      graph_path_exec(g, FSP_MODE, k2, g->cost, path2);
#endif
      c123 = c1 + c2 + c3;

      if (GT(path1[k2].dist + path1[k3].dist, c123)
       && GT(path1[k2].dist + path2[k3].dist, c123)
       && GT(path1[k3].dist + path2[k3].dist, c123))
         continue;

      elimins++;

      if (LT(path1[k2].dist, c1 + c2))
         graph_edge_del(g, e1);
      else
         redirect_edge(g, e1, k1, k2, c1 + c2);

      if (LT(path2[k3].dist, c2 + c3))
         graph_edge_del(g, e2);
      else
         redirect_edge(g, e2, k2, k3, c2 + c3);

      if (LT(path1[k3].dist, c1 + c3))
         graph_edge_del(g, e3);
      else
         redirect_edge(g, e3, k3, k1, c3 + c1);

      assert(g->grad[i] == 0);
   }
   free(heap);
   free(state);

   heap  = NULL;
   state = NULL;

   free(path1);
   free(path2);

   assert(graph_valid(g));

   SCIPdebugMessage("%d Knots deleted\n", elimins);

   return(elimins);
}


/* computes the shortest path from each terminal to every other vertex */
static void calculate_distances(
   const GRAPH* g,
   PATH** path)
{
   int i;

   SCIPdebug(fputc('C', stdout));
   SCIPdebug(fflush(stdout));

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]) && (g->grad[i] > 0))
      {
         if (path[i] == NULL)
            path[i] = malloc((size_t)g->knots * sizeof(PATH));

         assert(path[i] != NULL);

         graph_path_exec(g, FSP_MODE, i, g->cost, path[i]);
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

inline static double mst_cost(
   const GRAPH* g,
   const PATH*  mst)
{
   double cost = 0;
   int    i;
   int    e;

   for(i = 0; i < g->knots; i++)
      if ((e = mst[i].edge) >= 0)
         cost += g->cost[e];

   return(cost);
}

/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Nearest Special Vertex 3 Test
 */
int nsv_reduction(
   GRAPH*  g,
   double* fixed)
{
   PATH**  path;
   PATH*   mst1;
   PATH*   mst2;
   double* cost;
   int     i;
   int     e;
   int     k;
   int     j;
   double  min1;
   double  min2;
   double  cost1;
   double  cost2;
   int     elimins = 0;

   SCIPdebugMessage("NSV-Reduction: ");
   fflush(stdout);
/*
   graph_show(g);
*/
   path = malloc((size_t)g->knots * sizeof(PATH*));

   assert(path != NULL);

   mst1 = malloc((size_t)g->knots * sizeof(PATH));
   mst2 = malloc((size_t)g->knots * sizeof(PATH));

   assert(mst1 != NULL);
   assert(mst2 != NULL);

   cost = malloc((size_t)g->edges * sizeof(double));

   assert(cost != NULL);

   /* Check this cost setting. It may need to be changed for the directed case.
    */
   for(i = 0; i < g->edges; i++)
      cost[i] = g->cost[i];

   for(i = 0; i < g->knots; i++)
   {
      g->mark[i] = (g->grad[i] > 0);
      path[i] = NULL;
   }

   calculate_distances(g, path);

   for(i = 0; i < g->knots; i++)
   {
      if (!(i % 100))
      {
         SCIPdebug(fputc('.', stdout));
         SCIPdebug(fflush(stdout));
      }
      if (g->grad[i] < 3)
         continue;

      graph_path_exec(g, MST_MODE, i, g->cost, mst1);

      cost1 = mst_cost(g, mst1);

      for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
      {
         assert(g->tail[e] == i);

         if (mst1[g->head[e]].edge == e)
            break;
      }
      assert(e != EAT_LAST);

      cost[e]            = FARAWAY - 1.0;
      cost[Edge_anti(e)] = FARAWAY - 1.0;

      graph_path_exec(g, MST_MODE, i, cost, mst2);

      cost2 = mst_cost(g, mst2);

      cost[e]            = g->cost[e];
      cost[Edge_anti(e)] = g->cost[Edge_anti(e)];

      assert(GE(cost2, cost1));

      if (LE(cost2 - cost1, 2.0))
      {
/*
         SCIPdebugMessage("\t\te=%d i=%d k=%d cost1=%d cost2=%d\n",
            e, i, g->head[e], cost1, cost2);
*/
         continue;
      }
      k     = g->head[e];
      min1  = FARAWAY;
      min2  = FARAWAY;

      for(j = 0; j < g->knots; j++)
      {
         if (!Is_term(g->term[j]) || (g->grad[j] == 0))
            continue;

         assert(path[j] != NULL);

         if (LT(path[j][i].dist, min1) && LT(path[j][i].dist, path[j][k].dist))
            min1  = path[j][i].dist;

         if (LT(path[j][k].dist, min2) && LT(path[j][k].dist, path[j][i].dist))
            min2  = path[j][k].dist;
      }
      if (EQ(min1, FARAWAY) || EQ(min2, FARAWAY))
      {
/*
         SCIPdebugMessage("\te=%d i=%d k=%d min1=%d min2=%d cost1=%d cost2=%d\n",
            e, i, k, min1, min2, cost1, cost2);
*/
         continue;
      }
      if (LT(cost1 + min1 + min2, cost2))
      {
/*
         SCIPdebugMessage("e=%d i=%d k=%d min1=%d min2=%d cost1=%d cost2=%d\n",
            e, i, k, min1, min2, cost1, cost2);
 */
         *fixed += g->cost[e];

         graph_knot_contract(g, i, k);

         elimins++;

         calculate_distances(g, path);

         for(j = 0; j < g->edges; j++)
            cost[j] = g->cost[j];
      }
   }
   for(i = 0; i < g->knots; i++)
   {
      if (path[i] != NULL)
      {
         assert(Is_term(g->term[i]));

         free(path[i]);
      }
   }
   free(mst1);
   free(mst2);
   free(path);
   free(cost);

   assert(graph_valid(g));

   SCIPdebugMessage(" %d Knots deleted\n", elimins);

   return(elimins);
}



/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Nearest Special Vertex 3 Test
 *
 * Taken from:
 *
 * Maculan et. al.
 *
 * "An approach for the Steiner problem in directed graphs"
 *
 * Annals of Operations Research, Volume 33 (1991), 471-480
 *
 * Nearest Vertex Test (for optimal arcs)
 *
 * and
 *
 * T. Polzin
 *
 * "Algorithms for the Steiner problem in networks"
 *
 * Section 3.3.3 pp. 54-55
 *
 * This is only called for the directed Steiner tree problem
 */
int nv_reduction_optimal(
   GRAPH*  g,
   double* fixed)
{
   PATH**  path;
   PATH*   pathfromterm;
   PATH*   pathfromsource;
   int*    vregion;
   int*    heap;
   int*    state;
   int*    pred;
   int*    minArc1;
   int*    minArc2;
   int*    terms;
   int     termcount;
   int     i;
   int     j;
   int     k;
   int     e;
   double  min1;
   double  min2;
   int     shortarc1;
   int     elimins = 0;

   SCIPdebugMessage("NSV-Reduction: ");
   fflush(stdout);
/*
   graph_show(g);
*/
   path = malloc((size_t)g->knots * sizeof(PATH*));

   assert(path != NULL);

   pathfromterm = malloc((size_t)g->knots * sizeof(PATH));
   pathfromsource = malloc((size_t)g->knots * sizeof(PATH));

   assert(pathfromterm != NULL);
   assert(pathfromsource != NULL);

   vregion = malloc((size_t)g->knots * sizeof(int));

   assert(vregion != NULL);

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   assert(heap != NULL);
   assert(state != NULL);

   pred = malloc((size_t)g->knots * sizeof(int));
   minArc1 = malloc((size_t)g->knots * sizeof(int));
   minArc2 = malloc((size_t)g->knots * sizeof(int));
   terms = malloc((size_t)g->terms * sizeof(int));

   termcount = 0;
   for(i = 0; i < g->knots; i++)
   {
      if( Is_term(g->term[i]) )
      {
         terms[termcount] = i;
         termcount++;
      }
      g->mark[i] = (g->grad[i] > 0);
      minArc1[i] = -1;
      minArc2[i] = -1;
      path[i] = NULL;
   }

   assert(g->source[0] >= 0);

   /* computing the voronoi regions inward to a node */
   voronoi_term(g, g->cost, pathfromterm, vregion, heap, state, pred, 1);

   /* computing the shortest paths from the source node */
   graph_path_exec(g, FSP_MODE, g->source[0], g->cost, pathfromsource);

   /* computing the shortest paths from each terminal to every other node */
   //calculate_distances(g, path);

   for(i = 0; i < g->knots; i++)
   {
      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ((g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source[0] )
         continue;

      if (Is_term(g->term[i]) && g->grad[i] >= 3)
      {

         min1  = FARAWAY;
         min2  = FARAWAY;
         shortarc1 = -1;
         for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
         {
            if (g->cost[e] < min1)
            {
               shortarc1 = g->tail[e];

               min2 = min1;
               min1 = g->cost[e];
            }
         }

         if (LT(min1, FARAWAY) && LE(pathfromsource[shortarc1].dist + min1, min2))
         {
            *fixed += g->cost[shortarc1]; /* ????? this was e, but should probably be shortarc1 */

            graph_knot_contract(g, shortarc1, i);

            elimins++;

            /* computing the shortest paths from the source node */
            graph_path_exec(g, FSP_MODE, g->source[0], g->cost, pathfromsource);
            //calculate_distances(g, path);
         }
      }
      /* The knot is not a terminal so we can perform the short link test */
      else
      {
         for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
         {
            j = g->tail[e];
            if( vregion[i] != vregion[j] )
            {
               if( minArc1[vregion[i]] < 0 )
                  minArc1[vregion[i]] = e;
               else if( g->cost[e] < g->cost[minArc1[vregion[i]]] )
               {
                  minArc2[vregion[i]] = minArc1[vregion[i]];
                  minArc1[vregion[i]] = e;
               }
            }
         }
      }
   }

   for( k = 0; k < termcount; k++ )
   {
      assert(terms[k] >= 0 && terms[k] < g->knots);

      if( minArc1[terms[k]] >= 0 && minArc2[terms[k]] >= 0 && pathfromsource[g->tail[minArc1[terms[k]]]].dist
            + g->cost[minArc1[terms[k]]] + pathfromterm[g->head[minArc1[terms[k]]]].dist < g->cost[minArc2[terms[k]]] )
      {
         printf("Short link test %d\n", k);
         e = minArc1[terms[k]];
         i = g->head[e];
         j = g->tail[e];

         if( Is_term(g->term[i]) )
            *fixed += g->cost[e];

         graph_knot_contract(g, j, i);

         elimins++;
      }
   }

   free(terms);
   free(minArc2);
   free(minArc1);
   free(pred);
   free(vregion);
   free(pathfromsource);
   free(pathfromterm);
   free(path);

   assert(graph_valid(g));

   SCIPdebugMessage(" %d Knots deleted\n", elimins);

   return(elimins);
}
