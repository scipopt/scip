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

/**@file   reduce_alt.c
 * @brief  Altenative based reduction tests for Steiner problems
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 *
 * This file implements alternative-based reduction techniques for several Steiner problems.
 * All tests can be found in "Combining NP-Hard Reduction Techniques and Strong Heuristics in an Exact Algorithm for the
 * Maximum-Weight Connected Subgraph Problem" by Daniel Rehfeldt and Thorsten Koch,
 * or in "Reduction Techniques for the Prize-Collecting Steiner Tree Problem and the Maximum-Weight Connected Subgraph Problem"
 * by Daniel Rehfeldt et al.
 *
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "misc_stp.h"
#include "probdata_stp.h"
#include "scip/scip.h"

#define VERTEX_CONNECT      0
#define VERTEX_TEMPNEIGHBOR 1
#define VERTEX_NEIGHBOR     2
#define VERTEX_OTHER        3
#define STP_RED_CNSNN       25
#define STP_RED_ANSMAXCANDS 500
#define STP_RED_ANSEXMAXCANDS 50
#define STP_RED_ANSMAXNEIGHBORS 25


/* can edge be deleted in SD test in case of equality? If so, 'forbidden' array is adapted */
static
SCIP_Bool sddeltable(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 path,
   int*                  vbase,
   int*                  forbidden,
   int                   tpos,
   int                   hpos,
   int                   tail,
   int                   head,
   int                   edge,
   int                   nnodes
   )
{
   int e;

   assert(g != NULL);
   assert(path != NULL);
   assert(scip != NULL);
   assert(vbase != NULL);
   assert(forbidden != NULL);

   assert(tpos >= 0);
   assert(hpos >= 0);
   assert(tail >= 0);
   assert(head >= 0);
   assert(edge >= 0);
   assert(nnodes >= 0);

   /* edge between non-terminals */
   if( !Is_term(g->term[tail]) && !Is_term(g->term[head]) )
      return TRUE;

   /* has edge been marked as forbidden? */
   if( forbidden[edge] )
      return FALSE;

   /* edge between a terminal and a non terminal */
   if( !Is_term(g->term[tail]) || !Is_term(g->term[head]) )
   {

      /* check whether edge is used in shortest path */

      if( !Is_term(g->term[tail]) && Is_term(g->term[head]) )
      {
         e = path[tail + tpos * nnodes].edge;

         if( g->ieat[e] == EAT_FREE )
            return TRUE;

         assert(g->head[e] == tail);

         if( g->tail[e] == head )
            return FALSE;
      }
      else if( Is_term(g->term[tail]) && !Is_term(g->term[head]) )
      {
         e = path[head + hpos * nnodes].edge;

         if( g->ieat[e] == EAT_FREE )
            return TRUE;

         assert(g->head[e] == head);

         if( g->tail[e] == tail )
            return FALSE;
      }
   }

   /* update forbidden edges todo check bottleneck distance between terminals, and don't delete if distance is higher than ecost */

   if( Is_term(g->term[head]) )
   {
      SCIP_Real ecost = g->cost[edge];
      for( e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(e >= 0);

         if( SCIPisEQ(scip, g->cost[e], ecost) )
         {
            if( !(forbidden[e]) && Is_term(g->term[g->head[e]]) )
            {
               forbidden[e] = TRUE;
               forbidden[flipedge(e)] = TRUE;
            }
         }
      }
   }

   if( Is_term(g->term[tail]) )
   {
      SCIP_Real ecost = g->cost[edge];
      for( e = g->outbeg[head]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(e >= 0);

         if( SCIPisEQ(scip, g->cost[e], ecost) )
         {
            if( !(forbidden[e]) && Is_term(g->term[g->head[e]]) )
            {
               forbidden[e] = TRUE;
               forbidden[flipedge(e)] = TRUE;
            }
         }
      }
   }

   return TRUE;
}


static
int getcloseterms(
   SCIP*                 scip,
   const PATH*           vnoi,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   const int*            vbase,
   int*                  neighbterms,
   int                   i,
   int                   nnodes
   )
{
   int k;
   int pos = i;
   int nnterms = 0;

   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(vbase != NULL);
   assert(termdist != NULL);
   assert(neighbterms != NULL);

   for( k = 0; k < 4; k++ )
   {
      if( SCIPisLT(scip, vnoi[pos].dist, ecost) )
      {
         neighbterms[nnterms] = vbase[pos];
         termdist[nnterms++] = vnoi[pos].dist;
      }
      else
      {
         break;
      }
      pos += nnodes;
   }

   return nnterms;
}

static
int getcloseterms2term(
   SCIP*                 scip,
   const GRAPH*          g,
   const GRAPH*          netgraph,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   int*                  neighbterms,
   int*                  nodesid,
   int*                  nodesorg,
   int                   i
)
{
   int nnterms = 0;

   for( int k = 0; k < 4; k++ )
      neighbterms[k] = UNKNOWN;

   /* get three nearest terms */
   for( int ne = netgraph->outbeg[nodesid[i]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
   {
      if( SCIPisLT(scip, netgraph->cost[ne], ecost) )
      {
         const SCIP_Real necost = netgraph->cost[ne];
         int j = nodesorg[netgraph->head[ne]];

         assert(Is_term(g->term[j]));
         assert(j != i);

         if( nnterms < 4 )
            nnterms++;
         for( int k = 0; k < 4; k++ )
         {
            if( neighbterms[k] == UNKNOWN || SCIPisGT(scip, termdist[k], necost) )
            {
               for( int l = 3; l > k; l-- )
               {
                  neighbterms[l] = neighbterms[l - 1];
                  termdist[l] = termdist[l - 1];
               }
               neighbterms[k] = j;
               termdist[k] = necost;
               break;
            }
         }
      }
   }

   return nnterms;
}

static
int getlecloseterms(
   SCIP*                 scip,
   PATH*                 vnoi,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   int*                  vbase,
   int*                  neighbterms,
   int                   i,
   int                   nnodes
   )
{
   int k;
   int pos = i;
   int nnterms = 0;

   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(vbase != NULL);
   assert(termdist != NULL);
   assert(neighbterms != NULL);

   for( k = 0; k < 4; k++ )
   {
      if( SCIPisLE(scip, vnoi[pos].dist, ecost) )
      {
         neighbterms[nnterms] = vbase[pos];
         termdist[nnterms++] = vnoi[pos].dist;
      }
      else
      {
         break;
      }
      pos += nnodes;
   }

   return nnterms;
}

#if 0
static int issmaller(
   SCIP* scip,
   const double* pathdist,
   const double* pathtran,
   int          a,
   int          b)
{
   return (SCIPisLT(scip, pathdist[a], pathdist[b]) || (!SCIPisGT(scip, pathdist[a], pathdist[b]) && SCIPisLT(scip, pathtran[a], pathtran[b])));
}
static int islarger(
   SCIP* scip,
   const double* pathdist,
   const double* pathtran,
   int          a,
   int          b)
{
   return (SCIPisGT(scip, pathdist[a], pathdist[b]) || (!SCIPisLT(scip, pathdist[a], pathdist[b]) && SCIPisGT(scip, pathtran[a], pathtran[b])));
}
/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearest(
   SCIP*         scip,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   const double* pathdist,
   const double* pathtran)
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
      if (islarger(scip, pathdist, pathtran, heap[2], heap[3]))
         c++;

   while((c <= (*count)) && islarger(scip, pathdist, pathtran, heap[j], heap[c]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
         if (issmaller(scip, pathdist, pathtran, heap[c + 1], heap[c]))
            c++;
   }
   return(k);
}


/** insert respectively change element in heap */
inline static void correct(
   SCIP*         scip,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   double* pathdist,
   double* pathtran,
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

   while((j > 1) && (islarger(scip, pathdist, pathtran, heap[c], heap[j])))
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

/** Dijkstra's algorithm for shortest path or minimum spanning tree */
static
void compute_sd(
   SCIP* scip,
   const GRAPH*  p,
   int           start,
   const double* cost,
   const double* randarr,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   double* pathdist,
   double* pathtran,
   double* pathrand
   )
{
   int    k;
   int    m;
   int    i;
   int    done = 0;
   double tran;
   double dist;
   double temprand;

   assert(scip != NULL);
   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(pathdist   != NULL);
   assert(pathtran   != NULL);
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
      pathdist[i] = FARAWAY;
      pathtran[i] = FARAWAY;
   }
   /* Startknoten in den Heap
    */
   k            = start;
   pathdist[k] = 0.0;
   pathtran[k] = 0.0;
   pathrand[k] = 0.0;

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
         k = nearest(scip, heap, state, count, pathdist, pathtran);

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
               if( Is_term(p->term[m]) )
               {
                  tran = 0.0;
                  temprand = 0.0;
               }
               else
               {
                  tran = pathtran[k] + cost[i];
                  temprand = pathrand[k] + randarr[i];
               }

               if( SCIPisGE(scip, pathdist[k], pathtran[k] + cost[i]) )
                  dist = pathdist[k];
               else
                  dist = pathtran[k] + cost[i];

               if( SCIPisLT(scip, dist, pathdist[m])
                  || (SCIPisEQ(scip, dist, pathdist[m]) && SCIPisLT(scip, tran, pathtran[m])) )
               {
                  pathdist[m] = dist;
                  pathtran[m] = tran;
                  pathrand[m] = temprand;

                  correct(scip, heap, state, count, pathdist, pathtran, m);
               }
            }
         }
      }
   }
}
#endif



#if 0
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
SCIP_RETCODE reduce_nv_optimal(
   SCIP*   scip,
   GRAPH*  g,
   double* fixed,
   int* elimins,
   int runnum)
{
   PATH**  path;
   PATH*   pathfromterm;
   PATH*   pathfromsource;
   PATH*   pathhops;
   double* distance;
   double* radius;
   double* hopscost;
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
   int     e;
   double  min1;
   double  min2;
   int     minhops;
   int     shortarc;
   int     shortarctail;
   STP_Bool    antiedgeexists;
   int     knotoffset;

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

   pathhops = malloc((size_t)g->knots * sizeof(PATH));

   assert(pathhops != NULL);

   distance = malloc((size_t)g->knots * sizeof(double));
   radius = malloc((size_t)g->knots * sizeof(double));

   assert(distance != NULL);
   assert(radius != NULL);

   hopscost = malloc((size_t)g->edges * sizeof(double));

   assert(hopscost != NULL);

   vregion = malloc((size_t)g->knots * sizeof(int));

   assert(vregion != NULL);

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));

   assert(heap != NULL);
   assert(state != NULL);

   *elimins = 0;
   pred = malloc((size_t)g->edges * sizeof(int));
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

      for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
      {
         if( LT(g->cost[e], FARAWAY) )
            hopscost[e] = 1;
         else
            hopscost[e] = FARAWAY;

         if( LT(g->cost[Edge_anti(e)], FARAWAY) )
            hopscost[Edge_anti(e)] = 1;
         else
            hopscost[Edge_anti(e)] = FARAWAY;
      }
   }

   assert(g->source >= 0);

   /* computing the voronoi regions inward to a node */
   voronoi_term(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, 1);

   /* computing the shortest paths from the source node */
   graph_path_exec(scip, g, FSP_MODE, g->source, g->cost, pathfromsource);

   /* computing the shortest hops paths from the source node */
   graph_path_exec(scip, g, FSP_MODE, g->source, hopscost, pathhops);

   /* this is the offset used to minimise the number of knots to examine in large graphs. */
   srand(runnum*100);
   knotoffset = rand() % KNOTFREQ;

   for(i = 0; i < g->knots; i++)
   {
      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ((g->stp_type == STP_PCSPG || g->stp_type == STP_MWCSP) && i == g->source )
         continue;

      if( g->stp_type == STP_SAP && i != g->source )
         continue;


      if( g->knots > KNOTLIMIT && i % KNOTFREQ != knotoffset )
         continue;

      if (Is_term(g->term[i]) && g->grad[i] >= 3)
      {

         min1  = FARAWAY;
         min2  = FARAWAY;
         minhops = g->hoplimit;
         shortarctail = -1;
         shortarc = -1;
         antiedgeexists = FALSE;
         for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
         {
            if( g->stp_type == STP_SAP && i == g->source )
            {
               if( LT(g->cost[e], FARAWAY) )
               {
                  g->cost[e] = FARAWAY;
                  (*elimins)++;
               }

               continue;
            }

            if (g->cost[e] < min1)
            {
               shortarc = e;
               shortarctail = g->tail[e];

               min2 = min1;
               min1 = g->cost[e];
            }

            if( LT(pathfromsource[g->tail[e]].hops, minhops) )
               minhops = pathfromsource[g->tail[e]].hops;

            if( LT(g->cost[Edge_anti(e)], FARAWAY) )
               antiedgeexists = TRUE;

         }

         if( g->stp_type == STP_SAP )
            continue;

         if (LT(min1, FARAWAY) && LE(pathfromsource[shortarctail].dist + min1, min2))
         {
            assert(shortarc >= 0);
            assert(shortarctail >= 0);

            if ((g->stp_type == STP_PCSPG || g->stp_type == STP_MWCSP
                  || g->stp_type == STP_SAP) && shortarctail == g->source )
               continue;

            if( g->stp_type == STP_DHCSTP && GT(pathfromsource[shortarctail].hops, pathhops[i].dist - 1) )
               continue;

            if( antiedgeexists == TRUE )
            {
               if( LT(min2, FARAWAY) )
               {
                  for(e = g->inpbeg[i]; e != EAT_LAST; e = j)
                  {
                     j = g->ieat[e];

                     if( e != shortarc )
                     {
                        if( LT(g->cost[Edge_anti(e)], FARAWAY) )
                           g->cost[e] = FARAWAY;
                        else
                           graph_edge_del(scip, g, e, TRUE);
                     }
                  }
                  (*elimins)++;
               }
            }
            else
            {
               *fixed += min1;
               SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[shortarc], NULL) ); /* I think that this should be
                                                                                                         shortarc instead of shortarctail */
               SCIP_CALL( graph_knot_contract(scip, g, i, shortarctail) );

               if( g->stp_type == STP_DHCSTP )
               {
                  for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
                     g->cost[e] = FARAWAY;
               }

               (*elimins)++;
            }

            /* computing the shortest paths from the source node */
            graph_path_exec(scip, g, FSP_MODE, g->source, g->cost, pathfromsource);
         }
      }
      /* The knot is not a terminal so we can perform the short link test */
#if 0
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
#endif
   }

#if 0
   for( k = 0; k < termcount; k++ )
   {
      assert(terms[k] >= 0 && terms[k] < g->knots);

      if( minArc1[terms[k]] >= 0 && minArc2[terms[k]] >= 0 && pathfromsource[g->tail[minArc1[terms[k]]]].dist
         + g->cost[minArc1[terms[k]]] + pathfromterm[g->head[minArc1[terms[k]]]].dist < g->cost[minArc2[terms[k]]] )
      {
         e = minArc1[terms[k]];
         i = g->head[e];
         j = g->tail[e];

         if ((g->stp_type == STP_PCSPG || g->stp_type == STP_MWCSP) && (i == g->source || j == g->source) )
            continue;


         if( Is_term(g->term[i]) )
         {
            SCIPintListNodeAppendCopy(&(g->fixedges), g->ancestors[e], NULL);
            *fixed += g->cost[e];
         }
         graph_knot_contract(g, j, i);


         elimins++;
      }
   }
#endif

   free(terms);
   free(minArc2);
   free(minArc1);
   free(pred);
   free(state);
   free(heap);
   free(hopscost);
   free(radius);
   free(distance);
   free(vregion);
   free(pathhops);
   free(pathfromsource);
   free(pathfromterm);
   free(path);

   assert(graph_valid(g));
   SCIPdebugMessage(" %d Knots deleted\n", *elimins);

   return SCIP_OKAY;
}

* C. W. Duin and A. Volganant
 *
 * "An Edge Elimination Test for the Steiner Problem in Graphs"
 *
 * Operations Research Letters 8, (1989), 79-83
 *
 * Special Distance Test
 */
int reduce_sduction(
   SCIP* scip,
   GRAPH* g,
   double*  sddist,
   double*  sdtrans,
   double*  sdrand,
   double* cost,
   double* random,
   int*    heap,
   int*    state,
   int*    knotexamined,
   int     runnum
   )
{
   SCIP_Real redstarttime;
   SCIP_Real timelimit;
   SCIP_Real stalltime;
   int     count = 0;
   int     i;
   int     e;
   int     j;
   int     elimins = 0;
   int     knotoffset = 0;

   SCIPdebugMessage("SD-Reduktion: ");
   fflush(stdout);

   /*
     heap  = malloc((size_t)g->knots * sizeof(int));
     state = malloc((size_t)g->knots * sizeof(int));
   */
   assert(heap  != NULL);
   assert(state != NULL);
   /*
     sd = malloc((size_t)g->knots * sizeof(SDPTH));
   */
   assert(sddist != NULL);
   assert(sdtrans != NULL);
   /*
     cost  = malloc((size_t)g->edges * sizeof(double));
   */
   assert(cost != NULL);

   assert(knotexamined != NULL);

   redstarttime = SCIPgetTotalTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   stalltime = timelimit*0.1; /* this should be set as a parameter */

   for(i = 0; i < g->knots; i++)
   {
      g->mark[i] = (g->grad[i] > 0);
      for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
      {
         random[e] = (double)(rand() % 512);
         cost[e] = g->cost[e] * 1000.0 + random[e];
      }
   }

   /* this is the offset used to minimise the number of knots to examine in large graphs. */
   if( g->knots > KNOTLIMIT )
   {
      srand(runnum*100);
      i = 0;
      do
      {
         knotoffset = rand() % KNOTFREQ;
         i++;
      } while( g->knots > KNOTLIMIT && knotexamined[knotoffset] >= 0 && i < 50 );
      knotexamined[knotoffset]++;
   }


   for(i = 0; i < g->knots; i++)
   {
      if( i % 100 == 0 && elimins == 0 && SCIPgetTotalTime(scip) - redstarttime > stalltime)
         break;

      if (!(i % 100))
      {
         SCIPdebug(fputc('.', stdout));
         SCIPdebug(fflush(stdout));
      }
      if (g->grad[i] == 0)
         continue;


      if( g->knots > KNOTLIMIT && i % KNOTFREQ != knotoffset )
         continue;

      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ( (g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING
            || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source )
         continue;

      for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
      {
         assert(g->mark[g->head[e]] == 1);

         g->mark[g->head[e]] = 2;
      }

      compute_sd(g, i, cost, random, heap, state, &count, sddist, sdtrans, sdrand);

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

         if (LT(g->cost[e], FARAWAY) && LT(sddist[g->head[e]], cost[e])
            && LT(sddist[g->head[e]] - sdrand[g->head[e]], cost[e] - random[e]))
         {
       SCIPindexListNodeFree(&((g->ancestors)[e]));
       assert(g->ancestors[e] == NULL);
            graph_edge_del(g, e);
            elimins++;
         }
      }
   }
#if 0
   free(heap);
   free(state);

   heap  = NULL;
   state = NULL;

   free(sd);
   free(cost);
#endif
   assert(graph_valid(g));

   SCIPdebugMessage("%d Edges deleted\n", elimins * 2);
   /*printf("%d SD: Edges deleted\n", elimins * 2);*/
   return(elimins);
}
#endif

/** ans subtest */
static
void ansProcessCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count,              /**< pointer to number of reductions */
   SCIP_Real             min,                /**< value to not surpass */
   int                   candvertex          /**< candidate */
)
{
   int e2;
   int bridgeedge = -1;
   unsigned misses = 0;

   assert(g->mark[candvertex]);
   assert(candvertex != g->source);
   assert(!Is_pterm(g->term[candvertex]));

   for( e2 = g->outbeg[candvertex]; e2 != EAT_LAST; e2 = g->oeat[e2] )
   {
      if( !marked[g->head[e2]] )
      {
         misses++;
         if( misses >= 2 )
            break;
         bridgeedge = e2;
      }
   }

   /* neighbors of candvertex subset of those of k? */
   if( misses == 0 && SCIPisLE(scip, g->prize[candvertex], min) )
   {
      (*count) += g->grad[candvertex];
      while( g->outbeg[candvertex] != EAT_LAST )
         graph_edge_del(scip, g, g->outbeg[candvertex], TRUE);

      g->mark[candvertex] = FALSE;
      marked[candvertex] = FALSE;
   }
   else if( misses == 1 )
   {
      int e3;
      const int neighbor = g->head[bridgeedge];

      if( SCIPisGT(scip, g->prize[neighbor] + g->prize[candvertex], min) )
         return;

      for( e3 = g->outbeg[neighbor]; e3 != EAT_LAST; e3 = g->oeat[e3] )
      {
         const int head = g->head[e3];
         if( !marked[head] )
            break;
      }
      if( e3 == EAT_LAST )
      {
         // delete both vertices?
         if( SCIPisLE(scip, g->prize[neighbor], min) && SCIPisLE(scip, g->prize[candvertex], min) )
         {
            (*count) += g->grad[candvertex] + g->grad[neighbor] - 1;
            while( g->outbeg[candvertex] != EAT_LAST )
               graph_edge_del(scip, g, g->outbeg[candvertex], TRUE);
            while( g->outbeg[neighbor] != EAT_LAST )
               graph_edge_del(scip, g, g->outbeg[neighbor], TRUE);

            g->mark[candvertex] = FALSE;
            g->mark[neighbor] = FALSE;
            marked[candvertex] = FALSE;
            marked[neighbor] = FALSE;
         }
         else
         {
            graph_edge_del(scip, g, bridgeedge, TRUE);
         }
      }
   }
}

/** Special distance test */
SCIP_RETCODE reduce_sd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            edgepreds,          /**< array to store edge predecessors of auxiliary graph */
   SCIP_Real*            mstsdist,           /**< array to store mst distances in auxiliary graph */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  nodesid,            /**< array to map nodes in auxiliary graph to original ones */
   int*                  nodesorg,           /**< array to map terminals of original graph vertices of auxiliary graph */
   int*                  forbidden,          /**< array to mark whether an edge may be eliminated */
   int*                  nelims,             /**< point to store number of deleted edges */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int*                  edgestate           /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   )
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   SCIP_Real ecost;
   SCIP_Real dist;
   int neighbterms1[4];
   int neighbterms2[4];
   int e;
   int i;
   int j;
   int k;
   int l;
   int i2;
   int tj;
   int tk;
   int ne;
   int nj;
   int nk;
   int id1;
   int id2;
   int enext;
   int nnodes;
   int nterms;
   int nedges;
   int nnterms1;
   int nnterms2;
   int maxnedges;
   SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);
   assert(nodesid != NULL);
   assert(mstsdist != NULL);
   assert(nodesorg != NULL);
   assert(edgepreds != NULL);
   assert(forbidden != NULL);

   nnodes = g->knots;
   nterms = g->terms;
   nedges = g->edges;
   *nelims = 0;
   maxnedges = MIN(nedges, (nterms - 1) * nterms);

   /* only one terminal left? */
   if( nterms == 1 )
      return SCIP_OKAY;

   /* compute nearest four terminals to all non-terminals */
   graph_get4nextTerms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);

   /* construct auxiliary graph to compute paths between terminals */

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   j = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
         graph_knot_add(netgraph, -1);
         nodesid[k] = j;
         mstsdist[j] = -1.0;
         netgraph->mark[j] = TRUE;
         nodesorg[j++] = k;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   for( k = 0; k < nedges; k++ )
   {
      forbidden[k] = FALSE;
      edgepreds[k] = -1;
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   for( k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         continue;

      i = vbase[k];
      id1 = nodesid[i];

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( i != vbase[g->head[e]] )
         {
            i2 = vbase[g->head[e]];
            id2 = nodesid[i2];

            assert(id1 >= 0);
            assert(id2 >= 0);
            assert(Is_term(g->term[i]));
            assert(Is_term(g->term[i2]));

            for( ne = netgraph->outbeg[id1]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == id2 )
                  break;

            /* cost of the edge in the auxiliary graph */
            ecost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;

            /* does edge already exist? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->tail[ne] == id1);
               assert(netgraph->head[ne] == id2);

               /* is the new edge better than the existing one? */
               if( SCIPisGT(scip, netgraph->cost[ne], ecost) )
               {
                  netgraph->cost[ne]            = ecost;
                  netgraph->cost[flipedge(ne)]  = ecost;

                  edgepreds[ne] = e;
                  edgepreds[flipedge(ne)] = flipedge(e);

                  assert(ne <= maxnedges);
               }
            }
            else
            {
               edgepreds[netgraph->edges] = e;
               edgepreds[netgraph->edges + 1] = flipedge(e);

               graph_edge_add(scip, netgraph, id1, id2, ecost, ecost);

               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }

   /* compute MST on netgraph */
   graph_knot_chg(netgraph, 0, 0);
   netgraph->source = 0;

   SCIP_CALL( graph_path_init(scip, netgraph) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );

   graph_path_exec(scip, netgraph, MST_MODE, 0, netgraph->cost, mst);

   /* mark (original) edges of MST */
   for( k = 1; k < netgraph->knots; k++ )
   {
      assert(mst[k].edge != -1);
      assert((int) edgepreds[mst[k].edge] != -1);
      assert((int) edgepreds[flipedge(mst[k].edge)] != -1);

      e = (int) edgepreds[mst[k].edge];

      assert(vbase[g->tail[e]] == nodesorg[k] || vbase[g->head[e]] == nodesorg[k]);

      if( Is_term(g->tail[e]) && Is_term(g->head[e]) )
      {
         forbidden[e] = TRUE;
         forbidden[flipedge(e)] = TRUE;
      }
   }

   /* traverse all edges */
   for( i = 0; i < nnodes; i++ )
   {
      if( g->grad[i] <= 0 )
         continue;

      nnterms1 = 1;
      if( Is_term(g->term[i]) )
      {
         termdist1[0] = 0.0;
         neighbterms1[0] = i;
      }

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         e = enext;
         i2 = g->head[e];
         enext = g->oeat[e];

         if( i2 < i || !g->mark[i2] )
            continue;

         ecost = g->cost[e];

         /* is i a terminal? If not, get three closest terminals of distance <= ecost */
         if( !Is_term(g->term[i]) )
         {
            nnterms1 = getlecloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);

            if( nnterms1 == 0 )
               continue;
         }

         if( Is_term(g->term[i2]) )
         {
            nnterms2 = 1;
            termdist2[0] = 0.0;
            neighbterms2[0] = i2;
         }
         else
         {
            /* get closest terminals of distance <= ecost */
            nnterms2 = getlecloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);

            if( nnterms2 == 0 )
               continue;
         }

         for( j = 0; j < nnterms1; j++ )
         {
            /* has edge already been deleted? */
            if( g->oeat[e] == EAT_FREE )
               break;

            tj = neighbterms1[j];

            assert(tj >= 0);

            for( k = 0; k < nnterms2; k++ )
            {
               tk = neighbterms2[k];

               assert(tk >= 0);
               assert(Is_term(g->term[tk]));
               assert(Is_term(g->term[tj]));

               if( tj == tk )
               {
                  if( SCIPisGE(scip, termdist1[j], termdist2[k] ) )
                     dist = termdist1[j];
                  else
                     dist = termdist2[k];

                  assert(SCIPisGE(scip, ecost, dist));

                  if( SCIPisEQ(scip, dist, ecost) )
                     if( !sddeltable(scip, g, vnoi, vbase, forbidden, j, k, i, i2, e, nnodes ) )
                        continue;

                  if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
                     continue;

                  graph_edge_del(scip, g, e, TRUE);
                  (*nelims)++;
                  break;
               }
               else
               {
                  /* get sd between (terminals) tj and tk */
                  nj = nodesid[tj];
                  nk = nodesid[tk];

                  assert(nj != nk);

                  l = nj;
                  dist = 0.0;
                  mstsdist[l] = 0.0;

                  /* go down to the root */
                  while( l != 0 )
                  {
                     ne = mst[l].edge;

                     assert(netgraph->head[ne] == l);

                     l = netgraph->tail[ne];

                     if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                        dist = netgraph->cost[ne];

                     mstsdist[l] = dist;

                     if( l == nk )
                        break;
                  }

                  if( l == nk )
                  {
                     l = 0;
                  }
                  else
                  {
                     l = nk;
                     dist = 0.0;
                  }

                  while( l != 0 )
                  {
                     ne = mst[l].edge;
                     l = netgraph->tail[ne];

                     if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                        dist = netgraph->cost[ne];

                     if( mstsdist[l] >= 0 )
                     {
                        if( SCIPisGT(scip, mstsdist[l], dist) )
                           dist = mstsdist[l];
                        break;
                     }
                     if( l == 0 )
                        assert(nj == 0);
                  }

                  /* restore */
                  l = nj;
                  mstsdist[l] = -1.0;

                  while( l != 0 )
                  {
                     ne = mst[l].edge;
                     l = netgraph->tail[ne];
                     mstsdist[l] = -1.0;
                     if( l == nk )
                        break;
                  }

                  assert(SCIPisGT(scip, dist, 0.0));

                  if( SCIPisLT(scip, dist, termdist1[j]) )
                     dist = termdist1[j];

                  if( SCIPisLT(scip, dist, termdist2[k]) )
                     dist = termdist2[k];

                  if( SCIPisGE(scip, ecost, dist) )
                  {
                     if( SCIPisEQ(scip, ecost, dist) )
                        if( !(sddeltable(scip, g, vnoi, vbase, forbidden, j, k, i, i2, e, nnodes)) )
                           continue;

                     assert(SCIPisGE(scip, ecost, termdist1[j]));
                     assert(SCIPisGE(scip, ecost, termdist2[k]));

                     if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
                        continue;

                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               } /* tj != tk (else) */
            } /* k < nnterms2 */
         } /* j < nnterms1 */

      } /* while( enext != EAT_LAST ) */
   }

   if( nodereplacing )
   {
      SCIP_CALL( reduce_bdr(scip, g, netgraph, mst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, nelims) );
   }

   /* free memory*/
   SCIPfreeBufferArray(scip, &mst);
   graph_path_exit(scip, netgraph);
   graph_free(scip, &netgraph, TRUE);

   return SCIP_OKAY;
}


/** SD test for PC */
SCIP_RETCODE reduce_sdPc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nodesid,            /**< array */
   int*                  nodesorg,           /**< array */
   int*                  nelims              /**< pointer to store number of eliminated edges */
   )
{
   GRAPH* netgraph;
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   int neighbterms1[4];
   int neighbterms2[4];

   int j;
   int maxnedges;
   const int root = g->source;
   const int nnodes = g->knots;
   const int nterms = g->terms;
   const int nedges = g->edges;

   assert(g != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(nodesid != NULL);
   assert(nodesorg != NULL);

   *nelims = 0;

   if( nterms <= 1 )
      return SCIP_OKAY;
   else
   {
      const SCIP_Longint longedges = (SCIP_Longint) nedges;
      const SCIP_Longint longtermsq = (SCIP_Longint) (nterms - 1) * nterms;

      if( longedges <= longtermsq )
         maxnedges = nedges;
      else
         maxnedges = ((nterms - 1) * nterms);
   }

   /* compute nearest four terminals to each non-terminal */
   graph_get4nextTerms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);

   /*
    * construct auxiliary graph to compute paths between terminals
    */

   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   for( int k = 0; k < 4; k++ )
   {
      termdist1[k] = FARAWAY;
      termdist2[k] = FARAWAY;
   }

   j = 0;
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         assert(g->grad[k] > 0);
         graph_knot_add(netgraph, -1);
         nodesid[k] = j;
         nodesorg[j++] = k;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   /* insert Voronoi boundary paths as edges into netgraph */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         int i;
         if( !g->mark[g->head[e]] )
            continue;
         i = vbase[k];
         assert(i != UNKNOWN);
         if( i != vbase[g->head[e]] )
         {
            SCIP_Real ecost;
            int ne;
            const int i2 = vbase[g->head[e]];

            assert(i2 != UNKNOWN);
            assert(Is_term(g->term[i]));
            assert(Is_term(g->term[i2]));
            assert(nodesid[i] >= 0);
            assert(nodesid[i2] >= 0);

            for( ne = netgraph->outbeg[nodesid[i]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == nodesid[i2] )
                  break;

            ecost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;

            /* edge exists? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->head[ne] == nodesid[i2]);
               assert(netgraph->tail[ne] == nodesid[i]);

               if( SCIPisGT(scip, netgraph->cost[ne], ecost) )
               {
                  netgraph->cost[ne]            = ecost;
                  netgraph->cost[flipedge(ne)] = ecost;
                  assert(ne <= maxnedges);
               }
            }
            else
            {
               graph_edge_add(scip, netgraph, nodesid[i], nodesid[i2], ecost, ecost);
               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }

   /* compute four close terminals to each terminal */
   SCIP_CALL( graph_get4nextTTerms(scip, g, g->cost, vnoi, vbase, heap, state) );

   /* traverse all edges */
   for( int i = 0; i < nnodes; i++ )
   {
      int enext;
      if( !g->mark[i] )
         continue;

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         SCIP_Real ecost;
         int e = enext;
         int nnterms1;
         int nnterms2;
         const int i2 = g->head[e];
         enext = g->oeat[e];

         if( i2 < i || Is_term(g->term[i2]) || !g->mark[i2] )
            continue;
         ecost = g->cost[e];

         /* @todo: fix */
#if 1
         if( Is_term(g->term[i]) )
            nnterms1 = getcloseterms2term(scip, g, netgraph, termdist1, ecost, neighbterms1, nodesid, nodesorg, i);
         else
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
#else
         if( Is_term(g->term[i]) )
            nnterms1 = 0;
         else
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
#endif

         if( nnterms1 == 0 )
            continue;

         /* @todo: fix */
#if 1
         if( Is_term(g->term[i2]) )
            nnterms2 = getcloseterms2term(scip, g, netgraph, termdist2, ecost, neighbterms2, nodesid, nodesorg, i2);
         else
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#else
         if( Is_term(g->term[i2]) )
            nnterms2 = 0;
         else
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#endif

         if( nnterms2 == 0 )
            continue;

         /* todo: mark nearest terminals! */
         for( j = 0; j < nnterms1; j++ )
         {
            int tj;

            /* has edge already been deleted? */
            if( g->oeat[e] == EAT_FREE )
               break;

            tj = neighbterms1[j];

            assert(tj >= 0);
            assert(Is_term(g->term[tj]));

            for( int k = 0; k < nnterms2; k++ )
            {
               const int tk = neighbterms2[k];

               assert(tk >= 0);
               assert(Is_term(g->term[tk]));

               if( tj == tk )
               {
                  if( SCIPisGT(scip, ecost, termdist1[j] + termdist2[k] - g->prize[tj]) || tj == root )
                  {
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               }
               else
               {
                  SCIP_Real necost = FARAWAY;
                  int e2;
                  int pos;

                  /* get distance between terminals */
                  for( e2 = netgraph->outbeg[nodesid[tj]]; e2 != EAT_LAST; e2 = netgraph->oeat[e2] )
                  {
                     if( netgraph->head[e2] == nodesid[tk] )
                     {
                        necost = netgraph->cost[e2];
                        break;
                     }
                  }
                  pos = tj;
                  for( int l = 0; l < 4; l++ )
                  {
                     if( vbase[pos] == UNKNOWN )
                        break;
                     if( vbase[pos] == tk && SCIPisLT(scip, vnoi[pos].dist, necost) )
                     {
                        necost = vnoi[pos].dist;
                        break;
                     }
                     pos += nnodes;
                  }

                  if( SCIPisGT(scip, ecost, necost)
                     && SCIPisGT(scip, ecost, necost + termdist1[j] - g->prize[tj])
                     && SCIPisGT(scip, ecost, necost + termdist2[k] - g->prize[tk])
                     && SCIPisGT(scip, ecost, necost + termdist1[j] + termdist2[k] - g->prize[tj] - g->prize[tk]) )
                  {
                     SCIPdebugMessage("SDPC delete: %d %d (%d)\n", g->tail[e], g->head[e], e);
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               }
            }
         }
      }
   }

   SCIPdebugMessage("SDPC eliminations: %d \n", *nelims);
   graph_free(scip, &netgraph, TRUE);

   assert(graph_valid(g));

   return SCIP_OKAY;
}

/** get RSD to a single edge*/
static
SCIP_Real getRSD(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   GRAPH*                netgraph,           /**< auxiliary graph structure */
   PATH*                 mst,                /**< MST structure */
   PATH*                 vnoi,               /**< path structure */
   SCIP_Real*            mstsdist,           /**< MST distance in aux-graph */
   SCIP_Real*            termdist1,          /**< dist array */
   SCIP_Real*            termdist2,          /**< second dist array */
   SCIP_Real             sd_initial,         /**< initial sd or -1.0 */
   int*                  vbase,              /**< bases for nearest terminals */
   int*                  nodesid,            /**< nodes identification array */
   int*                  neighbterms1,       /**< neighbour terminals array */
   int*                  neighbterms2,       /**< second neighbour terminals array */
   int                   i,                  /**< first vertex */
   int                   i2,                 /**< second vertex */
   int                   limit               /**< limit for incident edges to consider */
   )
{
   SCIP_Real sd;
   SCIP_Real max;
   SCIP_Real dist;
   int l;
   int e;
   int j;
   int k;
   int ne;
   int tj;
   int tk;
   int nj;
   int nk;
   int nnodes;
   int nnterms1;
   int nnterms2;

   assert(scip != NULL);
   assert(g != NULL);
   assert(netgraph != NULL);
   assert(mst != NULL);
   assert(mstsdist != NULL);
   assert(termdist1 != NULL);
   assert(termdist2 != NULL);
   assert(neighbterms1 != NULL);
   assert(neighbterms2 != NULL);

   nnodes = g->knots;
   l = 0;

   if( sd_initial >= 0.0 )
      sd = sd_initial;
   else
      sd = FARAWAY;

   /* compare restricted sd with edge cost (if existing) */
   for( e = g->outbeg[i]; (l++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( g->cost[e] < sd )
            sd = g->cost[e];
         break;
      }
   }

   /* is i a terminal? If not, get three closest terminals of distance smaller sd */
   if( Is_term(g->term[i]) )
   {
      nnterms1 = 1;
      termdist1[0] = 0.0;
      neighbterms1[0] = i;
   }
   else
   {
      nnterms1 = getcloseterms(scip, vnoi, termdist1, sd, vbase, neighbterms1, i, nnodes);

      if( nnterms1 == 0 )
         return sd;
   }

   /* is i2 a terminal? If not, get three closest terminals of distance smaller sd */
   if( Is_term(g->term[i2]) )
   {
      nnterms2 = 1;
      termdist2[0] = 0.0;
      neighbterms2[0] = i2;
   }
   else
   {
      /* get closest terminals of distance smaller sd */
      nnterms2 = getcloseterms(scip, vnoi, termdist2, sd, vbase, neighbterms2, i2, nnodes);

      if( nnterms2 == 0 )
         return sd;
   }

   for( j = 0; j < nnterms1; j++ )
   {
      tj = neighbterms1[j];
      assert(tj >= 0);
      for( k = 0; k < nnterms2; k++ )
      {
         tk = neighbterms2[k];
         assert(Is_term(g->term[tk]));
         assert(Is_term(g->term[tj]));
         assert(tk >= 0);

         if( SCIPisGT(scip, termdist1[j], termdist2[k]) )
            max = termdist1[j];
         else
            max = termdist2[k];

         if( tj == tk )
         {
            if( SCIPisLT(scip, max, sd) )
               sd = max;
         }
         else
         {
            /* get sd between (terminals) tj and tk */
            nj = nodesid[tj];
            nk = nodesid[tk];
            assert(nj != nk);

            l = nj;
            dist = 0.0;
            mstsdist[l] = 0.0;

            while( l != 0 )
            {
               ne = mst[l].edge;

               assert(netgraph->head[ne] == l);
               l = netgraph->tail[ne];
               if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                  dist = netgraph->cost[ne];

               mstsdist[l] = dist;
               if( l == nk )
                  break;
            }

            if( l == nk )
            {
               l = 0;
            }
            else
            {
               l = nk;
               dist = 0.0;
            }
            while( l != 0 )
            {
               ne = mst[l].edge;
               l = netgraph->tail[ne];
               if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                  dist = netgraph->cost[ne];

               if( mstsdist[l] >= 0 )
               {
                  if( SCIPisGT(scip, mstsdist[l], dist) )
                     dist = mstsdist[l];
                  break;
               }
               if( l == 0 )
                  assert( nj == 0);
            }

            /* restore */
            l = nj;
            mstsdist[l] = -1.0;
            while( l != 0 )
            {
               ne = mst[l].edge;
               l = netgraph->tail[ne];
               mstsdist[l] = -1.0;
               if( l == nk )
                  break;
            }

            assert(SCIPisGT(scip, dist, 0.0));
            if( SCIPisGT(scip, dist, max) )
               max = dist;
            if( SCIPisLT(scip, max, sd) )
               sd = max;
         }
      } /* k < nnterms2 */
   } /* j < nnterms1 */
   return sd;
}



/** get SD to a single edge*/
SCIP_RETCODE reduce_getSd(
   SCIP* scip,
   GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   SCIP_Real*    sdist,
   SCIP_Real     distlimit,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int     i,
   int     i2,
   int     limit,
   SCIP_Bool pc,
   SCIP_Bool mw
   )
{
   SCIP_Real sd;
   SCIP_Real dist;
   int k;
   int l;
   int e;
   int nnodes;
   int nlbltail;
   int nlblhead;
   const SCIP_Bool pcmw = (pc || mw);

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(sdist != NULL);

   nnodes = g->knots;

   /* start from tail */
   graph_sdPaths(scip, g, pathtail, g->cost, distlimit, heap, statetail, memlbltail, &nlbltail, i, i2, limit);

   /* test whether edge e can be eliminated */
   graph_sdPaths(scip, g, pathhead, g->cost, distlimit, heap, statehead, memlblhead, &nlblhead, i2, i, limit);

   sd = FARAWAY;
#if 0
   if( statetail[i2] != UNKNOWN )
   {
      sd = pathtail[i2].dist;
      assert(SCIPisGT(scip, FARAWAY, sd));
   }
   if( statehead[i] != UNKNOWN && SCIPisGT(scip, sd, pathhead[i].dist) )
   {
      sd = pathhead[i].dist;
      assert(SCIPisGT(scip, FARAWAY, sd));
   }
#endif
   /* get restore state and path of tail and head */
   for( k = 0; k < nlbltail; k++ )
   {
      l = memlbltail[k];
      assert(statetail[l]     != UNKNOWN);
      if( statehead[l] != UNKNOWN )
      {
         assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
         assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));
         if( Is_term(g->term[l]) )
         {
            dist = 0.0;
            if( SCIPisLT(scip, dist, pathhead[l].dist) )
               dist = pathhead[l].dist;
            if( SCIPisLT(scip, dist, pathtail[l].dist) )
               dist = pathtail[l].dist;
            if( pcmw && SCIPisLT(scip, dist, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
               dist = pathhead[l].dist + pathtail[l].dist - g->prize[l];
            if( SCIPisGT(scip, sd, dist) )
               sd = dist;
         }
         else
         {
            if( mw && l != i && l != i2 )
               assert(SCIPisLE(scip, g->prize[l], 0.0));
            if( mw && SCIPisLT(scip, g->prize[l], 0.0) )
               dist = pathhead[l].dist + pathtail[l].dist + g->prize[l];
            else
               dist = pathhead[l].dist + pathtail[l].dist;
            if( SCIPisGT(scip, sd, dist) )
               sd = dist;
         }
      }

      statetail[l]     = UNKNOWN;
      pathtail[l].dist = FARAWAY;
      pathtail[l].edge = UNKNOWN;
   }
   /* restore state and path of tail and head */
   for( k = 0; k < nlblhead; k++ )
   {
      l = memlblhead[k];
      statehead[l]     = UNKNOWN;
      pathhead[l].dist = FARAWAY;
      pathhead[l].edge = UNKNOWN;
   }


   for( k = 0; k < nnodes; k++ )
   {
      assert(statetail[k]     == UNKNOWN);
      assert(pathtail[k].dist == FARAWAY);
      assert(pathtail[k].edge == UNKNOWN);
      assert(statehead[k]     == UNKNOWN);
      assert(pathhead[k].dist == FARAWAY);
      assert(pathhead[k].edge == UNKNOWN);
   }


   l = 0;
   /* compare restricted sd with edge cost (if existing) */
   for( e = g->outbeg[i]; (l++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( mw )
            sd = 0.0;
         else if( SCIPisGT(scip, sd, g->cost[e]) )
            sd = g->cost[e];
         break;
      }
   }

   *sdist = sd;
   return SCIP_OKAY;
}


/** get SD to a single edge*/
SCIP_RETCODE reduce_getSdPcMw(
   SCIP* scip,
   const GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   SCIP_Real*    sdist,
   SCIP_Real     distlimit,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int*    pathmaxnodetail,
   int*    pathmaxnodehead,
   int     i,
   int     i2,
   int     limit
   )
{
   SCIP_Real sd;
   int nlbltail;
   int nlblhead;
   const int nnodes = g->knots;
   const SCIP_Bool mw = g->stp_type == STP_MWCSP;

   assert((g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG) || mw);
   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(pathmaxnodetail != NULL);
   assert(pathmaxnodehead != NULL);
   assert(sdist != NULL);

   graph_path_PcMwSd(scip, g, pathtail, g->cost, distlimit, pathmaxnodetail, heap, statetail, NULL, memlbltail, &nlbltail, i, i2, limit);
   graph_path_PcMwSd(scip, g, pathhead, g->cost, distlimit, pathmaxnodehead, heap, statehead, statetail, memlblhead, &nlblhead, i2, i, limit);

   sd = FARAWAY;

   /* get restore state and path of tail and head */
   for( int k = 0; k < nlbltail; k++ )
   {
      const int l = memlbltail[k];
      assert(statetail[l] != UNKNOWN);

      if( statehead[l] != UNKNOWN )
      {
         SCIP_Real dist = FARAWAY;
         const int tailmaxterm = pathmaxnodetail[l];
         const int headmaxterm = pathmaxnodehead[l];

         assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
         assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));
         assert(tailmaxterm != i && headmaxterm != i);
         assert(tailmaxterm != i2 && headmaxterm != i2);

         /* any terminal on the path? */
         if( tailmaxterm >= 0 || headmaxterm >= 0 )
         {
            if( tailmaxterm == headmaxterm )
            {
               assert(tailmaxterm == l);
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

               dist = misc_stp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      pathtail[tailmaxterm].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[l]
                     }, 3);
               SCIPdebugMessage("sd1 %f \n", dist);
            }
            else if( tailmaxterm >= 0 && headmaxterm >= 0 )
            {
               const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
               const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

               assert(tailmaxterm != headmaxterm);
               assert(!SCIPisNegative(scip, distl2tailmax));
               assert(!SCIPisNegative(scip, distl2headmax));
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]) && SCIPisPositive(scip, g->prize[headmaxterm]));

               dist = misc_stp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      pathtail[tailmaxterm].dist,
                      distl2tailmax + distl2headmax,
                      distl2tailmax + pathhead[l].dist - g->prize[headmaxterm],
                      distl2headmax + pathtail[l].dist - g->prize[tailmaxterm],
                      pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm] - g->prize[headmaxterm]
                     }, 6);
               SCIPdebugMessage("sd2 %f \n", dist);
            }
            else if( tailmaxterm >= 0 )
            {
               const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;

               assert(headmaxterm < 0);
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

               dist = misc_stp_maxReal((SCIP_Real []) {
                      pathtail[tailmaxterm].dist,
                      distl2tailmax + pathhead[l].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm]
                     }, 3);
               SCIPdebugMessage("sd3 %f \n", dist);
            }
            else if( headmaxterm >= 0 )
            {
               const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

               assert(tailmaxterm < 0);
               assert(SCIPisPositive(scip, g->prize[headmaxterm]));

               dist = misc_stp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      distl2headmax + pathtail[l].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[headmaxterm]
                     }, 3);
               SCIPdebugMessage("sd4 %f \n", dist);
            }
         }
         else
         {
            dist = pathhead[l].dist + pathtail[l].dist;
         }

         if( dist < sd )
            sd = dist;

         if( Is_term(g->term[l]) )
         {
            dist = misc_stp_maxReal((SCIP_Real []) {
                   pathhead[l].dist,
                   pathtail[l].dist,
                   pathhead[l].dist + pathtail[l].dist - g->prize[l]
                  }, 3);
            if( dist < sd )
               sd = dist;
         }
      }
   }

   /* restore state and path of tail and head */

   for( int k = 0; k < nlbltail; k++ )
   {
      const int l = memlbltail[k];
      statetail[l] = UNKNOWN;
      pathtail[l].dist = FARAWAY;
      pathtail[l].edge = UNKNOWN;
      pathmaxnodetail[l] = -1;
   }

   for( int k = 0; k < nlblhead; k++ )
   {
      const int l = memlblhead[k];
      statehead[l] = UNKNOWN;
      pathhead[l].dist = FARAWAY;
      pathhead[l].edge = UNKNOWN;
      pathmaxnodehead[l] = -1;
   }

   for( int k = 0; k < nnodes; k++ )
   {
      assert(statetail[k]     == UNKNOWN);
      assert(pathtail[k].dist == FARAWAY);
      assert(pathtail[k].edge == UNKNOWN);
      assert(statehead[k]     == UNKNOWN);
      assert(pathhead[k].dist == FARAWAY);
      assert(pathhead[k].edge == UNKNOWN);
      assert(pathmaxnodehead[k] == -1);
      assert(pathmaxnodetail[k] == -1);
   }

   /* compare restricted sd with edge cost (if existing) */
   for( int e = g->outbeg[i], count = 0; (count++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( mw )
            sd = 0.0;
         else if( sd > g->cost[e] )
            sd = g->cost[e];
         break;
      }
   }

   *sdist = sd;

   return SCIP_OKAY;
}


/** SDC test for the SAP using a limited version of Dijkstra's algorithm from both endpoints of an arc */
SCIP_RETCODE reduce_sdspSap(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit
   )
{
   SCIP_Real sdist;
   SCIP_Real* costrev;
   int i;
   int k;
   int l;
   int e;
   int i2;
   int enext;
   int nnodes;
   int nedges;
   int nlbltail;
   int nlblhead;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   nedges = g->edges;
   nnodes = g->knots;
   *nelims = 0;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i] = (g->grad[i] > 0);
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   for( e = 0; e < nedges; e++ )
      costrev[e] = g->cost[flipedge(e)];

   /* iterate through all nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;
      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {

         enext = g->oeat[e];
         i2 = g->head[e];

         assert(g->mark[i2]);

         /* start limited dijkstra from i, marking all reached vertices */
         graph_sdPaths(scip, g, pathtail, g->cost, g->cost[e], heap, statetail, memlbltail, &nlbltail, i, i2, limit);

         /* start limited dijkstra from i2, marking all reached vertices */
         graph_sdPaths(scip, g, pathhead, costrev, g->cost[e], heap, statehead, memlblhead, &nlblhead, i2, i, limit);

         sdist = FARAWAY;
#if 0
         if( statetail[i2] != UNKNOWN )
         {
            sdist = pathtail[i2].dist;
            assert(SCIPisGT(scip, FARAWAY, sdist));
         }
         if( statehead[i] != UNKNOWN && SCIPisGT(scip, sdist, pathhead[i].dist) )
         {
            sdist = pathhead[i].dist;
            assert(SCIPisGT(scip, FARAWAY, sdist));
         }
#endif
         /* get restore state and path of tail and head */
         for( k = 0; k < nlbltail; k++ )
         {
            l = memlbltail[k];
            assert(g->mark[l]);
            assert(statetail[l]     != UNKNOWN);
            if( statehead[l] != UNKNOWN )
            {
               assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
               assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));

               if( SCIPisGT(scip, sdist, pathhead[l].dist + pathtail[l].dist) )
                  sdist = pathhead[l].dist + pathtail[l].dist;
            }

            statetail[l]     = UNKNOWN;
            pathtail[l].dist = FARAWAY;
            pathtail[l].edge = UNKNOWN;
         }
         /* restore state and path of tail and head */
         for( k = 0; k < nlblhead; k++ )
         {
            l = memlblhead[k];
            statehead[l]     = UNKNOWN;
            pathhead[l].dist = FARAWAY;
            pathhead[l].edge = UNKNOWN;
         }

#if 1
         for( k = 0; k < nnodes; k++ )
         {
            assert(statetail[k]     == UNKNOWN);
            assert(pathtail[k].dist == FARAWAY);
            assert(pathtail[k].edge == UNKNOWN);
            assert(statehead[k]     == UNKNOWN);
            assert(pathhead[k].dist == FARAWAY);
            assert(pathhead[k].edge == UNKNOWN);
         }
#endif
         /* can edge be deleted? */
         if( SCIPisGE(scip, g->cost[e], sdist) )
         {
            if( SCIPisGE(scip, costrev[e], FARAWAY) )
               graph_edge_del(scip, g, e, TRUE);
            else
               g->cost[e] = FARAWAY;
            (*nelims)++;
         }

         e = enext;
      }
   }

   SCIPfreeBufferArray(scip, &costrev);

   return SCIP_OKAY;
}

/** SD test using only limited dijkstra from both endpoints of an edge */
SCIP_RETCODE reduce_sdsp(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit,
   int*                  edgestate

   )
{
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   const int nnodes = g->knots;
   const SCIP_Bool pc = (g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG);
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   *nelims = 0;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );

      for( int i = 0; i < nnodes; i++ )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }
   else
   {
      for( int i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   /* iterate through all nodes */
   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         const SCIP_Real ecost = g->cost[e];
         const int i2 = g->head[e];
         const int enext = g->oeat[e];
         int nlbltail;
         int nlblhead;
         SCIP_Bool deletable;

         /* avoid double checking */
         if( i2 < i || !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         /* execute limited Dijkstra from both sides */

         if( pc )
         {
            graph_path_PcMwSd(scip, g, pathtail, g->cost, ecost, pathmaxnodetail, heap, statetail, NULL, memlbltail, &nlbltail, i, i2, limit);
            graph_path_PcMwSd(scip, g, pathhead, g->cost, ecost, pathmaxnodehead, heap, statehead, statetail, memlblhead, &nlblhead, i2, i, limit);
         }
         else
         {
            graph_sdPaths(scip, g, pathtail, g->cost, ecost, heap, statetail, memlbltail, &nlbltail, i, i2, limit);
            graph_sdPaths(scip, g, pathhead, g->cost, ecost, heap, statehead, memlblhead, &nlblhead, i2, i, limit);
         }

         deletable = FALSE;

         /* check whether edge e can be deleted and restore data structures */
         for( int k = 0; k < nlbltail && !deletable; k++ )
         {
            const int l = memlbltail[k];

            assert(g->mark[l]);
            assert(statetail[l] != UNKNOWN);

            if( statehead[l] != UNKNOWN )
            {
               assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
               assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));

               if( pc )
               {
                  const int tailmaxterm = pathmaxnodetail[l];
                  const int headmaxterm = pathmaxnodehead[l];

                  assert(tailmaxterm != i && headmaxterm != i);
                  assert(tailmaxterm != i2 && headmaxterm != i2);

                  /* any terminal on the path? */
                  if( tailmaxterm >= 0 || headmaxterm >= 0 )
                  {
                     if( tailmaxterm == headmaxterm )
                     {
                        assert(tailmaxterm == l);
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]));
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist) && SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));

                        if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("delete1Term \n");
                        }
                     }
                     else if( tailmaxterm >= 0 && headmaxterm >= 0 )
                     {
                        const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
                        const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

                        assert(tailmaxterm != headmaxterm);
                        assert(!SCIPisNegative(scip, distl2tailmax));
                        assert(!SCIPisNegative(scip, distl2headmax));
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]) && SCIPisPositive(scip, g->prize[headmaxterm]));
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist) && SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));

                        if( SCIPisGE(scip, ecost, distl2tailmax + distl2headmax)
                              && SCIPisGE(scip, ecost, distl2tailmax + pathhead[l].dist - g->prize[headmaxterm])
                              && SCIPisGE(scip, ecost, distl2headmax + pathtail[l].dist - g->prize[tailmaxterm])
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm] - g->prize[headmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("delete2Term \n");
                        }
                     }
                     else if( tailmaxterm >= 0 )
                     {
                        const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
                        // todo consider l == term?
                        assert(headmaxterm < 0);
                        assert(SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

                        if( SCIPisGE(scip, ecost, distl2tailmax + pathhead[l].dist)
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("deleteHalfTerm1 \n");
                        }
                     }
                     else if( headmaxterm >= 0 )
                     {
                        const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;
                        // todo consider l == term?
                        assert(tailmaxterm < 0);
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist));
                        assert(SCIPisPositive(scip, g->prize[headmaxterm]));

                        if( SCIPisGE(scip, ecost, distl2headmax + pathtail[l].dist)
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[headmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("deleteHalfTerm2 \n");
                        }
                     }
                  }
                  else if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist) )
                  {
                     deletable = TRUE;
                  }

                  if( Is_term(g->term[l]) )
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist) && SCIPisGE(scip, ecost, pathtail[l].dist)
                           && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                        deletable = TRUE;
                  }
               }
               else
               {
                  if( Is_term(g->term[l]) )
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist) && SCIPisGE(scip, ecost, pathtail[l].dist) )
                     {
                        deletable = TRUE;
#if 0
                        if( pc && SCIPisLT(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                           deletable = FALSE;
#endif
                     }
                  }
                  else
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist) )
                        deletable = TRUE;
                  }
               }
            }
         }

         /* restore data */

         for( int k = 0; k < nlbltail; k++ )
         {
            const int l = memlbltail[k];
            statetail[l] = UNKNOWN;
            pathtail[l].dist = FARAWAY;
            pathtail[l].edge = UNKNOWN;
            if( pc )
               pathmaxnodetail[l] = -1;
         }

         for( int k = 0; k < nlblhead; k++ )
         {
            const int l = memlblhead[k];
            statehead[l]     = UNKNOWN;
            pathhead[l].dist = FARAWAY;
            pathhead[l].edge = UNKNOWN;
            if( pc )
               pathmaxnodehead[l] = -1;
         }

#ifndef NDEBUG
         for( int k = 0; k < nnodes; k++ )
         {
            assert(statetail[k]     == UNKNOWN);
            assert(pathtail[k].dist == FARAWAY);
            assert(pathtail[k].edge == UNKNOWN);
            assert(statehead[k]     == UNKNOWN);
            assert(pathhead[k].dist == FARAWAY);
            assert(pathhead[k].edge == UNKNOWN);
            if( pc )
            {
               assert(pathmaxnodetail[k] == -1);
               assert(pathmaxnodehead[k] == -1);
            }
         }
#endif
         /* can edge be deleted? */
         if( deletable )
         {
            if( !checkstate || (edgestate[e] != EDGE_BLOCKED) )
            {
               graph_edge_del(scip, g, e, TRUE);
               (*nelims)++;
            }
         }

         e = enext;
      }
   }

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   assert(graph_valid(g));

   return SCIP_OKAY;
}


#define STP_BDR_MAXDEGREE 4
#define STP_BDR_MAXDNEDGES 6

/** bd_k test for given Steiner bottleneck distances todo bd5 */
SCIP_RETCODE reduce_bdr(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   GRAPH*                netgraph,           /**< auxiliary graph structure */
   PATH*                 netmst,             /**< MST structure */
   PATH*                 vnoi,               /**< path structure */
   SCIP_Real*            mstsdist,           /**< MST distance in aux-graph */
   SCIP_Real*            termdist1,          /**< dist array */
   SCIP_Real*            termdist2,          /**< second dist array */
   int*                  vbase,              /**< bases for nearest terminals */
   int*                  nodesid,            /**< nodes identification array */
   int*                  neighbterms1,       /**< neighbour terminals array */
   int*                  neighbterms2,       /**< second neighbour terminals array */
   int*                  nelims              /**< number of eliminations */
   )
{
   SCIP_Real cutoffs[STP_BDR_MAXDNEDGES];
   SCIP_Real sd[STP_BDR_MAXDEGREE];
   SCIP_Real ecost[STP_BDR_MAXDEGREE];
   int adjvert[STP_BDR_MAXDEGREE];
   GRAPH* auxg;
   PATH* mst;
   SCIP_Real csum;
   SCIP_Real mstcost;
   const int nnodes = g->knots;
   const SCIP_Bool pc = g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG;

   assert(g != NULL);
   assert(netgraph  != NULL);
   assert(netmst != NULL);

   /* initialize mst struct and new graph for bd4 test */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, STP_BDR_MAXDEGREE) );
   SCIP_CALL( graph_init(scip, &auxg, STP_BDR_MAXDEGREE, 2 * STP_BDR_MAXDNEDGES, 1) );

   for( int k = 0; k < STP_BDR_MAXDEGREE; k++ )
      graph_knot_add(auxg, -1);

   for( int k = 0; k < STP_BDR_MAXDEGREE - 1; k++ )
      for( int k2 = STP_BDR_MAXDEGREE - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   assert(auxg->edges == 2 * STP_BDR_MAXDNEDGES);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   SCIPdebugMessage("BD3-R Reduction: ");

   for( int i = 0; i < STP_BDR_MAXDEGREE; i++ )
      sd[i] = 0.0;

   if( !pc )
      for( int i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

   for( int degree = 3; degree <= STP_BDR_MAXDEGREE; degree ++ )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         int k;

         if( Is_term(g->term[i]) || g->grad[i] != degree )
            continue;

         assert(g->mark[i]);

         k = 0;
         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            ecost[k] = g->cost[e];
            adjvert[k++] = g->head[e];
         }

         assert(k == degree);

         /* vertex of degree 3? */
         if( degree == 3 )
         {
            csum = ecost[0] + ecost[1] + ecost[2];

            sd[0] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, ecost[0] + ecost[1], vbase, nodesid, neighbterms1, neighbterms2, adjvert[0],
                  adjvert[1], 300);
            sd[1] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, ecost[1] + ecost[2], vbase, nodesid, neighbterms1, neighbterms2, adjvert[1],
                  adjvert[2], 300);
            sd[2] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, ecost[2] + ecost[0], vbase, nodesid, neighbterms1, neighbterms2, adjvert[2],
                  adjvert[0], 300);

            if( sd[0] + sd[1] <= csum || sd[0] + sd[2] <= csum || sd[1] + sd[2] <= csum )
            {
               SCIP_Bool success;

               cutoffs[0] = sd[0];
               cutoffs[1] = sd[2];
               cutoffs[2] = sd[1];

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               assert(success);
               assert(g->grad[i] == 0);

               SCIPdebugMessage("BD3-R Reduction: %f %f %f csum: %f\n ", sd[0], sd[1], sd[2], csum);
               (*nelims)++;
            }
         }
         /* vertex of degree 4? */
         else if( degree == 4 )
         {
            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = TRUE;
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     auxg->cost[e] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, ecost[k] + ecost[k2], vbase, nodesid, neighbterms1,
                           neighbterms2, adjvert[k], adjvert[k2], 200);
                     auxg->cost[flipedge(e)] = auxg->cost[e];
                  }
               }
            }

            for( int l = 0; l < 4; l++ )
               mst[l].dist = UNKNOWN;

            /* compute mst on all neighbours */
            graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
            mstcost = 0.0;

#ifndef NDEBUG
            for( int l = 1; l < 4; l++ )
               assert(mst[l].dist != UNKNOWN);
#endif

            for( int l = 1; l < 4; l++ )
               mstcost += mst[l].dist;

            k = UNKNOWN;
            csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

            if( csum >= mstcost )
            {
               /* compute mst on all 3-subsets of all neigbours */
               for( k = 0; k < 4; k++ )
               {
                  auxg->mark[k] = FALSE;
                  mstcost = 0.0;

                  if( k != 0 )
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                     for( int l = 1; l < 4; l++ )
                        if( auxg->mark[l] )
                           mstcost += mst[l].dist;
                  }
                  else
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                     for( int l = 2; l < 4; l++ )
                        mstcost += mst[l].dist;
                  }

                  auxg->mark[k] = TRUE;
                  csum -= ecost[k];

                  if( csum < mstcost )
                     break;

                  csum += ecost[k];
               }
            }

            if( k == 4 )
            {

               int edgecount = 0;
               SCIP_Bool success;
               for( k = 0; k < 3; k++ )
               {
                  for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
                  {
                     const int k2 = auxg->head[e];
                     if( k2 > k )
                        cutoffs[edgecount++] = auxg->cost[e];
                  }
               }

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               if( success )
                  (*nelims)++;
            }
         }
      }
   }

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


#define STP_BD_MAXDEGREE 4
#define STP_BD_MAXDNEDGES 6


/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Bottleneck Degree 3,4 Test
 */
SCIP_RETCODE reduce_bd34(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   PATH*                 pathtail,           /**< array for internal use */
   PATH*                 pathhead,           /**< array for internal use */
   int*                  heap,               /**< array for internal use */
   int*                  statetail,          /**< array for internal use */
   int*                  statehead,          /**< array for internal use */
   int*                  memlbltail,         /**< array for internal use */
   int*                  memlblhead,         /**< array for internal use */
   int*                  nelims,             /**< point to return number of eliminations */
   int                   limit               /**< limit for edges to consider for each vertex */
   )
{
   SCIP_Real cutoffs[STP_BD_MAXDNEDGES];
   PATH mst[STP_BD_MAXDEGREE];
   SCIP_Real sd[STP_BD_MAXDEGREE];
   SCIP_Real ecost[STP_BD_MAXDEGREE];
   GRAPH* auxg;
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   int adjvert[STP_BD_MAXDEGREE];
   const int nnodes = g->knots;
   const SCIP_Bool pc = g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG;

   SCIPdebugMessage("BD34-Reduction: ");

   assert(scip != NULL);
   assert(g != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);

   /* initialize new graph for bd4 tests */
   SCIP_CALL( graph_init(scip, &auxg, STP_BD_MAXDEGREE, 2 * STP_BD_MAXDNEDGES, 1) );

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
      graph_knot_add(auxg, -1);

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
      for( int k2 = STP_BD_MAXDEGREE - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   *nelims = 0;

   for( int i = 0; i < STP_BD_MAXDEGREE; i++ )
      sd[i] = 0.0;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );
   }
   else
   {
      for( int  i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
      if( pc )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }

   for( int degree = 3; degree <= STP_BD_MAXDEGREE; degree++ )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         int edgecount;
         if( Is_term(g->term[i]) || g->grad[i] != degree )
            continue;

         edgecount = 0;
         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            ecost[edgecount] = g->cost[e];
            adjvert[edgecount++] = g->head[e];
         }

         /* vertex of degree 3? */
         if( degree == 3 )
         {
            const SCIP_Real csum = ecost[0] + ecost[1] + ecost[2];

            assert(edgecount == 3);

            if( pc )
            {
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[0]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[0], adjvert[1], limit));
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[1]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[1], adjvert[2], limit));
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[2]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[2], adjvert[0], limit));
            }
            else
            {
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[0]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[0], adjvert[1], limit, pc, FALSE));
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[1]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[1], adjvert[2], limit, pc, FALSE));
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[2]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[2], adjvert[0], limit, pc, FALSE));
            }

            if( SCIPisLE(scip, sd[0] + sd[1], csum) || SCIPisLE(scip, sd[0] + sd[2], csum) || SCIPisLE(scip, sd[1] + sd[2], csum) )
            {
               SCIP_Bool success;

               cutoffs[0] = sd[0];
               cutoffs[1] = sd[2];
               cutoffs[2] = sd[1];

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));
               assert(success);
               assert(g->grad[i] == 0);

               SCIPdebugMessage("BD3 Reduction: %f %f %f csum: %f\n ", sd[0], sd[1], sd[2], csum);
               (*nelims)++;
            }
         }
         /* vertex of degree 4? */
         else if( degree == 4 )
         {
            int k;
            SCIP_Real csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];
            SCIP_Real mstcost;

            assert(edgecount == 4);

            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = TRUE;
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     SCIP_Real s1 = -1.0;
                     if( pc )
                        SCIP_CALL(
                              reduce_getSdPcMw(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                                    pathmaxnodetail, pathmaxnodehead, adjvert[k], adjvert[k2], limit));
                     else
                        SCIP_CALL(
                              reduce_getSd(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[k], adjvert[k2], limit, pc, FALSE));
                     assert(s1 >= 0);
                     auxg->cost[e] = s1;
                     auxg->cost[flipedge(e)] = s1;
                  }
               }
            }

            for( int l = 0; l < 4; l++ )
               mst[l].dist = UNKNOWN;

            /* compute mst on all neighbours */
            graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
            mstcost = 0.0;

            /* sum cost of (root == 0) */
            for( int l = 1; l < 4; l++ )
            {
               assert(mst[l].dist != UNKNOWN);
               mstcost += mst[l].dist;
            }

            k = 0;
            if( SCIPisGE(scip, csum, mstcost) )
            {
               /* compute mst on all 3-subsets of all neighbors */
               for( k = 0; k < 4; k++ )
               {
                  auxg->mark[k] = FALSE;
                  mstcost = 0.0;

                  if( k != 0 )
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                     for( int l = 1; l < 4; l++ )
                        if( auxg->mark[l] )
                           mstcost += mst[l].dist;
                  }
                  else
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                     for( int l = 2; l < 4; l++ )
                        mstcost += mst[l].dist;
                  }

                  auxg->mark[k] = TRUE;

                  if( SCIPisLT(scip, csum - ecost[k], mstcost) )
                     break;
               }
            }

            if( k == 4 )
            {
               SCIP_Bool success;
               edgecount = 0;

               for( k = 0; k < 3; k++ )
               {
                  for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
                  {
                     const int k2 = auxg->head[e];
                     if( k2 > k )
                        cutoffs[edgecount++] = auxg->cost[e];
                  }
               }

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               if( success )
                  (*nelims)++;
            }
         }
      }
   }

   /* free memory */

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   SCIPdebugMessage("bd34: %d nodes deleted\n", *nelims);

   assert(graph_valid(g));

   return SCIP_OKAY;
}


/** Non-Terminal-Set test*/
SCIP_RETCODE reduce_nts(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   PATH*                 pathtail,           /**< array for internal use */
   PATH*                 pathhead,           /**< array for internal use */
   int*                  heap,               /**< array for internal use */
   int*                  statetail,          /**< array for internal use */
   int*                  statehead,          /**< array for internal use */
   int*                  memlbltail,         /**< array for internal use */
   int*                  memlblhead,         /**< array for internal use */
   int*                  nelims,             /**< point to return number of eliminations */
   int                   limit               /**< limit for edges to consider for each vertex */
   )
{
   PATH mst[STP_BD_MAXDEGREE];
   SCIP_Real ecost[STP_BD_MAXDEGREE];
   int outedge[STP_BD_MAXDEGREE];
   int adjvert[STP_BD_MAXDEGREE];
   GRAPH* auxg;
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool pc = g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG;

   SCIPdebugMessage("NTS-Reduction: ");

   assert(scip != NULL);
   assert(g != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);
   assert(limit > 0);

   // todo extra method, also called from bd34

   /* initialize new graph for nts tests */
   SCIP_CALL( graph_init(scip, &auxg, STP_BD_MAXDEGREE, 2 * STP_BD_MAXDNEDGES, 1) );

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
   {
      graph_knot_add(auxg, -1);
      auxg->mark[k] = TRUE;
   }

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
      for( int k2 = STP_BD_MAXDEGREE - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   *nelims = 0;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );
   }
   else
   {
      for( int  i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
      if( pc )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }

   for( int degree = 4; degree <= STP_BD_MAXDEGREE; degree++ )
   {
      for( int edge = 0; edge < nedges; edge += 2 )
      {
         /* edge not deleted? */
         if( g->oeat[edge] != EAT_FREE )
         {
            SCIP_Real csum;
            SCIP_Real mstcost;
            const SCIP_Real edgecost = g->cost[edge];
            int edgecount;
            const int tail = g->tail[edge];
            const int head = g->head[edge];
            SCIP_Bool success = TRUE;
            SCIP_Bool duplicate;

            assert(edge >= 0);
            assert(tail >= 0 && head >= 0);

            // todo
            if( g->grad[tail] != 3 || g->grad[head] != 3 )
               continue;

            if( Is_term(g->term[tail]) || Is_term(g->term[head]) )
               continue;

            edgecount = 0;
            for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == head )
                  continue;
               outedge[edgecount] = e;
               ecost[edgecount] = g->cost[e];
               adjvert[edgecount++] = g->head[e];
            }

            for( int e = g->outbeg[head]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == tail )
                  continue;
               outedge[edgecount] = e;
               ecost[edgecount] = g->cost[e];
               adjvert[edgecount++] = g->head[e];
            }

            assert(edgecount == degree);

            duplicate = FALSE;

            /* check for shared neighbor */
            for( int i = 0; i < degree - 1; i++ )
               for( int j = i + 1; j < degree; j++ )
                  if( adjvert[i] == adjvert[j] && adjvert[j] >= 0 )
                  {
                     adjvert[j] = -adjvert[j] - 1;
                     duplicate = TRUE;
                  }

            // todo reshuffle
            if( duplicate )
            {
               continue;
            }

            assert(g->mark[head] && g->mark[tail]);

            g->mark[head] = FALSE;
            g->mark[tail] = FALSE;

            csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

            /* compute sd values and check 2-subsets of neighbors */
            for( int k = 0; k < degree - 1 && success; k++ )
            {
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     SCIP_Real s1 = -1.0;

                     if( pc )
                        SCIP_CALL( reduce_getSdPcMw(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                                    pathmaxnodetail, pathmaxnodehead, adjvert[k], adjvert[k2], limit) );
                     else
                        SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                              adjvert[k], adjvert[k2], limit, pc, FALSE) );

                     assert(s1 >= 0);

                     if( g->tail[outedge[k]] != g->tail[outedge[k2]] )
                     {
                        const SCIP_Real innercost = ecost[k] + ecost[k2] + edgecost;

                        if( SCIPisGT(scip, s1, innercost) )
                        {
          //                 success = FALSE;
            //               break;
                        }
                     }

                     auxg->cost[e] = s1;
                     auxg->cost[flipedge(e)] = s1;
                  }
               }
            }

            g->mark[head] = TRUE;
            g->mark[tail] = TRUE;

            for( int l = 0; l < 4; l++ )
               mst[l].dist = UNKNOWN;

            /* compute MST on all neighbors */
            graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
            mstcost = 0.0;

            /* sum cost of (root == 0) */
            for( int l = 1; l < 4; l++ )
            {
               assert(mst[l].dist != UNKNOWN);
               mstcost += mst[l].dist;
            }

            success = (success && SCIPisGE(scip, csum, mstcost));


            if( success && 0)
            {
               /* compute MST on all 3-subsets of all neighbors */
               for( int k = 0; k < degree; k++ )
               {
                  assert(auxg->mark[k]);

                  auxg->mark[k] = FALSE;
                  mstcost = 0.0;

                  if( k != 0 )
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                     for( int l = 1; l < 4; l++ )
                        if( auxg->mark[l] )
                           mstcost += mst[l].dist;
                  }
                  else
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                     for( int l = 2; l < 4; l++ )
                        mstcost += mst[l].dist;
                  }

                  auxg->mark[k] = TRUE;

                  if( SCIPisLT(scip, csum - ecost[k], mstcost) )
                  {
                     success = FALSE;
                     break;
                  }
               }
            }

            if( success )
            {
             //  graph_edge_del(scip, g, edge, TRUE);

               (*nelims)++;
            }
         } /* check single edge */
      } /* for all edges */
   } /* for all degrees */


   /* free memory */

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   assert(graph_valid(g));

   return SCIP_OKAY;
}


/* shortest link reduction */
SCIP_RETCODE reduce_sl(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   double*               fixed,              /**< offset pointer */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< shortest path array */
   int*                  vbase,              /**< Voronoi/shortest path bases array */
   int*                  vrnodes,            /**< Voronoi/shortest path array  */
   STP_Bool*             visited,            /**< Voronoi/shortest path array */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                  or NULL */
   int*                  nelims              /**< pointer to store number of eliminations */
   )
{
   SCIP_QUEUE* queue;
   SCIP_Real cost;
   SCIP_Real mincost2;
   SCIP_Real mincost3;
   int     i;
   int     k;
   int     e;
   int     j;
   int     t;
   int     old;
   int     head;
   int     tail;
   int     root;
   int     nnodes;
   int     vrcount;
   int     minedge;
   int*    qnode;
   STP_Bool    contract;
   STP_Bool    foundterms;
   STP_Bool*   forbidden;
   STP_Bool*   newterm;
   SCIP_Bool pc;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(vrnodes != NULL);
   assert(visited != NULL);

   *nelims = 0;
   foundterms = FALSE;
   nnodes = g->knots;
   root = g->source;
   pc = (g->stp_type == STP_PCSPG) || (g->stp_type == STP_RPCSPG);

   if( g->terms <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &forbidden, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newterm, nnodes) );

   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

   graph_voronoiTerms(scip, g, g->cost, vnoi, vbase, heap, state);

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );
   for( j = 0; j < nnodes; j++ )
   {
      newterm[j] = FALSE;
      forbidden[j] = FALSE;
      visited[j] = FALSE;
   }
   for( i = 0; i < nnodes; i++ )
   {
      /* is i terminal and not disabled? */
      if( Is_term(g->term[i]) && g->mark[i] && !forbidden[i] )
      {
         /* traverse voronoi-region of (terminal) i */
         assert(SCIPqueueIsEmpty(queue));
         t = i;
         SCIP_CALL( SCIPqueueInsert(queue, &t) );
         vrcount = 1;
         vrnodes[0] = i;
         visited[i] = TRUE;
         minedge = UNKNOWN;
         mincost2 = FARAWAY;
         mincost3 = FARAWAY;

         while( !SCIPqueueIsEmpty(queue) )
         {
            qnode = (SCIPqueueRemove(queue));
            /* traverse all adjacent edges */
            for( e = g->outbeg[*qnode]; e != EAT_LAST; e = g->oeat[e] )
            {
               j = g->head[e];

               if( !g->mark[j] )
                  continue;

               k = vbase[j];
               assert(k != UNKNOWN);
               if( !visited[j] && k == i )
               {
                  visited[j] = TRUE;
                  vrnodes[vrcount++] = j;
                  SCIP_CALL( SCIPqueueInsert(queue, &(g->head[e])) );
               }
               else if( k != i )
                  /* update shortest and second shortest edge (cost) leaving the voronoi region */
               {
                  cost = g->cost[e];
                  if( minedge == UNKNOWN )
                  {
                     minedge = e;
                  }
                  else if( SCIPisLE(scip, cost, g->cost[minedge]) )
                  {
                     mincost3 = mincost2;
                     mincost2 = g->cost[minedge];
                     minedge = e;
                  }
                  else if( SCIPisLE(scip, cost, mincost2) )
                  {
                     mincost3 = mincost2;
                     mincost2 = g->cost[e];
                  }
                  else if( SCIPisLT(scip, cost, mincost3) )
                  {
                     mincost3 = g->cost[e];
                  }
               }
            }
         }
         for( j = 0; j < vrcount; j++ )
            visited[vrnodes[j]] = FALSE;
         if( minedge == UNKNOWN )
            continue;
         e = minedge;
         tail = g->tail[e];
         head = g->head[e];
         assert(vbase[tail] == i);

         contract = FALSE;
         cost = vnoi[tail].dist + g->cost[e] + vnoi[head].dist;
         if( SCIPisGE(scip, mincost2, cost) )
         {
            contract = TRUE;
         }

         /* check whether minedge can be removed */
         if( contract )
         {
            if( pc )
            {
               if( root != vbase[head] && !SCIPisLE(scip, g->cost[e] + vnoi[tail].dist + vnoi[head].dist, g->prize[vbase[head]]) )
                  continue;
               if( i == tail )
               {
                  if( root != i && !SCIPisLE(scip, vnoi[tail].dist + g->cost[e], g->prize[i]) )
                     continue;
               }
               else
               {
                  if( root != i && !SCIPisLT(scip, vnoi[tail].dist + g->cost[e], g->prize[i]) )
                     continue;
               }
               if( Is_term(g->term[head]) && Is_term(g->term[tail]) )
                  continue;
            }

            *fixed += g->cost[e];
            assert(g->mark[tail] && g->mark[head]);
            assert(!Is_pterm(g->term[tail]) && !Is_pterm(g->term[head]));

            if( Is_term(g->term[head]) )
            {
               j = head;
               k = tail;
            }
            else
            {
               j = tail;
               k = head;
            }

            old = g->grad[j] + g->grad[k] - 1;

            if( pc )
            {
               SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, j, k, i) );
            }
            else
            {
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e], NULL) );
               SCIP_CALL( graph_knot_contract(scip, g, solnode, j, k) );

               assert(g->grad[k] == 0 && g->grad[j] >= 0);

               if( !Is_term(g->term[j]) )
               {
                  newterm[j] = TRUE;
                  foundterms = TRUE;
               }
            }

            assert(old - g->grad[j] - g->grad[k] > 0);
            (*nelims) += old - g->grad[j] - g->grad[k];
            forbidden[vbase[j]] = TRUE;
            forbidden[vbase[k]] = TRUE;
         }
      }
   }

   for( i = 0; i < nnodes && foundterms; i++ )
      if( newterm[i] && !Is_term(g->term[i]) && g->grad[i] > 0 )
         graph_knot_chg(g, i, 0);


   /* free memory */
   SCIPqueueFree(&queue);
   SCIPfreeBufferArray(scip, &newterm);
   SCIPfreeBufferArray(scip, &forbidden);

   return SCIP_OKAY;
}


/* NV reduction from T. Polzin's "Algorithms for the Steiner problem in networks" */
SCIP_RETCODE reduce_nv(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   double*               fixed,              /**< offset pointer */
   int*                  edgearrint,         /**< edge int array for internal computations */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array for internal computations */
   int*                  vbase,              /**< array for internal computations */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                  or NULL */
   int* nelims                               /**< pointer to store number of eliminations */
   )
{
   SCIP_Real* distance;
   SCIP_Real  min1;
   SCIP_Real  min2;
   SCIP_Real  pi;
   SCIP_Real  pt;
   SCIP_Real  ttdist;
   int* term;
   int*    minedge1;
   int*    distnode;
   int     i;
   int     l;
   int     k;
   int     e;
   int     t;
   int     old;
   int     edge1;
   int     nnodes;
   int     nterms;
   int     mingrad;
   int     termcount;
   SCIP_Bool pc;
   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   t = 0;
   termcount = 0;
   *nelims = 0;
   pi = 0;
   pt = 0;

   nnodes = g->knots;
   nterms = g->terms;
   pc = (g->stp_type == STP_PCSPG) || (g->stp_type == STP_RPCSPG);
   SCIP_CALL( SCIPallocBufferArray(scip, &term, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distance, nnodes) );

   /* minimal grad of a vertex to be scrutinized */
   if( pc )
   {
      if( g->stp_type == STP_RPCSPG )
         mingrad = 3;
      else
         mingrad = 4;

      SCIP_CALL( SCIPallocBufferArray(scip, &distnode, nnodes) );
   }
   else
   {
      mingrad = 2;
      distnode = NULL;
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && g->mark[i] && g->grad[i] > 0 )
      {
         /* compute shortest incident edge */
         edge1 = UNKNOWN;
         if( g->grad[i] >= 1 )
         {
            min1  = FARAWAY;

            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( !g->mark[g->head[e]] )
                  continue;
               if( SCIPisLE(scip, g->cost[e], min1) )
               {
                  edge1 = e;
                  min1 = g->cost[e];
               }
            }
         }
         minedge1[termcount] = edge1;
         term[termcount++] = i;
      }
   }

   /* compute Voronoi regions and distances */
   SCIP_CALL( graph_voronoiWithDist(scip, g, g->cost, distance, edgearrint, vbase, minedge1, heap, state, distnode, vnoi) );

   for( l = 0; l < termcount; l++ )
   {
      /* get l'th terminal */
      i = term[l];

      if( g->grad[i] < mingrad )
         continue;

      assert(minedge1[l] != UNKNOWN);
      /* get shortest two edges */
      edge1 = UNKNOWN;
      min2 = FARAWAY;
      min1 = FARAWAY;
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( !g->mark[g->head[e]] )
            continue;
         if( SCIPisLE(scip, g->cost[e], min1) )
         {
            edge1 = e;
            min2 = min1;
            min1 = g->cost[e];
         }
         else if( SCIPisLE(scip, g->cost[e], min2) )
         {
            min2 = g->cost[e];
         }
      }

      assert(edge1 != UNKNOWN);
      assert(i == g->tail[edge1]);
      k = g->head[edge1];

      /* covered in degree test */
      if( Is_term(g->term[k]) )
         continue;

      if( vbase[k] != i )
      {
         if( pc )
            t = vbase[k];
         ttdist = g->cost[edge1] + vnoi[k].dist;
      }
      else
      {
         if( distnode != NULL )
            t = distnode[i];
         ttdist = distance[i];
      }
      if( pc )
      {
         if( i != g->source )
            pi = g->prize[i];
         else
            pi = FARAWAY;

         if( t != g->source )
            pt = g->prize[t];
         else
            pt = FARAWAY;
      }

      if( SCIPisGE(scip, min2, ttdist)
         && (!pc || (SCIPisLE(scip, g->cost[edge1], pi) && SCIPisLE(scip, ttdist, pt))) )
      {
         old = g->grad[i] + g->grad[k] - 1;
         *fixed += g->cost[edge1];

         if( pc )
         {
            SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i, k, i) );
         }
         else
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1], NULL) );
            SCIP_CALL( graph_knot_contract(scip, g, solnode, i, k) );
         }
         assert(old - g->grad[i] - g->grad[k] > 0);
         (*nelims) += old - g->grad[i] - g->grad[k];
      }
   }

   SCIPfreeBufferArrayNull(scip, &distnode);
   SCIPfreeBufferArray(scip, &distance);
   SCIPfreeBufferArray(scip, &minedge1);
   SCIPfreeBufferArray(scip, &term);

   return SCIP_OKAY;
}


/* advanced NV reduction */
SCIP_RETCODE reduce_nvAdv(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            distance,           /**< nodes-sized distance array */
   double*               fixed,              /**< offset pointer */
   int*                  edgearrint,         /**< edges-sized array */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< Voronoi base array  */
   int*                  neighb,             /**< nodes-sized neighborhood array  */
   int*                  distnode,           /**< nodes-sized distance array  */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                  or NULL */
   int*                  nelims              /**< pointer to store number of eliminations */
   )
{
   SCIP_Real  min1;
   SCIP_Real  min2;
   SCIP_Real  min3;
   SCIP_Real  pi;
   SCIP_Real  pt;
   SCIP_Real  ttdist;
   int* term;
   int*    minedge1;
   int     i;
   int     l;
   int     k;
   int     e;
   int     t;
   int     edge1;
   int     edge2;
   int     nnodes;
   int     nterms;
   int     mingrad;
   int     termcount;
   SCIP_Bool pc;
   SCIP_Bool contract;

   assert(g != NULL);
   assert(neighb != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   t = 0;
   pi = 0;
   pt = 0;
   pc = (g->stp_type == STP_PCSPG) || (g->stp_type == STP_RPCSPG);
   nnodes = g->knots;
   nterms = g->terms;
   *nelims = 0;
   termcount = 0;

   if( nterms <= 1 )
      return SCIP_OKAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &term, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nterms) );

   /* set minimal grad of a vertex to be scrutinized */
   if( pc )
   {
      if( g->stp_type == STP_RPCSPG )
         mingrad = 3;
      else
         mingrad = 4;

      assert(distnode != NULL);
   }
   else
   {
      mingrad = 2;
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   /* compute shortest incident edge to each terminal */
   for( i = 0; i < nnodes; i++ )
   {
      neighb[i] = FALSE;
      if( Is_term(g->term[i]) && g->mark[i] && g->grad[i] > 0 )
      {
         /* compute shortest incident edge */
         edge1 = UNKNOWN;
         if( g->grad[i] >= 1 )
         {
            min1  = FARAWAY;

            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->mark[g->head[e]] && SCIPisLE(scip, g->cost[e], min1) )
               {
                  edge1 = e;
                  min1 = g->cost[e];
               }
            }
         }

         minedge1[termcount] = edge1;
         term[termcount++] = i;
      }
   }

   /* compute Voronoi regions and distances */
   SCIP_CALL( graph_voronoiWithDist(scip, g, g->cost, distance, edgearrint, vbase, minedge1, heap, state, distnode, vnoi) );

   /* main loop: try to contract (shortest) edges into terminals */
   for( l = 0; l < termcount; l++ )
   {
      /* get l'th terminal */
      i = term[l];

      if( g->grad[i] < mingrad )
         continue;

      assert(minedge1[l] != UNKNOWN);

      /* get shortest two edges */

      min3 = FARAWAY;
      min2 = FARAWAY;
      min1 = FARAWAY;
      edge1 = UNKNOWN;
      edge2 = UNKNOWN;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( !g->mark[g->head[e]] )
            continue;
         neighb[g->head[e]] = TRUE;

         if( SCIPisLE(scip, g->cost[e], min1) )
         {
            edge2 = edge1;
            edge1 = e;

            min3 = min2;
            min2 = min1;
            min1 = g->cost[e];
         }
         else if( SCIPisLE(scip, g->cost[e], min2) )
         {
            edge2 = e;

            min3 = min2;
            min2 = g->cost[e];
         }
         else if( SCIPisLE(scip, g->cost[e], min3) )
         {
            min3 = g->cost[e];
         }
      }

      assert(edge1 != UNKNOWN);
      assert(i == g->tail[edge1]);

      k = g->head[edge1];

      /* covered in degree test */
      if( Is_term(g->term[k]) )
         continue;

      if( vbase[k] != i )
      {
         if( pc )
            t = vbase[k];
         ttdist = g->cost[edge1] + vnoi[k].dist;
      }
      else
      {
         if( pc )
         {
            assert(distnode != NULL);
            t = distnode[i];
         }
         ttdist = distance[i];
      }
      if( pc )
      {
         if( i != g->source )
            pi = g->prize[i];
         else
            pi = FARAWAY;

         if( t == UNKNOWN )
            pt = -1.0;
         else if( t != g->source )
            pt = g->prize[t];
         else
            pt = FARAWAY;
      }

      contract = FALSE;

      if( SCIPisGE(scip, min2, ttdist) )
      {
         contract = TRUE;
      }
      else if( edge2 != UNKNOWN && !Is_term(g->term[g->head[edge2]]) && SCIPisGE(scip, min3, ttdist) )
      {
         t = g->head[edge2];
         for( e = g->outbeg[t]; e != EAT_LAST; e = g->oeat[e] )
            if( e != flipedge(edge2) && SCIPisLT(scip, g->cost[e], ttdist) )/*&& !neighb[g->head[e]] ) */
               break;

         if( e == EAT_LAST )
            contract = TRUE;
      }

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         neighb[g->head[e]] = FALSE;

      if( contract && (!pc || (SCIPisLE(scip, g->cost[edge1], pi) && SCIPisLE(scip, ttdist, pt))) )
      {
         (*nelims)++;
         *fixed += g->cost[edge1];

         if( pc )
         {
            SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i, k, i) );
         }
         else
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1], NULL) );
            SCIP_CALL( graph_knot_contract(scip, g, solnode, i, k) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &minedge1);
   SCIPfreeBufferArray(scip, &term);

   return SCIP_OKAY;
}


/*  longest edge reduction test from T. Polzin's "Algorithms for the Steiner problem in networks" (Lemma 20) */
SCIP_RETCODE reduce_ledge(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int* edgestate
)
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real cost;
   SCIP_Real maxcost;
   int v1;
   int v2;
   int k;
   int e;
   int ne;
   int nedges;
   int nnodes;
   int nterms;
   int maxnedges;
   int netnnodes;
   int* nodesid;
   int* edgeorg;
   SCIP_Bool checkstate = (edgestate != NULL);
   STP_Bool* blocked;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   *nelims = 0;
   nedges = g->edges;
   nnodes = g->knots;
   assert(graph_valid(g));

   nterms = 0;
   for( k = 0; k < nnodes; k++ )
   {
      g->mark[k] = (g->grad[k] > 0);
      if( Is_term(g->term[k]) && g->mark[k] )
         nterms++;
   }

   if( nterms <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &blocked, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgeorg, nedges / 2) );

   graph_voronoiTerms(scip, g, g->cost, vnoi, vbase, heap, state);

   if( nedges >= (nterms - 1) * nterms )
      maxnedges = (nterms - 1) * nterms;
   else
      maxnedges = nedges;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   e = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
         if( e == 0 )
            graph_knot_add(netgraph, 0);
         else
            graph_knot_add(netgraph, -1);
         nodesid[k] = e++;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   netnnodes = netgraph->knots;
   assert(netnnodes == e);
   assert(netnnodes == nterms);

   for( e = 0; e < nedges / 2; e++ )
   {
      blocked[e] = FALSE;
      edgeorg[e] = UNKNOWN;
   }

   for( k = 0; k < nnodes; k++ )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         v1 = vbase[k];
         assert(k == g->tail[e]);

         if( v1 != vbase[g->head[e]] )
         {
            v2 = vbase[g->head[e]];
            assert(Is_term(g->term[v1]));
            assert(Is_term(g->term[v2]));
            assert(nodesid[v1] >= 0);
            assert(nodesid[v2] >= 0);

            for( ne = netgraph->outbeg[nodesid[v1]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == nodesid[v2] )
                  break;

            cost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;
            /* edge exists? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->head[ne] == nodesid[v2]);
               assert(netgraph->tail[ne] == nodesid[v1]);
               if( SCIPisGT(scip, netgraph->cost[ne], cost) )
               {
                  netgraph->cost[ne]            = cost;
                  netgraph->cost[Edge_anti(ne)] = cost;
                  edgeorg[ne / 2] = e;
                  assert(ne <= maxnedges);
               }
            }
            else
            {
               edgeorg[netgraph->edges / 2] = e;
               graph_edge_add(scip, netgraph, nodesid[v1], nodesid[v2], cost, cost);
               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }
   netgraph->source = 0;

   assert(graph_valid(netgraph));

   for( k = 0; k < netnnodes; k++ )
      netgraph->mark[k] = TRUE;

   /* compute a MST on netgraph */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, netnnodes) );
   SCIP_CALL( graph_path_init(scip, netgraph) );
   graph_path_exec(scip, netgraph, MST_MODE, 0, netgraph->cost, mst);

   maxcost = -1;
   assert(mst[0].edge == -1);

   for( k = 1; k < netnnodes; k++ )
   {
      assert(netgraph->path_state[k] == CONNECT);
      e = mst[k].edge;
      assert(e >= 0);
      cost = netgraph->cost[e];
      if( SCIPisGT(scip, cost, maxcost) )
         maxcost = cost;

      ne = edgeorg[e / 2];
      blocked[ne / 2] = TRUE;
      for( v1 = g->head[ne]; v1 != vbase[v1]; v1 = g->tail[vnoi[v1].edge] )
         blocked[vnoi[v1].edge / 2] = TRUE;

      for( v1 = g->tail[ne]; v1 != vbase[v1]; v1 = g->tail[vnoi[v1].edge] )
         blocked[vnoi[v1].edge / 2] = TRUE;
      assert(e != EAT_LAST);
   }

   for( k = 0; k < nnodes; k++ )
   {
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
         assert(e >= 0);
         if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
         {
            e = g->oeat[e];
            continue;
         }

         if( SCIPisGE(scip, g->cost[e], maxcost) && !blocked[e / 2] )
         {
            (*nelims)++;
            v1 = g->oeat[e];
            graph_edge_del(scip, g, e, TRUE);
            e = v1;
         }
         else
         {
            e = g->oeat[e];
         }
      }
   }

   /* graph might have become disconnected */
   if( *nelims > 0 )
   {
      SCIP_CALL( level0(scip, g) );
   }

   /* free netgraph and  MST data structure */
   graph_path_exit(scip, netgraph);
   graph_free(scip, &netgraph, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &nodesid);
   SCIPfreeBufferArray(scip, &edgeorg);
   SCIPfreeBufferArray(scip, &blocked);

   assert(graph_valid(g));
   return SCIP_OKAY;
}


#if 0
/** domination vertex reduction for the SPG */
void reduce_alt_dv(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   STP_Bool*             marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   const int nnodes = g->knots;
   int nreds = 0;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_SPG);

   BMSclearMemoryArray(marked, nnodes);

   /* main loop */
   for( int k = 0; k < nnodes; k++ )
   {
      const SCIP_Real maxgrad = g->grad[k];

      if( maxgrad < 3 )
         continue;

      assert(g->mark[k]);

      /* mark adjacent vertices and k*/
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = TRUE;

      marked[k] = TRUE;

      /* check all other nodes */
      for( int i = 0; i < nnodes; i++ )
      {
         /* valid candidate? */
         if( !Is_term(g->term[i]) && g->grad[i] <= maxgrad && g->mark[i] && k != i )
         {
            SCIP_Real min;
            int e2;
            for( e2 = g->outbeg[i]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               if( !marked[g->head[e2]] )
                  break;

            /* neighbors of j subset of those of k? */
            if( e2 == EAT_LAST )
            {
#if 0
               int maxe = g->outbeg[i];
               for( e2 = g->outbeg[i]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  if( g->cost[e] > g->cost[maxe] )
                     maxe = e;

               min = 0.0;




               (*count) += g->grad[i];
               while( g->outbeg[i] != EAT_LAST )
                  graph_edge_del(scip, g, g->outbeg[j] TRUE);

               g->mark[i] = FALSE;
               marked[i] = FALSE;
#endif
            }
         }
      }

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      marked[k] = FALSE;

      for( int j = 0; j < nnodes; j++ )
         assert(marked[j] == FALSE);
   }

   *count += nreds;
}

#endif

/** adjacent neighbourhood reduction for the MWCSP */
void reduce_ans(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   const int nnodes = g->knots;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MWCSP);
   assert(graph_valid(g));

   *count = 0;

   /* unmark all nodes */
   for( int k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighborhood of all nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real min;
      int e;
      int neighborcount = 0;

      if( !g->mark[k] )
         continue;

      /* mark adjacent vertices and k*/
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = TRUE;

      marked[k] = TRUE;

      if( SCIPisLT(scip, g->prize[k], 0.0) )
         min = g->prize[k];
      else
         min = 0.0;

      /* check all neighbors of k */
      e = g->outbeg[k];
      while( e != EAT_LAST && neighborcount++ < STP_RED_ANSMAXNEIGHBORS )
      {
         const int j = g->head[e];
         e = g->oeat[e];

         /* valid candidate? */
         if( g->grad[j] <= g->grad[k] && !Is_term(g->term[j]) && g->mark[j] )
            ansProcessCandidate(scip, g, marked, count, min, j);
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      marked[k] = FALSE;

      for( int j = 0; j < nnodes; j++ )
         assert(marked[j] == FALSE);
   }

   assert(graph_valid(g));
}

/** advanced adjacent neighbourhood reduction for the MWCSP */
void reduce_ansAdv(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count,              /**< pointer to number of performed reductions */
   SCIP_Bool             extneigbhood        /**< use extended neighbour hood */
   )
{
   int candidates[MAX(STP_RED_ANSMAXCANDS, STP_RED_ANSEXMAXCANDS)];
   int neighbarr[STP_RED_CNSNN];
   const int nnodes = g->knots;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MWCSP);

   *count = 0;

   /* unmark all nodes */
   for( int k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighborhood of all nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real min;
      int nn;
      int ncands;
      int maxgrad;

      if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
         continue;

      nn = 0;

      /* mark adjacent vertices and k */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int j = g->head[e];
         marked[j] = TRUE;
         if( SCIPisGT(scip, g->prize[j], 0.0) && nn < STP_RED_CNSNN )
            neighbarr[nn++] = j;
      }

      marked[k] = TRUE;
      maxgrad = g->grad[k];
      ncands = 0;

      for( int l = 0; l < nn; l++ )
      {
         for( int e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
         {
            marked[g->head[e]] = TRUE;
         }
         maxgrad += g->grad[neighbarr[l]];
      }

      assert(SCIPisLE(scip, g->prize[k], 0.0));

      min = g->prize[k];

      if( extneigbhood )
      {
         assert(0 && "implement me");
      }
      else
      {
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int j = g->head[e];
            if( g->grad[j] <= maxgrad && g->mark[j] && !Is_term(g->term[j]) )
            {
               candidates[ncands++] = j;
               if( ncands >= STP_RED_ANSMAXCANDS )
               {
                  SCIPdebugMessage("REACHED ANS LIMIT %d \n", ncands);
                  break;
               }
            }
         }
      }

      /* check all neighbors of k */
      for( int l = 0; l < ncands; l++ )
         ansProcessCandidate(scip, g, marked, count, min, candidates[l]);

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      for( int l = 0; l < nn; l++ )
         for( int e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

      marked[k] = FALSE;

      for( int k2 = 0; k2 < nnodes; k2++ )
         assert(marked[k2] == FALSE);
   }

   assert(graph_valid(g));
}


/** alternative advanced adjacent neighbourhood reduction for the MWCSP */
void reduce_ansAdv2(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int neighbarr[STP_RED_CNSNN + 1];
   SCIP_Real min;
   const int nnodes = g->knots;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MWCSP);

   *count = 0;

   /* unmark all nodes */
   for( int k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

      /* check neighbourhood of all nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real maxprize;
      int neighbor0;

      if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
         continue;

      maxprize = g->prize[k];
      neighbor0 = -1;

      for( int run = 0; run < 2; run++ )
      {
         int e;
         int nn = 0;
         int k2 = UNKNOWN;
         int maxgrad = g->grad[k];

         /* mark adjacent vertices and k */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int j = g->head[e];
            marked[j] = TRUE;
            if( SCIPisGE(scip, g->prize[j], 0.0) && nn < STP_RED_CNSNN - 1 )
            {
               neighbarr[nn++] = j;
            }
            else if( SCIPisGT(scip, g->prize[j], maxprize) && g->mark[j] && j != neighbor0 )
            {
               maxprize = g->prize[j];
               k2 = j;
            }
         }

         marked[k] = TRUE;

         if( run == 0 && k2 != UNKNOWN )
            neighbarr[nn++] = k2;

         for( int l = 0; l < nn; l++ )
         {
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            {
               const int j = g->head[e];
               if( run == 1 && g->mark[j] && !Is_term(g->term[j]) && SCIPisGT(scip, g->prize[j], maxprize) && j != neighbor0 )
               {
                  maxprize = g->prize[j];
                  k2 = j;
               }
               marked[j] = TRUE;
            }
            maxgrad += g->grad[neighbarr[l]];
         }

         if( run == 1 && k2 != UNKNOWN )
         {
            maxgrad += g->grad[k2];
            neighbarr[nn++] = k2;

            for( e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = TRUE;
         }

         assert(SCIPisLE(scip, g->prize[k], 0.0));

         min = g->prize[k];
         if( k2 != UNKNOWN )
         {
            neighbor0 = k2;
            min += g->prize[k2];
         }

         /* check all neighbours of k */
         e = g->outbeg[k];
         while( e != EAT_LAST )
         {
            const int j = g->head[e];
            e = g->oeat[e];

            /* valid candidate? */
            if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) && k2 != j )
               ansProcessCandidate(scip, g, marked, count, min, j);
         }

         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

         for( int l = 0; l < nn; l++ )
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = FALSE;

         marked[k] = FALSE;

         for( k2 = 0; k2 < nnodes; k2++ )
            assert(marked[k2] == FALSE);
      }
   }

   assert(graph_valid(g));
}


/** advanced connected neighborhood subset reduction test for the MWCSP */
SCIP_RETCODE reduce_cnsAdv(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array for internal use */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   SCIP_Real kprize;
   int neighbarr[STP_RED_CNSNN + 1];
   int neighbarr2[STP_RED_CNSNN + 1];
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int nn;
   int nn2;
   int k2grad;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MWCSP);

   k2grad = 0;
   *count = 0;
   nnodes = g->knots;

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = VERTEX_OTHER;

   /* first run: consider node plus adjacent terminals */

   /* check neighborhood of all nodes */
   for( k = 0; k < nnodes; k++ )
   {
      if( !(g->mark[k]) || (g->grad[k] < 2) )
         continue;

      nn = 0;
      k2 = UNKNOWN;
      nn2 = 0;
      kprize = g->prize[k];
      maxgrad = g->grad[k];

      /* mark adjacent vertices and k */
      for (e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e])
      {
         j = g->head[e];

         if( !g->mark[j] )
            continue;

         if( SCIPisGE(scip, g->prize[j], 0.0) && nn < STP_RED_CNSNN - 1 )
         {
            neighbarr[nn++] = j;
            marked[j] = VERTEX_CONNECT;
         }
         else
         {
            marked[j] = VERTEX_NEIGHBOR;
         }
      }

      marked[k] = VERTEX_CONNECT;

      /* traverse all connected non-negative nodes and mark their neighbors */
      for (l = 0; l < nn; l++)
      {
         for (e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e])
         {
            j = g->head[e];
            if( !g->mark[j] )
               continue;

            if( marked[j] == VERTEX_OTHER )
               marked[j] = VERTEX_NEIGHBOR;
         }
         maxgrad += g->grad[neighbarr[l]] - 1;
      }

      if( Is_term(g->term[k]) )
         min = 0.0;
      else
         min = g->prize[k];

      /* traverse all vertices (main loop) */
      for (j = 0; j < nnodes; j++)
      {
         /* vertex part of the current connected subset? Or terminal? Or belonging to the extension of the graph? */
         if( marked[j] != VERTEX_CONNECT && g->mark[j] && !Is_term(g->term[j])
               &&
               /* valid candidate? */
               g->grad[j] <= maxgrad && SCIPisLE(scip, g->prize[j], min) )
         {
            for (e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2])
               if( marked[g->head[e2]] == VERTEX_OTHER )
                  break;

            /* neighbors of j subset of those of k? */
            if( e2 == EAT_LAST )
            {
               /* yes, delete vertex */
               while (g->outbeg[j] != EAT_LAST)
               {
                  e2 = g->outbeg[j];
                  (*count)++;
                  graph_edge_del(scip, g, e2, TRUE);
               }
               g->mark[j] = FALSE;
               marked[j] = VERTEX_OTHER;
            }
         }
      }

      for (e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e])
         marked[g->head[e]] = VERTEX_OTHER;

      for (l = 0; l < nn; l++)
         for (e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e])
            marked[g->head[e]] = VERTEX_OTHER;

      marked[k] = VERTEX_OTHER;

#ifdef DEBUG
      for( l = 0; l < nnodes; l++ )
      assert(marked[l] == VERTEX_OTHER);
#endif

   }
    /* second run: consider the same plus an additional (non-positive) vertex  */

   for (k = 0; k < nnodes; k++)
   {
      if( !(g->mark[k]) || g->grad[k] < 2 || Is_term(g->term[k]) )
         continue;

      nn = 0;
      k2 = UNKNOWN;
      nn2 = 0;
      kprize = g->prize[k];
      maxgrad = g->grad[k];

      /* mark adjacent vertices and k */
      for (e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e])
      {
         j = g->head[e];

         if( !g->mark[j] )
            continue;

         if( SCIPisGE(scip, g->prize[j], 0.0) && nn < STP_RED_CNSNN - 1 )
         {
            neighbarr[nn++] = j;
            marked[j] = VERTEX_CONNECT;
         }
         else if( (SCIPisGT(scip, g->prize[j], kprize) && nn2 < STP_RED_CNSNN)
               || (SCIPisGE(scip, g->prize[j], kprize) && j > k && nn2 < 3) )
         {
            neighbarr2[nn2++] = j;
            marked[j] = VERTEX_NEIGHBOR;
         }
         else
         {
            marked[j] = VERTEX_NEIGHBOR;
         }
      }

      marked[k] = VERTEX_CONNECT;

      /* traverse all connected non-negative nodes and mark their neighbors */
      for (l = 0; l < nn; l++)
      {
         for (e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e])
         {
            j = g->head[e];
            if( !g->mark[j] )
               continue;

            if( marked[j] == VERTEX_OTHER && nn2 < STP_RED_CNSNN
                  && (SCIPisGT(scip, g->prize[j], kprize)
                        || (SCIPisGE(scip, g->prize[j], kprize) && j > k
                              && nn2 < 3)) )
            {
               neighbarr2[nn2++] = j;
               marked[j] = VERTEX_NEIGHBOR;
            }
            else if( marked[j] == VERTEX_OTHER )
            {
               marked[j] = VERTEX_NEIGHBOR;
            }
         }
         maxgrad += g->grad[neighbarr[l]] - 1;
      }

      if( Is_term(g->term[k]) )
         min = 0.0;
      else
         min = g->prize[k];

      for (l = 0; l < nn2; l++)
      {
         k2 = neighbarr2[l];

         if( !g->mark[k2] )
            continue;

         for (e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e])
            if( marked[g->head[e]] == VERTEX_OTHER && g->mark[g->head[e]] )
               marked[g->head[e]] = VERTEX_TEMPNEIGHBOR;
         min += g->prize[k2];
         k2grad = g->grad[k2];
         maxgrad += k2grad - 1;
         assert(SCIPisLE(scip, g->prize[k2], 0.0));

         /* traverse all vertices (main loop) */
         for (j = 0; j < nnodes; j++)
         {
            /* vertex part of the current connected subset? Or terminal? Or belonging to the extension of the graph? */
            if( marked[j] != VERTEX_CONNECT && g->mark[j]
                  && !Is_term(g->term[j]) &&
                  /* valid candidate? */
                  g->grad[j] <= maxgrad && SCIPisLE(scip, g->prize[j], min) )
            {
               for (e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2])
                  if( marked[g->head[e2]] == VERTEX_OTHER )
                     break;

               /* neighbors of j subset of those of k? */
               if( e2 == EAT_LAST )
               {
                  /* yes, delete vertex */
                  while (g->outbeg[j] != EAT_LAST)
                  {
                     e2 = g->outbeg[j];
                     (*count)++;
                     graph_edge_del(scip, g, e2, TRUE);
                  }
                  g->mark[j] = FALSE;
                  marked[j] = VERTEX_OTHER;
               }
            }
         }
         for (e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e])
            if( marked[g->head[e]] == VERTEX_TEMPNEIGHBOR
                  && g->mark[g->head[e]] )
               marked[g->head[e]] = VERTEX_OTHER;
         min -= g->prize[k2];
         maxgrad -= k2grad - 1;

      }
      for (e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e])
         marked[g->head[e]] = VERTEX_OTHER;

      for (l = 0; l < nn; l++)
         for (e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e])
            marked[g->head[e]] = VERTEX_OTHER;

      marked[k] = VERTEX_OTHER;
#ifdef DEBUG
      for( l = 0; l < nnodes; l++ )
      assert(marked[l] == VERTEX_OTHER);
#endif
   }

   return SCIP_OKAY;
}


/** non-positive vertex reduction for the MWCSP */
SCIP_RETCODE reduce_npv(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit
   )
{
   GRAPH* auxg;
   PATH mst[5];
   int adjverts[5];

   SCIP_Real prize;
   SCIP_Real sdist0;
   SCIP_Real sdist1;
   SCIP_Real sdist2;

   const int nnodes = g->knots;;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   *nelims = 0;

   /* initialize arrays */
   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }


   /* --- NPV3 test --- */

   /* try to eliminate non-positive vertices of degree 3 */
   for( int i = 0; i < nnodes; i++ )
   {
      int k;

      assert(g->grad[i] >= 0);

      /* only non-positive vertices of degree 3 */
      if( !g->mark[i] || g->grad[i] != 3 || Is_term(g->term[i]) )
         continue;

      k = 0;
      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->head[e] != g->source);
         assert(k < 3);

         adjverts[k++] = g->head[e];
      }

      assert(k == 3);

      g->mark[i] = FALSE;

      prize = g->prize[i];
      SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &sdist0, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[0], adjverts[1], limit, FALSE, TRUE) );
      SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &sdist1, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[1], adjverts[2], limit, FALSE, TRUE) );
      SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &sdist2, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[2], adjverts[0], limit, FALSE, TRUE) );

      /* can vertex be deleted? */
      if( (SCIPisGE(scip, -sdist0 - sdist1, prize) && SCIPisGE(scip, -sdist2, prize))
         || (SCIPisGE(scip, -sdist1 - sdist2, prize) && SCIPisGE(scip, -sdist0, prize))
         || (SCIPisGE(scip, -sdist2 - sdist0, prize) && SCIPisGE(scip, -sdist1, prize))
         )
      {
         SCIPdebugMessage("npv3Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
         (*nelims) +=  g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }

   /* --- NPV4 test --- */

   /* initialize mst struct and new graph for further tests */
   SCIP_CALL( graph_init(scip, &auxg, 5, 40, 1) );

   for( int k = 0; k < 4; k++ )
      graph_knot_add(auxg, -1);
   for( int k = 0; k < 4; k++ )
      for( int k2 = k + 1; k2 < 4; k2++ )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 4 */
   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      int k;

      /* only non-positive vertices of degree 4 */
      if( !g->mark[i] || g->grad[i] != 4 || Is_term(g->term[i]) )
         continue;
      k = 0;

      /* store neighbours */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjverts[k++] = g->head[e];

      assert(k == 4);

      g->mark[i] = FALSE;
      prize = g->prize[i];

      /* compute mw bottleneck distance to each pair of neighbours */
      for( k = 0; k < 4; k++ )
      {
         auxg->mark[k] = TRUE;
         for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
         {
            const int k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &(sdist0), -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
               auxg->cost[e] = sdist0;
               if( SCIPisGT(scip, prize, -auxg->cost[e]) )
                  break;

               auxg->cost[flipedge(e)] = auxg->cost[e];
            }
         }
         if( e != EAT_LAST )
            break;
      }

      k = UNKNOWN;
      if( e == EAT_LAST )
      {
         /* compute mst on all neighbours */
         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);

         /* calculate mst cost */
         sdist0 = 0.0;
         for( int l = 1; l < 4; l++ )
            sdist0 += mst[l].dist;

         if( SCIPisLE(scip, prize, -sdist0) )
         {
            /* compute subset msts on all neighbours */
            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = FALSE;
               sdist0 = 0.0;
               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( int l = 1; l < 4; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( int l = 2; l < 4; l++ )
                     sdist0 += mst[l].dist;
               }
               auxg->mark[k] = TRUE;
               if( SCIPisGT(scip, prize, -sdist0) )
                  break;
            }
         }
      }

      /* can node be eliminated? */
      if( k == 4 )
      {
         SCIPdebugMessage("npv4Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
         (*nelims) += g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }

   /* --- NPV5 test --- */

   /* enlarge graph for NPV5 test*/
   graph_knot_add(auxg, -1);
   for( int k = 0; k < 4; k++ )
      graph_edge_add(scip, auxg, k, 4, 1.0, 1.0);
   graph_path_exit(scip, auxg);
   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 5 */
   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      int k;

      /* only non-positive vertices of degree 5 */
      if( !g->mark[i] || g->grad[i] != 5 || Is_term(g->term[i]) )
         continue;
      k = 0;

      /* store neighbours */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjverts[k++] = g->head[e];

      assert(k == 5);

      g->mark[i] = FALSE;
      prize = g->prize[i];

      for( k = 0; k < 5; k++ )
      {
         auxg->mark[k] = TRUE;
         for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
         {
            const int k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &(sdist0), -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
               auxg->cost[e] = sdist0;
               if( SCIPisGT(scip, prize, -auxg->cost[e]) )
                  break;

               auxg->cost[flipedge(e)] = auxg->cost[e];
            }
         }
         if( e != EAT_LAST )
            break;
      }
      k = UNKNOWN;
      if( e == EAT_LAST )
      {
         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
         sdist0 = 0.0;
         for( int l = 1; l < 5; l++ )
            sdist0 += mst[l].dist;

         if( SCIPisLE(scip, prize, -sdist0) )
         {
            for( k = 0; k < 5; k++ )
            {
               int k2 = UNKNOWN;
               auxg->mark[k] = FALSE;
               sdist0 = 0.0;
               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( int l = 1; l < 5; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( int l = 2; l < 5; l++ )
                     sdist0 += mst[l].dist;
               }

               if( SCIPisLE(scip, prize, -sdist0) )
               {
                  for( k2 = k + 1; k2 < 5; k2++ )
                  {
                     if( k2 == k )
                        continue;
                     auxg->mark[k2] = FALSE;
                     sdist0 = 0.0;
                     if( k2 != 0 && k != 0)
                     {
                        graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                        for( int l = 1; l < 5; l++ )
                           if( auxg->mark[l] )
                              sdist0 += mst[l].dist;
                     }
                     else
                     {
                        int s;
                        if( k != 1 && k2 != 1 )
                           s = 1;
                        else
                           s = 2;
                        graph_path_exec(scip, auxg, MST_MODE, s, auxg->cost, mst);
                        for( int l = 0; l < 5; l++ )
                           if( auxg->mark[l] && l != s  )
                              sdist0 += mst[l].dist;
                     }
                     auxg->mark[k2] = TRUE;
                     if( SCIPisGT(scip, prize, -sdist0) )
                        break;
                  }
               }
               auxg->mark[k] = TRUE;
               if( k2 != 5 )
                  break;
            }
         }
      }

      if( k == 5 )
      {
         SCIPdebugMessage(" \n npv5Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
         (*nelims) += g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);

      }
      else
      {
         g->mark[i] = TRUE;
      }
   }

   /* free memory*/
   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   return SCIP_OKAY;
}


/** chain reduction (NPV_2) for the MWCSP */
SCIP_RETCODE reduce_chain2(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit
   )
{
   SCIP_Real sdist;
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   *nelims = 0;
   nnodes = g->knots;

   /* initialize arrays */
   for( i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   for( i = 0; i < nnodes; i++ )
   {
      assert(g->grad[i] >= 0);
      if( !g->mark[i] || g->grad[i] == 0 || Is_term(g->term[i]) || g->grad[i] != 2 )
         continue;

      /* non-positive chains */
      e1 = g->outbeg[i];
      e2 = g->oeat[e1];
      i1 = g->head[e1];
      i2 = g->head[e2];

      assert(e1 >= 0);
      assert(e2 >= 0);
      assert(g->mark[i1]);
      assert(g->mark[i2]);
      g->mark[i] = FALSE;
      SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &sdist, -(g->prize[i]), heap, statetail, statehead, memlbltail, memlblhead, i1, i2, limit, FALSE, TRUE) );
      if( SCIPisGE(scip, -sdist, g->prize[i]) )
      {
         SCIPdebugMessage("delete : %d prize: %f sd: %f \n", i,  g->prize[i], -sdist );
         graph_edge_del(scip, g, e1, TRUE);
         graph_edge_del(scip, g, e2, TRUE);
         (*nelims) += 2;
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }
   return SCIP_OKAY;
}



/** non-negative path reduction for the MWCSP */
void reduce_nnp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int localcount = 0;
   const int nnodes = g->knots;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MWCSP);

   /* unmark all nodes */
   for( int k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighborhood of terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] || g->prize[k] < 0.0 )
         continue;

      /* mark adjacent vertices of k */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         if( g->mark[g->head[e]] )
            marked[g->head[e]] = TRUE;

      /* ... and traverse them */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int j = g->head[e];

         if( marked[j] )
         {
            int e2 = g->outbeg[j];

            while( e2 != EAT_LAST )
            {
               const int candedge = e2;
               e2 = g->oeat[e2];
               if( marked[g->head[candedge]] )
               {
                  graph_edge_del(scip, g, candedge, TRUE);
                  localcount++;
               }
            }
         }
      }

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      for( int j = 0; j < nnodes; j++ )
         assert(marked[j] == FALSE);
   }

   *count = localcount;

   assert(graph_valid(g));
}
