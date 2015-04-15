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

#define KNOTFREQ 100
#define KNOTLIMIT 1e+20

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

#if 0
static int compare(
   const double* pathdist,
   const double* pathtran,
   int          a,
   int          b)
{
   if (NE(pathdist[a], pathdist[b]))
      return(LT(pathdist[a], pathdist[b]) ? -1 : 1);

   if (EQ(pathtran[a], pathtran[b]))
      return(0);

   return(LT(pathtran[a], pathtran[b]) ? -1 : 1);
}
#endif


#if 0
/** for debug purposes only */
static
SCIP_RETCODE printGraph(
   SCIP* scip,
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   int*                  result
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;
   char* stnodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, graph->knots ) );

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->knots; e++ )
   {
      stnodes[e] = FALSE;
   }
   for( e = 0; e < graph->edges; e++ )
   {
      if( 1 )
      {
	 stnodes[graph->tail[e]] = TRUE;
	 stnodes[graph->head[e]] = TRUE;
      }
   }

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
      if( stnodes[n] )
      {
         if( n == graph->source[0] )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
            m = 1;
         }
         else if( graph->term[n] == 0 )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
            e += 1;
         }
         else
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
         }

      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e ++ )
   {
      if( 1 )
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);
	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      }
   }
   SCIPfreeBufferArray(scip, &stnodes);
   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

#endif

static int issmaller(
   const double* pathdist,
   const double* pathtran,
   int          a,
   int          b)
{
   return (LT(pathdist[a], pathdist[b]) || (!GT(pathdist[a], pathdist[b]) && LT(pathtran[a], pathtran[b])));
}
static int islarger(
   const double* pathdist,
   const double* pathtran,
   int          a,
   int          b)
{
   return (GT(pathdist[a], pathdist[b]) || (!LT(pathdist[a], pathdist[b]) && GT(pathtran[a], pathtran[b])));
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
      if (islarger(pathdist, pathtran, heap[2], heap[3]))
         c++;

   while((c <= (*count)) && islarger(pathdist, pathtran, heap[j], heap[c]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
         if (issmaller(pathdist, pathtran, heap[c + 1], heap[c]))
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

   while((j > 1) && (islarger(pathdist, pathtran, heap[c], heap[j])))
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
   const double* random,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   double* pathdist,
   double* pathtran,
   double* pathrand)
{
   int    k;
   int    m;
   int    i;
   int    done = 0;
   double tran;
   double dist;
   double temprand;

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
         k = nearest(heap, state, count, pathdist, pathtran);

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
                  temprand = pathrand[k] + random[i];
               }


               dist = Max(pathdist[k], pathtran[k] + cost[i]);

               if (LT(dist, pathdist[m])
                  || (EQ(dist, pathdist[m]) && LT(tran, pathtran[m])))
               {
                  pathdist[m] = dist;
                  pathtran[m] = tran;
                  pathrand[m] = temprand;

                  correct(heap, state, count, pathdist, pathtran, m);
               }
            }
         }
      }
   }
}

/* Function to compute the inward or outward special distance for directed graphs */
static void compute_sd_dir(
   const GRAPH*  p,
   int           start,
   const double* cost,
   int*          heap,
   int*          state,
   int*          count, /* pointer to store number of elements of heap */
   double*       pathdist,
   double*       pathtran,
   char          inward)
{
   int    k;
   int    m;
   int    i;
   int    curr_edge;
   double tran;
   double dist;

   assert(p          != NULL);
   assert(start      >= 0);
   assert(start      <  p->knots);
   assert(heap       != NULL);
   assert(state      != NULL);
   assert(pathdist   != NULL);
   assert(pathtran   != NULL);
   assert(cost       != NULL);
   assert(count      != NULL);
   assert(*count     >= 0);

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
         k = nearest(heap, state, count, pathdist, pathtran);

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
            if( inward )
               curr_edge = Edge_anti(i);
            else
               curr_edge = i;

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
               tran = Is_term(p->term[m]) ? 0.0 : pathtran[k] + cost[curr_edge];
               dist = Max(pathdist[k], pathtran[k] + cost[curr_edge]);

               if (LT(dist, pathdist[m])
                  || (EQ(dist, pathdist[m]) && LT(tran, pathtran[m])))
               {
                  pathdist[m] = dist;
                  pathtran[m] = tran;

                  correct(heap, state, count, pathdist, pathtran, m);
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
SCIP_RETCODE sd_reduction(
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
   int*    elimins,
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

   *elimins = 0;
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
      if( i % 100 == 0 && *elimins == 0 && SCIPgetTotalTime(scip) - redstarttime > stalltime)
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
            || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source[0] )
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
	    SCIPindexListNodeFree(scip, &((g->ancestors)[e]));
	    SCIPindexListNodeFree(scip, &((g->ancestors)[flipedge(e)]));
	    assert(g->ancestors[e] == NULL);
            graph_edge_del(g, e);
            (*elimins)++;
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

   SCIPdebugMessage("%d Edges deleted\n", *elimins * 2);
   /*printf("%d SD: Edges deleted\n", elimins * 2);*/
   return SCIP_OKAY;
}

/* C. W. Duin and A. Volganant
 *
 * "An Edge Elimination Test for the Steiner Problem in Graphs"
 *
 * Operations Research Letters 8, (1989), 79-83
 *
 * Special Distance Test for directed graphs
 */
SCIP_RETCODE sd_reduction_dir(
   SCIP*    scip,
   GRAPH*   g,
   double** sd_indist,
   double** sd_intran,
   double** sd_outdist,
   double** sd_outtran,
   double*  cost,
   int*     heap,
   int*     state,
   int*     outterms,
   int*     elimins
   )
{
   int*    sourceadj;
   int     outtermcount = 0;
   int     count = 0;
   int     i;
   int     e;
   int     j;
   int     k;
   int     l;
   double  tempsd;
   double  specialdist;

   assert(sd_indist  != NULL);
   assert(sd_intran  != NULL);
   assert(sd_outdist != NULL);
   assert(sd_outtran != NULL);
   assert(cost       != NULL);
   assert(heap       != NULL);
   assert(state      != NULL);
   assert(elimins     != NULL);


   SCIPdebugMessage("SD-Reduktion: ");
   fflush(stdout);

   assert(outterms != NULL);

   sourceadj = malloc((size_t)g->knots * sizeof(int));

   for(i = 0; i < g->knots; i++)
   {
      assert(sd_indist[i]  != NULL);
      assert(sd_intran[i]  != NULL);
      assert(sd_outdist[i] != NULL);
      assert(sd_outtran[i] != NULL);

      if( Is_term(g->term[i]) )
      {
         for(e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
         {
            if( Is_term(g->term[i]) && LT(g->cost[e], FARAWAY) )
            {
               //printf("Outgoing edge for terminal %d: %d\n", i, e);
               outterms[outtermcount] = i;
               outtermcount++;

               break;
            }
         }
      }
      g->mark[i] = (g->grad[i] > 0);
      sourceadj[i] = -1;
   }

   /* getting the knots that are adjacent to the source */
   for( e = g->outbeg[g->source[0]]; e != EAT_LAST; e = g->oeat[e] )
   {
      l = g->head[e];
      sourceadj[l] = l;
   }

   for(i = 0; i < g->edges; i++)
      cost[i] = g->cost[i] * 1000.0 + (double)(rand() % 512);

   for( i = 0; i < outtermcount; i++ )
   {
      compute_sd_dir(g, outterms[i], cost, heap, state, &count, sd_indist[i], sd_intran[i], TRUE);
      compute_sd_dir(g, outterms[i], cost, heap, state, &count, sd_outdist[i], sd_outtran[i], FALSE);
   }

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

      /* for the hop constrained problems we only want to examine the nodes adjacent to the source. */
      if( g->stp_type == STP_HOP_CONS && sourceadj[i] < 0 )
         continue;


      for(e = g->outbeg[i]; e != EAT_LAST; e = j)
      {
         assert(g->tail[e] == i);
         l = g->head[e];

         j = g->oeat[e];

         if ( (g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) && l == g->source[0] )
            continue;

         /* for the hop constrained problems we only want to examine the nodes adjacent to the source. */
         if( g->stp_type == STP_HOP_CONS && sourceadj[l] < 0 )
            continue;

         specialdist = FARAWAY;
         for( k = 0; k < outtermcount; k++ )
         {
            assert(l >= 0 && l < g->knots);
            tempsd = FARAWAY;
            if( outterms[k] == g->source[0] )
            {
               if( (LT(sd_indist[k][i], FARAWAY) || LT(sd_outdist[k][i], FARAWAY)) && LT(sd_outdist[k][l], FARAWAY) )
               {
                  if( !LT(sd_indist[k][i], FARAWAY) )
                     tempsd = sd_outdist[k][i];
                  else if( !LT(sd_outdist[k][i], FARAWAY) )
                     tempsd = sd_indist[k][i];
                  else if( GT(sd_indist[k][i], sd_outdist[k][i]) )
                     tempsd = sd_indist[k][i];
                  else
                     tempsd = sd_outdist[k][i];


                  if( GT(sd_outdist[k][l], tempsd) )
                     tempsd = sd_outdist[k][l];
               }
            }
            else
            {
               if( LT(sd_indist[k][i], FARAWAY) && LT(sd_outdist[k][l], FARAWAY) )
               {
                  tempsd = sd_indist[k][i];

                  if( GT(sd_outdist[k][l], tempsd) )
                     tempsd = sd_outdist[k][l];
               }
            }

            if( LT(tempsd, specialdist) )
               specialdist = tempsd;
         }

         if (LT(cost[e], FARAWAY) && LT(specialdist, cost[e]))
         {
            if( LT(g->cost[Edge_anti(e)], FARAWAY) )
            {
               g->cost[e] = FARAWAY;
               cost[e] = FARAWAY;
            }
            else
               graph_edge_del(g, e);

            (*elimins)++;

         }
      }
   }

   assert(graph_valid(g));

   SCIPdebugMessage("%d Edges deleted\n", *elimins * 2);

   free(sourceadj);

   return SCIP_OKAY;
}

static int redirect_edge(
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
      else
      {
	 e = -1;
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
      return eki;
   }
   return e;
}

/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Bottleneck Degree 3 Test
 */
SCIP_RETCODE bd3_reduction(
   SCIP* scip,
   GRAPH* g,
   double* pathdist1,
   double* pathtran1,
   int* heap,
   int*  state,
   int*  nelims
   )
{
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real* pathdist2;
   SCIP_Real* pathtran2;
   SCIP_Real* pathrand;

   int     count = 0;
   int    i;
   int    k;
   int    n1;
   int    k1;
   int    k2;
   int    k3;
   int    e1;
   int    e2;
   int    e3;
   SCIP_Real c1;
   SCIP_Real c1a;
   SCIP_Real c2;
   SCIP_Real c2a;
   SCIP_Real c3;
   SCIP_Real c3a;
   SCIP_Real c123;

   SCIPdebugMessage("BD3-Reduction: ");
   fflush(stdout);

   assert(heap  != NULL);
   assert(state != NULL);
   assert(nelims != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist2, g->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathtran2, g->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathrand, g->knots) );

   /*
     ancestors = malloc((size_t)(3) * sizeof(IDX*));
     revancestors = malloc((size_t)(3) * sizeof(IDX*));
     pathdist2 = malloc((size_t)g->knots * sizeof(double));
     pathtran2 = malloc((size_t)g->knots * sizeof(double));
     pathrand  = malloc((size_t)g->knots * sizeof(double));
   */
   *nelims = 0;

   for( i = 0; i < 3; i++ )
   {
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }
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
      c1a = g->cost[Edge_anti(e1)];
      e2 = g->oeat[e1];

      assert(e2 != EAT_LAST);

      k2 = g->head[e2];
      c2 = g->cost[e2];
      c2a = g->cost[Edge_anti(e2)];
      e3 = g->oeat[e2];

      assert(e3 != EAT_LAST);

      k3 = g->head[e3];
      c3 = g->cost[e3];
      c3a = g->cost[Edge_anti(e3)];

      assert(g->oeat[e3] == EAT_LAST);

      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ( (g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) &&
         (g->head[e1] == g->source[0] || g->head[e2] == g->source[0] || g->head[e3] == g->source[0]))
         continue;

      if( !(EQ(c1, c1a) && EQ(c2, c2a) && EQ(c3, c3a)) )
         continue;

      if( !(LT(c1, FARAWAY) && LT(c2, FARAWAY) && LT(c3, FARAWAY) &&
            LT(c1a, FARAWAY) && LT(c2a, FARAWAY) && LT(c3a, FARAWAY)) )
         continue;

      for(k = 0; k < 3; k++)
      {
	 SCIPindexListNodeFree(scip, &(ancestors[k]));
	 SCIPindexListNodeFree(scip, &(revancestors[k]));
         assert(ancestors[k] == NULL);
         assert(revancestors[k] == NULL);
      }

      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(ancestors[0]), g->ancestors[e1]) );
      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(revancestors[0]), g->ancestors[Edge_anti(e1)]) );

      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(ancestors[1]), g->ancestors[e2]) );
      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(revancestors[1]), g->ancestors[Edge_anti(e2)]) );

      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(ancestors[2]), g->ancestors[e3]) );
      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(revancestors[2]), g->ancestors[Edge_anti(e3)]) );

      compute_sd(g, k1, g->cost, g->cost, heap, state, &count, pathdist1, pathtran1, pathrand);
      compute_sd(g, k2, g->cost, g->cost, heap, state, &count, pathdist2, pathtran2, pathrand);

      c123 = c1 + c2 + c3;

      if (GT(pathdist1[k2] + pathdist1[k3], c123)
         && GT(pathdist1[k2] + pathdist2[k3], c123)
         && GT(pathdist1[k3] + pathdist2[k3], c123))
         continue;

      (*nelims)++;

      if (LT(pathdist1[k2], c1 + c2))
      {
	 SCIPindexListNodeFree(scip, &((g->ancestors)[e1]));
	 SCIPindexListNodeFree(scip, &((g->ancestors)[Edge_anti(e1)]));
         graph_edge_del(g, e1);
      }
      else
      {
	 n1 = redirect_edge(g, e1, k1, k2, c1 + c2);
         if( n1 >= 0 )
	 {
            SCIPindexListNodeFree(scip, &(g->ancestors[n1]));
            SCIPindexListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors[0]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors[1]) );

            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors[0]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors[1]) );
	 }
      }

      if (LT(pathdist2[k3], c2 + c3))
      {
	 SCIPindexListNodeFree(scip, &((g->ancestors)[e2]));
	 SCIPindexListNodeFree(scip, &((g->ancestors)[Edge_anti(e2)]));
         graph_edge_del(g, e2);
      }
      else
      {
	 n1 = redirect_edge(g, e2, k2, k3, c2 + c3);
	 if( n1 >= 0 )
	 {
            SCIPindexListNodeFree(scip, &(g->ancestors[n1]));
            SCIPindexListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors[1]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors[2]) );

            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors[1]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors[2]) );
	 }
      }

      if (LT(pathdist1[k3], c1 + c3))
      {
	 SCIPindexListNodeFree(scip, &((g->ancestors)[e3]));
	 SCIPindexListNodeFree(scip, &((g->ancestors)[Edge_anti(e3)]));
         graph_edge_del(g, e3);
      }
      else
      {
         n1 = redirect_edge(g, e3, k3, k1, c3 + c1);
         if( n1 >= 0 )
	 {
            SCIPindexListNodeFree(scip, &(g->ancestors[n1]));
            SCIPindexListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors[2]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors[0]) );

            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors[2]) );
            SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors[0]) );
	 }
      }

      assert(g->grad[i] == 0);
   }

   for(k = 0; k < 3; k++)
   {
      SCIPindexListNodeFree(scip, &(ancestors[k]));
      SCIPindexListNodeFree(scip, &(revancestors[k]));
      assert(ancestors[k] == NULL);
      assert(revancestors[k] == NULL);
   }

   SCIPfreeBufferArray(scip, &pathrand);
   SCIPfreeBufferArray(scip, &pathtran2);
   SCIPfreeBufferArray(scip, &pathdist2);
   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);

   assert(graph_valid(g));
   /*
     printf(" Knots deleted %d\n", nelims);
   */
   SCIPdebugMessage("bd3: %d Knots deleted\n", *nelims);

   return SCIP_OKAY;
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
SCIP_RETCODE nsv_reduction(
   SCIP*   scip,
   GRAPH*  g,
   double* cost,
   double* fixed,
   int* nelims
   )
{
   SCIP_Real redstarttime;
   SCIP_Real timelimit;
   SCIP_Real stalltime;
   PATH**  path;
   PATH*   mst1;
   PATH*   mst2;
   int     i;
   int     e;
   int     k;
   int     j;
   double  min1;
   double  min2;
   double  cost1;
   double  cost2;

   SCIPdebugMessage("NSV-Reduction: ");
   fflush(stdout);
   /*
     graph_show(g);
   */
   *nelims = 0;

   path = malloc((size_t)g->knots * sizeof(PATH*));

   assert(path != NULL);

   mst1 = malloc((size_t)g->knots * sizeof(PATH));
   mst2 = malloc((size_t)g->knots * sizeof(PATH));

   assert(mst1 != NULL);
   assert(mst2 != NULL);
   assert(cost != NULL);

   redstarttime = SCIPgetTotalTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   stalltime = timelimit*0.1; /* this should be set as a parameter */

   /* Check this cost setting. It may need to be changed for the directed case.
    */
   for(i = 0; i < g->edges; i++)
      cost[i] = g->cost[i];

   for(i = 0; i < g->knots; i++)
   {
      g->mark[i] = (g->grad[i] > 0);
      path[i] = NULL;
   }

   calculate_distances(g, path, g->cost, FSP_MODE);

   for(i = 0; i < g->knots; i++)
   {
      if( i % 100 == 0 && SCIPgetTotalTime(scip) > timelimit )
         break;

      if( i % 100 == 0 && (*nelims) == 0 && SCIPgetTotalTime(scip) - redstarttime > stalltime)
         break;

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
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]) );
         SCIP_CALL( graph_knot_contract(scip, g, i, k) );

         (*nelims)++;

         calculate_distances(g, path, g->cost, FSP_MODE);

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

   assert(graph_valid(g));

   SCIPdebugMessage(" %d Knots deleted\n", *nelims);
   /*printf("nsv_reduction: %d Knots deleted\n", elimins);*/

   return SCIP_OKAY;
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
SCIP_RETCODE nv_reduction_optimal(
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
   char    antiedgeexists;
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

   assert(g->source[0] >= 0);

   /* computing the voronoi regions inward to a node */
   voronoi_term(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, 1);

   /* computing the shortest paths from the source node */
   graph_path_exec(g, FSP_MODE, g->source[0], g->cost, pathfromsource);

   /* computing the shortest hops paths from the source node */
   graph_path_exec(g, FSP_MODE, g->source[0], hopscost, pathhops);

   /* this is the offset used to minimise the number of knots to examine in large graphs. */
   srand(runnum*100);
   knotoffset = rand() % KNOTFREQ;

   for(i = 0; i < g->knots; i++)
   {
      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ((g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source[0] )
         continue;

      if( g->stp_type == STP_DIRECTED && i != g->source[0] )
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
            if( g->stp_type == STP_DIRECTED && i == g->source[0] )
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

         if( g->stp_type == STP_DIRECTED )
            continue;

         if (LT(min1, FARAWAY) && LE(pathfromsource[shortarctail].dist + min1, min2))
         {
            assert(shortarc >= 0);
            assert(shortarctail >= 0);

            if ((g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT
                  || g->stp_type == STP_DIRECTED) && shortarctail == g->source[0] )
               continue;

            if( g->stp_type == STP_HOP_CONS && GT(pathfromsource[shortarctail].hops, pathhops[i].dist - 1) )
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
                           graph_edge_del(g, e);

                     }
                  }
                  (*elimins)++;
               }
            }
            else
            {
               *fixed += min1;
               SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[shortarc]) ); /* I think that this should be
                                                                                                           shortarc instead of shortarctail */
               SCIP_CALL( graph_knot_contract(scip, g, i, shortarctail) );

               if( g->stp_type == STP_HOP_CONS )
               {
                  for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
                     g->cost[e] = FARAWAY;
               }

               (*elimins)++;
            }

            /* computing the shortest paths from the source node */
            graph_path_exec(g, FSP_MODE, g->source[0], g->cost, pathfromsource);
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

         if ((g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_MAX_NODE_WEIGHT) && (i == g->source[0] || j == g->source[0]) )
            continue;


         if( Is_term(g->term[i]) )
         {
            SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e]);
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



SCIP_RETCODE sl_reduction(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   double* fixed,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
   )
{
   double* mincost2;
   int*    minedge1;
   int*    minedgehead;
   int*    minedgetail;
   int     i;
   int     k;
   int     e;
   int     j;
   int     min1;
   int     head;
   int     tail;
   int     nnodes;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   *nelims = 0;
   nnodes = g->knots;

   /* TODO ID */
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincost2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedgehead, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedgetail, nnodes) );

   assert(vbase != NULL);
   assert(vnoi != NULL);
   assert(minedge1 != NULL);
   assert(mincost2 != NULL);
   assert(minedgehead != NULL);
   assert(minedgetail != NULL);

   for( k = 0; k < nnodes; k++ )
   {
      minedge1[k] = UNKNOWN;
      mincost2[k] = FARAWAY;
   }

   voronoi_pres(g, g->cost, vnoi, vbase, heap, state);

   for( i = 0; i < nnodes; i++ )
   {
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( vbase[g->tail[e]] != vbase[g->head[e]] )
         {
            k = vbase[g->tail[e]];
            assert(k != UNKNOWN);

            min1 = minedge1[k];

            if( min1 == UNKNOWN )
            {
               minedge1[k] = e;
            }
            else if( SCIPisLT(scip, g->cost[e], g->cost[min1]) )
	    {
	       mincost2[k] = g->cost[min1];
	       minedge1[k] = e;
	    }
            else if( SCIPisLT(scip, g->cost[e], mincost2[k]) )
            {
               mincost2[k] = g->cost[e];
            }
         }
      }
   }

   for( k = 0; k < nnodes; k++ )
   {
      if( !SCIPisEQ(scip, FARAWAY, mincost2[k]) )
      {
	 e = minedge1[k];
	 assert(e >= 0);
	 assert(LT(mincost2[k], FARAWAY));
	 minedgetail[k] = g->tail[e];
	 minedgehead[k] = g->head[e];
      }
      else
      {
	 minedgehead[k] = UNKNOWN;
         minedgetail[k] = UNKNOWN;
      }
   }

   for( i = 0; i < nnodes; i++ )
   {
      if( !SCIPisEQ(scip, FARAWAY, mincost2[i]) )
      {
	 e = minedge1[i];
	 assert(e >= 0);
	 assert(LT(mincost2[i], FARAWAY));

	 tail = g->tail[e];
	 head = g->head[e];
	 if( tail == minedgetail[i] && head == minedgehead[i] )
	 {
            if( SCIPisGT(scip, mincost2[i], vnoi[tail].dist + g->cost[e] + vnoi[head].dist) )
            {
               for( j = g->outbeg[tail]; j != EAT_LAST; j = g->oeat[j] )
                  if( j == e )
                     break;
               if( j == EAT_LAST )
                  continue;

               for( j = g->inpbeg[head]; j != EAT_LAST; j = g->ieat[j] )
                  if( j == e )
                     break;
               if( j == EAT_LAST )
                  continue;
               (*nelims)++;
               *fixed += g->cost[e];
               //printf("contrSL: %d-%d \n", tail, head);
               SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]) );
	       SCIP_CALL( graph_knot_contract(scip, g, tail, head) );
            }
	 }
      }
   }

   SCIPfreeBufferArray(scip, &minedgetail);
   SCIPfreeBufferArray(scip, &minedgehead);
   SCIPfreeBufferArray(scip, &mincost2);
   SCIPfreeBufferArray(scip, &minedge1);
   //printf("sl: nelims: %d \n", nelims);
   return SCIP_OKAY;
}


/* NV reduction from T. Polzin's "Algorithms for the Steiner problem in networks" */
int nv_reduction(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   double* fixed,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
   )
{
   double* distance;
   double* mincost2;
   double  min1;
   double  min2;
   int*    minedge1;
   int*    min1head;
   int*    min1tail;
   int     edge1;
   int     nnodes;
   int     nterms;
   int     termcount;
   int     i;
   int     j;
   int     k;
   int     e;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   termcount = 0;
   *nelims = 0;
   nnodes = g->knots;
   nterms = g->terms;
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincost2, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &min1head, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &min1tail, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distance, nnodes) );
   /*
     minedge1 = malloc((size_t) * sizeof(int));
     mincost2 = malloc((size_t)g->terms * sizeof(double));
     min1head = malloc((size_t)g->terms * sizeof(int));
     min1tail = malloc((size_t)g->terms * sizeof(int));
     distance = malloc((size_t)nnodes * sizeof(double));
   */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         edge1 = UNKNOWN;
	 min2  = FARAWAY;
	 if( g->grad[i] >= 2 )
         {
            min1  = FARAWAY;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
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
         }
         minedge1[termcount] = edge1;

	 mincost2[termcount] = min2;
	 if( SCIPisLT(scip, min2, FARAWAY) )
	 {
	    assert(edge1 != UNKNOWN);
	    min1head[termcount] = g->head[edge1];
	    min1tail[termcount] = g->tail[edge1];
	 }
	 else
	 {
	    min1head[termcount] = UNKNOWN;
	    min1tail[termcount] = UNKNOWN;
	 }
	 termcount++;
      }
   }
   assert(termcount == g->terms);
   voronoi_dist(g, g->cost, distance, vbase, minedge1, heap, state, vnoi);
   for( i = 0; i < termcount; i++ )
   {
      min2 = mincost2[i];
      if( SCIPisEQ(scip, min2, FARAWAY) )
	 continue;
      edge1 = minedge1[i];
      assert(LT(min2, FARAWAY));
      assert(edge1 != UNKNOWN);
      j = min1tail[i];
      k = min1head[i];
      if( j == g->tail[edge1] && k == g->head[edge1] )
      {
         if( vbase[k] != j )
         {
            if( SCIPisGE(scip, min2, g->cost[edge1] + vnoi[k].dist) )
            {
               for( e = g->outbeg[j]; e != EAT_LAST; e = g->oeat[e] )
                  if( e == edge1 )
                     break;
               if( e == EAT_LAST )
                  continue;

               for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
                  if( e == edge1 )
                     break;
               if( e == EAT_LAST )
                  continue;

               (*nelims)++;
               *fixed += g->cost[edge1];
               SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1]) );
               SCIP_CALL( graph_knot_contract(scip, g, j, k) );
            }
         }
         else
         {
            if( SCIPisGE(scip, min2, distance[j]) )
            {
               for( e = g->outbeg[j]; e != EAT_LAST; e = g->oeat[e] )
                  if( e == edge1 )
                     break;
               if( e == EAT_LAST )
                  continue;

               for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
                  if( e == edge1 )
                     break;
               if( e == EAT_LAST )
                  continue;
               *fixed += g->cost[edge1];
               SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1]) );
               SCIP_CALL( graph_knot_contract(scip, g, j, k) );
               (*nelims)++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &distance);
   SCIPfreeBufferArray(scip, &min1tail);
   SCIPfreeBufferArray(scip, &min1head);
   SCIPfreeBufferArray(scip, &mincost2);
   SCIPfreeBufferArray(scip, &minedge1);

   assert(graph_valid(g));
   return SCIP_OKAY;
}

/*  longest edge reduction test from T. Polzin's "Algorithms for the Steiner problem in networks" (Lemma 20) */
SCIP_RETCODE ledge_reduction(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
   )
{
   GRAPH* netgraph;
  // GRAPH* dg;
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
   //int ledge;
   int* nodesid;
   ////int* edgerev;
  // char* blocked;
  // SCIP_CALL( SCIPallocBufferArray(scip, &blocked, g->edges) );
  // SCIP_CALL( SCIPallocBufferArray(scip, &edgerev, g->edges) );
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
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
         nterms++;

   if( nterms <= 1 )
      return SCIP_OKAY;

   voronoi_pres(g, g->cost, vnoi, vbase, heap, state);

   if( SCIPisGT(scip, nedges, (nterms - 1) * nterms) )
      maxnedges = (nterms - 1) * nterms;
   else
      maxnedges = nedges;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize the new graph */
   netgraph = graph_init(nterms, maxnedges, 1, 0);

   //dg = graph_init(nnodes, 2 * nedges, 1, 0);

   e = 0;
   for( k = 0; k < nnodes; k++ )
   {

      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
	 if( e == 0 )
            graph_knot_add(netgraph, 0);
         else
            graph_knot_add(netgraph, -1);
	 netgraph->mark[e] = TRUE;
	 nodesid[k] = e++;
      }
      else
      {
	 nodesid[k] = UNKNOWN;
      }
    //  graph_knot_add(dg, 0);
   }
   //printf("nterms: %d\n\n", e);

  // dg->source[0] = g->source[0];
  // assert(Is_term(dg->term[g->source[0]]));

   netnnodes = netgraph->knots;
   assert(netnnodes == e);
   assert(netnnodes == nterms);
   //for( e = 0; e < nedges; e++ )
  // {
    //  blocked[e] = FALSE;
     // edgerev[e] = UNKNOWN;
   //}
   for( k = 0; k < nnodes; k++ )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         /* TODO */
	 assert(g->grad[k] > 0 );
	 v1 = vbase[k];
         assert(k == g->tail[e]);

	 if( v1 != vbase[g->head[e]] )
	 {
            v2 = vbase[g->head[e]];
	    assert(Is_term(g->term[v1]));
	    assert(Is_term(g->term[v2]));
	    assert(nodesid[v1] >= 0);
	    assert(nodesid[v2] >= 0);
	    ne = EAT_LAST;
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
		//  edgerev[ne] = e;
		 // edgerev[flipedge(ne)] = flipedge(e);
		  assert(ne <= maxnedges);
               }
	    }
	    else
	    {
	    //   edgerev[netgraph->edges] = e;
             //  edgerev[flipedge(netgraph->edges)] = flipedge(e);
	       graph_edge_add(netgraph, nodesid[v1], nodesid[v2], cost, cost);
	       assert(netgraph->edges <= maxnedges);
	    }
	 }
      }
   }
   netgraph->source[0] = 0;

   assert( graph_valid(netgraph ) );

   /* compute a MST on netgraph */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, netnnodes) );
   graph_path_init(netgraph);
   graph_path_exec(netgraph, MST_MODE, 0, netgraph->cost, mst);

   maxcost = -1;
   assert(mst[0].edge == -1);

   for( k = 1; k < netnnodes; k++ )
   {
      assert(netgraph->path_state[k] == CONNECT);
      e = mst[k].edge;
      assert(e >= 0);
      cost = netgraph->cost[e];
      if( SCIPisGT(scip, cost, maxcost) )
      {
         maxcost = cost;
      }
   }
#if 0
   ledge = -1;
   for( k = 1; k < netnnodes; k++ )
   {

      ledge = edgerev[mst[k].edge];
      assert(ledge >= 0);
      v1 = g->tail[ledge];
      v2 = g->head[ledge];
      // printf("e: %d->%d \n", v1,v2);
      blocked[ledge] = TRUE;

      graph_edge_add(dg, v1, v2, 1, 1);
      int k1 = v1;
      int k2 = v2;

      SCIP_Real obj = g->cost[ledge];

      // printf("k1: %d, %d\n",k1,k2);
      assert(Is_term(g->term[vbase[k1]]));
      assert(Is_term(g->term[vbase[k2]]));
      assert(vbase[k1] != vbase[k2]);
      while( k1 != vbase[k1] )
      {
         //  printf("et: %d->%d \n", k1, g->tail[vnoi[k1].edge] );
         blocked[vnoi[k1].edge] = TRUE;
         obj += g->cost[vnoi[k1].edge];
         graph_edge_add(dg, k1, g->tail[vnoi[k1].edge], 1, 1);
         k1 = g->tail[vnoi[k1].edge];
      }
      assert(k1 == vbase[k1]);

      while( k2 != vbase[k2] )
      {
         // printf("eh: %d->%d \n", k2, g->tail[vnoi[k2].edge] );
         blocked[vnoi[k2].edge] = TRUE;
         obj += g->cost[vnoi[k2].edge];
         graph_edge_add(dg, k2, g->tail[vnoi[k2].edge], 1, 1);
         k2 = g->tail[vnoi[k2].edge];
      }

      assert(k1 != k2);
      nodesid[k1] = UNKNOWN;
      nodesid[k2] = UNKNOWN;
      assert(SCIPisEQ(scip, obj, netgraph->cost[mst[k].edge]));
   }

   assert(graph_valid(dg));
#endif
   // printf("maxcost: %f \n", maxcost);
   for( k = 0; k < nnodes; k++ )
   {
      //assert(nodesid[k] == UNKNOWN);
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
         assert(e >= 0);
         if( SCIPisGT(scip, g->cost[e], maxcost) )
         {
            (*nelims)++;
            SCIPindexListNodeFree(scip, &((g->ancestors)[e]));
            SCIPindexListNodeFree(scip, &((g->ancestors)[Edge_anti(e)]));
	    v1 = g->oeat[e];
            // printf("LE: elim: %d->%d (%d) %f \n", g->tail[e], g->head[e], e, g->cost[e]);
	 //   if( blocked[e] || blocked[flipedge(e)] )
          //     assert(0);
            graph_edge_del(g, e);
            e = v1;
            //	    assert(graph_valid(g));
         }
         else
	 {
            e = g->oeat[e];
	 }
      }
   }

   if( *nelims > 0 )
   {
      int count = 0;
      for(k = 0; k < g->knots; k++)
         g->mark[k] = FALSE;

      graph_trail(g, g->source[0]);

      for( k = 0; k < g->knots; k++ )
      {
         if( !g->mark[k] && (g->grad[k] > 0) )
	 {
	    assert(!Is_term(g->term[k]));
	    e = g->inpbeg[k];
	    while( e != EAT_LAST )
            {
               graph_edge_del(g, e);
               SCIPindexListNodeFree(scip, &(g->ancestors[e]));
               SCIPindexListNodeFree(scip, &(g->ancestors[Edge_anti(e)]));
               (count)++;
               e = g->inpbeg[k];
            }
	 }
      }
      //printf("ndeletes: %d \n", count);
   }
   /* free netgraph and  MST data structure */
   graph_path_exit(netgraph);
   graph_free(scip, netgraph, TRUE);
   //graph_free(scip, dg, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &nodesid);
  // SCIPfreeBufferArray(scip, &blocked);
   //SCIPfreeBufferArray(scip, &edgerev);
   //printf("LE elims: %d \n", *nelims);
   assert(graph_valid(g));
   return SCIP_OKAY;
}
#if 0

/* T. Polzin
 *
 * "Algorithms for the Steiner problem in networks"
 *
 * Section 3.3.3 pp. 54-55
 *
 * This is undirected nearest vertex test
 */
int nv_reduction(
   GRAPH*  g,
   double* fixed)
{
   PATH**  path;
   PATH*   pathfromterm;
   PATH*   pathfromsource;
   double* distance;
   double* radius;
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
   double  mindist;
   double  minshortarcdist;
   int     shortarc;
   int     shortarctail;
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

   distance = malloc((size_t)g->knots * sizeof(double));
   radius = malloc((size_t)g->knots * sizeof(double));

   assert(distance != NULL);
   assert(radius != NULL);

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
   voronoi_term(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, 1);

   /* computing the shortest paths from each terminal to every other node */
   calculate_distances(g, path, g->cost, FSP_MODE);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]) && g->grad[i] >= 2)
      {
         min1  = FARAWAY;
         min2  = FARAWAY;
         shortarctail = -1;
         shortarc = -1;
         for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
         {
            assert(g->cost[e] == g->cost[Edge_anti(e)]);
            if ( LE(g->cost[e], min1) )
            {
               shortarc = e;
               shortarctail = g->tail[e];

               min2 = min1;
               min1 = g->cost[e];
            }
            else if( LE(g->cost[e], min2) )
               min2 = g->cost[e];

            if( Is_term(g->term[g->tail[e]]) )
            {
               min1 = FARAWAY;
               break;
            }
         }

         // NOTE: if min1 and min2 are equal, the distance array may not be updated correctly.
         // What can occur is distance is calculated using min2, the contraction occurs with min1. Hence, the wrong
         // reduction occurs.
         if( LT(min1, FARAWAY) && !EQ(min1, min2) )
         {
            //printf("i: %d, min1: %f, min2: %f, distance: %f, pathfromterm: %f, path: %f\n", i, min1, min2, distance[i],
            //      pathfromterm[shortarctail].dist, path[vregion[shortarctail]][shortarctail].dist);
            mindist = FARAWAY;
            minshortarcdist = FARAWAY;
            k = -1;
            if( vregion[shortarctail] == vregion[i] )
            {
               for( j = 0; j < g->knots; j++ )
               {
                  if( i != j && Is_term(g->term[j]) && g->grad[j] > 0 )
                  {
                     if( LT(path[j][shortarctail].dist, minshortarcdist) )
                        minshortarcdist = path[j][shortarctail].dist;

                     if( LT(path[j][i].dist, mindist) )
                     {
                        //printf("Terminal: %d\n", j);
                        mindist = path[j][i].dist;
                        k = j;
                     }
                  }
               }
               //assert(LE(mindist + min1, distance[i]));
               assert(LE(mindist, min1 + minshortarcdist));
               assert(LE(mindist, distance[i]));
               if( path[k][i].edge != shortarc && path[k][i].edge != Edge_anti(shortarc) )
                  continue;
               //printf("Minimum Distance: %f %f %f\n", mindist, min1 + minshortarcdist, min1 + distance[i]);
               assert(EQ(mindist, min1 + minshortarcdist));
               assert(LE(min1 + minshortarcdist, min1 + distance[i]));
            }

            //continue;

            // distance is only given as a upper bound on the length of the path from i to the nearest terminal using
            // shortarc.
            if( (vregion[shortarctail] == vregion[i] && GE(min2, min1 + distance[i])) ||
               (vregion[shortarctail] != vregion[i] && GE(min2, min1 + pathfromterm[shortarctail].dist)) )
            {

               assert(vregion[i] == i);
               assert(vregion[shortarctail] != vregion[i] || distance[i] < FARAWAY);
               assert(min2 < FARAWAY);

               *fixed += min1;
               SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[shortarc]); /* I think that this should be
                                                                                                         shortarc instead of shortarctail */

                  graph_knot_contract(scip, g, shortarctail, i);
                  //graph_knot_contract(scip, g, i, shortarctail);

                  elimins++;

                  //printf("i: %d, shortarctail: %d, isterm: %d, radius[i]: %f, radius[shortarctail] %f\n", i,
                  //   shortarctail, g->term[i], radius[vregion[i]], radius[vregion[shortarctail]]);

                  voronoi_term(g, g->cost, distance, radius, pathfromterm, vregion, heap, state, pred, 1);
                  calculate_distances(g, path, g->cost, FSP_MODE);
                  }
            }
         }
#if 0
         /* The knot is not a terminal so we can perform the short link test */
         else if ( !Is_term(g->term[i]) )
         {
            for(e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
            {
               j = g->tail[e];
               if( vregion[i] != vregion[j] )
               {
                  if( minArc1[vregion[i]] < 0 )
                     minArc1[vregion[i]] = e;
                  else if( LE(g->cost[e], g->cost[minArc1[vregion[i]]]) )
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

         if( minArc1[terms[k]] >= 0 && minArc2[terms[k]] >= 0 && GE(g->cost[minArc2[terms[k]]],
               pathfromterm[g->tail[minArc1[terms[k]]]].dist + g->cost[minArc1[terms[k]]]
               + pathfromterm[g->head[minArc1[terms[k]]]].dist) )
         {
            e = minArc1[terms[k]];
            i = g->head[e];
            j = g->tail[e];

            if( !Is_term(g->term[i]) && !Is_term(g->term[j]) )
            {
               SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]);
                  *fixed += g->cost[e];
                  graph_knot_contract(scip, g, j, i);

                  elimins++;
                  }
            }
         }
#endif

         for( i = g->knots - 1; i >= 0; i-- )
         {
            if( path[i] != NULL )
               free(path[i]);
         }

         free(terms);
         free(minArc2);
         free(minArc1);
         free(pred);
         free(state);
         free(heap);
         free(vregion);
         free(radius);
         free(distance);
         free(pathfromsource);
         free(pathfromterm);
         free(path);

         assert(graph_valid(g));
         printf("nv_reduction: %d Knots deleted\n", elimins);
         SCIPdebugMessage("nv_reduction: %d Knots deleted\n", elimins);
         /*printf("nv_reduction: %d Knots deleted\n", elimins);*/

         return(elimins);
      }
#endif
