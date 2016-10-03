/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_alt.c
 * @brief  Altenative based reduction tests for Steiner problems
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 *
 * This file implements alternative-based reduction techniques for several Steiner problems.
 * All tests can be found in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
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
#include "probdata_stp.h"
#include "scip/scip.h"

#define CNSNN     25
#define KNOTFREQ  100
#define KNOTLIMIT 1e+20


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
#if 0
   int e;
#endif

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

   return FALSE;

#if 0 /* @todo */
   /* edge between non-terminals */
   if( !Is_term(g->term[tail]) && !Is_term(g->term[head]) )
   {
      printf("deletable! \n");
      return TRUE;
   }

   /* has edge been marked as forbidden? */
   if( forbidden[edge] )
      return FALSE;

   /* edge between a terminal and a non terminal */
   if( !Is_term(g->term[tail]) || !Is_term(g->term[head]) )
   {
#if 0
      int k;
      int base;
      int shift;
      int antiedge;

      antiedge = flipedge(edge);
#endif

      /* check whether edge is used in shortest path */

      if( !Is_term(g->term[tail]) && Is_term(g->term[head]) )
      {
         e = path[tail + tpos * nnodes].edge;

         assert(g->ieat[e] != EAT_FREE);
         assert(g->head[e] == tail);

         if( g->tail[e] == head )
            return FALSE;
#if 0
         k = tail;
         base = vbase[tail + tpos * nnodes];
         shift = tpos * nnodes;
         while( k != base )
         {
            e = path[k + shift].edge;

            assert(g->ieat[e] != EAT_FREE);

            if( e == edge || e == antiedge || g->ieat[e] == EAT_FREE )
               return FALSE;

            k = g->tail[e];
         }
#endif
      }
      else if( Is_term(g->term[tail]) && !Is_term(g->term[head]) )
      {
         e = path[head + hpos * nnodes].edge;

         assert(g->ieat[e] != EAT_FREE);
         assert(g->head[e] == head);

         if( g->tail[e] == tail )
            return FALSE;
      }
   }

   /* update forbidden edges */

   if( Is_term(g->term[head]) )
   {

      SCIP_Real ecost = g->cost[edge];
      for( e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(e >= 0);

         if( SCIPisEQ(scip, g->cost[e], ecost) )
         {
            if( forbidden[e] == FALSE && Is_term(g->term[g->head[e]]) )
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
            if( forbidden[e] == FALSE  && Is_term(g->term[g->head[e]]) )
            {
               forbidden[e] = TRUE;
               forbidden[flipedge(e)] = TRUE;
            }
         }
      }
   }
   printf("deletable! \n");
   return TRUE;
#endif
}


static
int getcloseterms(
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

/** Special distance test */
SCIP_RETCODE sd_red(
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
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real* termdist1;
   SCIP_Real* termdist2;
   SCIP_Real ecost;
   SCIP_Real dist;

   int* neighbterms1;
   int* neighbterms2;
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

   SCIP_CALL( SCIPallocBufferArray(scip, &termdist1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termdist2, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms2, 4) );

   /* compute nearest four terminals to all non-terminals */
   getnext4terms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);
#if 0
   int tpos;
   int base;
   int shift;
   int x = 0;
   for( tpos = 0; tpos < 4; tpos++ )
   {
      for( k = 0; k < nnodes; k++ )
      {
         if(  Is_term(g->term[k]) || !g->mark[k] )
            continue;


         base = vbase[k + tpos * nnodes];
         if( base == -1 )
            continue;
         shift = tpos * nnodes;
         x = 0;
         printf("in %d\n", base);
         i = k;
         while( i != base )
         {
            if( x++ >= 10000 )
            {
               printf("FAIL k: %d, tpos: %d, base: %d\n",k, tpos, base);
               assert(0);
            }
            e = vnoi[i + shift].edge;

            assert(g->ieat[e] != EAT_FREE);

            i = g->tail[e];
         }
      }
   }
#endif
#if 0
   if( 1 )
   {
      SCIP_Real** pathdist = NULL;
      int** pathedge = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
      BMSclearMemoryArray(pathdist, nnodes);
      BMSclearMemoryArray(pathedge, nnodes);

      for( k = 0; k < nnodes; k++ )
      {
         g->mark[k] = (g->grad[k] > 0);

         if( 1 || !Is_term(g->term[k]) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathdist[k]), nnodes) ); /*lint !e866*/
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathedge[k]), nnodes) ); /*lint !e866*/
         }
      }
      SCIP_CALL( SCIPprobdataPrintGraph2(g, "exgraph.gml", NULL) );
      for( k = 0; k < nnodes; k++ )
      {
         if( (1 || !Is_term(g->term[k])) && g->mark[k] )
         {
            assert(pathdist[k] != NULL);
            assert(pathedge[k] != NULL);
            graph_path_execX(scip, g, k, g->cost,  pathdist[k], pathedge[k]);
            if( !Is_term(g->term[k]) )
               continue;
            j = vbase[k];
            i = j;

            printf("x %f, %f ", pathdist[k][i], vnoi[k].dist);
            printf("node terminal %d, %d \n\n", k, i);
            assert(Is_term(g->term[i]));
            if( !SCIPisEQ(scip, pathdist[k][i], vnoi[k].dist ) )
            {
               printf("FAILS0 %f, %f ", pathdist[k][i], vnoi[k].dist);
               printf("node terminal %d, %d \n\n", k, i);

               assert(0);
            }

            i = vbase[k + nnodes];
            printf("x2 %f, %f ", pathdist[k][i], vnoi[k + nnodes].dist);
            printf("node terminal %d, %d \n\n", k, i);
            assert(Is_term(g->term[i]));
            if( i != -1 && !SCIPisEQ(scip, pathdist[k][i], vnoi[k + nnodes].dist ) )
            {
               printf("FAILS1 %f, %f ", pathdist[k][i], vnoi[k + nnodes].dist);
               printf("node terminal %d, %d \n\n", k, i);
               assert(0);
            }
            i = vbase[k + 2 * nnodes];
            printf("x3 %f, %f ", pathdist[k][i], vnoi[k + 2 *  nnodes].dist);
            printf("node terminal %d, %d \n\n", k, i);
            assert(Is_term(g->term[i]));
            if( i != -1 && !SCIPisEQ(scip, pathdist[k][i], vnoi[k + 2 *  nnodes].dist ) )
            {
               printf("FAILS2 %f, %f ", pathdist[k][i], vnoi[k + 2 *  nnodes].dist);
               printf("node terminal %d, %d \n\n", k, i);
               assert(0);
            }

            i = vbase[k + 3 * nnodes];
            printf("x4 %f, %f ", pathdist[k][i], vnoi[k + 3 *  nnodes].dist);
            printf("node terminal %d, %d \n\n", k, i);
            assert(i == -1 || Is_term(g->term[i]));
            if( i != -1 && SCIPisLT(scip, pathdist[k][i], vnoi[k + 3 *  nnodes].dist ) )
            {
               printf("FAILS3 %f, %f ", pathdist[k][i], vnoi[k + 3 *  nnodes].dist);
               printf("node terminal %d, %d \n\n", k, i);
            }

         }
      }

      for( k = nnodes - 1; k >= 0; k-- )
      {
         SCIPfreeBufferArrayNull(scip, &(pathedge[k]));
         SCIPfreeBufferArrayNull(scip, &(pathdist[k]));
      }

      SCIPfreeBufferArray(scip, &pathedge);
      SCIPfreeBufferArray(scip, &pathdist);

   }
#endif
   /* construct auxiliary graph to compute paths between terminals */

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1, 0) );

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
      forbidden[k] = FALSE;
      edgepreds[k] = -1.0;
   }

   for( k = nnodes - 1; k < nedges; k++ )
   {
      forbidden[k] = FALSE;
      edgepreds[k] = -1.0;
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   for( k = 0; k < nnodes; k++ )
   {
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
                  netgraph->cost[Edge_anti(ne)] = ecost;
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
   netgraph->source[0] = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
   SCIP_CALL( graph_path_init(scip, netgraph) );

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
#if 0
      i1 = g->tail[e];
      i2 = g->head[e];
      while( i1 != vbase[i1] )
      {
         e = vnoi[i1].edge;
         forbidden[e] = TRUE;
         forbidden[flipedge(e)] = TRUE;
         i1 = g->tail[e];
      }

      while( i2 != vbase[i2] )
      {
         e = vnoi[i2].edge;
         forbidden[e] = TRUE;
         forbidden[flipedge(e)] = TRUE;
         i2 = g->tail[e];
      }
#endif
   }

   /* traverse all edges */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
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

         /* is i a terminal? If not, get three closest terminals of distance smaller ecost */
         if( !Is_term(g->term[i]) )
         {
#if 0
            if( Is_term(g->term[i2]) )
               nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
            else
               nnterms2 = getlecloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#endif
            nnterms1 = getlecloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);

            if( nnterms1 == 0 )
               continue;
         }

         nnterms2 = 0;

         if( Is_term(g->term[i2]) )
         {
            nnterms2 = 1;
            termdist2[0] = 0.0;
            neighbterms2[0] = i2;
         }
         else
         {
            /* get closest terminals of distance smaller ecost */
#if 0
            if(  Is_term(g->term[i]) )
               nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
            else
               nnterms2 = getlecloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#endif
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

                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               } /* tj != tk (else) */
            } /* k < nnterms2 */
         } /* j < nnterms1 */

      } /* while( enext != EAT_LAST ) */
   }

   SCIP_CALL( bdr_reduction(scip, g, netgraph, mst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, nelims) );

   /* free memory*/
   graph_path_exit(scip, netgraph);
   graph_free(scip, netgraph, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &neighbterms2);
   SCIPfreeBufferArray(scip, &neighbterms1);
   SCIPfreeBufferArray(scip, &termdist2);
   SCIPfreeBufferArray(scip, &termdist1);

   return SCIP_OKAY;
}


/** SD test for PC */
SCIP_RETCODE sdpc_reduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            boundedges,         /**< array to store bound edges */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nodesid,            /**< array */
   int*                  nodesorg,           /**< array */
   int*                  nelims              /**< pointer to store number of eliminated edges */
   )
{
   GRAPH* netgraph;
   SCIP_Real* termdist1;
   SCIP_Real* termdist2;
   SCIP_Real ecost;
   SCIP_Real necost;
   int* neighbterms1;
   int* neighbterms2;
   int e;
   int e2;
   int i;
   int i2;
   int j;
   int k;
   int l;
   int tj;
   int tk;
   int ne;
   int pos;
   int root;
   int enext;
   int nnodes;
   int nterms;
   int nedges;
   int nnterms1;
   int nnterms2;
   int maxnedges;

   assert(g != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(nodesid != NULL);
   assert(nodesorg != NULL);
   assert(boundedges != NULL);

   root = g->source[0];
   nnodes = g->knots;
   nterms = g->terms;
   nedges = g->edges;
   *nelims = 0;
   maxnedges = MIN(nedges, (nterms - 1) * nterms);

   SCIP_CALL( SCIPallocBufferArray(scip, &termdist1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termdist2, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms2, 4) );

   /* compute nearest four terminals to each non-terminal */
   getnext4terms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);

   /* construct auxiliary graph to compute paths between terminals */

   /* initialize new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1, 0) );

   for( k = 0; k < 4; k++ )
   {
      termdist1[k] = FARAWAY;
      termdist2[k] = FARAWAY;
   }

   j = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
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
   for( k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( !g->mark[g->head[e]] )
            continue;
         i = vbase[k];
         assert(i != UNKNOWN);
         if( i != vbase[g->head[e]] )
         {
            i2 = vbase[g->head[e]];
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
                  netgraph->cost[Edge_anti(ne)] = ecost;
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
   getnext4tterms(scip, g, g->cost, boundedges, vnoi, vbase, heap, state);

   /* traverse all edges */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         e = enext;
         enext = g->oeat[e];
         i2 = g->head[e];
         if( i2 < i || Is_term(g->term[i2]) || !g->mark[i2] )
            continue;
         ecost = g->cost[e];

         /* @todo: fix */
#if 0
         if( Is_term(g->term[i]) )
         {
            for( k = 0; k < 4; k++ )
               neighbterms1[k] = UNKNOWN;
            nnterms1 = 0;
            /* get three nearest terms */
            for( ne = netgraph->outbeg[nodesid[i]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
            {
               if( SCIPisLT(scip, netgraph->cost[ne], ecost) )
               {
                  j = nodesorg[netgraph->head[ne]];
                  necost = netgraph->cost[ne];
                  assert(Is_term(g->term[j]));
                  assert(j != i2);
                  if( nnterms1 < 4 )
                     nnterms1++;
                  for( k = 0; k < 4; k++ )
                  {
                     if( neighbterms1[k] == UNKNOWN || SCIPisGT(scip, termdist1[k], necost) )
                     {
                        for( l = 3; l > k; l-- )
                        {
                           neighbterms1[l] = neighbterms1[l - 1];
                           termdist1[l] = termdist1[l - 1];
                        }
                        neighbterms1[k] = j;
                        termdist1[k] = necost;
                        break;
                     }
                  }
               }
            }
         }
         else
         {
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
         }
#endif
         nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);

         if( nnterms1 == 0 )
            continue;

         /* @todo: fix */
#if 0
         if( Is_term(g->term[i2]) )
         {
            for( k = 0; k < 4; k++ )
               neighbterms2[k] = UNKNOWN;
            nnterms2 = 0;
            /* get three nearest terms */
            for( ne = netgraph->outbeg[nodesid[i2]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
            {
               if( SCIPisLT(scip, netgraph->cost[ne], ecost) )
               {
                  j = nodesorg[netgraph->head[ne]];
                  necost = netgraph->cost[ne];
                  assert(Is_term(g->term[j]));
                  assert(j != i);
                  if( nnterms2 < 4 )
                     nnterms2++;
                  for( k = 0; k < 4; k++ )
                  {
                     if( neighbterms2[k] == UNKNOWN || SCIPisGT(scip, termdist2[k], necost) )
                     {
                        for( l = 3; l > k; l-- )
                        {
                           neighbterms2[l] = neighbterms2[l - 1];
                           termdist2[l] = termdist2[l - 1];
                        }
                        neighbterms2[k] = j;
                        termdist2[k] = necost;
                        break;
                     }
                  }
               }
            }
         }
         else
         {
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
         }
#endif
         nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);

         if( nnterms2 == 0 )
            continue;

         /* @todo: mark nearest terminals!!!! */
         for( j = 0; j < nnterms1; j++ )
         {
            /* has edge already been deleted? */
            if( g->oeat[e] == EAT_FREE )
               break;
            tj = neighbterms1[j];
            assert(tj >= 0);
            assert(Is_term(g->term[tj]));
            for( k = 0; k < nnterms2; k++ )
            {
               tk = neighbterms2[k];
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
                  necost = FARAWAY;

                  for( e2 = netgraph->outbeg[nodesid[tj]]; e2 != EAT_LAST; e2 = netgraph->oeat[e2] )
                  {
                     if( netgraph->head[e2] == nodesid[tk] )
                     {
                        necost = netgraph->cost[e2];

                        break;
                     }
                  }
                  pos = tj;
                  for( l = 0; l < 4; l++ )
                  {
                     if( vbase[pos] == UNKNOWN )
                        break;
                     if( vbase[pos] == tk && SCIPisLT(scip, vnoi[pos].dist, necost) )
                     {
                        necost = vnoi[pos].dist;
                        e2 = 0;
                        break;
                     }
                     pos += nnodes;
                  }

#if 1
                  if( e2 != EAT_LAST && SCIPisGT(scip, ecost, necost)
                     && SCIPisGT(scip, ecost, necost + termdist1[j] - g->prize[tj])
                     && SCIPisGT(scip, ecost, necost + termdist2[k] - g->prize[tk])
                     && SCIPisGT(scip, ecost, necost + termdist1[j] + termdist2[k] - g->prize[tj] - g->prize[tk]) )
                  {
                     SCIPdebugMessage("SDSP delete: %d %d (%d)\n", g->tail[e], g->head[e], e);
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
#endif
               }
            }
         }
      }
   }

   graph_free(scip, netgraph, TRUE);
   SCIPfreeBufferArray(scip, &neighbterms2);
   SCIPfreeBufferArray(scip, &neighbterms1);
   SCIPfreeBufferArray(scip, &termdist2);
   SCIPfreeBufferArray(scip, &termdist1);
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
   sd = FARAWAY;
   /* compare restriced sd with edge cost (if existing) */
   for( e = g->outbeg[i]; (l++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
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
SCIP_RETCODE getSD(
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
   SCIP_Bool pcmw;

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

   pcmw = pc || mw;
   nnodes = g->knots;

   /* start from tail */
   sdpaths(scip, g, pathtail, g->cost, distlimit, heap, statetail, memlbltail, &nlbltail, i, i2, limit);

   /* test whether edge e can be eliminated */
   sdpaths(scip, g, pathhead, g->cost, distlimit, heap, statehead, memlblhead, &nlblhead, i2, i, limit);

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
   /* compare restriced sd with edge cost (if existing) */
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


/** SD test using only two edges */
SCIP_RETCODE sd2_reduction(
   SCIP* scip,
   GRAPH* g,
   SCIP_Real* nodecost,
   int*    nelims,
   int*    adjacent
   )
{
   SCIP_Real cost;
   int i;
   int e;
   int i2;
   int i3;
   int e2;
   int enext;
   int nnodes;
   SCIP_Bool pc;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(adjacent != NULL);
   assert(nodecost != NULL);

   nnodes = g->knots;
   *nelims = 0;
   pc = g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   for( i = 0; i < nnodes; i++ )
      adjacent[i] = FALSE;

   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;

      /* mark neighbours */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         i2 = g->head[e];
         if( g->mark[i2] )
         {
            adjacent[i2] = TRUE;
            nodecost[i2] = g->cost[e];
         }
      }
      assert(!adjacent[i]);

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         enext = g->oeat[e];
         i2 = g->head[e];
         /* avoid double check */
         if( i2 < i || !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         /* test whether edge e can be eliminated */
         cost = g->cost[e];

         for( e2 = g->outbeg[i2]; e2 != EAT_LAST; e2 = g->oeat[e2] )
         {
            i3 = g->head[e2];
            if( !adjacent[i3] )
               continue;

            assert(g->mark[i3]);

            if( (Is_term(g->term[i3]) && SCIPisGE(scip, cost, g->cost[e2]) && SCIPisGE(scip, cost, nodecost[i3])) ||
               (!Is_term(g->term[i3]) && SCIPisGE(scip, cost, g->cost[e2] + nodecost[i3])) )
            {
               if( pc && Is_term(g->term[i3]) && !SCIPisGE(scip, cost, g->cost[e2] + nodecost[i3] - g->prize[i3]) )
                  continue;

               adjacent[i2] = FALSE;
               graph_edge_del(scip, g, e, TRUE);
               (*nelims)++;
               break;
            }
         }
         e = enext;
      }

      /* restore neigbour adjacency array */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjacent[g->head[e]] = FALSE;
   }

   return SCIP_OKAY;
}



/** SDC test for the SAP using a limited version of Dijkstra's algorithm from both endpoints of an arc */
SCIP_RETCODE sdsp_sap_reduction(
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
         sdpaths(scip, g, pathtail, g->cost, g->cost[e], heap, statetail, memlbltail, &nlbltail, i, i2, limit);

         /* start limited dijkstra from i2, marking all reached vertices */
         sdpaths(scip, g, pathhead, costrev, g->cost[e], heap, statehead, memlblhead, &nlblhead, i2, i, limit);

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
               printf("ELIM edge \n");
            else
               printf("counteredge \n");

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
SCIP_RETCODE sdsp_reduction(
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
   SCIP_Real dist;
   SCIP_Real sdist;
   int i;
   int k;
   int l;
   int e;
   int i2;
   int enext;
   int nnodes;
   int nlbltail;
   int nlblhead;
   SCIP_Bool pc;

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

   nnodes = g->knots;
   *nelims = 0;
   pc = g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

   for( i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

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

         /* avoid double checking */
         if( i2 < i || !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         /* start limited dijkstra from i, marking all reached vertices */
         sdpaths(scip, g, pathtail, g->cost, g->cost[e], heap, statetail, memlbltail, &nlbltail, i, i2, limit);

         /* start limited dijkstra from i2, marking all reached vertices */
         sdpaths(scip, g, pathhead, g->cost, g->cost[e], heap, statehead, memlblhead, &nlblhead, i2, i, limit);

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
               if( Is_term(g->term[l]) )
               {
                  dist = 0.0;
                  if( SCIPisLT(scip, dist, pathhead[l].dist) )
                     dist = pathhead[l].dist;
                  if( SCIPisLT(scip, dist, pathtail[l].dist) )
                     dist = pathtail[l].dist;
                  if( pc && SCIPisLT(scip, dist, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                     dist = pathhead[l].dist + pathtail[l].dist - g->prize[l];
                  if( SCIPisGT(scip, sdist, dist) )
                     sdist = dist;
               }
               else
               {
                  if( SCIPisGT(scip, sdist, pathhead[l].dist + pathtail[l].dist) )
                     sdist = pathhead[l].dist + pathtail[l].dist;
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
            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;
         }

         e = enext;
      }
   }
   return SCIP_OKAY;
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
   SCIP_Real*  sddist,
   SCIP_Real*  sdtrans,
   SCIP_Real*  sdrand,
   SCIP_Real* cost,
   SCIP_Real* randarr,
   int*    heap,
   int*    state,
   int*    knotexamined,
   int*    elimins,
   int     runnum,
   unsigned int* seed
   )
{
   SCIP_Real redstarttime;
   SCIP_Real timelimit;
   SCIP_Real stalltime;
   int     count = 0;
   int     i;
   int     e;
   int     j;
   int     nnodes;
   int     knotoffset = 0;

   SCIPdebugMessage("SD-Reduktion: ");

   assert(g != NULL);
   assert(scip != NULL);
   assert(heap  != NULL);
   assert(state != NULL);
   assert(sddist != NULL);
   assert(sdtrans != NULL);
   assert(cost != NULL);
   assert(knotexamined != NULL);

   nnodes = g->knots;
   *elimins = 0;
   redstarttime = SCIPgetTotalTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   stalltime = timelimit * 0.1; /* @todo this should be set as a parameter */

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i] = (g->grad[i] > 0);
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         randarr[e] = SCIPgetRandomReal(0.0, (g->cost[e]), seed);/* @todo: org (double)(rand() % 512); */
         cost[e] = g->cost[e] * 1.0 + randarr[e];
      }
   }

   /* this is the offset used to minimise the number of knots to examine in large graphs. */
   if( nnodes > KNOTLIMIT )
   {
      srand(((unsigned int)runnum * 100));
      i = 0;
      do
      {
         knotoffset = rand() % KNOTFREQ;
         i++;
      } while( nnodes > KNOTLIMIT && knotexamined[knotoffset] >= 0 && i < 50 );
      knotexamined[knotoffset]++;
   }


   for( i = 0; i < nnodes; i++ )
   {
      if( i % 100 == 0 && *elimins == 0 && SCIPgetTotalTime(scip) - redstarttime > stalltime)
         break;

      if( !g->mark[i] || g->grad[i] == 0 )
         continue;

      if( nnodes > KNOTLIMIT && i % KNOTFREQ != knotoffset )
         continue;

#if 0
      /* For the prize collecting variants all edges from the "dummy" root node must be retained. */
      if ( (g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING
            || g->stp_type == STP_MAX_NODE_WEIGHT) && i == g->source[0] )
         continue;
#endif

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->mark[g->head[e]] == 1);
         g->mark[g->head[e]] = 2;
      }

      compute_sd(scip, g, i, cost, randarr, heap, state, &count, sddist, sdtrans, sdrand);

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->mark[g->head[e]] == 2);
         g->mark[g->head[e]] = 1;
      }

      for( e = g->outbeg[i]; e != EAT_LAST; e = j )
      {
         assert(g->tail[e] == i);

         j = g->oeat[e];

         if( SCIPisLT(scip, g->cost[e], FARAWAY) && SCIPisLT(scip, sddist[g->head[e]], cost[e])
            && SCIPisLT(scip, sddist[g->head[e]] - sdrand[g->head[e]], cost[e] - randarr[e]) )
         {
            graph_edge_del(scip, g, e, TRUE);
            (*elimins)++;
         }
      }
   }

   assert(graph_valid(g));

   SCIPdebugMessage("%d Edges deleted\n", *elimins * 2);
   return SCIP_OKAY;
}

#if 0
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
               graph_edge_del(scip, g, e, TRUE);

            (*elimins)++;

         }
      }
   }

   assert(graph_valid(g));

   SCIPdebugMessage("%d Edges deleted\n", *elimins * 2);

   free(sourceadj);

   return SCIP_OKAY;
}
#endif



SCIP_RETCODE bdr_reduction(
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
   GRAPH* auxg;
   PATH* mst;
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real csum;
   SCIP_Real mstcost;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   int    i;
   int    k;
   int    e;
   int    l;
   int    k2;
   int    tail;
   int    head;
   int    nnodes;
   int*   adjvert;
   int*   incedge;
   int*   reinsert;
   SCIP_Bool pc;

   assert(g != NULL);
   assert(netgraph  != NULL);
   assert(netmst != NULL);

   /* initialize mst struct and new graph for bd4, bd5 tests */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, 5) );
   SCIP_CALL( graph_init(scip, &auxg, 5, 40, 1, 0) );

   for( k = 0; k < 4; k++ )
      graph_knot_add(auxg, -1);
   for( k = 0; k < 4; k++ )
      for( k2 = k + 1; k2 < 4; k2++ )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sd, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ecost, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adjvert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &reinsert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &incedge, 4) );

   SCIPdebugMessage("BD3-R Reduction: ");

   nnodes = g->knots;
   pc = g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   for( i = 0; i < 4; i++ )
   {
      sd[i] = 0.0;
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }
   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) || (g->grad[i] != 3 && g->grad[i] != 4) )
         continue;

      k = 0;
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         incedge[k] = e;
         ecost[k] = g->cost[e];
         adjvert[k++] = g->head[e];
         assert(k <= 4);
      }

      /* vertex of degree 3? */
      if( g->grad[i] == 3 )
      {
         assert(k == 3);

         csum = ecost[0] + ecost[1] + ecost[2];

         sd[0] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, adjvert[0], adjvert[1], 300);
         sd[1] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, adjvert[1], adjvert[2], 300);
         sd[2] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, adjvert[2], adjvert[0], 300);

         if( SCIPisLE(scip, sd[0] + sd[1], csum) || SCIPisLE(scip, sd[0] + sd[2], csum) || SCIPisLE(scip, sd[1] + sd[2], csum) )
         {
            SCIPdebugMessage("BD3-R Reduction: %f %f %f csum: %f\n ", sd[0], sd[1], sd[2], csum);
            (*nelims)++;
            /* save ancestors */
            for( k = 0; k < 3; k++ )
            {
               SCIPintListNodeFree(scip, &(ancestors[k]));
               SCIPintListNodeFree(scip, &(revancestors[k]));
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[k]), g->ancestors[incedge[k]]) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[k]), g->ancestors[Edge_anti(incedge[k])]) );
            }

            for( k = 0; k < 3; k++ )
            {
               k2 = (k + 1) % 3;

               if( SCIPisLT(scip, sd[k], ecost[k] + ecost[k2]) )
               {
                  graph_edge_del(scip, g, incedge[k], TRUE);
               }
               else
                  SCIP_CALL( graph_edge_reinsert(scip, g, incedge[k], adjvert[k], adjvert[k2], ecost[k] + ecost[k2], ancestors[k], ancestors[k2], revancestors[k], revancestors[k2]) );
            }

            assert(g->grad[i] == 0);
         }
      }
      /* vertex of degree 4? */
      else if( g->grad[i] == 4 )
      {
         assert(k == 4);

         for( k = 0; k < 4; k++ )
         {
            auxg->mark[k] = TRUE;
            for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
            {
               k2 = auxg->head[e];
               if( k2 > k )
               {
                  auxg->cost[e] = getRSD(scip, g, netgraph, netmst, vnoi, mstsdist, termdist1, termdist2, vbase, nodesid, neighbterms1, neighbterms2, adjvert[k], adjvert[k2], 200);
                  auxg->cost[flipedge(e)] = auxg->cost[e];
               }
            }
         }

         for( l = 0; l < 4; l++ )
            mst[l].dist = UNKNOWN;
         /* compute mst on all neighbours */
         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
         mstcost = 0.0;
         for( l = 1; l < 4; l++ )
            assert(mst[l].dist != UNKNOWN);
         for( l = 1; l < 4; l++ )
            mstcost += mst[l].dist;

         k = UNKNOWN;
         csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

         if( SCIPisGE(scip, csum, mstcost) )
         {
            /* compute mst on all 3-subsets of all neigbours */
            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = FALSE;
               mstcost = 0.0;

               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( l = 1; l < 4; l++ )
                     if( auxg->mark[l] )
                        mstcost += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 4; l++ )
                     mstcost += mst[l].dist;
               }

               auxg->mark[k] = TRUE;
               csum -= ecost[k];

               if( SCIPisLT(scip, csum, mstcost) )
                  break;

               csum += ecost[k];
            }
         }

         if( k == 4 )
         {
            l = 0;
            for( k = 0; k < 4; k++ )
            {
               for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     if( SCIPisGE(scip, auxg->cost[e], ecost[k] + ecost[k2]) )
                     {
                        if( l >= 4 )
                           break;
                        reinsert[l++] = e;
                     }
                  }
               }
            }

            if( k == 4 )
            {
               (*nelims) += g->grad[i] - l;

               /* save ancestors */
               for( k = 0; k < 4; k++ )
               {
                  SCIPintListNodeFree(scip, &(ancestors[k]));
                  SCIPintListNodeFree(scip, &(revancestors[k]));
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[k]), g->ancestors[incedge[k]]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[k]), g->ancestors[Edge_anti(incedge[k])]) );
               }

               for( k = 0; k < l; k++ )
               {
                  e = reinsert[k];
                  tail = auxg->tail[e];
                  head = auxg->head[e];
                  SCIP_CALL( graph_edge_reinsert(scip, g, incedge[k], adjvert[tail], adjvert[head], ecost[tail] + ecost[head], ancestors[tail], ancestors[head], revancestors[tail], revancestors[head]) );
               }

               /* delete remaining edges */
               while( g->outbeg[i] >= 0 )
                  graph_edge_del(scip, g, g->outbeg[i], TRUE);
            }
            assert(g->grad[i] == 0);
         }

      }
   }

   for( k = 0; k < 4; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &incedge);
   SCIPfreeBufferArray(scip, &reinsert);
   SCIPfreeBufferArray(scip, &adjvert);
   SCIPfreeBufferArray(scip, &ecost);
   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);
   SCIPfreeBufferArray(scip, &sd);

   graph_path_exit(scip, auxg);
   graph_free(scip, auxg, TRUE);
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}




/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Bottleneck Degree 3,4 Test
 */
SCIP_RETCODE bd3_reduction(
   SCIP* scip,
   GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int*    nelims,
   int     limit
   )
{
   GRAPH* auxg;
   PATH* mst;
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real s1;
   SCIP_Real csum;
   SCIP_Real mstcost;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   int    i;
   int    k;
   int    e;
   int    l;
   int    k2;
   int    tail;
   int    head;
   int    nnodes;
   int*   adjvert;
   int*   incedge;
   int*   reinsert;

   SCIP_Bool pc;

   SCIPdebugMessage("BD3-Reduction: ");

   assert(g != NULL);
   assert(heap  != NULL);
   assert(nelims != NULL);

   /* initialize mst struct and new graph for bd4, bd5 tests */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, 5) );
   SCIP_CALL( graph_init(scip, &auxg, 5, 40, 1, 0) );

   for( k = 0; k < 4; k++ )
      graph_knot_add(auxg, -1);
   for( k = 0; k < 4; k++ )
      for( k2 = k + 1; k2 < 4; k2++ )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sd, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ecost, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adjvert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &reinsert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &incedge, 4) );

   *nelims = 0;
   nnodes = g->knots;
   pc = g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   for( i = 0; i < 4; i++ )
   {
      sd[i] = 0.0;
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }
   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
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
      if( Is_term(g->term[i]) || (g->grad[i] != 3 && g->grad[i] != 4) )
         continue;

      k = 0;
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         incedge[k] = e;
         ecost[k] = g->cost[e];
         adjvert[k++] = g->head[e];
         assert(k <= 4);
      }

      /* vertex of degree 3? */
      if( g->grad[i] == 3 )
      {
         assert(k == 3);

         csum = ecost[0] + ecost[1] + ecost[2];

         SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sd[0]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[0], adjvert[1], limit, pc, FALSE) );
         SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sd[1]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[1], adjvert[2], limit, pc, FALSE) );
         SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sd[2]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[2], adjvert[0], limit, pc, FALSE) );

         if( SCIPisLE(scip, sd[0] + sd[1], csum) || SCIPisLE(scip, sd[0] + sd[2], csum) || SCIPisLE(scip, sd[1] + sd[2], csum) )
         {
            (*nelims)++;

            /* save ancestors */
            for( k = 0; k < 3; k++ )
            {
               SCIPintListNodeFree(scip, &(ancestors[k]));
               SCIPintListNodeFree(scip, &(revancestors[k]));
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[k]), g->ancestors[incedge[k]]) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[k]), g->ancestors[Edge_anti(incedge[k])]) );
            }

            for( k = 0; k < 3; k++ )
            {
               k2 = (k + 1) % 3;

               if( SCIPisLT(scip, sd[k], ecost[k] + ecost[k2]) )
               {
                  graph_edge_del(scip, g, incedge[k], TRUE);
               }
               else
                  SCIP_CALL( graph_edge_reinsert(scip, g, incedge[k], adjvert[k], adjvert[k2], ecost[k] + ecost[k2], ancestors[k], ancestors[k2], revancestors[k], revancestors[k2]) );
            }
#if 0
            if( SCIPisLE(scip, s1, ecost[0] + ecost[1]) )
               graph_edge_del(scip, g, incedge[0], TRUE);
            else
               SCIP_CALL( graph_edge_reinsert(scip, g, incedge[0], adjvert[0], adjvert[1], ecost[0] + ecost[1], ancestors[0], ancestors[1], revancestors[0], revancestors[1]) );

            if( SCIPisLE(scip, s2, ecost[1] + ecost[2]) )
               graph_edge_del(scip, g, incedge[1], TRUE);
            else
               SCIP_CALL( graph_edge_reinsert(scip, g, incedge[1], adjvert[1], adjvert[2], ecost[1] + ecost[2], ancestors[1], ancestors[2], revancestors[1], revancestors[2]) );

            if( SCIPisLE(scip, s3, ecost[0] + ecost[2]) )
               graph_edge_del(scip, g, incedge[2], TRUE);
            else
               SCIP_CALL( graph_edge_reinsert(scip, g, incedge[2], adjvert[2], adjvert[0], ecost[2] + ecost[0], ancestors[2], ancestors[0], revancestors[2], revancestors[0]) );
#endif
#if 0
            {
               n1 = graph_edge_redirect(g, incedge[2], adjvert[2], adjvert[0], ecost[2] + ecost[0]);
               if( n1 >= 0 )
               {
                  SCIPintListNodeFree(scip, &(g->ancestors[n1]));
                  SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
                  SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors[2]) );
                  SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors[0]) );

                  SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors[2]) );
                  SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors[0]) );
               }
            }
#endif
            assert(g->grad[i] == 0);
         }

      }
      /* vertex of degree 4? */
      else if( g->grad[i] == 4 )
      {
         assert(k == 4);
         csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

         for( k = 0; k < 4; k++ )
         {
            auxg->mark[k] = TRUE;
            for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
            {
               k2 = auxg->head[e];
               if( k2 > k )
               {
                  SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[k], adjvert[k2], limit, pc, FALSE) );
                  auxg->cost[e] = s1;
                  auxg->cost[flipedge(e)] = s1;
               }
            }
         }

         for( l = 0; l < 4; l++ )
            mst[l].dist = UNKNOWN;
         /* compute mst on all neighbours */
         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
         mstcost = 0.0;
         for( l = 1; l < 4; l++ )
            assert(mst[l].dist != UNKNOWN);
         for( l = 1; l < 4; l++ )
            mstcost += mst[l].dist;

         k = UNKNOWN;
         if( SCIPisGE(scip, csum, mstcost) )
         {
            /* compute mst on all 3-subsets of all neigbours */
            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = FALSE;
               mstcost = 0.0;

               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( l = 1; l < 4; l++ )
                     if( auxg->mark[l] )
                        mstcost += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 4; l++ )
                     mstcost += mst[l].dist;
               }

               auxg->mark[k] = TRUE;
               csum -= ecost[k];

               if( SCIPisLT(scip, csum, mstcost) )
                  break;

               csum += ecost[k];
            }
         }

         if( k == 4 )
         {
            l = 0;
            for( k = 0; k < 4; k++ )
            {
               for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     if( SCIPisGE(scip, auxg->cost[e], ecost[k] + ecost[k2]) )
                     {
                        if( l >= 4 )
                           break;

                        reinsert[l++] = e;
                     }
                  }
               }
            }

            if( k == 4 )
            {
               SCIPdebugMessage("npv4Reduction delete: %d\n", i);
               (*nelims) += g->grad[i] - l;

               /* save ancestors */
               for( k = 0; k < 4; k++ )
               {
                  SCIPintListNodeFree(scip, &(ancestors[k]));
                  SCIPintListNodeFree(scip, &(revancestors[k]));
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[k]), g->ancestors[incedge[k]]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[k]), g->ancestors[Edge_anti(incedge[k])]) );
               }

               for( k = 0; k < l; k++ )
               {
                  e = reinsert[k];
                  tail = auxg->tail[e];
                  head = auxg->head[e];
                  SCIP_CALL( graph_edge_reinsert(scip, g, incedge[k], adjvert[tail], adjvert[head], ecost[tail] + ecost[head], ancestors[tail], ancestors[head], revancestors[tail], revancestors[head]) );
               }

               /* delete remaining edges */
               while( g->outbeg[i] >= 0 )
                  graph_edge_del(scip, g, g->outbeg[i], TRUE);
            }
            assert(g->grad[i] == 0);
         }

      }
   }

   for( k = 0; k < 4; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &incedge);
   SCIPfreeBufferArray(scip, &reinsert);
   SCIPfreeBufferArray(scip, &adjvert);
   SCIPfreeBufferArray(scip, &ecost);
   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);
   SCIPfreeBufferArray(scip, &sd);

   graph_path_exit(scip, auxg);
   graph_free(scip, auxg, TRUE);
   SCIPfreeBufferArray(scip, &mst);

   SCIPdebugMessage("bd3: %d Knots deleted\n", *nelims);

   return SCIP_OKAY;
}

#if 0
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

   calculate_distances(scip, g, path, g->cost, FSP_MODE);

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

      graph_path_exec(scip, g, MST_MODE, i, g->cost, mst1);

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

      graph_path_exec(scip, g, MST_MODE, i, cost, mst2);

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
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]) );
         SCIP_CALL( graph_knot_contract(scip, g, i, k) );

         (*nelims)++;

         calculate_distances(scip, g, path, g->cost, FSP_MODE);

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
   graph_path_exec(scip, g, FSP_MODE, g->source[0], g->cost, pathfromsource);

   /* computing the shortest hops paths from the source node */
   graph_path_exec(scip, g, FSP_MODE, g->source[0], hopscost, pathhops);

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
                           graph_edge_del(scip, g, e, TRUE);
                     }
                  }
                  (*elimins)++;
               }
            }
            else
            {
               *fixed += min1;
               SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[shortarc]) ); /* I think that this should be
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
            graph_path_exec(scip, g, FSP_MODE, g->source[0], g->cost, pathfromsource);
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
            SCIPintListNodeAppendCopy(&(g->fixedges), g->ancestors[e]);
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

#endif
/* shortest link reduction */
SCIP_RETCODE sl_reduction(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   double* fixed,
   int* heap,
   int* state,
   int* vbase,
   int*  vrnodes,
   char*  visited,
   int* nelims
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
#if 0
   int     e2;
#endif
   int     old;
   int     head;
   int     tail;
   int     root;
   int     nnodes;
   int     vrcount;
   int     minedge;
#if 0
   int     minedge2;
#endif
   int*    qnode;
   char    contract;
   char*   forbidden;
   SCIP_Bool pc;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(vrnodes != NULL);
   assert(visited != NULL);

   *nelims = 0;
   nnodes = g->knots;
   root = g->source[0];
   pc = (g->stp_type == STP_PRIZE_COLLECTING) || (g->stp_type == STP_ROOTED_PRIZE_COLLECTING);

   SCIP_CALL( SCIPallocBufferArray(scip, &forbidden, nnodes) );

   if( !pc )
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

   voronoi_terms(scip, g, g->cost, vnoi, vbase, heap, state);

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );
   for( j = 0; j < nnodes; j++ )
   {
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
#if 0
         minedge2 = UNKNOWN;
#endif
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
#if 0
                     minedge2 = minedge;
#endif
                     minedge = e;
                  }
                  else if( SCIPisLE(scip, cost, mincost2) )
                  {
#if 0
                     minedge2 = e;
#endif
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
         /* @todo fix */
#if 0
         else if( minedge2 != UNKNOWN && !Is_term(g->term[g->head[minedge2]]) && SCIPisGE(scip, mincost3, cost) )
         {
            t = g->head[minedge2];
            for( e2 = g->outbeg[t]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               if( e2 != flipedge(minedge2) && SCIPisLT(scip, g->cost[e2], cost) && vbase[g->head[e2]] != i )
                  break;

            if( e2 == EAT_LAST )
               contract = TRUE;
            if( contract )
               printf("sladv \n");

         }
#endif

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
               SCIP_CALL( graph_knot_contractpc(scip, g, j, k, i) );
            }
            else
            {
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]) );
               SCIP_CALL( graph_knot_contract(scip, g, j, k) );
            }

            assert(old - g->grad[j] - g->grad[k] > 0);
            (*nelims) += old - g->grad[j] - g->grad[k];
            forbidden[vbase[j]] = TRUE;
            forbidden[vbase[k]] = TRUE;
         }
      }
   }

   /* free memory */
   SCIPqueueFree(&queue);
   SCIPfreeBufferArray(scip, &forbidden);

   return SCIP_OKAY;
}


/* NV reduction from T. Polzin's "Algorithms for the Steiner problem in networks" */
SCIP_RETCODE nv_reduction(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   double* fixed,
   int* edgearrint,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
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
   pc = (g->stp_type == STP_PRIZE_COLLECTING) || (g->stp_type == STP_ROOTED_PRIZE_COLLECTING);
   SCIP_CALL( SCIPallocBufferArray(scip, &term, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distance, nnodes) );

   /* minimal grad of a vertex to be scrutinized */
   if( pc )
   {
      if( g->stp_type == STP_ROOTED_PRIZE_COLLECTING )
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
   SCIP_CALL( voronoi_dist(scip, g, g->cost, distance, edgearrint, vbase, minedge1, heap, state, distnode, vnoi) );

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
         if( i != g->source[0] )
            pi = g->prize[i];
         else
            pi = FARAWAY;

         if( t != g->source[0] )
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
            SCIP_CALL( graph_knot_contractpc(scip, g, i, k, i) );
         }
         else
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1]) );
            SCIP_CALL( graph_knot_contract(scip, g, i, k) );
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
SCIP_RETCODE nv_reductionAdv(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 vnoi,
   SCIP_Real*            distance,
   double*               fixed,
   int*                  edgearrint,
   int*                  heap,
   int*                  state,
   int*                  vbase,
   int*                  neighb,
   int*                  distnode,
   int*                  nelims
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
   pc = (g->stp_type == STP_PRIZE_COLLECTING) || (g->stp_type == STP_ROOTED_PRIZE_COLLECTING);
   nnodes = g->knots;
   nterms = g->terms;
   *nelims = 0;
   termcount = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &term, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minedge1, nterms) );

   /* set minimal grad of a vertex to be scrutinized */
   if( pc )
   {
      if( g->stp_type == STP_ROOTED_PRIZE_COLLECTING )
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
   SCIP_CALL( voronoi_dist(scip, g, g->cost, distance, edgearrint, vbase, minedge1, heap, state, distnode, vnoi) );

#if 0
   printf("rootmarkterm?: %d \n", Is_term(g->term[g->source[0]]));
   printf("rootmarked?: %d \n", g->mark[g->source[0]]);
   printf("vbase 187?: %d \n", vbase[187]);

   int a;
   int h;
   int* pnode;
   int v = 187;
   SCIP_QUEUE* queue;

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes + 1, 2.0) );
   for( k = 0; k < nnodes; k++ )
      g->mark[k] = FALSE;
   g->mark[v] = TRUE;
   SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[g->outbeg[v]])) );

   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));

      /* traverse incoming arcs */
      for( a = g->outbeg[*pnode]; a != EAT_LAST; a = g->oeat[a] )
      {
         h = g->head[a];

         if( !g->mark[h] && h != 0 )
         {
            printf("goto: %d->%d \n", *pnode, h);
            /* if an active vertex has been hit, break */
            if( Is_term(g->term[h]) )
            {
               printf("hit: %d \n", h);
               assert(0);
            }

            g->mark[h] = TRUE;

            SCIP_CALL( SCIPqueueInsert(queue, &(g->head[a])) );
         }

      }
   }

   assert(0);
#endif
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
#if 0
      printf("root: %d \n", g->source[0]);
      printf("edge: %d->%d\n", g->tail[edge1], g->head[edge1]);
      printf("i: %d distance: %f \n", i, distance[0]);
      assert(0);
#endif
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
         if( i != g->source[0] )
            pi = g->prize[i];
         else
            pi = FARAWAY;

         if( t == UNKNOWN )
            pt = -1.0;
         else if( t != g->source[0] )
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
#if 0
         contract = FALSE;
#endif
      }

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         neighb[g->head[e]] = FALSE;

      if( contract && (!pc || (SCIPisLE(scip, g->cost[edge1], pi) && SCIPisLE(scip, ttdist, pt))) )
      {
         (*nelims)++;
         *fixed += g->cost[edge1];

         if( pc )
         {
            SCIP_CALL( graph_knot_contractpc(scip, g, i, k, i) );
         }
         else
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[edge1]) );
            SCIP_CALL( graph_knot_contract(scip, g, i, k) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &minedge1);
   SCIPfreeBufferArray(scip, &term);

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
   char* blocked;

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

   voronoi_terms(scip, g, g->cost, vnoi, vbase, heap, state);

   if( nedges >= (nterms - 1) * nterms )
      maxnedges = (nterms - 1) * nterms;
   else
      maxnedges = nedges;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1, 0) );

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
   netgraph->source[0] = 0;

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
      level0(scip, g);

   /* free netgraph and  MST data structure */
   graph_path_exit(scip, netgraph);
   graph_free(scip, netgraph, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &nodesid);
   SCIPfreeBufferArray(scip, &edgeorg);
   SCIPfreeBufferArray(scip, &blocked);

   assert(graph_valid(g));
   return SCIP_OKAY;
}


/** adjacent neighbourhood reduction for the MWCSP */
SCIP_RETCODE ansReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   int k;
   int j;
   int e;
   int e2;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighbourhood of all nodes */
   for( k = 0; k < nnodes; k++ )
   {
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

      maxgrad = g->grad[k];

      /* check all neighbours of k */
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
#if 0
         if( !g->mark[k] || k == k2 )
            continue;
         j = k2;
#else
         j = g->head[e];
         e = g->oeat[e];
#endif
         /* valid candidate? */
         if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) )
         {
            for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               if( !marked[g->head[e2]] )
                  break;

            /* neighbours of j subset of those of k? */
            if( e2 == EAT_LAST )
            {
               while( g->outbeg[j] != EAT_LAST )
               {
                  e2 = g->outbeg[j];
                  (*count)++;
                  graph_edge_del(scip, g, e2, TRUE);
               }
               g->mark[j] = FALSE;
               marked[j] = FALSE;
            }
         }
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;
      marked[k] = FALSE;

      for( j = 0; j < nnodes; j++ )
         assert(marked[j] == FALSE);
   }

   return SCIP_OKAY;
}

#if 1

/** advanced connected neighbourhood subset reduction test for the MWCSP */
SCIP_RETCODE cnsAdvReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   SCIP_Real kprize;
   int* neighbarr;
   int* neighbarr2;
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int run;
   int nn;
   int nn2;
   int k2grad;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   k2grad = 0;
   *count = 0;
   kprize = 0.0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr, CNSNN + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr2, CNSNN + 1) );

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = 0;

   /* first run: consider node plus adjacent terminals; second run: consider the same plus an additional (non-positive) vertex  */
   for( run = 0; run < 2; run++ )
   {
      /* check neighbourhood of all nodes */
      for( k = 0; k < nnodes; k++ )
      {
         if( (!(g->mark[k])) || (g->grad[k] < 2) )
            continue;

         if( run == 1 && Is_term(g->term[k]) )
            continue;
         nn = 0;
         k2 = UNKNOWN;
         nn2 = 0;
         kprize = g->prize[k];
         maxgrad = g->grad[k];

         /* mark adjacent vertices and k */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];

            if( !g->mark[j] )
               continue;

            if( SCIPisGE(scip, g->prize[j], 0.0) && nn < CNSNN - 1 )
            {
               neighbarr[nn++] = j;
               marked[j] = 3;
            }
            else if( (run == 1) &&
               ((SCIPisGT(scip, g->prize[j], kprize) && nn2 < CNSNN) || (SCIPisGE(scip, g->prize[j], kprize) && j > k && nn2 < 3))  )
            {
               neighbarr2[nn2++] = j;
               marked[j] = 3;
            }
            else
            {
               marked[j] = 2;
            }
         }


         for( l = 0; l < nn; l++ )
         {
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            {
               j = g->head[e];
               if( !g->mark[j] )
                  continue;
               if( run == 1 && SCIPisGT(scip, g->prize[j], kprize) && nn2 < CNSNN )
               {
                  neighbarr2[nn2++] = j;
                  marked[j] = 3;
               }
               else if( marked[j] == 0 )
               {
                  marked[j] = 2;
               }
            }
            maxgrad += g->grad[neighbarr[l]];
         }


         marked[k] = 3;

         if( Is_term(g->term[k]) )
            min = 0.0;
         else
            min = g->prize[k];

         for( l = 0; l < nn2 + 1 - run; l++ )
         {
            if( run == 1 )
            {
               k2 = neighbarr2[l];
               for( e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e] )
                  if( marked[g->head[e]] < 2 && g->mark[g->head[e]] )
                     marked[g->head[e]] = 1;
               min += g->prize[k2];
               k2grad = g->grad[k2];
               maxgrad += k2grad;
               assert(SCIPisLE(scip, g->prize[k2], 0.0));

            }

            /* traverse all vertices */
            for( j = 0; j < nnodes; j++ )
            {
               /* vertex part of the current connected subset? Or terminal? Or belonging to the extension of the graph? */
               if( marked[j] == 3 || Is_term(g->term[j]) || !(g->mark[j]) )
                  continue;

               /* valid candidate? */
               if( g->grad[j] <= maxgrad && SCIPisLE(scip, g->prize[j], min) )
               {
                  for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                     if( marked[g->head[e2]] == 0 )
                        break;

                  /* neighbours of j subset of those of k? */
                  if( e2 == EAT_LAST )
                  {
                     while( g->outbeg[j] != EAT_LAST )
                     {
                        e2 = g->outbeg[j];
                        (*count)++;
                        graph_edge_del(scip, g, e2, TRUE);
                     }
                     g->mark[j] = FALSE;
                     marked[j] = 0;

                  }
               }
            }
            if( run == 1 )
            {
               for( e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e] )
                  if( marked[g->head[e]] == 1 && g->mark[g->head[e]] )
                     marked[g->head[e]] = 0;
               min -= g->prize[k2];
               maxgrad -= k2grad;
            }
         }
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            marked[g->head[e]] = 0;

         }

         for( l = 0; l < nn; l++ )
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = 0;


         marked[k] = 0;

         for( l = 0; l < nnodes; l++ )
            assert(marked[l] == 0);

      }
   }

   SCIPfreeBufferArray(scip, &neighbarr2);
   SCIPfreeBufferArray(scip, &neighbarr);
   return SCIP_OKAY;
}
#endif

/** advanced adjacent neighbourhood reduction for the MWCSP */
SCIP_RETCODE ansadvReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   int* neighbarr;
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int nn;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr, CNSNN) );

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighbourhood of all nodes */
   for( k = 0; k < nnodes; k++ )
   {
      if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
         continue;

      maxgrad = 0;
      nn = 0;

      /* mark adjacent vertices and k */
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         marked[j] = TRUE;
         if( SCIPisGT(scip, g->prize[j], 0.0) && nn < CNSNN )
            neighbarr[nn++] = j;
      }

      marked[k] = TRUE;

      maxgrad = g->grad[k];
      for( l = 0; l < nn; l++ )
      {
         for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = TRUE;
         maxgrad += g->grad[neighbarr[l]];
      }

      assert(SCIPisLE(scip, g->prize[k], 0.0));

      min = g->prize[k];

      /* check all neighbours of k */
#if 1
      e = g->outbeg[k];
      while( e != EAT_LAST )
#else
         for( j = 0; j < nnodes; j++ )
#endif
         {
#if 0
            if( j == k )
               continue;
#else
            j = g->head[e];
            e = g->oeat[e];
#endif
            /* valid candidate? */
            if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) )
            {
               for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  if( !marked[g->head[e2]] )
                     break;

               /* neighbours of j subset of those of k? */
               if( e2 == EAT_LAST )
               {
                  while( g->outbeg[j] != EAT_LAST )
                  {
                     e2 = g->outbeg[j];
                     (*count)++;
                     graph_edge_del(scip, g, e2, TRUE);
                  }
                  g->mark[j] = FALSE;
                  marked[j] = FALSE;
               }
            }
         }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      for( l = 0; l < nn; l++ )
         for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

      marked[k] = FALSE;

      for( k2 = 0; k2 < nnodes; k2++ )
         assert(marked[k2] == FALSE);
   }

   SCIPfreeBufferArray(scip, &neighbarr);
   return SCIP_OKAY;
}


/** alternative advanced adjacent neighbourhood reduction for the MWCSP */
SCIP_RETCODE ansadv2Reduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   SCIP_Real maxprize;
   int* neighbarr;
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int run;
   int nn;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr, CNSNN + 1) );

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   for( run = 0; run < 2; run++ )
   {
      /* check neighbourhood of all nodes */
      for( k = 0; k < nnodes; k++ )
      {
         if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
            continue;

         nn = 0;
         k2 = UNKNOWN;
         maxgrad = g->grad[k];
         maxprize = g->prize[k];

         /* mark adjacent vertices and k */
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            marked[j] = TRUE;
            if( SCIPisGT(scip, g->prize[j], 0.0) && nn < CNSNN - 1 )
            {
               neighbarr[nn++] = j;
            }
            else if( SCIPisGT(scip, g->prize[j], maxprize) && g->mark[j] )
            {
               maxprize = g->prize[j];
               k2 = j;
            }
         }

         marked[k] = TRUE;

         if( run == 0 && k2 != UNKNOWN )
            neighbarr[nn++] = k2;

         for( l = 0; l < nn; l++ )
         {
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            {
               j = g->head[e];
               if( run == 1 && g->mark[j] && !Is_term(g->term[j]) && SCIPisGT(scip, g->prize[j], maxprize) )
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
            min += g->prize[k2];

         /* check all neighbours of k */
         e = g->outbeg[k];
         while( e != EAT_LAST )
         {
#if 0
            if( j == k )
               continue;
#else
            j = g->head[e];
            e = g->oeat[e];
#endif
            /* valid candidate? */
            if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) && k2 != j )
            {
               for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  if( !marked[g->head[e2]] )
                     break;

               /* neighbours of j subset of those of k? */
               if( e2 == EAT_LAST )
               {
                  while( g->outbeg[j] != EAT_LAST )
                  {
                     e2 = g->outbeg[j];
                     (*count)++;
                     graph_edge_del(scip, g, e2, TRUE);
                  }
                  g->mark[j] = FALSE;
                  marked[j] = FALSE;
               }
            }
         }

         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

         for( l = 0; l < nn; l++ )
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = FALSE;

         marked[k] = FALSE;

         for( k2 = 0; k2 < nnodes; k2++ )
            assert(marked[k2] == FALSE);
      }
   }

   SCIPfreeBufferArray(scip, &neighbarr);
   return SCIP_OKAY;
}


/** non-positive vertex reduction for the MWCSP */
SCIP_RETCODE npvReduction(
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
   PATH* mst;
   SCIP_Real prize;
   SCIP_Real sdist0;
   SCIP_Real sdist1;
   SCIP_Real sdist2;
   int* adjverts;
   int i;
   int k;
   int k2;
   int s;
   int l;
   int e;
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

   SCIP_CALL( SCIPallocBufferArray(scip, &adjverts, 5) );

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


   /* --- NPV3 test --- */

   /* try to eliminate non-positive vertices of degree 3 */
   for( i = 0; i < nnodes; i++ )
   {
      assert(g->grad[i] >= 0);
      /* only non-positive vertices of degree 3 */
      if( !g->mark[i] || g->grad[i] != 3 || Is_term(g->term[i]) )
         continue;

      k = 0;
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjverts[k++] = g->head[e];

      assert(k == 3);

      g->mark[i] = FALSE;
      prize = g->prize[i];
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist0, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[0], adjverts[1], limit, FALSE, TRUE) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist1, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[1], adjverts[2], limit, FALSE, TRUE) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist2, -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[2], adjverts[0], limit, FALSE, TRUE) );

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
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, 5) );
   SCIP_CALL( graph_init(scip, &auxg, 5, 40, 1, 0) );

   for( k = 0; k < 4; k++ )
      graph_knot_add(auxg, -1);
   for( k = 0; k < 4; k++ )
      for( k2 = k + 1; k2 < 4; k2++ )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 4 */
   for( i = 0; i < nnodes; i++ )
   {
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
            k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sdist0), -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
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
         for( l = 1; l < 4; l++ )
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
                  for( l = 1; l < 4; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 4; l++ )
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
   for( k = 0; k < 4; k++ )
      graph_edge_add(scip, auxg, k, 4, 1.0, 1.0);
   graph_path_exit(scip, auxg);
   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 5 */
   for( i = 0; i < nnodes; i++ )
   {
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
            k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sdist0), -prize, heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
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
         for( l = 1; l < 5; l++ )
            sdist0 += mst[l].dist;

         if( SCIPisLE(scip, prize, -sdist0) )
         {
            for( k = 0; k < 5; k++ )
            {
               auxg->mark[k] = FALSE;
               sdist0 = 0.0;
               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( l = 1; l < 5; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 5; l++ )
                     sdist0 += mst[l].dist;
               }
               k2 = UNKNOWN;
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
                        for( l = 1; l < 5; l++ )
                           if( auxg->mark[l] )
                              sdist0 += mst[l].dist;
                     }
                     else
                     {
                        if( k != 1 && k2 != 1 )
                           s = 1;
                        else
                           s = 2;
                        graph_path_exec(scip, auxg, MST_MODE, s, auxg->cost, mst);
                        for( l = 0; l < 5; l++ )
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
   graph_free(scip, auxg, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &adjverts);

   return SCIP_OKAY;
}


/** chain reduction (NPV_2) for the MWCSP */
SCIP_RETCODE chain2Reduction(
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
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist, -(g->prize[i]), heap, statetail, statehead, memlbltail, memlblhead, i1, i2, limit, FALSE, TRUE) );
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
SCIP_RETCODE nnpReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  mem,                /**< nodes array */
   int*                  marked,             /**< nodes array */
   int*                  visited,            /**< nodes array */
   int*                  count,              /**< pointer to number of reductions */
   int                   maxniter,           /**< max number of edges to check */
   char*                 positive            /**< nodes array */
   )
{
   SCIP_QUEUE* queue;
   int* pnode;
   int i;
   int j;
   int k;
   int e;
   int e2;
   int erev;
   int enext;
   int nnodes;
   int nvisited1;
   int nvisited2;
   int iterations;
   SCIP_Bool success;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(fixed     != NULL);
   assert(count     != NULL);
   assert(mem       != NULL);
   assert(positive  != NULL);
   assert(visited   != NULL);
   assert(marked    != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   /* unmark all nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && SCIPisGE(scip, g->prize[i], 0.0) )
         positive[i] = TRUE;
      else
         positive[i] = FALSE;
      marked[i] = FALSE;
      visited[i] = FALSE;
   }

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );

   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;

      /* mark adjacent vertices and i*/
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         j = g->head[e];
         enext = g->oeat[e];
         if( g->mark[j] )
         {
            nvisited1 = 0;
            iterations = 0;
            success = FALSE;

            assert(SCIPqueueIsEmpty(queue));

            SCIP_CALL( SCIPqueueInsert(queue, &g->head[g->inpbeg[i]]) );
            marked[i] = TRUE;
            mem[nvisited1++] = i;

            /* BFS */
            while( !SCIPqueueIsEmpty(queue) )
            {
               pnode = SCIPqueueRemove(queue);
               for( e2 = g->outbeg[*pnode]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               {
                  if( e2 == e )
                     continue;
                  k = g->head[e2];
                  if( iterations++ >= maxniter || k == j )
                  {
                     if( k == j )
                        success = TRUE;
                     SCIPqueueClear(queue);
                     break;
                  }

                  if( positive[k] && !marked[k] )
                  {
                     marked[k] = TRUE;
                     assert(nvisited1 < nnodes);
                     mem[nvisited1++] = k;
                     SCIP_CALL( SCIPqueueInsert(queue, &g->head[e2]) );
                  }
               }
            }

            nvisited2 = 0;
            /* vertex j not reached yet? */
            if( !success )
            {
               assert(SCIPqueueIsEmpty(queue));
               SCIP_CALL( SCIPqueueInsert(queue, &j) );
               visited[j] = TRUE;
               assert(nvisited2 + nvisited1 < nnodes);
               mem[nvisited1 + nvisited2++] = j;
               erev = flipedge(e);

               iterations = 0;
               while( !SCIPqueueIsEmpty(queue) )
               {
                  pnode = (SCIPqueueRemove(queue));
                  for( e2 = g->outbeg[*pnode]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  {
                     if( e2 == erev )
                        continue;
                     k = g->head[e2];

                     if( iterations++ >= maxniter || marked[k] )
                     {
                        if( marked[k] )
                           success = TRUE;
                        SCIPqueueClear(queue);
                        break;
                     }

                     if( positive[k] && !visited[k] )
                     {
                        visited[k] = TRUE;
                        assert(nvisited2 + nvisited1 < nnodes);
                        mem[nvisited1 + nvisited2++] = k;
                        SCIP_CALL( SCIPqueueInsert(queue, &g->head[e2]) );
                     }
                  }
               }
            }

            if( success )
            {
               (*count)++;
               graph_edge_del(scip, g, e, TRUE);
            }
            for( j = 0; j < nvisited1; j++ )
               marked[mem[j]] = FALSE;

            for( j = nvisited1; j < nvisited1 + nvisited2; j++ )
               visited[mem[j]] = FALSE;

            for( j = 0; j < nnodes; j++ )
            {
               assert(marked[j] == FALSE);
               assert(visited[j] == FALSE);
            }
         }
         e = enext;
      }
   }
   SCIPqueueFree(&queue);

   return SCIP_OKAY;
}
