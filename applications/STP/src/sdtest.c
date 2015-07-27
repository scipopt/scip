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

/**@file   sdtest.c
 * @brief  special distance based reduction tests for Steiner problems
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"

#define KNOTFREQ 100
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

static
int getcloseterms(
   SCIP* scip,
   PATH* vnoi,
   SCIP_Real* termdist,
   SCIP_Real ecost,
   int* vbase,
   int* neighbterms,
   int i,
   int nnodes,
   int dnnodes
   )
{
   int nnterms = 0;

   if( SCIPisLT(scip, vnoi[i].dist, ecost) )
   {
      neighbterms[nnterms] = vbase[i];
      termdist[nnterms++] = vnoi[i].dist;
      if( SCIPisLT(scip, vnoi[i + nnodes].dist, ecost) )
      {
         neighbterms[nnterms] = vbase[i + nnodes];
         termdist[nnterms++] = vnoi[i + nnodes].dist;
         if( SCIPisLT(scip, vnoi[i + dnnodes].dist, ecost) )
         {
            neighbterms[nnterms] = vbase[i + dnnodes];
            termdist[nnterms++] = vnoi[i + dnnodes].dist;
         }
      }
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
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  nodesid,            /**< array to map nodes in auxiliary graph to original ones */
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real* mstsdist;
   SCIP_Real* termdist1;
   SCIP_Real* termdist2;
   SCIP_Real ecost;
   SCIP_Real dist;
   int* nodesorg;
   int* neighbterms1;
   int* neighbterms2;
   int e;
   int i;
   int i2;
   int j;
   int k;
   int l;
   int tj;
   int tk;
   int ne;
   int nj;
   int nk;
   int enext;
   int nnodes;
   int nterms;
   int nedges;
   int dnnodes;
   int nnterms1;
   int nnterms2;
   int maxnedges;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nodesid != NULL);

   *nelims = 0;
   nnodes = g->knots;
   nterms = g->terms;
   nedges = g->edges;
   dnnodes = 2 * nnodes;
   maxnedges = MAX(nedges, (nterms - 1) * nterms);

   /* only one terminal left? */
   if( nterms == 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &mstsdist, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodesorg, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termdist1, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termdist2, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms1, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms2, 3) );

   /* compute nearest three terminals to all non-terminals */
   getnext3terms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);

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

         if( !Is_term(g->term[k]) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathdist[k]), nnodes) ); /*lint !e866*/
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathedge[k]), nnodes) ); /*lint !e866*/
         }
      }
      for( k = 0; k < nnodes; k++ )
      {
         if( !Is_term(g->term[k]) )
         {
            assert(pathdist[k] != NULL);
            assert(pathedge[k] != NULL);
            graph_path_execX(scip, g, k, g->cost,  pathdist[k], pathedge[k]);

            j = vbase[k];
            i = j;
            assert(Is_term(g->term[i]));
            if( !SCIPisEQ(scip, pathdist[k][i], vnoi[k].dist ) )
            {
               printf("FAILS0 %f, %f ", pathdist[k][i], vnoi[k].dist);
               printf("node terminal %d, %d \n\n", k, i);
               assert(0);
            }
            i = vbase[k + nnodes];
            assert(Is_term(g->term[i]));
            if( i != -1 && !SCIPisEQ(scip, pathdist[k][i], vnoi[k + nnodes].dist ) )
            {
               printf("FAILS1 %f, %f ", pathdist[k][i], vnoi[k + nnodes].dist);
               printf("node terminal %d, %d \n\n", k, i);
               assert(0);
            }
            i = vbase[k + 2 * nnodes];
            assert(Is_term(g->term[i]));
            if( i != -1 && !SCIPisEQ(scip, pathdist[k][i], vnoi[k + 2 *  nnodes].dist ) )
            {
               printf("FAILS2 %f, %f ", pathdist[k][i], vnoi[k + 2 *  nnodes].dist);
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
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   for( k = 0; k < nnodes; k++ )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         i = vbase[k];

         if( i != vbase[g->head[e]] )
         {
            i2 = vbase[g->head[e]];
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

   /* compute MST on netgraph */
   graph_knot_chg(netgraph, 0, 0);
   netgraph->source[0] = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
   SCIP_CALL( graph_path_init(scip, netgraph) );
   graph_path_exec(scip, netgraph, MST_MODE, 0, netgraph->cost, mst);

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
         enext = g->oeat[e];
         i2 = g->head[e];
         if( i2 < i || !g->mark[i2] )
            continue;
         ecost = g->cost[e];

         /* is i a terminal? If not, get three closest terminals of distance smaller ecost */
         if( !Is_term(g->term[i]) )
         {
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes, dnnodes);
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
            /* get three closest terminals of distance smaller ecost */
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes, dnnodes);
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
	       assert(Is_term(g->term[tk]));
	       assert(Is_term(g->term[tj]));
               assert(tk >= 0);
               if( tj == tk )
               {
                  graph_edge_del(scip, g, e, TRUE);
                  (*nelims)++;
                  break;

               }
               else
               {
                  /* get sd between (terminals) tj and tk */
                  /* TODO start from t far away*/
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
                  if( SCIPisGT(scip, ecost, dist) )
                  {
                     assert(SCIPisGT(scip, ecost, termdist1[j]));
                     assert(SCIPisGT(scip, ecost, termdist2[k]));
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               }
            } /* k < nnterms2 */
         } /* j < nnterms1 */

      } /* while( enext != EAT_LAST ) */
   }

   /* free memory*/
   graph_path_exit(scip, netgraph);
   graph_free(scip, netgraph, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &neighbterms2);
   SCIPfreeBufferArray(scip, &neighbterms1);
   SCIPfreeBufferArray(scip, &termdist2);
   SCIPfreeBufferArray(scip, &termdist1);
   SCIPfreeBufferArray(scip, &nodesorg);
   SCIPfreeBufferArray(scip, &mstsdist);

   return SCIP_OKAY;
}


/* SD test for PC */
SCIP_RETCODE sdpc_reduction(
   SCIP* scip,
   GRAPH* g,
   PATH*   vnoi,
   int* heap,
   int* state,
   int* vbase,
   int* nodesid,
   int* nodesorg,
   int* nelims
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
   int root;
   int enext;
   int nnodes;
   int nterms;
   int nedges;
   int dnnodes;
   int nnterms1;
   int nnterms2;
   int maxnedges;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nodesid != NULL);
   assert(nodesorg != NULL);

   *nelims = 0;
   nnodes = g->knots;
   nterms = g->terms;
   nedges = g->edges;
   dnnodes = 2 * nnodes;
   root = g->source[0];
   maxnedges = MAX(nedges, (nterms - 1) * nterms);

   SCIP_CALL( SCIPallocBufferArray(scip, &termdist1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termdist2, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms1, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbterms2, 4) );

   /* compute nearest three terminals to all non-terminals */
   getnext3terms(scip, g, g->cost, g->cost, vnoi, vbase, heap, state);

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
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes, dnnodes);
	 }
         if( nnterms1 == 0 )
	    continue;


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
#if 0
	    /* get three shortest terms */
	    for( ne = netgraph->outbeg[nodesid[i2]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
	    {
               if( SCIPisLT(scip, netgraph->cost[ne], ecost ) )
               {
                  j = nodesorg[netgraph->head[ne]];
                  assert(Is_term(g->term[j]));
                  assert(j != i);
                  neighbterms2[nnterms2] = j;
                  termdist2[nnterms2++] = netgraph->cost[ne];
                  if( nnterms2 >= 3 )
                     break;
               }
	    }
#endif
	 }
	 else
	 {
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes, dnnodes);

	 }
         if( nnterms2 == 0 )
	    continue;
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
		  for( e2 = netgraph->outbeg[nodesid[tj]]; e2 != EAT_LAST; e2 = netgraph->oeat[e2] )
                     if( netgraph->head[e2] == nodesid[tk] )
                        break;
#if 1
                  if( e2 != EAT_LAST && SCIPisGT(scip, ecost, netgraph->cost[e2])
                     && SCIPisGT(scip, ecost, netgraph->cost[e2] + termdist1[j] - g->prize[tj])
                     && SCIPisGT(scip, ecost, netgraph->cost[e2] + termdist2[k] - g->prize[tk])
                     && SCIPisGT(scip, ecost, netgraph->cost[e2] + termdist1[j] + termdist2[k] - g->prize[tj] - g->prize[tk]) )
                  {
		     printf("SDSP del; %d %d (%d)\n", g->tail[e], g->head[e], e);
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


/* get SD to a single edge*/
SCIP_RETCODE getSD(
   SCIP* scip,
   GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   SCIP_Real*    sdist,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int     i,
   int     i2,
   int     limit,
   SCIP_Bool pc
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
   sdpaths(scip, g, pathtail, g->cost, heap, statetail, memlbltail, &nlbltail, i, i2, limit);

   /* test whether edge e can be eliminated */
   sdpaths(scip, g, pathhead, g->cost, heap, statehead, memlblhead, &nlblhead, i2, i, limit);

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
            if( pc && SCIPisLT(scip, dist, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
               dist = pathhead[l].dist + pathtail[l].dist - g->prize[l];
            if( SCIPisGT(scip, sd, dist) )
               sd = dist;
         }
         else
         {
            if( SCIPisGT(scip, sd, pathhead[l].dist + pathtail[l].dist) )
               sd = pathhead[l].dist + pathtail[l].dist;
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


   /* compare restriced sd with edge cost (if existing) */
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( g->head[e] == i2 && SCIPisGT(scip, sd, g->cost[e]) )
         sd = g->cost[e];

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
         sdpaths(scip, g, pathtail, g->cost, heap, statetail, memlbltail, &nlbltail, i, i2, limit);

         /* start limited dijkstra from i2, marking all reached vertices */
	 sdpaths(scip, g, pathhead, g->cost, heap, statehead, memlblhead, &nlblhead, i2, i, limit);

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
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real c1;
   SCIP_Real c2;
   SCIP_Real c3;
   SCIP_Real s1;
   SCIP_Real s2;
   SCIP_Real s3;
   SCIP_Real c123;
   int    i;
   int    k;
   int    k1;
   int    k2;
   int    k3;
   int    e1;
   int    e2;
   int    e3;
   int    nnodes;

   SCIP_Bool pc;

   SCIPdebugMessage("BD3-Reduction: ");

   assert(g != NULL);
   assert(heap  != NULL);
   assert(nelims != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 3) );

   *nelims = 0;
   nnodes = g->knots;
   pc = g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   for( i = 0; i < 3; i++ )
   {
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
      if( g->grad[i] != 3 || Is_term(g->term[i]) )
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

      c123 = c1 + c2 + c3;

      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &s1, heap, statetail, statehead, memlbltail, memlblhead, k1, k2, limit, pc) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &s2, heap, statetail, statehead, memlbltail, memlblhead, k2, k3, limit, pc) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &s3, heap, statetail, statehead, memlbltail, memlblhead, k3, k1, limit, pc) );

      if( SCIPisGT(scip, s1 + s2, c123)
         &&  SCIPisGT(scip, s1 + s3, c123)
         &&  SCIPisGT(scip, s2 + s3, c123) )
         continue;

      /* save ancestors */
      for( k = 0; k < 3; k++ )
      {
	 SCIPintListNodeFree(scip, &(ancestors[k]));
	 SCIPintListNodeFree(scip, &(revancestors[k]));
      }

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[0]), g->ancestors[e1]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[0]), g->ancestors[Edge_anti(e1)]) );

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[1]), g->ancestors[e2]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[1]), g->ancestors[Edge_anti(e2)]) );

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[2]), g->ancestors[e3]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[2]), g->ancestors[Edge_anti(e3)]) );

      (*nelims)++;

      if( SCIPisLT(scip, s1, c1 + c2) )
         graph_edge_del(scip, g, e1, TRUE);
      else
	 SCIP_CALL( graph_edge_reinsert(scip, g, e1, k1, k2, c1 + c2, ancestors[0], ancestors[1], revancestors[0], revancestors[1]) );

      if( SCIPisLT(scip, s2, c2 + c3) )
         graph_edge_del(scip, g, e2, TRUE);
      else
         SCIP_CALL( graph_edge_reinsert(scip, g, e2, k2, k3, c2 + c3, ancestors[1], ancestors[2], revancestors[1], revancestors[2]) );
#if 0
      {
	 n1 = graph_edge_redirect(g, e2, k2, k3, c2 + c3);
	 if( n1 >= 0 )
	 {
            SCIPintListNodeFree(scip, &(g->ancestors[n1]));
            SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
            SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors[1]) );
            SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors[2]) );

            SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors[1]) );
            SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors[2]) );
	 }
      }
#endif
      if( SCIPisLT(scip, s3, c1 + c3) )
         graph_edge_del(scip, g, e3, TRUE);
      else
	 SCIP_CALL( graph_edge_reinsert(scip, g, e3, k3, k1, c3 + c1, ancestors[2], ancestors[0], revancestors[2], revancestors[0]) );
#if 0
      {
         n1 = graph_edge_redirect(g, e3, k3, k1, c3 + c1);
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

   for( k = 0; k < 3; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);
   /*
     if( !pc )
     assert(graph_valid(g));
   */
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
   int* nelims
   )
{
   SCIP_QUEUE* queue;
   SCIP_Real mincost2;
   int     i;
   int     k;
   int     e;
   int     j;
   int     t;
   int     head;
   int     tail;
   int     root;
   int     nnodes;
   int     vrcount;
   int     minedge;
   int*    qnode;
   int*    vrnodes;
   char*   visited;
   char*   forbidden;
   SCIP_Bool pc;
   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   *nelims = 0;
   nnodes = g->knots;
   root = g->source[0];
   pc = (g->stp_type == STP_PRIZE_COLLECTING) || (g->stp_type == STP_ROOTED_PRIZE_COLLECTING);

   SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vrnodes, nnodes) );
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
	 mincost2 = FARAWAY;

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
                  if( minedge == UNKNOWN )
                  {
                     minedge = e;
                  }
                  else if( SCIPisLE(scip, g->cost[e], g->cost[minedge]) )
                  {
                     mincost2 = g->cost[minedge];
                     minedge = e;
                  }
                  else if( SCIPisLT(scip, g->cost[e], mincost2) )
                  {
                     mincost2 = g->cost[e];
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

	 /* check whether minedge can be removed */
         if( SCIPisGE(scip, mincost2, vnoi[tail].dist + g->cost[e] + vnoi[head].dist) )
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
            (*nelims)++;
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

            if( pc )
	    {
               SCIP_CALL( graph_knot_contractpc(scip, g, j, k, i) );
	    }
            else
	    {
	       SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e]) );
               SCIP_CALL( graph_knot_contract(scip, g, j, k) );
	    }
            forbidden[vbase[j]] = TRUE;
	    forbidden[vbase[k]] = TRUE;
         }
      }
   }

   /* free memory */
   SCIPqueueFree(&queue);
   SCIPfreeBufferArray(scip, &forbidden);
   SCIPfreeBufferArray(scip, &vrnodes);
   SCIPfreeBufferArray(scip, &visited);

   return SCIP_OKAY;
}


/* NV reduction from T. Polzin's "Algorithms for the Steiner problem in networks" */
SCIP_RETCODE nv_reduction(
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

   /* compute voronoi regions and distances */
   SCIP_CALL( voronoi_dist(scip, g, g->cost, distance, vbase, minedge1, heap, state, distnode, vnoi) );

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

   SCIPfreeBufferArrayNull(scip, &distnode);
   SCIPfreeBufferArray(scip, &distance);
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

   assert( graph_valid(netgraph ) );

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
