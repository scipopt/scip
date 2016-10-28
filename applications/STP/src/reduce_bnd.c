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

/**@file   reduce_bnd.c
 * @brief  bound based reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements bound-based reduction techniques for several Steiner problems.
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
#include "scip/scip.h"
#include "grph.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "heur_local.h"
#include "misc_stp.h"
#include "prop_stp.h"
#include "probdata_stp.h"

#define DEFAULT_HEURRUNS 100                 /**< number of runs of constructive heuristic */
#define DEFAULT_DARUNS     5                 /**< number of runs for dual ascent heuristic */

#if 0
/** print graph (in undirected form) in GML format */
static
SCIP_RETCODE probdataPrintGraph(
   GRAPH*                graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   SCIP_Bool*            edgemark            /**< Array of (undirected) edges to highlight */
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "stpgraph.gml", "w");

#ifndef NDEBUG
   for( e = 0; e < graph->edges; e += 2 )
   {
      assert(graph->tail[e] == graph->head[e + 1]);
      assert(graph->tail[e + 1] == graph->head[e]);
   }
#endif

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
#if 1
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
#else
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "");
      if( n == graph->source[0] )
      {

	 SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
	 m = 1;
      }
      else if( graph->term[n] == 0 )
      {

	 SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
	 e += 1;
      }
      else
      {
         SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#666666", NULL);
      }
#endif

   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e += 2 )
   {
#if 1
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);
#endif
      if( edgemark != NULL && edgemark[e / 2] )
	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      else
         SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, NULL);
   }

   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}
#endif

/** compute starting points for constructive heuristics */
static
void compTMstarts(
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< starting points array */
   unsigned int*         seed,               /**< random seed */
   int                   runs                /**< number of runs */
   )
{
   int r;
   int k;
   int l;
   int e;
   int root;
   int nnodes;
   int nterms;
   int randval;

   assert(graph != NULL);
   assert(starts != NULL);

   root = graph->source[0];
   nnodes = graph->knots;
   nterms = graph->terms;

   r = 0;
   if( graph->mark[root] )
      starts[r++] = root;
   randval = SCIPgetRandomInt(0, nnodes - 1, seed);

   /* use non-isolated terminals as starting points for TM heuristic */
   for( k = 0; k < nnodes; k++ )
   {
      if( r >= runs || r >= nterms )
         break;

      l = (k + randval) % nnodes;
      if( Is_term(graph->term[l]) && graph->mark[l] && l != root )
         starts[r++] = l;
   }

   /* still empty slots in start array? */

   /* fill empty slots with terminal neighbours */
   for( k = 0; k < r && r < runs; k++ )
   {
      for( e = graph->outbeg[starts[k]]; e != EAT_LAST && r < runs; e = graph->oeat[e] )
      {
         l = graph->head[e];
         if( !Is_term(graph->term[l]) && graph->mark[l] )
            starts[r++] = l;
      }
   }

   /* fill empty slots randomly */
   for( k = 0; k < nnodes && r < runs; k++ )
   {
      l = (k + randval) % nnodes;
      if( !Is_term(graph->term[l]) && graph->mark[l] )
         starts[r++] = l;
   }
}

/** dual ascent based reductions */
SCIP_RETCODE da_reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   int*                  edgearrint,         /**< int edges array for internal computations or NULL */
   int*                  vbase,              /**< array for Voronoi bases */
   int*                  pathedge,           /**< array for predecessor edge on a path */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  heursources,        /**< array to store starting points of TM heuristic */
   char*                 nodearrchar,        /**< char node array for internal computations or NULL */
   int*                  nelims              /**< pointer to store number of reduced edges */
   )
{
   SCIP_HEURDATA* tmheurdata;
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real maxcost;
   SCIP_Real lpobjval;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   SCIP_Bool rpc;
   SCIP_Bool success;
   SCIP_Real hopfactor;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   int* terms;
   int* result;
   int* starts;
   int* adjvert;
   int* incedge;
   int* reinsert;
   int* termdegs;
   int i;
   int k;
   int e;
   int run;
   int etmp;
   int runs;
   int root;
   int nruns;
   int nterms;
   int nedges;
   int nnodes;
   int nfixed;
   int best_start;
   unsigned int seed;
   char* marked;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);

   seed = 0;
   root = graph->source[0];
   rpc = (graph->stp_type == STP_ROOTED_PRIZE_COLLECTING);
   nfixed = 0;
   nedges = graph->edges;
   nnodes = graph->knots;
   maxcost = 0.0;
   hopfactor = DEFAULT_HOPFACTOR;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, nedges) );

   if( !rpc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &terms, graph->terms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &termdegs, graph->terms) );
   }
   else
   {
      terms = NULL;
      termdegs = NULL;
   }

   /* allocate length-4 buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sd, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ecost, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adjvert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &reinsert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &incedge, 4) );

   for( i = 0; i < 4; i++ )
   {
      sd[i] = 0.0;
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }

   /* 1. step: compute upper bound */

   /* initialize */
   k = 0;
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( !rpc )
         graph->mark[i] = (graph->grad[i] > 0);
      if( graph->mark[i] )
      {
         k++;
         if( Is_term(graph->term[i]) )
            nterms++;
      }
   }

   /* not more than two terminals? */
   if( nterms <= 2 )
      goto TERMINATE;

   /* number of runs should not exceed number of connected vertices */
   runs = MIN(k, DEFAULT_HEURRUNS);

   /* neither PC, MW, RPC nor HC? */
   if( !rpc && heursources != NULL )
   {
      /* choose starting points for TM heuristic */
      starts = heursources;
      compTMstarts(graph, starts, &seed, runs);
   }
   else
   {
      starts = NULL;
   }

   /* get TM heuristic data */
   assert(SCIPfindHeur(scip, "TM") != NULL);
   tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   for( e = 0; e < nedges; e++ )
   {
      cost[e] = graph->cost[e];
      result[e] = UNKNOWN;
      costrev[e] = graph->cost[flipedge(e)];

      if( graph->stp_type == STP_HOP_CONS && SCIPisLT(scip, graph->cost[e], BLOCKED) && SCIPisGT(scip, graph->cost[e], maxcost) )
         maxcost = graph->cost[e];
   }

   /* RPC? Then restore transformed graph */
   if( rpc )
      SCIP_CALL( pcgraphtrans(scip, graph) );

   SCIP_CALL( SCIPheurComputeSteinerTree(scip, tmheurdata, graph, starts, &best_start, result, runs, root, cost, costrev, &hopfactor, NULL, maxcost, &success) );

   /* RPC? Then restore original graph */
   if( rpc )
   {
      SCIP_CALL( extendSteinerTreePcMw(scip, graph, vnoi, costrev, vbase, result, nodearrchar, &e) );
      SCIP_CALL( pcgraphorg(scip, graph) );
   }

   /* no feasible solution found? */
   if( !success )
      goto TERMINATE;

   /* calculate objective value of solution */
   upperbound = 0.0;
   for( e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         upperbound += graph->cost[e];

   /* 2. step: repeatedly compute lower bound and reduced costs and try to eliminate */
   for( i = 0; i < nnodes; i++ )
      if( Is_term(graph->term[i]) )
         assert(graph->grad[i] > 0);

   /* if not RPC, select roots for dual ascent */
   if( !rpc )
   {
      assert(terms != NULL);
      assert(termdegs != NULL);
      k = 0;
      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(graph->term[i]) && (graph->grad[i] > 0) )
         {
            termdegs[k] = graph->grad[i];
            terms[k++] = i;
         }
      }
      nruns = MIN(nterms, DEFAULT_DARUNS);
      SCIPsortIntInt(termdegs, terms, nterms);
   }
   else
   {
      nruns = 1;
      SCIP_CALL( pcgraphtrans(scip, graph) );
   }

   if( graph->stp_type == STP_HOP_CONS )
      nruns = 1;

   for( run = 0; run < nruns; run++ )
   {
      /* graph vanished? */
      if( graph->grad[graph->source[0]] == 0 )
	 break;

      if( !rpc && graph->stp_type != STP_HOP_CONS )
      {
         assert(terms != NULL);
         root = terms[nterms - run - 1];
      }

      SCIP_CALL( SCIPdualAscentStp(scip, graph, cost, &lpobjval, FALSE, gnodearr, edgearrint, state, root, 1, NULL, nodearrchar) );

      /* the required reduced path cost to be surpassed */
      minpathcost = upperbound - lpobjval;
      SCIPdebugMessage("da: upperbound %f, lpobjval %f\n", upperbound, lpobjval);

      for( e = 0; e < nedges; e++ )
      {
	 marked[e] = FALSE;
         costrev[e] = cost[flipedge(e)];
      }

      for( k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

      /* distance from root to all nodes */
      graph_path_execX(scip, graph, root, cost, pathdist, pathedge);

      /* no paths should go back to the root */
      for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         costrev[e] = FARAWAY;

      /* build voronoi diagram */
      getnext4terms(scip, graph, costrev, costrev, vnoi, vbase, graph->path_heap, state);

      for( k = 0; k < nnodes; k++ )
         if( !Is_term(graph->term[k]) )
            assert(vbase[k + nnodes] != root );

      /* RPC? If yes, restore original graph */
      if( rpc )
      {
         SCIP_CALL( pcgraphorg(scip, graph) );
	 graph->mark[root] = FALSE;
      }

      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
         {
            while( graph->outbeg[k] != EAT_LAST )
            {
               graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
               nfixed++;
            }
         }
         else
         {
            e = graph->outbeg[k];
            while( e != EAT_LAST )
            {
               etmp = graph->oeat[e];

	       /* for rpc not artificial terminal arcs should be deleted */
	       if( rpc && !graph->mark[graph->head[e]] )
	       {
		  e = etmp;
		  continue;
	       }

               if( SCIPisGT(scip, pathdist[k] + cost[e] + vnoi[graph->head[e]].dist, minpathcost) )
               {
                  if( marked[flipedge(e)] )
                  {
                     graph_edge_del(scip, graph, e, TRUE);
                     nfixed++;
                  }
                  else
                  {
                     marked[e] = TRUE;
                  }
               }
               e = etmp;
            }
         }
      }

      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] || Is_term(graph->term[k]) )
            continue;
         if( graph->grad[k] == 3 )
	 {

	    if( SCIPisGT(scip,pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist, minpathcost) )
	    {
               int i2;

               SCIPdebugMessage("DA eliminates 3-Vertex: %d %f min: %f \n", k, pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist,  minpathcost);
               nfixed++;
               i = 0;

               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
                  incedge[i] = e;
                  ecost[i] = graph->cost[e];
                  adjvert[i++] = graph->head[e];
                  assert(i <= 4);
               }

               /* save ancestors */
               for( i = 0; i < 3; i++ )
               {
                  SCIPintListNodeFree(scip, &(ancestors[i]));
                  SCIPintListNodeFree(scip, &(revancestors[i]));
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[i]), graph->ancestors[incedge[i]]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[i]), graph->ancestors[Edge_anti(incedge[i])]) );
               }

               for( i = 0; i < 3; i++ )
               {
                  i2 = (i + 1) % 3;
                  SCIP_CALL( graph_edge_reinsert(scip, graph, incedge[i], adjvert[i], adjvert[i2], ecost[i] + ecost[i2], ancestors[i], ancestors[i2], revancestors[i], revancestors[i2]) );
               }
	    }
	 }
	 /* @todo */
         if( graph->grad[k] == 4 )
	 {

	    if( SCIPisGT(scip,pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist, minpathcost) )
	    {
               SCIPdebugMessage("DA could eliminate 4-Vertex: %d %f min: %f \n", k, pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist,  minpathcost);
	    }
	 }
      }
      SCIPdebugMessage("deleted by da: %d \n", nfixed );

      if( rpc )
         graph->mark[root] = TRUE;

   } /* root loop */

 TERMINATE:
   *nelims = nfixed;

   /* free memory */
   for( k = 0; k < 4; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   SCIPfreeBufferArray(scip, &incedge);
   SCIPfreeBufferArray(scip, &reinsert);
   SCIPfreeBufferArray(scip, &adjvert);
   SCIPfreeBufferArray(scip, &ecost);
   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);
   SCIPfreeBufferArray(scip, &sd);
   SCIPfreeBufferArrayNull(scip, &termdegs);
   SCIPfreeBufferArrayNull(scip, &terms);
   SCIPfreeBufferArray(scip, &marked);
   SCIPfreeBufferArray(scip, &result);

   return SCIP_OKAY;
}


/** dual ascent based reductions for PCSPG and MWCSP */
SCIP_RETCODE daPc_reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure array */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   int*                  vbase,              /**< Voronoi base array */
   int*                  pathedge,           /**< shorest path incoming edge array for shortest path calculations */
   int*                  edgearrint,         /**< int edges array for internal computations or NULL */
   int*                  state,              /**< int 4 * vertices array  */
   char*                 nodearrchar,        /**< char node array for internal computations */
   int*                  nelims              /**< pointer to store number of reduced edges */
   )
{
   SCIP_HEURDATA* tmheurdata;
   IDX** ancestors;
   IDX** revancestors;
   GRAPH* transgraph;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   SCIP_Real lpobjval;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   SCIP_Bool success;
   SCIP_Real offset;
   SCIP_Real hopfactor;
   int* result;
   int* adjvert;
   int* incedge;
   int* reinsert;
   int i;
   int k;
   int e;
   int run;
   int etmp;
   int runs;
   int root;
   int nruns;
   int nterms;
   int nedges;
   int nnodes;
   int nfixed;
   int best_start;
   int transnnodes;
   int transnedges;
   char* marked;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);
   assert(nodearrchar != NULL);

   root = graph->source[0];
   nfixed = 0;
   nnodes = graph->knots;
   nedges = graph->edges;
   hopfactor = DEFAULT_HOPFACTOR;
   upperbound = 0.0;

   /* allocate memory */
   if( edgearrint == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   }
   else
   {
      result = edgearrint;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, nedges + 2 * (graph->terms - 1)) );

   /* allocate length-4 buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sd, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ecost, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &adjvert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &reinsert, 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &incedge, 4) );

   for( i = 0; i < 4; i++ )
   {
      sd[i] = 0.0;
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }

   /* 1. step: compute upper bound */

   /* initialize */
   k = 0;
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( graph->mark[i] )
      {
         k++;
         if( Is_term(graph->term[i]) )
            nterms++;
      }
   }

   /* not more than two terminals? */
   if( nterms <= 2 )
      goto TERMINATE;

   /* number of runs should not exceed number of connected vertices */
   runs = MIN(k, DEFAULT_HEURRUNS);

   /* get TM heuristic data */
   assert(SCIPfindHeur(scip, "TM") != NULL);
   tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   for( e = 0; e < nedges; e++ )
   {
      result[e] = UNKNOWN;
      cost[e] = graph->cost[e];
      costrev[e] = graph->cost[flipedge(e)];
   }

   /* restore transformed graph */
   SCIP_CALL( pcgraphtrans(scip, graph) );

   /* compute Steiner tree to obtain upper bound */
   SCIP_CALL( SCIPheurComputeSteinerTree(scip, tmheurdata, graph, NULL, &best_start, result, runs, root, cost, costrev, &hopfactor, NULL, 0.0, &success) );

   SCIP_CALL( extendSteinerTreePcMw(scip, graph, vnoi, costrev, vbase, result, nodearrchar, &e) );

   /* restore oringinal graph */
   SCIP_CALL( pcgraphorg(scip, graph) );

   /* no feasible solution found? */
   if( !success )
      goto TERMINATE;

   /* calculate objective value of solution */
   for( e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         upperbound += graph->cost[e];

   /* 2. step: compute lower bound and reduced costs and try to eliminate */

   SCIP_CALL( pcgraphtrans(scip, graph) );

   /* @todo: vary */
   nruns = 1;

   for( run = 0; run < nruns; run++ )
   {
      offset = 0.0;
      SCIP_CALL( graph_PcSapCopy(scip, graph, &transgraph, &offset) );
      transnnodes = transgraph->knots;
      transnedges = transgraph->edges;

      SCIP_CALL( SCIPdualAscentStp(scip, transgraph, cost, &lpobjval, FALSE, gnodearr, edgearrint, state, root, 1, marked, nodearrchar) );

      lpobjval += offset;

      /* the required reduced path cost to be surpassed */
      minpathcost = upperbound - lpobjval;
      SCIPdebugMessage("xupperbound %f, lpobjval %f\n", upperbound, lpobjval);

      for( e = 0; e < transnedges; e++ )
      {
         costrev[e] = cost[flipedge(e)];
         marked[e] = FALSE;
      }

      for( k = 0; k < transnnodes; k++ )
         transgraph->mark[k] = (transgraph->grad[k] > 0);

      /* init data structures for shortest paths and history */
      SCIP_CALL( graph_path_init(scip, transgraph) );
      SCIP_CALL( graph_init_history(scip, transgraph ) );

      /* distance from root to all nodes */
      graph_path_execX(scip, transgraph, root, cost, pathdist, pathedge);

      for( i = 0; i < transnnodes; i++ )
         if( Is_term(transgraph->term[i]) )
            assert(SCIPisEQ(scip, pathdist[i], 0.0));

      /* no paths should go back to the root */
      for( e = transgraph->outbeg[root]; e != EAT_LAST; e = transgraph->oeat[e] )
         costrev[e] = FARAWAY;

      /* build voronoi diagram */
      voronoi_terms(scip, transgraph, costrev, vnoi, vbase, transgraph->path_heap, transgraph->path_state);

      /* restore original transgraph */
      SCIP_CALL( pcgraphorg(scip, graph) );

      /* try to eliminate vertices and edges */
      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] || Is_gterm(graph->term[k]) )
            continue;

         if( SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
         {
            while( transgraph->outbeg[k] != EAT_LAST )
            {
	       e = transgraph->outbeg[k];
	       graph_edge_del(scip, transgraph, e, FALSE);
	       graph_edge_del(scip, graph, e, TRUE);
               nfixed++;
            }
            assert(graph->outbeg[k] == EAT_LAST);
         }
         else
         {
            e = transgraph->outbeg[k];
            while( e != EAT_LAST )
            {
               etmp = transgraph->oeat[e];

               if( SCIPisGT(scip, pathdist[k] + cost[e] + vnoi[transgraph->head[e]].dist, minpathcost) )
               {
                  if( marked[flipedge(e)] )
                  {
		     graph_edge_del(scip, graph, e, TRUE);
                     graph_edge_del(scip, transgraph, e, FALSE);
                     nfixed++;
                  }
                  else
                  {
                     marked[e] = TRUE;
                  }
               }
               e = etmp;
            }
         }
      }

      SCIPdebugMessage("fixed by da: %d \n", nfixed );

      graph_path_exit(scip, transgraph);
      graph_free(scip, transgraph, TRUE);
   } /* root loop */

 TERMINATE:
   *nelims = nfixed;

   /* free memory */
   for( k = 0; k < 4; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   SCIPfreeBufferArray(scip, &incedge);
   SCIPfreeBufferArray(scip, &reinsert);
   SCIPfreeBufferArray(scip, &adjvert);
   SCIPfreeBufferArray(scip, &ecost);
   SCIPfreeBufferArray(scip, &revancestors);
   SCIPfreeBufferArray(scip, &ancestors);
   SCIPfreeBufferArray(scip, &sd);
   SCIPfreeBufferArray(scip, &marked);

   if( edgearrint == NULL )
      SCIPfreeBufferArray(scip, &result);

   return SCIP_OKAY;
}


/** bound-based reductions for the (R)PCSTP, the MWCSP and the STP */
SCIP_RETCODE bound_reduce(
   SCIP*  scip,
   GRAPH* graph,
   PATH* vnoi,
   SCIP_Real* cost,
   SCIP_Real* prize,
   SCIP_Real* radius,
   SCIP_Real* costrev,
   SCIP_Real* offset,
   SCIP_Real* upperbound,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
   )
{
   SCIP_HEURDATA* tmheurdata;
   GRAPH* adjgraph;
   PATH* mst;
   SCIP_Real  r;
   SCIP_Real  obj;
   SCIP_Real  max;
   SCIP_Real  radii;
   SCIP_Real  bound;
   SCIP_Real  tmpcost;
   SCIP_Real  mstobj;
   SCIP_Real  mstobj2;
   SCIP_Real  maxcost;
   SCIP_Real  radiim2;
#if 1
   SCIP_Real* cost3;
   SCIP_Real  radiim3;
   IDX** ancestors;
   IDX** revancestors;
   int* edges3;
   int* nodes3;
#endif
   int* perm;
   int* result;
   int* starts;
   int e;
   int k;
   int l;
   int head;
   int tail;
   int runs;
   int root;
   int etemp;
   int nterms;
   int nnodes;
   int nedges;
   int best_start;
   unsigned int seed;
   char* stnode;
   SCIP_Bool ub;
   SCIP_Bool pc;
   SCIP_Bool mw;
   SCIP_Bool pcmw;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vnoi != NULL);
   assert(cost != NULL);
   assert(radius != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);
   assert(graph->source[0] >= 0);
   assert(upperbound != NULL);

   mst = NULL;
   obj = DEFAULT_HOPFACTOR;
   seed = 0;
   perm = NULL;
   root = graph->source[0];
   nedges = graph->edges;
   nnodes = graph->knots;
   mstobj = 0.0;
   *nelims = 0;
   mstobj2 = 0.0;
   best_start = 0;
   ub = SCIPisGT(scip, *upperbound, 0.0);
   mw = (graph->stp_type == STP_MAX_NODE_WEIGHT);
   pc = (graph->stp_type == STP_ROOTED_PRIZE_COLLECTING) || (graph->stp_type == STP_PRIZE_COLLECTING);
   pcmw = pc || mw;
   success = TRUE;
#if 1
   cost3 = NULL;
   edges3 = NULL;
   nodes3 = NULL;
   ancestors = NULL;
   revancestors = NULL;
#endif

   /* no upper bound provided? */
   if( !ub )
   {
      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &stnode, nnodes) );
   }
   else
   {
      result = NULL;
      stnode = NULL;
   }

   /* initialize */
   e = 0;
   nterms = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( !ub )
      {
	 assert(stnode != NULL);
         stnode[k] = FALSE;
      }
      if( !pcmw )
         graph->mark[k] = (graph->grad[k] > 0);
      if( graph->mark[k] )
      {
         e++;
         if( Is_term(graph->term[k]) )
            nterms++;
      }
   }

   /* not more than two terminals? */
   if( nterms <= 2 )
   {
      /* free memory and return */
      SCIPfreeBufferArrayNull(scip, &stnode);
      SCIPfreeBufferArrayNull(scip, &result);
      return SCIP_OKAY;
   }

   assert(nterms == (graph->terms - ((graph->stp_type == STP_PRIZE_COLLECTING || mw)? 1 : 0)));

   runs = MIN(e, 100);

   /* neither PC, MW, RPC nor HC? */
   if( !pcmw && graph->stp_type != STP_HOP_CONS )
   {
      /* choose starting points for TM heuristic */

      SCIP_CALL( SCIPallocBufferArray(scip, &starts, nnodes) );

      compTMstarts(graph, starts, &seed, runs);
#if 0
      int randval;

      r = 0;
      if( graph->mark[root] )
         starts[r++] = root;
      randval = SCIPgetRandomInt(0, nnodes - 1, &seed);

      /* use non-isolated terminals as starting points for TM heuristic */
      for( k = 0; k < nnodes; k++ )
      {
         if( r >= runs || r >= nterms )
            break;

         l = (k + randval) % nnodes;
         if( Is_term(graph->term[l]) && graph->mark[l] && l != root )
            starts[r++] = l;
      }

      /* still empty slots in start array? */

      /* fill empty slots with terminal neighbours */
      for( k = 0; k < r && r < runs; k++ )
      {
         for( e = graph->outbeg[starts[k]]; e != EAT_LAST && r < runs; e = graph->oeat[e] )
         {
            l = graph->head[e];
            if( !Is_term(graph->term[l]) && graph->mark[l] )
               starts[r++] = l;
         }
      }

      /* fill empty slots randomly */
      for( k = 0; k < nnodes && r < runs; k++ )
      {
         l = (k + randval) % nnodes;
         if( !Is_term(graph->term[l]) && graph->mark[l] )
            starts[r++] = l;
      }
#endif
   }
   else
   {
      starts = NULL;
   }

   /* initialise cost and costrev array */
   maxcost = 0.0;
   for( e = 0; e < nedges; e++ )
   {
      if( !ub )
      {
	 assert(result != NULL);
         result[e] = UNKNOWN;
      }
      cost[e] = graph->cost[e];
      costrev[e] = graph->cost[flipedge(e)];

      assert(SCIPisGE(scip, cost[e], 0.0));

      if( graph->stp_type == STP_HOP_CONS && SCIPisLT(scip, graph->cost[e], BLOCKED) && SCIPisGT(scip, graph->cost[e], maxcost) )
         maxcost = graph->cost[e];
   }

   /* init auxiliary graph */
   if( !mw )
   {
      SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1, 0) );
   }
   else
   {
      adjgraph = NULL;
   }

   /* build voronoi regions, concomitantly building adjgraph and computing radii values*/
   SCIP_CALL( voronoi_radius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

   /* get 2nd next terminals to all non-terminal nodes */
   get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* get 3th next terminals to all non-terminal nodes */
   get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* for (rooted) prize collecting get 4th next terminals to all non-terminal nodes */
   if( pc )
      get4next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* no MWCS problem? */
   if( !mw )
   {
      assert(adjgraph != NULL);
      graph_knot_chg(adjgraph, 0, 0);
      adjgraph->source[0] = 0;

      /* compute MST on adjgraph */
      SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
      SCIP_CALL( graph_path_init(scip, adjgraph) );
      graph_path_exec(scip, adjgraph, MST_MODE, 0, adjgraph->cost, mst);

      max = -1.0;
      r = -1.0;
      for( k = 1; k < nterms; k++ )
      {
         assert(adjgraph->path_state[k] == CONNECT);
         e = mst[k].edge;
         assert(e >= 0);
         tmpcost = adjgraph->cost[e];
         mstobj += tmpcost;
         if( SCIPisGT(scip, tmpcost, max) )
            max = tmpcost;
	 else if( SCIPisGT(scip, tmpcost, r) )
	    r = tmpcost;
      }
      mstobj -= max;
      mstobj2 = mstobj - r;
   }

   /* for (rooted) prize collecting an maximum weight problems: correct radius values */
   if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
   {
      assert(graph->mark[graph->source[0]]);
      for( k = 0; k < nnodes; k++ )
      {
         if( !Is_term(graph->term[k]) || !graph->mark[k] )
            continue;

         if( SCIPisGT(scip, radius[k], prize[k]) && k != graph->source[0] )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_PRIZE_COLLECTING )
   {
      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( Is_term(graph->term[k]) && SCIPisGT(scip, radius[k], prize[k])  )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_MAX_NODE_WEIGHT )
   {
      max = 0.0;
      for( k = 0; k < nnodes; k++ )
      {
         if( !Is_term(graph->term[k]) || !graph->mark[k] )
            continue;

         assert(SCIPisGE(scip, prize[k], 0.0));
         if( SCIPisGT(scip, prize[k], max) )
            max = prize[k];

         if( SCIPisGE(scip, radius[k], 0.0 )  )
	    radius[k] = 0.0;
	 else
	    radius[k] = -radius[k];
      }
   }

   /* sort radius values */
   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );
      for( k = 0; k < nnodes; k++ )
         perm[k] = k;
      SCIPsortRealInt(radius, perm, nnodes);
   }
   else
   {
      SCIPsortReal(radius, nnodes);
   }

   radiim2 = 0.0;

   /* sum all but two radius values of highest/lowest value */
   if( mw )
   {
      for( k = 2; k < nterms; k++ )
      {
         assert( SCIPisGT(scip, FARAWAY, radius[k]) );
         radiim2 += radius[k];
      }
   }
   else
   {
      for( k = 0; k < nterms - 2; k++ )
      {
         assert( SCIPisGT(scip, FARAWAY, radius[k]) );
         radiim2 += radius[k];
      }
   }
   radii = radiim2 + radius[nterms - 2] + radius[nterms - 1];
#if 1
   if( nterms >= 3 )
      radiim3 = radiim2 - radius[nterms - 3];
   else
      radiim3 = 0;
#endif
   /* get TM heuristic data */
   assert(SCIPfindHeur(scip, "TM") != NULL);
   tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   /* no upper bound available? */
   if( !ub )
   {
      assert(stnode != NULL);
      assert(result != NULL);

      /* PC or RPC? Then restore transformed graph */
      if( pcmw )
         SCIP_CALL( pcgraphtrans(scip, graph) );

      SCIP_CALL( SCIPheurComputeSteinerTree(scip, tmheurdata, graph, starts, &best_start, result, runs, root, cost, costrev, &obj, NULL, maxcost, &success) );

      /* PC or RPC? Then restore oringinal graph */
      if( pcmw )
         SCIP_CALL( pcgraphorg(scip, graph) );

      /* no feasible solution found? */
      if( !success )
      {
         /* free memory and return */
         if( !mw )
         {
            graph_path_exit(scip, adjgraph);
            graph_free(scip, adjgraph, TRUE);
         }
         SCIPfreeBufferArrayNull(scip, &mst);
         SCIPfreeBufferArrayNull(scip, &starts);
         SCIPfreeBufferArray(scip, &result);
         SCIPfreeBufferArray(scip, &stnode);
         return SCIP_OKAY;
      }

      /* calculate objective value of solution */
      obj = 0.0;
      for( e = 0; e < nedges; e++ )
      {
         if( result[e] == CONNECT )
         {
            head = graph->head[e];
            if( mw )
            {
               if( graph->mark[head] )
               {
                  assert(stnode[head] == FALSE);
                  obj += prize[head];
               }
            }
            else
            {
               obj += graph->cost[e];
               stnode[head] = TRUE;
               stnode[graph->tail[e]] = TRUE;
            }
         }
      }
      *upperbound = obj + *offset;
   }
   else
   {
      obj = *upperbound - *offset;
      assert(SCIPisGE(scip, obj, 0.0));
   }
#if 0
   printf("radiim2: %f \n", radiim2);
   printf("mstobj:  %f \n", mstobj);
   printf("totalobj: %f \n", obj);
#endif
   if( SCIPisGT(scip, radiim2, mstobj) )
      bound = radiim2;
   else
      bound = mstobj;

   /* traverse all node, try to eliminate each node or incident edges */
   for( k = 0; k < nnodes; k++ )
   {
      if( (!graph->mark[k] && (pcmw)) || graph->grad[k] == 0 )
         continue;

      if( mw )
      {
         tmpcost = -vnoi[k].dist - vnoi[k + nnodes].dist + bound - graph->prize[k];

	 if( !Is_term(graph->term[k]) && (SCIPisLT(scip, tmpcost, obj)) )
         {
	    while( graph->outbeg[k] != EAT_LAST )
	    {
               (*nelims)++;
               graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
	    }
         }
         else
         {
            e = graph->outbeg[k];
            while( e != EAT_LAST )
            {
               etemp = graph->oeat[e];
               tail = graph->tail[e];
               head = graph->head[e];
               if( !graph->mark[head] )
               {
                  e = etemp;
                  continue;
               }
               tmpcost = bound - graph->prize[k];
               if( vbase[tail] != vbase[head] )
               {
                  tmpcost -= vnoi[head].dist + vnoi[tail].dist;
               }
               else
               {
                  if( SCIPisGT(scip, vnoi[tail].dist + vnoi[head + nnodes].dist, vnoi[tail + nnodes].dist + vnoi[head].dist) )
                     tmpcost -= vnoi[tail + nnodes].dist + vnoi[head].dist;
                  else
                     tmpcost -= vnoi[tail].dist + vnoi[head + nnodes].dist;
               }
               if( (SCIPisLT(scip, tmpcost, obj)) )
               {
                  (*nelims)++;
                  graph_edge_del(scip, graph, e, TRUE);
               }
               e = etemp;
            }
         }
	 continue;
      }

      tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + bound;

      /* can node k be deleted? */
      if( !Is_term(graph->term[k]) && (SCIPisGT(scip, tmpcost, obj)
            || (((stnode != NULL) ? !stnode[k] : 0) && SCIPisGE(scip, tmpcost, obj))) )
      {
         SCIPdebugMessage("delete vertex: %d of degree: %d\n", k, graph->grad[k]);
         /* delete all incident edges */
         while( graph->outbeg[k] != EAT_LAST )
         {
            e = graph->outbeg[k];
            (*nelims)++;
            assert(!pc || graph->tail[e] != root);
            assert(!pc || graph->mark[graph->head[e]]);
            assert(!Is_pterm(graph->term[graph->head[e]]));
            assert(!Is_pterm(graph->term[graph->tail[e]]));

            graph_edge_del(scip, graph, e, TRUE);
         }
      }
      else
      {
         e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            etemp = graph->oeat[e];
            tail = graph->tail[e];
            head = graph->head[e];
            tmpcost = graph->cost[e] + bound;

            if( vbase[tail] != vbase[head] )
            {
               tmpcost += vnoi[head].dist + vnoi[tail].dist;
            }
            else
            {
               if( SCIPisGT(scip, vnoi[tail].dist + vnoi[head + nnodes].dist, vnoi[tail + nnodes].dist + vnoi[head].dist) )
                  tmpcost += vnoi[tail + nnodes].dist + vnoi[head].dist;
               else
                  tmpcost += vnoi[tail].dist + vnoi[head + nnodes].dist;
               assert(SCIPisGE(scip, tmpcost, vnoi[head].dist + vnoi[tail].dist + graph->cost[e] + bound));
            }
            /* can edge e or arc e be deleted? */
            if( (SCIPisGT(scip, tmpcost, obj) || (((result != NULL) ? (result[e] != CONNECT) : 0) && result[flipedge(e)] != CONNECT && SCIPisGE(scip, tmpcost, obj)))
               && SCIPisLT(scip, graph->cost[e], FARAWAY) && (!(pc || mw) || graph->mark[head]) )
            {
               SCIPdebugMessage("delete edge: %d->%d \n", graph->tail[e], graph->head[e]);
               if( graph->stp_type == STP_HOP_CONS && SCIPisGT(scip, graph->cost[e], graph->cost[flipedge(e)]) )
               {
                  graph->cost[e] = FARAWAY;
                  (*nelims)++;
               }
               else
               {
                  assert(!Is_pterm(graph->term[head]));
                  assert(!Is_pterm(graph->term[tail]));
                  graph_edge_del(scip, graph, e, TRUE);
                  (*nelims)++;
               }
            }
            e = etemp;
         }
      }
   }
#if 1
   /* traverse all node, try to eliminate 3 degree nodes */
   if( !mw )
   {
      for( k = 0; k < nnodes; k++ )
      {
         if( (!graph->mark[k] && pc) || graph->grad[k] == 0 )
            continue;

         if( graph->grad[k] == 3 && !Is_term(graph->term[k]) )
         {
            tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist + radiim3;
            if( SCIPisGT(scip, tmpcost, obj) )
            {
               /* first 3-node elimination? */
               if( ancestors == NULL )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &cost3, 3) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &edges3, 3) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &nodes3, 3) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 3) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 3) );
                  for( l = 0; l < 3; l++ )
                  {
                     ancestors[l] = NULL;
                     revancestors[l] = NULL;
                  }
               }

               assert(cost3 != NULL);
               assert(edges3 != NULL);
               assert(nodes3 != NULL);
               assert(ancestors != NULL);
               assert(revancestors != NULL);
               SCIPdebugMessage("eliminated 3 knot %d\n", k);
               /* get incident edges, cost and adjacent nodes */
               l = 0;
               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
                  assert(l < 3);
                  edges3[l] = e;
                  nodes3[l] = graph->head[e];
                  cost3[l++] = graph->cost[e];
               }

               /* clear */
               for( l = 0; l < 3; l++ )
               {
                  SCIPintListNodeFree(scip, &(ancestors[l]));
                  SCIPintListNodeFree(scip, &(revancestors[l]));
               }

               /* store ancestors of incident edges */
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[0]), graph->ancestors[edges3[0]]) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[0]), graph->ancestors[Edge_anti(edges3[0])]) );

               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[1]), graph->ancestors[edges3[1]]) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[1]), graph->ancestors[Edge_anti(edges3[1])]) );

               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[2]), graph->ancestors[edges3[2]]) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[2]), graph->ancestors[Edge_anti(edges3[2])]) );

               SCIP_CALL( graph_edge_reinsert(scip, graph, edges3[0], nodes3[0], nodes3[1], cost3[0] + cost3[1], ancestors[0], ancestors[1], revancestors[0], revancestors[1]) );
               SCIP_CALL( graph_edge_reinsert(scip, graph, edges3[1], nodes3[1], nodes3[2], cost3[1] + cost3[2], ancestors[1], ancestors[2], revancestors[1], revancestors[2]) );
               SCIP_CALL( graph_edge_reinsert(scip, graph, edges3[2], nodes3[2], nodes3[0], cost3[2] + cost3[0], ancestors[2], ancestors[0], revancestors[2], revancestors[0]) );

               assert(graph->grad[k] == 0);
            }
         }
      }

      if( ancestors != NULL )
      {
         assert(revancestors != NULL);
         for( k = 0; k < 3; k++ )
         {
            SCIPintListNodeFree(scip, &(ancestors[k]));
            SCIPintListNodeFree(scip, &(revancestors[k]));
         }
         SCIPfreeBufferArray(scip, &revancestors);
         SCIPfreeBufferArray(scip, &ancestors);
         SCIPfreeBufferArray(scip, &nodes3);
         SCIPfreeBufferArray(scip, &edges3);
         SCIPfreeBufferArray(scip, &cost3);
      }
   }
#endif

   /* for(R)PC: try to eliminate terminals */
   if( pc )
   {
      getnext4tterms(scip, graph, cost, costrev, vnoi, vbase, heap, state);

      for( k = 0; k < nnodes; k++ )
      {
         /* is k a terminal other than the root? */
         if( Is_term(graph->term[k]) && graph->mark[k] && graph->grad[k] > 2 && k != graph->source[0] )
         {
            assert(perm != NULL);
            for( l = 0; l < nterms; l++ )
               if( perm[l] == k )
                  break;

            tmpcost = vnoi[k].dist + radii - radius[l];

            if( l == nterms - 1 )
               tmpcost -= radius[nterms - 2];
            else
               tmpcost -= radius[nterms - 1];


            if( SCIPisGT(scip, tmpcost, obj) )
            {
	       SCIPdebugMessage("alternative bnd elimination!!! \n\n");
               (*nelims) += deleteterm(scip, graph, k);
	       (*offset) += graph->prize[k];
            }
            else
	    {
               tmpcost += vnoi[k + nnodes].dist;
               if( l >= nterms - 2 )
                  tmpcost -= radius[nterms - 3];
               else
                  tmpcost -= radius[nterms - 2];
               if( SCIPisGT(scip, tmpcost, obj)  || SCIPisGT(scip, mstobj2 + vnoi[k].dist + vnoi[k + nnodes].dist, obj) )
               {
                  for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                     if( graph->mark[graph->head[e]] && SCIPisLT(scip, graph->cost[e], graph->prize[k]) )
                        break;

                  if( e == EAT_LAST )
                  {
                     SCIPdebugMessage("second elimination!!! prize: %f \n\n", graph->prize[k]);
                     (*offset) += graph->prize[k];
                     (*nelims) += deleteterm(scip, graph, k);
                  }
               }
	    }
         }
      }
   }

   SCIPdebugMessage("nelims (edges) in bound reduce: %d,\n", *nelims);

   /* free adjgraph */
   if( !mw )
   {
      graph_path_exit(scip, adjgraph);
      graph_free(scip, adjgraph, TRUE);
   }

   /* free memory*/
   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArrayNull(scip, &mst);
   SCIPfreeBufferArrayNull(scip, &starts);
   SCIPfreeBufferArrayNull(scip, &stnode);
   SCIPfreeBufferArrayNull(scip, &result);

   return SCIP_OKAY;
}


/** bound-based reduction test for the HCDSTP */
SCIP_RETCODE hopbound_reduce(
   SCIP*  scip,
   GRAPH* graph,
   PATH* vnoi,
   SCIP_Real* cost,
   SCIP_Real* radius,
   SCIP_Real* costrev,
   int* heap,
   int* state,
   int* vbase,
   int* nelims
   )
{
   SCIP_Real  max;
   SCIP_Real  tmpcost;
   SCIP_Real  bound;
   SCIP_Real  mstobj;
   SCIP_Real  radiim2;

   GRAPH* adjgraph;
   PATH* mst;
   int e;
   int k;
   int tail;
   int head;
   int etemp;
   int nnodes;
   int nedges;
   int nterms;
   SCIP_Real hoplimit;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vnoi != NULL);
   assert(cost != NULL);
   assert(radius != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);

   *nelims = 0;
   nterms = 0;
   nedges = graph->edges;
   nnodes = graph->knots;

   for( k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);
      if( graph->mark[k] && Is_term(graph->term[k]) )
         nterms++;
   }

   for( e = 0; e < nedges; e++ )
   {
      if( SCIPisLT(scip, graph->cost[e], FARAWAY) )
         cost[e] =  1.0;
      else
         cost[e] =  FARAWAY;
      if( SCIPisLT(scip, graph->cost[flipedge(e)], FARAWAY) )
         costrev[e] =  1.0;
      else
         costrev[e] =  FARAWAY;
   }

   /* init auxiliary graph */
   SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1, 0) );

   SCIP_CALL( voronoi_radius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

   /* get 2nd next terminals to all nodes */
   get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* compute MST on adjgraph */
   graph_knot_chg(adjgraph, 0, 0);
   adjgraph->source[0] = 0;
   assert(graph_valid(adjgraph));
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
   SCIP_CALL( graph_path_init(scip, adjgraph) );
   graph_path_exec(scip, adjgraph, MST_MODE, 0, adjgraph->cost, mst);

   max = -1;
   assert(mst[0].edge == -1);
   mstobj = 0.0;

   /* compute MST cost ...*/
   for( k = 1; k < nterms; k++ )
   {
      e = mst[k].edge;
      assert(adjgraph->path_state[k] == CONNECT);
      assert(e >= 0);
      tmpcost = adjgraph->cost[e];
      mstobj += tmpcost;
      if( SCIPisGT(scip, tmpcost, max) )
         max = tmpcost;
   }
   /* ...minus longest edge */
   mstobj -= max;

   SCIPsortReal(radius, nnodes);
   radiim2 = 0.0;

   for( e = 0; e < nterms - 2; e++ )
   {
      assert( SCIPisGT(scip, FARAWAY, radius[e]) );
      radiim2 += radius[e];
   }

   hoplimit = (SCIP_Real) graph->hoplimit;
#if 0
   printf("radiim2: %f \n", radiim2);
   printf("mstobj: %f \n", mstobj);
   printf("hoplimit: %f \n", hoplimit);
#endif
   if( SCIPisGT(scip, radiim2, mstobj) )
      bound = radiim2;
   else
      bound = mstobj;

   /* traverse all node, try to eliminate first the node and then all incident edges */
   for( k = 0; k < nnodes; k++ )
   {
      /* can node k be deleted? */
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, vnoi[k].dist + vnoi[k + nnodes].dist + bound, hoplimit) )
      {
         e = graph->outbeg[k];

         /* delete incident edges */
         while( e != EAT_LAST )
         {
            assert(e >= 0);
            (*nelims)++;
            etemp = graph->oeat[e];
            graph_edge_del(scip, graph, e, TRUE);
            e = etemp;
         }
      }
      else
      {
         e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            assert(e >= 0);
            tail = graph->tail[e];
            head = graph->head[e];
            tmpcost = 1.0 + bound;
            if( vbase[tail] != vbase[head] )
            {
               tmpcost += vnoi[head].dist + vnoi[tail].dist;
            }
            else
            {
               if( SCIPisGT(scip, vnoi[tail].dist + vnoi[head + nnodes].dist, vnoi[tail + nnodes].dist + vnoi[head].dist) )
                  tmpcost += vnoi[tail + nnodes].dist + vnoi[head].dist;
               else
                  tmpcost += vnoi[tail].dist + vnoi[head + nnodes].dist;
               assert(SCIPisGE(scip, tmpcost, vnoi[head].dist + vnoi[tail].dist + 1.0 + bound));
            }

            /* can edge e (i.e. both arc e and its reverse arc) or arc e be deleted? */
            if( SCIPisGT(scip, tmpcost, hoplimit) && SCIPisLT(scip, graph->cost[e], FARAWAY) )
            {
               etemp = graph->oeat[e];
               if( SCIPisGT(scip, graph->cost[e], graph->cost[flipedge(e)]) )
               {
                  graph->cost[e] = FARAWAY;
                  (*nelims)++;
               }
               else
               {
                  graph_edge_del(scip, graph, e, TRUE);
                  (*nelims)++;
               }
               e = etemp;
            }
            else
            {
               e = graph->oeat[e];
            }
         }
      }
   }

   SCIPdebugMessage("nelimsX (edges) in hop bound reduce: %d,\n", *nelims);

   /* free adjgraph */
   graph_path_exit(scip, adjgraph);
   graph_free(scip, adjgraph, TRUE);

   /* free memory*/
   SCIPfreeBufferArray(scip, &mst);
   assert(graph_valid(graph));

   return SCIP_OKAY;
}

/** hop bound-based reduction test for the HCDSTP */
SCIP_RETCODE hcrbound_reduce(
   SCIP*  scip,
   GRAPH* graph,
   PATH* vnoi,
   SCIP_Real* cost,
   SCIP_Real* costrev,
   SCIP_Real* pathdist,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int* pathedge
   )
{
   SCIP_Real tmpcost;
   int e;
   int k;
   int root;
   int head;
   int etemp;
   int bound;
   int nnodes;
   int nedges;
   int nterms;
   int hoplimit;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vnoi != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);

   *nelims = 0;
   nterms = 0;
   root = graph->source[0];
   nedges = graph->edges;
   nnodes = graph->knots;
   hoplimit = graph->hoplimit;

   for( k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);
      if( graph->mark[k] && Is_term(graph->term[k]) )
         nterms++;
   }
   bound = nterms - 2;
   for( e = 0; e < nedges; e++ )
   {
      if( SCIPisLT(scip, graph->cost[e], FARAWAY) )
         cost[e] = 1.0;
      else
         cost[e] = graph->cost[e];
      if( SCIPisLT(scip, graph->cost[flipedge(e)], FARAWAY) )
         costrev[e] = 1.0;
      else
         costrev[e] = graph->cost[flipedge(e)];
   }

   /* distance from root to all nodes */
   graph_path_execX(scip, graph, root, cost, pathdist, pathedge);

   /* no paths should go back to the root */
   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build voronoi diagram */
   voronoi_terms(scip, graph, costrev, vnoi, vbase, heap, state);

   /* traverse all node, try to eliminate first the node and then all incident edges */
   for( k = 0; k < nnodes; k++ )
   {
      /* can node k be deleted? */
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, vnoi[k].dist + pathdist[k] + (double) bound, (double) hoplimit) )
      {
         e = graph->outbeg[k];

         /* delete incident edges */
         while( e != EAT_LAST )
         {
            assert(e >= 0);
            (*nelims)++;
            etemp = graph->oeat[e];
            graph_edge_del(scip, graph, e, TRUE);
            e = etemp;
         }
      }
      else
      {
         e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            assert(e >= 0);
            head = graph->head[e];
            tmpcost = pathdist[k] + 1.0 + vnoi[head].dist + bound;

            etemp = graph->oeat[e];
            /* can edge e (i.e. both arc e and its reverse arc) or arc e be deleted? */
            if( SCIPisGT(scip, tmpcost, (double) hoplimit) && SCIPisLT(scip, graph->cost[e], FARAWAY) )
            {

               if( SCIPisGT(scip, FARAWAY, graph->cost[flipedge(e)]) )
               {
                  graph->cost[e] = FARAWAY;
                  (*nelims)++;
               }
               else
               {
                  graph_edge_del(scip, graph, e, TRUE);
                  (*nelims)++;
               }
            }
            e = etemp;
         }
      }
   }

   SCIPdebugMessage("eliminated (edges) in hcr bound reduce: %d,\n", *nelims);

   assert(graph_valid(graph));

   return SCIP_OKAY;
}

/* reduction method for HCSTP */
SCIP_RETCODE hcrcbound_reduce(
   SCIP*  scip,
   GRAPH* graph,
   PATH* vnoi,
   SCIP_Real* cost,
   SCIP_Real* costrev,
   SCIP_Real* pathdist,
   SCIP_Real objval,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int* pathedge,
   SCIP_Bool fix
   )
{
   SCIP_VAR** vars;
   SCIP_HEURDATA* tmheurdata;
   SCIP_Real min;
   SCIP_Real bound;
   SCIP_Real maxmin;
   SCIP_Real maxcost;
   SCIP_Real tmpcost;
   SCIP_Real hopfactor;
   int* result;
   int e;
   int k;
   int root;
   int head;
   int etemp;
   int nnodes;
   int nedges;
   int best_start;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vnoi != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);

   hopfactor = DEFAULT_HOPFACTOR;
   bound = 0.0;
   *nelims = 0;
   best_start = 0;
   success = TRUE;
   vars = NULL;
   root = graph->source[0];
   nedges = graph->edges;
   nnodes = graph->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   maxcost = 0.0;
   if( fix )
   {
      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);
      for( e = 0; e < nedges; e += 2 )
      {
         result[e] = UNKNOWN;
         result[e + 1] = UNKNOWN;

         if( SCIPvarGetUbLocal(vars[e + 1]) < 0.5 )
         {
            costrev[e] = BLOCKED;
         }
         else
         {
            costrev[e] = graph->cost[e + 1];

            if( SCIPisGT(scip, costrev[e], maxcost) && SCIPisLT(scip, costrev[e], BLOCKED) )
               maxcost = costrev[e];
         }
         cost[e + 1] = costrev[e];
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
         {
            costrev[e + 1] = BLOCKED;
         }
         else
         {
            costrev[e + 1] = graph->cost[e];

            if( SCIPisGT(scip, graph->cost[e], maxcost) && SCIPisLT(scip, costrev[e + 1], BLOCKED) )
               maxcost = graph->cost[e];
         }
         cost[e] = costrev[e + 1];
      }
   }
   else
   {
      for( e = 0; e < nedges; e++ )
      {
         result[e] = UNKNOWN;
         cost[e] = graph->cost[e];
         costrev[e] = graph->cost[flipedge(e)];
         if( SCIPisGT(scip, graph->cost[e], maxcost) )
            maxcost = graph->cost[e];
      }
   }

   maxmin = -1.0;
   for( k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);
      if( graph->mark[k] && Is_term(graph->term[k]) )
      {
         if( k != root )
         {
            min = FARAWAY;
            for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
               if( SCIPisLT(scip, cost[e], min) )
                  min = cost[e];
            assert(SCIPisGT(scip, BLOCKED, min));
            if( SCIPisGT(scip, min, maxmin) )
               maxmin = min;
            bound += min;
         }
      }
   }
   bound -= maxmin;


   /* distance from root to all nodes */
   graph_path_execX(scip, graph, root, cost, pathdist, pathedge);

   /* no paths should go back to the root */
   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build voronoi diagram */
   voronoi_terms(scip, graph, costrev, vnoi, vbase, heap, state);

   if( SCIPisLT(scip, objval, 0.0) )
   {
      /* get TM heuristic data */
      assert(SCIPfindHeur(scip, "TM") != NULL);
      tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

      /* compute UB */
      SCIP_CALL( SCIPheurComputeSteinerTree(scip, tmheurdata, graph, NULL, &best_start, result, 50, root, cost, costrev, &hopfactor, NULL, maxcost, &success) );

      objval = 0.0;
      for( e = 0; e < nedges; e++ )
         if( result[e] == CONNECT )
            objval += graph->cost[e];
   }
   else
   {
      /* objval = objval - fixed; */
      objval = SCIPgetCutoffbound(scip);
      assert(SCIPisGT(scip, objval, 0.0));
   }

   /* traverse all node, try to eliminate first the node and then all incident edges */
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) )
         continue;
      /* can node k be deleted? */
      if( SCIPisGT(scip, vnoi[k].dist + pathdist[k] + bound, objval) )
      {
         e = graph->outbeg[k];

         /* delete incident edges */
         while( e != EAT_LAST )
         {
            assert(e >= 0);

            etemp = graph->oeat[e];
            if( fix )
            {
               assert(vars != NULL);
               /* try to fix edge */
               SCIP_CALL( fixedgevar(scip, vars[e], nelims) );

               /* try to fix reversed edge */
               SCIP_CALL( fixedgevar(scip, vars[flipedge(e)], nelims) );
            }
            else
            {
               graph_edge_del(scip, graph, e, TRUE);
               (*nelims)++;
            }
            e = etemp;
         }
      }
      else
      {
         e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            assert(e >= 0);
            head = graph->head[e];
            tmpcost = pathdist[k] + graph->cost[e] + vnoi[head].dist + bound;

            etemp = graph->oeat[e];
            /* can edge e (i.e. both arc e and its reverse arc) or arc e be deleted? */
            if( SCIPisGT(scip, tmpcost, objval) && SCIPisLT(scip, graph->cost[e], FARAWAY) )
            {
               if( fix )
               {
                  assert(vars != NULL);

                  /* try to fix edge */
                  SCIP_CALL( fixedgevar(scip, vars[e], nelims) );
               }
               else
               {
                  if( SCIPisGT(scip, FARAWAY, graph->cost[flipedge(e)]) )
                  {
                     graph->cost[e] = FARAWAY;
                     (*nelims)++;
                  }
                  else
                  {
                     graph_edge_del(scip, graph, e, TRUE);
                     (*nelims)++;
                  }
               }
            }
            e = etemp;
         }
      }
   }

   SCIPdebugMessage("CCC eliminated (edges) in hcrc bound reduce: %d,\n", *nelims);
   /* free memory */
   SCIPfreeBufferArray(scip, &result);

   assert(graph_valid(graph));

   return SCIP_OKAY;
}
