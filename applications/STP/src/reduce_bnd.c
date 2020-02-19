
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_bnd.c
 * @brief  bound based reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements bound-based reduction techniques for several Steiner problems.
 * Most tests can be found in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "heur_ascendprune.h"
#include "misc_stp.h"
#include "probdata_stp.h"
#include "prop_stp.h"


#define BND_TMHEUR_NRUNS 100                  /**< number of runs of constructive heuristic */


/** bound-based reductions for the (R)PCSTP, the MWCSP and the STP */
static
SCIP_RETCODE computeSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real*            cost,               /**< edge cost array                    */
   SCIP_Real*            costrev,            /**< reversed edge cost array           */
   int*                  result,             /**< result edges */
   STP_Bool*             stnode,             /**< result nodes */
   SCIP_Real*            upperbound,         /**< pointer to an upper bound          */
   SCIP_Bool*            success             /**< success? */
   )
{
   int* starts = NULL;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   int runs;
   SCIP_Real obj;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   SCIPdebugMessage("compute Steiner tree \n");

   assert(!graph_pc_isMw(graph));

   runs = 0;
   *upperbound = -FARAWAY;

   for( int k = 0; k < nnodes; ++k )
   {
      stnode[k] = FALSE;

      if( graph->mark[k] )
         runs++;
   }

   for( int e = 0; e < nedges; ++e )
      result[e] = UNKNOWN;

   runs = MIN(runs, BND_TMHEUR_NRUNS);

   if( !pcmw )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &starts, nnodes) );
      SCIPStpHeurTMCompStarts(graph, starts, &runs);
   }

   if( pcmw )
   {
      graph_pc_2trans(scip, graph);

      graph_get_edgeCosts(graph, cost, costrev);
   }

   SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, graph, starts, NULL, result, runs, graph->source, cost, costrev, &obj, NULL, success, FALSE));

   if( pcmw )
   {
      obj = graph_pc_solGetObj(scip, graph, result, 0.0);

      graph_pc_2org(scip, graph);

      assert(SCIPisEQ(scip, obj, graph_pc_solGetObj(scip, graph, result, 0.0)));

      graph_get_edgeCosts(graph, cost, costrev);
   }
   else
   {
      obj = graph_solGetObj(graph, result, 0.0, nedges);
   }

   if( !(*success) )
      return SCIP_OKAY;


   if( graph_pc_isPc(graph) && !graph->extended )

   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         const int head = graph->head[e];
         const int tail = graph->tail[e];

         stnode[head] = TRUE;
         stnode[tail] = TRUE;
      }
   }

   *upperbound = obj;

   return SCIP_OKAY;
}


/** bound-based reductions for the (R)PCSTP, the MWCSP and the STP */
SCIP_RETCODE reduce_bound(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            radius,             /**< radius array                       */
   SCIP_Real*            offset,             /**< pointer to the offset              */
   SCIP_Real*            upperbound,         /**< pointer to an upper bound          */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nelims              /**< pointer to store number of eliminated edges */
   )
{
   GRAPH* adjgraph = NULL;
   PATH* mst = NULL;
   SCIP_Real* prize = NULL;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real  r;
   SCIP_Real  obj;
   SCIP_Real  max;
   SCIP_Real  radii;
   SCIP_Real  bound;
   SCIP_Real  tmpcost;
   SCIP_Real  mstobj;
   SCIP_Real  mstobj2;
   SCIP_Real  radiim2;
   SCIP_Real  radiim3;
   int* perm = NULL;
   int* result = NULL;
   int e;
   int nterms;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   STP_Bool* stnode = NULL;
   SCIP_Bool ub;
   const SCIP_Bool pc = graph_pc_isPc(graph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   SCIP_Bool success = TRUE;

   assert(scip && nelims);
   assert(graph->source >= 0);
   assert(upperbound != NULL);
   assert(graph->stp_type != STP_MWCSP);
   assert(graph->stp_type != STP_RMWCSP);

   obj = FARAWAY;
   mstobj = 0.0;
   *nelims = 0;
   mstobj2 = 0.0;
   ub = SCIPisGT(scip, *upperbound, 0.0);

   graph_mark(graph);

   if( pcmw )
      prize = graph->prize;

   nterms = 0;
   for( int k = 0; k < nnodes; k++ )
      if( graph->mark[k] && Is_term(graph->term[k]) )
         nterms++;

   /* not more than two terminals? */
   if( nterms <= 2 || (pcmw && nterms <= 3) )
      return SCIP_OKAY;

   assert(nterms == (graph->terms - ((graph->stp_type == STP_PCSPG)? 1 : 0)));

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   graph_get_edgeCosts(graph, cost, costrev);

   /* init auxiliary graph */
   SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1) );

   /* build voronoi regions, concomitantly building adjgraph and computing radii values*/
   SCIP_CALL( graph_voronoiWithRadius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );
   graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
   graph_get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   if( pc )
      graph_get4next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   assert(adjgraph != NULL);
   graph_knot_chg(adjgraph, 0, 0);
   adjgraph->source = 0;

   /* compute MST on adjgraph */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
   SCIP_CALL( graph_path_init(scip, adjgraph) );
   graph_path_exec(scip, adjgraph, MST_MODE, 0, adjgraph->cost, mst);

   max = -1.0;
   r = -1.0;
   for( int k = 1; k < nterms; k++ )
   {
      assert(adjgraph->path_state[k] == CONNECT); /* basically asserts that adjgraph is connected */

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

   /* for (rooted) prize collecting an maximum weight problems: correct radius values */
   if( graph->stp_type == STP_RPCSPG )
   {
      assert(graph->mark[graph->source]);
      for( int k = 0; k < nnodes; k++ )
      {
         if( !Is_term(graph->term[k]) || !graph->mark[k] )
            continue;

         if( SCIPisGT(scip, radius[k], prize[k]) && k != graph->source )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_PCSPG )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( Is_term(graph->term[k]) && SCIPisGT(scip, radius[k], prize[k])  )
            radius[k] = prize[k];
      }
   }

   /* sort radius values */
   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );
      for( int k = 0; k < nnodes; k++ )
         perm[k] = k;

      SCIPsortRealInt(radius, perm, nnodes);
   }
   else
   {
      SCIPsortReal(radius, nnodes);
   }

   radiim2 = 0.0;

   /* sum all but two radius values of highest/lowest value */
   for( int k = 0; k < nterms - 2; k++ )
   {
      assert( SCIPisGT(scip, FARAWAY, radius[k]) );
      radiim2 += radius[k];
   }

   radii = radiim2 + radius[nterms - 2] + radius[nterms - 1];

   if( nterms >= 3 )
      radiim3 = radiim2 - radius[nterms - 3];
   else
      radiim3 = 0;

   /* no upper bound available? */
   if( !ub )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &stnode, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

      SCIP_CALL( computeSteinerTree(scip, graph, cost, costrev, result, stnode, &obj, &success) );

      if( !success )
      {
         assert(graph->stp_type == STP_DHCSTP);

         /* free memory and return */
         graph_path_exit(scip, adjgraph);
         graph_free(scip, &adjgraph, TRUE);
         SCIPfreeBufferArrayNull(scip, &mst);
         SCIPfreeBufferArrayNull(scip, &result);
         SCIPfreeBufferArrayNull(scip, &stnode);
         SCIPfreeBufferArray(scip, &costrev);
         SCIPfreeBufferArray(scip, &cost);


         return SCIP_OKAY;
      }

      *upperbound = obj;
   }
   else
   {
      obj = *upperbound;
      assert(SCIPisGE(scip, obj, 0.0));
   }

   assert(SCIPisLT(scip, obj, FARAWAY));

   if( SCIPisGT(scip, radiim2, mstobj) )
   {
      SCIPdebugMessage("select radii bound \n");

      bound = radiim2;
   }
   else
   {
      SCIPdebugMessage("select MST bound \n");

      bound = mstobj;
   }

   SCIPdebugMessage("bound=%f obj=%f \n", bound, obj);

   assert(SCIPisLE(scip, bound, obj));

   /* traverse all node, try to eliminate each node or incident edges */
   for( int k = 0; k < nnodes; k++ )
   {
      if( (!graph->mark[k] && pcmw) || graph->grad[k] == 0 )
         continue;

      tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + bound;

      /* can node k be deleted? */
      if( !Is_term(graph->term[k]) && (SCIPisGT(scip, tmpcost, obj)
            || (((stnode != NULL) ? !stnode[k] : 0) && SCIPisGE(scip, tmpcost, obj))) )
      {
         /* delete all incident edges */
         while( graph->outbeg[k] != EAT_LAST )
         {
            e = graph->outbeg[k];
            (*nelims)++;
            assert(!pc || graph->tail[e] != graph->source);
            assert(!pc || graph->mark[graph->head[e]]);
            assert(!Is_pseudoTerm(graph->term[graph->head[e]]));
            assert(!Is_pseudoTerm(graph->term[graph->tail[e]]));

            SCIPdebugMessage("delete non-terminal edge \n");

            graph_edge_del(scip, graph, e, TRUE);
         }
      }
      else
      {
         e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            const int etemp = graph->oeat[e];
            const int tail = graph->tail[e];
            const int head = graph->head[e];
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
               && SCIPisLT(scip, graph->cost[e], FARAWAY) && (!pc || graph->mark[head]) )
            {
               if( graph->stp_type == STP_DHCSTP && SCIPisGT(scip, graph->cost[e], graph->cost[flipedge(e)]) )
               {
                  graph->cost[e] = FARAWAY;
                  (*nelims)++;
               }
               else
               {
                  assert(!Is_pseudoTerm(graph->term[head]));
                  assert(!Is_pseudoTerm(graph->term[tail]));

                  SCIPdebugMessage("delete terminal edge \n");

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
   for( int k = 0; k < nnodes; k++ )
   {
      if( (!graph->mark[k] && pc) || graph->grad[k] == 0 )
         continue;

      if( (graph->grad[k] == 3 || graph->grad[k] == 4) && !Is_term(graph->term[k]) )
      {
         tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist + radiim3;
         if( SCIPisGT(scip, tmpcost, obj) )
         {
#if 0
            printf("pseudo-delete non-terminal node \n");
            graph_knot_printInfo(graph, k);
            graph_printInfo(graph);
#endif

            SCIP_CALL( graph_knot_delPseudo(scip, graph, graph->cost, NULL, NULL, k, &success) );

            assert(graph->grad[k] == 0 || (graph->grad[k] == 4 && !success));
         }
      }
   }
#endif

   /* for(R)PC: try to eliminate terminals */
   if( pc )
   {
      assert(!graph->extended);

      SCIP_CALL( graph_get4nextTTerms(scip, graph, cost, vnoi, vbase, heap, state) );

      for( int k = 0; k < nnodes; k++ )
      {
         /* is k a terminal other than the root? */
         if( Is_term(graph->term[k]) && graph->mark[k] && graph->grad[k] > 2 && !graph_pc_knotIsFixedTerm(graph, k) )
         {
            int l;
            assert(graph->source != k);

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
               SCIPdebugMessage("alternative bnd elimination! \n\n");

               (*nelims) += graph_pc_deleteTerm(scip, graph, k, offset);
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
                     SCIPdebugMessage("second elimination! prize: %f \n\n", graph->prize[k]);
                     (*nelims) += graph_pc_deleteTerm(scip, graph, k, offset);
                  }
               }
            }
         }
      }
   }

   SCIPdebugMessage("nelims (edges) in bound reduce: %d,\n", *nelims);

   /* free adjgraph */
   graph_path_exit(scip, adjgraph);
   graph_free(scip, &adjgraph, TRUE);

   /* free memory*/
   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArrayNull(scip, &mst);
   SCIPfreeBufferArrayNull(scip, &stnode);
   SCIPfreeBufferArrayNull(scip, &result);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}




/** Bound-based reduction method for the MWCSP .
 * The reduction method tries to eliminate nodes and vertices
 * by checking whether an upper bound for each solution that contains them
 * is smaller than the best known solution value.
 * The essence of the approach is a decomposition of the graph such that this upper bound
 * is minimized.
 * */
SCIP_RETCODE reduce_boundMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure (size 3 * nnodes) */
   PATH*                 path,               /**< shortest path data structure (size nnodes) */
   SCIP_Real*            cost,               /**< edge cost array                    */
   SCIP_Real*            radius,             /**< radius array                       */
   SCIP_Real*            costrev,            /**< reversed edge cost array           */
   SCIP_Real*            offset,             /**< pointer to the offset              */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  result,             /**< solution array or NULL */
   int*                  nelims              /**< pointer to store number of eliminated edges */
   )
{
   PATH* mst;
   SCIP_Real* prize;
   SCIP_Real  obj;
   SCIP_Real  bound;
   SCIP_Real  tmpcost;
   SCIP_Real  radiim2;
   int e;
   int k;
   int head;
   int nterms;
   int nnodes;
   int nedges;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(vnoi != NULL);
   assert(path != NULL);
   assert(cost != NULL);
   assert(radius != NULL);
   assert(costrev != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);
   assert(graph->source >= 0);

   mst = NULL;
   prize = graph->prize;
   nedges = graph->edges;
   nnodes = graph->knots;
   nterms = graph->terms - 1;
   *nelims = 0;

   assert(prize != NULL);

   /* not more than two nodes of positive weight? */
   if( nterms <= 2 )
      /* return */
      return SCIP_OKAY;

   /* initialize cost and costrev array */
   for( e = 0; e < nedges; e++ )
   {
      cost[e] = graph->cost[e];
      costrev[e] = graph->cost[flipedge(e)];

      assert(SCIPisGE(scip, cost[e], 0.0));
   }

   /* compute decomposition of graph and radius values */
   graph_voronoiWithRadiusMw(scip, graph, path, cost, radius, vbase, heap, state);

   /* sum all radius values, exclude two radius values of lowest value */
   for( k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) || !graph->mark[k] )
         continue;

      assert(vbase[k] == k);
      assert(SCIPisGE(scip, prize[k], 0.0));

      if( SCIPisGE(scip, radius[k], FARAWAY) )
         radius[k] = 0.0;
      else
      {
         if( SCIPisGE(scip, radius[k], prize[k] ) )
            radius[k] = 0.0;
         else
            radius[k] = prize[k] - radius[k];
      }
   }

   for( k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] || graph->grad[k] == 0 || Is_term(graph->term[k]) )
         continue;
   }

   /* build Voronoi regions */
   graph_voronoiMw(scip, graph, costrev, vnoi, vbase, heap, state);

   /* get 2nd next positive node to all non-positive nodes */
   graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   for( k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) || !graph->mark[k] )
         continue;

      assert(vbase[k] == k);
   }

   SCIPsortReal(radius, nnodes);

   radiim2 = 0.0;

   for( k = 2; k < nterms; k++ )
   {
      assert( SCIPisGT(scip, FARAWAY, radius[k]) );
      radiim2 += radius[k];
   }

   /* solution available? */
   if( result != NULL)
   {
      /* calculate objective value of solution */
      obj = 0.0;
      for( e = 0; e < nedges; e++ )
      {
         if( result[e] == CONNECT )
         {
            head = graph->head[e];

            if( graph->mark[head] )
               obj += prize[head];
         }
      }
   }
   else
   {
      obj = 0.0;
   }

   bound = radiim2;

   /* traverse all nodes, try to eliminate each non-positive node */
   for( k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] || graph->grad[k] == 0 || Is_term(graph->term[k]) )
         continue;

      assert(SCIPisLE(scip, graph->prize[k], 0.0));

      tmpcost = -vnoi[k].dist - vnoi[k + nnodes].dist + bound + graph->prize[k];

      if( (SCIPisLT(scip, tmpcost, obj)) )
      {
         while( graph->outbeg[k] != EAT_LAST )
         {
            (*nelims)++;
            graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
         }
      }
   }

   SCIPdebugMessage("nelims (edges) in MWCSP bound reduce: %d,\n", *nelims);

   /* free memory*/
   SCIPfreeBufferArrayNull(scip, &mst);

   return SCIP_OKAY;
}



/** bound-based reductions for the (R)PCSTP, the MWCSP and the STP; used by prune heuristic */
SCIP_RETCODE reduce_boundPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            cost,               /**< edge cost array                    */
   SCIP_Real*            radius,             /**< radius array                       */
   SCIP_Real*            costrev,            /**< reversed edge cost array           */
   SCIP_Real*            offset,             /**< pointer to the offset              */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation */
   int*                  vbase,              /**< Voronoi base to each node */
   const int*            solnode,            /**< array of nodes of current solution that is not to be destroyed */
   const int*            soledge,            /**< array of edges of current solution that is not to be destroyed */
   int*                  nelims,             /**< pointer to store number of eliminated edges */
   int                   minelims            /**< minimum number of edges to be eliminated */
   )
{
   SCIP_Real* const prize = graph->prize;
   SCIP_Real obj;
   SCIP_Real max;
   SCIP_Real bound;
   SCIP_Real tmpcost;
   SCIP_Real mstobj;
   SCIP_Real radiim2;
   SCIP_Real radiim3;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const SCIP_Bool mw = (graph->stp_type == STP_MWCSP);
   const SCIP_Bool pc = (graph->stp_type == STP_RPCSPG) || (graph->stp_type == STP_PCSPG);
   const SCIP_Bool pcmw = (pc || mw);
   const int nterms = graph->terms - ( (mw || graph->stp_type == STP_PCSPG) ? 1 : 0 );
   SCIP_Bool eliminate;

   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(cost != NULL);
   assert(heap != NULL);
   assert(graph != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);
   assert(radius != NULL);
   assert(costrev != NULL);
   assert(solnode != NULL);
   assert(soledge != NULL);
   assert(graph->source >= 0);
   assert(graph_valid(scip, graph));
   assert(!graph->extended);
   assert( graph->stp_type != STP_RPCSPG || Is_term(graph->term[root]) );

   mstobj = 0.0;
   *nelims = 0;

   if( nterms <= 2 )
      return SCIP_OKAY;

   /* initialize cost and costrev array */
   for( int e = 0; e < nedges; e++ )
   {
      cost[e] = graph->cost[e];
      costrev[e] = graph->cost[flipedge(e)];

      assert(SCIPisGE(scip, cost[e], 0.0));
   }

   if( !pcmw )
   {
      GRAPH* adjgraph;
      PATH* mst;

      for( int k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

      SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1) );

      /* build Voronoi regions, concomitantly building adjgraph and computing radii values*/
      SCIP_CALL( graph_voronoiWithRadius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

      graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
      graph_get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

      assert(adjgraph != NULL);
      graph_knot_chg(adjgraph, 0, 0);
      adjgraph->source = 0;
      assert(graph_valid(scip, adjgraph));

      /* compute MST on adjgraph */
      SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );
      SCIP_CALL( graph_path_init(scip, adjgraph) );
      graph_path_exec(scip, adjgraph, MST_MODE, 0, adjgraph->cost, mst);

      max = -1.0;
      for( int k = 1; k < nterms; k++ )
      {
         const int e = mst[k].edge;
         assert(adjgraph->path_state[k] == CONNECT);
         assert(e >= 0);
         tmpcost = adjgraph->cost[e];
         mstobj += tmpcost;
         if( SCIPisGT(scip, tmpcost, max) )
            max = tmpcost;
      }
      mstobj -= max;

      /* free memory*/
      SCIPfreeBufferArray(scip, &mst);
      graph_path_exit(scip, adjgraph);
      graph_free(scip, &adjgraph, TRUE);

   }
   else
   {
      if( mw )
      {
         graph_voronoiMw(scip, graph, costrev, vnoi, vbase, heap, state);
         graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
      }
      else
      {
         SCIP_CALL( graph_voronoiWithRadius(scip, graph, NULL, vnoi, radius, cost, costrev, vbase, heap, state) );
         graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
         graph_get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
      }
   }

   /* for (rooted) prize collecting an maximum weight problems: correct radius values */
   if( graph->stp_type == STP_RPCSPG )
   {
      assert(graph->mark[graph->source]);
      for( int k = 0; k < nnodes; k++ )
      {
         if( !Is_term(graph->term[k]) || !graph->mark[k] )
            continue;

         if( SCIPisGT(scip, radius[k], prize[k]) && k != graph->source )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_PCSPG )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( Is_term(graph->term[k]) && SCIPisGT(scip, radius[k], prize[k])  )
            radius[k] = prize[k];
      }
   }

   /* sort radius values */
   if( !mw )
      SCIPsortReal(radius, nnodes);

   radiim2 = 0.0;

   if( mw )
   {
      for( int e = 0; e < nedges; e++ )
         costrev[e] = FARAWAY;
      bound = 0.0;
      radiim3 = 0.0;
   }
   else
   {
      for( int e = 0; e < nedges; e++ )
         costrev[e] = -1.0;

      /* sum all but two radius values of highest/lowest value */
      for( int k = 0; k < nterms - 2; k++ )
      {
         assert( SCIPisGT(scip, FARAWAY, radius[k]) );
         radiim2 += radius[k];
      }
      if( nterms >= 3 )
         radiim3 = radiim2 - radius[nterms - 3];
      else
         radiim3 = 0.0;

      if( SCIPisGT(scip, radiim2, mstobj) )
         bound = radiim2;
      else
      {
         assert(!pcmw);
         bound = mstobj;
      }
   }

   for( int redrounds = 0; redrounds < 3; redrounds++ )
   {
      int nrealelims = MIN(2 * minelims, nedges - 1);

      if( redrounds == 0 )
      {
         eliminate = FALSE;
         obj = FARAWAY;
      }
      else if( redrounds == 1 )
      {
         assert(minelims > 0);
         assert(2 * minelims < nedges);
         eliminate = TRUE;
         SCIPsortReal(costrev, nedges);

         if( mw )
         {
            obj = costrev[nrealelims];
         }
         else
         {
            obj = costrev[nedges - nrealelims];

            if( SCIPisLT(scip, obj, 0.0) )
               obj = 0.0;
         }
      }
      else
      {
         if( mw )
         {
            obj = costrev[nrealelims] + 2 * SCIPepsilon(scip);
         }
         else
         {
            obj = costrev[nedges - nrealelims] - 2 * SCIPepsilon(scip);

            if( SCIPisLT(scip, obj, 0.0) )
               obj = 0.0;
         }
         eliminate = TRUE;
      }

      if( mw )
      {
         for( int k = 0; k < nnodes; k++ )
         {
            if( (*nelims) >= minelims )
               break;

            if( root == k )
               continue;

            if( !graph->mark[k] || graph->grad[k] == 0 || Is_anyTerm(graph->term[k] ) )
               continue;

            tmpcost = -vnoi[k].dist - vnoi[k + nnodes].dist + graph->prize[k];

            if( (!eliminate || SCIPisLT(scip, tmpcost, obj)) && solnode[k] != CONNECT )
            {
               if( eliminate )
               {
                  while (graph->outbeg[k] != EAT_LAST)
                  {
                     (*nelims)++;
                     graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
                  }
               }
               else
               {
                  for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     if( SCIPisLT(scip, tmpcost, costrev[e]) )
                     {
                        costrev[e] = tmpcost;
                        costrev[flipedge(e)] = tmpcost;
                     }
                  }
               }
            }
            else if( solnode[k] == CONNECT )
            {
               int e = graph->outbeg[k];

               while( e != EAT_LAST )
               {
                  const int etemp = graph->oeat[e];
                  const int tail = graph->tail[e];
                  const int head = graph->head[e];

                  if( tail > head || soledge[e] == CONNECT || soledge[flipedge(e)] == CONNECT )
                  {
                     e = etemp;
                     continue;
                  }

                  tmpcost = graph->prize[head] + graph->prize[tail];

                  if( vbase[tail] != vbase[head] )
                  {
                     tmpcost -= vnoi[head].dist + vnoi[tail].dist;
                  }
                  else
                  {
                     if( SCIPisGT(scip, -vnoi[tail].dist -vnoi[head + nnodes].dist, -vnoi[tail + nnodes].dist -vnoi[head].dist) )
                        tmpcost -= vnoi[tail].dist + vnoi[head + nnodes].dist;
                     else
                        tmpcost -= vnoi[tail + nnodes].dist + vnoi[head].dist;
                  }
                  /* can edge e or arc e be deleted? */
                  if( (!eliminate || SCIPisLT(scip, tmpcost, obj))
                     && SCIPisLT(scip, graph->cost[e], FARAWAY) && (graph->mark[head]) )
                  {
                     assert(!Is_pseudoTerm(graph->term[head]));
                     assert(!Is_pseudoTerm(graph->term[tail]));

                     if( eliminate )
                     {
                        graph_edge_del(scip, graph, e, TRUE);
                        (*nelims)++;
                     }
                     else if( SCIPisLT(scip, tmpcost, costrev[e]) )
                     {
                        costrev[e] = tmpcost;
                        costrev[flipedge(e)] = tmpcost;
                     }

                  }
                  e = etemp;
               }
            }
         }
      }
      /* no MWCSP */
      else
      {
         /* traverse all nodes, try to eliminate each node or incident edges */
         for( int k = 0; k < nnodes; k++ )
         {
            if( (*nelims) >= minelims )
               break;
            if( root == k )
               continue;

            if( (!graph->mark[k] && (pcmw)) || graph->grad[k] == 0 )
               continue;

            tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + bound;

            /* can node k be deleted? */
            if( !Is_term(graph->term[k]) && (!eliminate || SCIPisGT(scip, tmpcost, obj)) && solnode[k] != CONNECT )
            {
               /* delete all incident edges */
               if( eliminate )
               {
                  while( graph->outbeg[k] != EAT_LAST )
                  {
                     const int e = graph->outbeg[k];
                     (*nelims)++;

                     assert(!pc || graph->tail[e] != graph->source);
                     assert(!pc || graph->mark[graph->head[e]]);
                     assert(!Is_pseudoTerm(graph->term[graph->head[e]]));
                     assert(!Is_pseudoTerm(graph->term[graph->tail[e]]));

                     graph_edge_del(scip, graph, e, TRUE);
                  }
               }
               else
               {
                  for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     if( SCIPisGT(scip, tmpcost, costrev[e]) )
                     {
                        costrev[e] = tmpcost;
                        costrev[flipedge(e)] = tmpcost;
                     }
                  }
               }
            }
            else
            {
               int e = graph->outbeg[k];
               while( e != EAT_LAST )
               {
                  const int etemp = graph->oeat[e];
                  const int tail = graph->tail[e];
                  const int head = graph->head[e];

                  if( tail > head )
                  {
                     e = etemp;
                     continue;
                  }

                  if( soledge[e] == CONNECT || soledge[flipedge(e)] == CONNECT )
                  {
                     e = etemp;
                     continue;
                  }

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
                  if( (!eliminate || SCIPisGT(scip, tmpcost, obj))
                     && SCIPisLT(scip, graph->cost[e], FARAWAY) && (!(pc) || graph->mark[head]) )
                  {
                     assert(!Is_pseudoTerm(graph->term[head]));
                     assert(!Is_pseudoTerm(graph->term[tail]));

                     if( eliminate )
                     {
                        graph_edge_del(scip, graph, e, TRUE);
                        (*nelims)++;
                     }
                     else if( SCIPisGT(scip, tmpcost, costrev[e]) )
                     {
                        costrev[e] = tmpcost;
                        costrev[flipedge(e)] = tmpcost;
                     }
                  }
                  e = etemp;
               }
            }
         }
#if 1
         /* traverse all nodes, try to eliminate 3,4  degree nodes */
         for( int k = 0; k < nnodes; k++ )
         {
            if( (*nelims) >= minelims )
               break;

            if( (!graph->mark[k] && pc) || graph->grad[k] <= 0 )
               continue;

            if( solnode[k] == CONNECT )
               continue;

            if( !eliminate )
            {
               if( graph->grad[k] >= 3 && graph->grad[k] <= 4 && !Is_term(graph->term[k]) )
               {
                  tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist + radiim3;
                  for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     if( SCIPisGT(scip, tmpcost, costrev[e]) )
                     {
                        costrev[e] = tmpcost;
                        costrev[flipedge(e)] = tmpcost;
                     }
                  }
               }
               continue;
            }

            if( graph->grad[k] >= 3 && graph->grad[k] <= 4 && !Is_term(graph->term[k]) )
            {
               tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist + radiim3;
               if( SCIPisGT(scip, tmpcost, obj) )
               {
                  SCIP_Bool success;

                  SCIP_CALL( graph_knot_delPseudo(scip, graph, graph->cost, NULL, NULL, k, &success) );

                  assert(graph->grad[k] == 0 || (graph->grad[k] == 4 && !success));

                  if( success )
                     (*nelims)++;
               }
            }
         }
      } /* no MWCS */
#endif
   } /* redrounds */

   SCIPdebugMessage("nelims (edges) in bound reduce: %d,\n", *nelims);

   return SCIP_OKAY;
}



/** bound-based reduction test for the HCDSTP */
SCIP_RETCODE reduce_boundHop(
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
   SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1) );

   SCIP_CALL( graph_voronoiWithRadius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

   /* get 2nd next terminals to all nodes */
   graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* compute MST on adjgraph */
   graph_knot_chg(adjgraph, 0, 0);
   adjgraph->source = 0;
   assert(graph_valid(scip, adjgraph));
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
   graph_free(scip, &adjgraph, TRUE);

   /* free memory*/
   SCIPfreeBufferArray(scip, &mst);
   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}

/** hop bound-based reduction test for the HCDSTP */
SCIP_RETCODE reduce_boundHopR(
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
   root = graph->source;
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
   graph_voronoiTerms(scip, graph, costrev, vnoi, vbase, heap, state);

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

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}

/* reduction method for HCSTP */
SCIP_RETCODE reduce_boundHopRc(
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
   success = TRUE;
   vars = NULL;
   root = graph->source;
   nedges = graph->edges;
   nnodes = graph->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

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
         }
         cost[e + 1] = costrev[e];
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
         {
            costrev[e + 1] = BLOCKED;
         }
         else
         {
            costrev[e + 1] = graph->cost[e];
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
   graph_voronoiTerms(scip, graph, costrev, vnoi, vbase, heap, state);

   if( SCIPisLT(scip, objval, 0.0) )
   {
      /* get TM heuristic data */
      assert(SCIPfindHeur(scip, "TM") != NULL);
      tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

      /* compute UB */
      SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, graph, NULL, NULL, result, 50, root, cost, costrev, &hopfactor, NULL, &success, FALSE) );

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
               SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nelims) );

               /* try to fix reversed edge */
               SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[flipedge(e)], nelims) );
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
                  SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nelims) );
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

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}
