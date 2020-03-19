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

/**@file   reduce_da.c
 * @brief  dual-ascent based reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements dual-ascent based techniques for several Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "extreduce.h"
#include "heur_tm.h"
#include "heur_ascendprune.h"
#include "heur_local.h"
#include "heur_rec.h"
#include "solpool.h"
#include "dualascent.h"
#include "probdata_stp.h"

#define BND_TMHEUR_NRUNS 100                  /**< number of runs of constructive heuristic */
#define DEFAULT_DARUNS     7                  /**< number of runs for dual ascent heuristic */
#define DEFAULT_NMAXROOTS  8                  /**< max number of roots to use for new graph in dual ascent heuristic */
#define PERTUBATION_RATIO   0.05              /**< pertubation ratio for dual-ascent primal bound computation */
#define PERTUBATION_RATIO_PC   0.005          /**< pertubation ratio for dual-ascent primal bound computation */
#define SOLPOOL_SIZE 20                       /**< size of presolving solution pool */
#define STP_RED_MINBNDTERMS   750
#define STP_DABD_MAXDEGREE 5
#define STP_DABD_MAXDNEDGES 10
#define DAMAXDEVIATION_RANDOM_LOWER 0.15  /**< random upper bound for max deviation for dual ascent */
#define DAMAXDEVIATION_RANDOM_UPPER 0.30  /**< random upper bound for max deviation for dual ascent */
#define DAMAXDEVIATION_FAST         0.75


/** returns solution value for given edge-solution array (CONNECT/UNKNOWN) and offset, takes prizes into account! */
static
SCIP_Real getSolObj(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge             /**< solution */
)
{
   SCIP_Real obj;
   if( graph_pc_isPc(g) )
      obj = graph_pc_solGetObj(scip, g, soledge, 0.0);
   else
      obj = graph_solGetObj(g, soledge, 0.0, g->edges);

   return obj;
}


/** returns maximum allowed deviation for dual-ascent*/
static
SCIP_Real getDaMaxDeviation(
   const RPDA*           paramsda,           /**< parameters */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
)
{
   SCIP_Real damaxdeviation;

   assert(paramsda && randnumgen);
   assert(paramsda->prevrounds >= 0);

   if( paramsda->prevrounds > 0 )
      damaxdeviation = SCIPrandomGetReal(randnumgen, DAMAXDEVIATION_RANDOM_LOWER, DAMAXDEVIATION_RANDOM_UPPER);
   else
      damaxdeviation = -1.0;

   return damaxdeviation;
}


/** computes dual solution with dual-ascent and guided solution (and possibly reroots given solution) */
static
SCIP_RETCODE computeDualSolutionGuided(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real             damaxdeviation,     /**< maximum deviation for DA */
   REDCOST*              redcostdata,        /**< reduced cost data */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  result              /**< solution array */
)
{
   const int daroot = redcostdata->redCostRoot;
   SCIP_Real* const redcost = redcostdata->redEdgeCost;
   SCIP_Real* dualobjval = &(redcostdata->dualBound);

   /* solution might not be valid anymore */
   if( !graph_solIsValid(scip, graph, result) )
   {
      SCIPdebugMessage("solution not valid; run normal dual-ascent \n");
      SCIP_CALL(SCIPStpDualAscent(scip, graph, redcost, dualobjval, FALSE, FALSE, NULL, NULL, NULL, state, daroot, FALSE, damaxdeviation));
   }
   else
   {
      SCIPdebugMessage("reroot solution and run guided dual-ascent \n");
      SCIP_CALL(graph_solReroot(scip, graph, result, daroot));

      SCIP_CALL(SCIPStpDualAscent(scip, graph, redcost, dualobjval, FALSE, FALSE, NULL, result, NULL, state, daroot, FALSE, damaxdeviation));
   }

   if( STP_RPCSPG == graph->stp_type )
   {
      SCIPdebugMessage("RPC: add %f to dual objective \n", graph_pc_getNonLeafTermOffset(scip, graph));

      *dualobjval += graph_pc_getNonLeafTermOffset(scip, graph);
   }

   return SCIP_OKAY;
}


/** computes dual solution with dual-ascent */
static
SCIP_RETCODE computeDualSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real             damaxdeviation,     /**< maximum deviation for DA */
   REDCOST*              redcostdata,        /**< reduced cost data */
   int*                  state               /**< int 4 * nnodes array for internal computations */
)
{
   const int daroot = redcostdata->redCostRoot;
   SCIP_Real* const redcost = redcostdata->redEdgeCost;
   SCIP_Real* dualobjval = &(redcostdata->dualBound);

   SCIPdebugMessage("no rerooting, run normal dual-ascent \n");
   SCIP_CALL( SCIPStpDualAscent(scip, graph, redcost, dualobjval, FALSE, FALSE, NULL, NULL, NULL, state, daroot, FALSE, damaxdeviation));

   if( STP_RPCSPG == graph->stp_type )
   {
      SCIPdebugMessage("RPC: add %f to dual objective \n", graph_pc_getNonLeafTermOffset(scip, graph));

      *dualobjval += graph_pc_getNonLeafTermOffset(scip, graph);
   }

   return SCIP_OKAY;
}


/** computes TM solution */
static
SCIP_RETCODE computeSteinerTreeTM(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  result,             /**< solution array */
   SCIP_Real*            bestobjval          /**< pointer to the objective value */
)
{
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real obj;
   int* startstm = NULL;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   SCIP_Bool success = FALSE;
   const SCIP_Bool directed = (graph->stp_type == STP_SAP || graph->stp_type == STP_NWSPG);

   /* number of runs should not exceed number of connected vertices */
   int runstm = BND_TMHEUR_NRUNS / (directed ? 1 : 5);

   assert(graph->stp_type != STP_RPCSPG || !graph->extended);

   SCIP_CALL( SCIPallocBufferArray(scip, &startstm, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   graph_get_edgeCosts(graph, cost, costrev);

   SCIPStpHeurTMCompStarts(graph, startstm, &runstm);

   SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fromheurdata,
      graph, startstm, NULL, result, runstm, graph->source, cost, costrev, NULL, NULL, &success) );
   assert(success);

   obj = getSolObj(scip, graph, result);

   if( obj < *bestobjval )
      *bestobjval = obj;

   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &startstm);

   return SCIP_OKAY;
}


/** computes solution from reduced costs */
static
SCIP_RETCODE computeSteinerTreeRedCosts(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   SCIP_Bool             userec,             /**< use recombination? */
   STPSOLPOOL*           pool,               /**< solution pool */
   int*                  result,             /**< result array */
   SCIP_Bool*            bestimproved,       /**< could best solution be improved ? */
   SCIP_Real*            bestobjval          /**< pointer to the objective value */
)
{
   const SCIP_Real* redcosts = redcostdata->redEdgeCost;
   int* nodearrint;
   STP_Bool* nodearrchar;
   const int daroot = redcostdata->redCostRoot;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   SCIP_Bool success;
   SCIP_Bool soladded;
   SCIP_Real objval;

   *bestimproved = FALSE;
   soladded = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, redcosts, result, nodearrint, daroot, nodearrchar, &success, FALSE));
   assert(success && graph_solIsValid(scip, graph, result));

   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &nodearrchar);

   SCIP_CALL(SCIPStpHeurLocalRun(scip, graph, result));
   assert(graph_solIsValid(scip, graph, result));

   objval = getSolObj(scip, graph, result);

   if( userec )
      SCIP_CALL(solpool_addSol(scip, objval, result, pool, &soladded));

   /* should we try recombination? */
   if( userec && soladded && pool->size >= 2 && objval < *bestobjval )
   {
      /* get index of just added solution */
      int solindex = pool->maxindex;
      SCIP_Bool solfound;

      SCIPdebugMessage("POOLSIZE %d \n", pool->size);

      SCIP_CALL(SCIPStpHeurRecRun(scip, pool, NULL, NULL, graph, NULL, &solindex, 1, pool->size, FALSE, &solfound));

      if( solfound )
      {
         const STPSOL* const sol = solpool_solFromIndex(pool, solindex);
         SCIP_Real solobjval;

         assert(sol != NULL);

         solobjval = sol->obj;

         if( graph_pc_isRootedPcMw(graph) )
            solobjval += graph_pc_getNonLeafTermOffset(scip, graph);

         assert(SCIPisEQ(scip, getSolObj(scip, graph, sol->soledges), solobjval));

         SCIPdebugMessage("DA: rec found better solution with obj %f vs %f \n", sol->obj, objval);

         if( SCIPisLT(scip, solobjval, objval) )
         {
            assert(SCIPisLT(scip, getSolObj(scip, graph, sol->soledges), getSolObj(scip, graph, result)));

            BMScopyMemoryArray(result, sol->soledges, nedges);

            SCIP_CALL(SCIPStpHeurLocalRun(scip, graph, result));
            objval = getSolObj(scip, graph, result);

            assert(SCIPisLE(scip, objval, solobjval));

            if( objval < solobjval )
               SCIP_CALL(solpool_addSol(scip, objval, result, pool, &solfound));
         }
      }
   }

   if( SCIPisLE(scip, objval, *bestobjval) )
   {
      *bestimproved = TRUE;
      *bestobjval = objval;
   }

   return SCIP_OKAY;
}


/** compute primal solution during dual-ascent routine for PCSTP or MWCSP based on reduced costs */
static
SCIP_RETCODE computeSteinerTreeRedCostsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   STPSOLPOOL*           pool,               /**< solution pool */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   SCIP_Real*            upperbound,         /**< upperbound pointer */
   int*                  result1,            /**< sol int array corresponding to upper bound */
   int*                  result2,            /**< sol int array corresponding to best new solution (might be worse than upper bound) */
   int*                  vbase,              /**< int array */
   int*                  pathedge,           /**< int array */
   STP_Bool*             nodearrchar,        /**< node array storing solution vertices */
   SCIP_Bool*            apsol               /**< ascend-prune sol? */
)
{
   SCIP_Real ub2;
   const int nedges = graph->edges;
   SCIP_Bool success;

   assert(graph_pc_isPcMw(graph));

   /* compute new solution and store it in result2 */

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, result2, vbase, -1, nodearrchar, &success, FALSE) );
   assert(success);

   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, result2) );

   assert(graph_solIsValid(scip, graph, result2));

   ub2 = getSolObj(scip, graph, result2);
   SCIPdebugMessage("DA: first new sol value in computeSteinerTreeRedCostsPcMw: %f ... old value: %f \n", ub2, *upperbound);

   /* try recombination? */
   if( pool != NULL )
   {
      SCIPdebugMessage("ub %f vs best sol %f\n", ub2, pool->sols[0]->obj);
      SCIP_CALL( solpool_addSol(scip, ub2, result2, pool, &success) );

#ifdef SCIP_DEBUG
      for( int i = 0; i < pool->size; i++ )
         printf(" %f ", pool->sols[i]->obj);
      printf("\n ");
#endif

      if( success && pool->size >= 2 )
      {
         /* get index of just added solution */
         int solindex = pool->maxindex;

         SCIP_Bool solfound;

         SCIPdebugMessage("POOLSIZE %d \n", pool->size);

         SCIP_CALL( SCIPStpHeurRecRun(scip, pool, NULL, NULL, graph, NULL, &solindex, 3, pool->size, FALSE, &solfound) );

         if( solfound )
         {
            STPSOL* sol = solpool_solFromIndex(pool, solindex);
            SCIP_Real solobjval;

            assert(sol != NULL);

            solobjval = sol->obj + graph_pc_getNonLeafTermOffset(scip, graph);

            assert(SCIPisEQ(scip, getSolObj(scip, graph, sol->soledges), solobjval));

            SCIPdebugMessage("DA: rec found better solution with obj %f vs %f \n", sol->obj, ub2);

            if( SCIPisLT(scip, solobjval, ub2) )
            {
               assert(SCIPisLT(scip, getSolObj(scip, graph, sol->soledges), getSolObj(scip, graph, result2)));

               BMScopyMemoryArray(result2, sol->soledges, nedges);

               SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, result2) );

               assert(SCIPisLT(scip, getSolObj(scip, graph, result2), ub2));

               ub2 = getSolObj(scip, graph, result2);

               if( SCIPisLT(scip, ub2, sol->obj) )
                  SCIP_CALL( solpool_addSol(scip, ub2, result2, pool, &success) );
            }
         }
      }
   }

   if( SCIPisLE(scip, ub2, *upperbound) )
   {
      SCIPdebugMessage("DA: improved incumbent %f vs %f, return \n", ub2, *upperbound);

      *apsol = TRUE;
      *upperbound = ub2;
      BMScopyMemoryArray(result1, result2, nedges);
   }

   if( graph->stp_type != STP_MWCSP || !(*apsol) )
     return SCIP_OKAY;

#if 1
   SCIP_CALL( SCIPStpHeurRecExclude(scip, graph, result1, result2, pathedge, nodearrchar, &success) );

   if( success )
   {
      BMScopyMemoryArray(result1, result2, nedges);
      *upperbound = getSolObj(scip, graph, result1);
      SCIPdebugMessage("DA: afterLastExclusion %f \n", *upperbound);
   }
#endif

   return SCIP_OKAY;
}



/** collected terminals (fixed ones for RPC) */
static
void collectFixedTerminals(
   const GRAPH*          graph,              /**< graph data structure */
   int*                  terminals,          /**< terminals array (of size graph->terms) */
   int*                  nterms              /**< number of terminals (might differ for RPC) */
)
{
   int n = 0;
   const int nnodes = graph->knots;
   const SCIP_Bool rpc = graph_pc_isRootedPcMw(graph);

   assert(graph->stp_type != STP_PCSPG && graph->stp_type != STP_MWCSP);
   assert(!rpc || !graph->extended);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(graph->term[i]) )
      {
         assert(graph->mark[i]);
         if( !rpc || graph_pc_knotIsFixedTerm(graph, i) )
            terminals[n++] = i;
      }
   }

   assert(rpc || graph->terms == n);
   *nterms = n;
}

/* initialize dual-ascent distances */
static
SCIP_RETCODE daInitializeDistances(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata,        /**< reduced cost data */
   int*                  pathedge,           /**< path edge */
   int*                  state               /**< state */
   )
{
   const int daroot = redcostdata->redCostRoot;
   const SCIP_Real* const redcosts = redcostdata->redEdgeCost;
   PATH* const vnoi = redcostdata->nodeTo3TermsPaths;
   SCIP_Real* const pathdist = redcostdata->rootToNodeDist;
   int* const vbase = redcostdata->nodeTo3TermsBases;
   SCIP_Real* costrev = NULL;
   const int nedges = g->edges;
   const SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   const SCIP_Bool directed = (g->stp_type == STP_SAP || g->stp_type == STP_NWSPG);

   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   /* distance from root to all nodes */
   graph_path_execX(scip, g, daroot, redcosts, pathdist, pathedge);

   for( int e = 0; e < nedges; e++ )
      costrev[e] = redcosts[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = g->outbeg[daroot]; e != EAT_LAST; e = g->oeat[e] )
      costrev[e] = FARAWAY;

   if( rpc )
      graph_pc_2trans(scip, g);

   /* build Voronoi diagram */
   if( directed )
   {
      assert(!rpc);
      graph_voronoiTerms(scip, g, costrev, vnoi, vbase, g->path_heap, state);
   }
   else
   {
      graph_get4nextTerms(scip, g, costrev, costrev, vnoi, vbase, g->path_heap, state);

#ifndef NDEBUG
      {
         const int nnodes = g->knots;

         for( int i = 0; i < nnodes; i++ )
         {
            if( !g->mark[i] )
               continue;

            if( !Is_term(g->term[i]) )
            {
               assert(vbase[i] != daroot || vnoi[i].dist >= FARAWAY);
               assert(vbase[i + nnodes] != daroot || vnoi[i + nnodes].dist >= FARAWAY);
            }
            else
               assert(vbase[i] == i);
         }
      }
#endif
   }

   if( rpc )
      graph_pc_2org(scip, g);

   SCIPfreeBufferArray(scip, &costrev);

   return SCIP_OKAY;
}

/** updates node bounds for reduced cost fixings */
static
void updateNodeFixingBounds(
   SCIP_Real*            fixingbounds,       /**< fixing bounds */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   SCIP_Real             lpobjval,           /**< LP objective  */
   SCIP_Bool             initialize          /**< initialize fixing bounds? */
)
{
   const int nnodes = graph->knots;

   assert(graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || !graph->extended);

   if( initialize )
      for( int k = 0; k < nnodes; k++ )
         fixingbounds[k] = -FARAWAY;

   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real fixbnd;

      if( !graph->mark[k] )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      fixbnd = pathdist[k] + vnoi[k].dist + lpobjval;

      if( fixbnd > fixingbounds[k] )
      {
         fixingbounds[k] = fixbnd;
      }
   }
}

/** updates node bounds for reduced cost replacement */
static
SCIP_RETCODE updateNodeReplaceBounds(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real*            replacebounds,      /**< replacement bounds */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   const int*            vbase,              /**< bases to Voronoi paths */
   SCIP_Real             lpobjval,           /**< LP objective  */
   SCIP_Real             upperbound,         /**< upper bound */
   int                   root,               /**< DA root */
   SCIP_Bool             initialize,         /**< initialize fixing bounds? */
   SCIP_Bool             extendedsearch      /**< perform extended searching? */
)
{
   unsigned int* eqstack = NULL;
   SCIP_Bool* eqmark = NULL;
   int outedges[STP_DABD_MAXDEGREE];
   const int nnodes = graph->knots;
   const SCIP_Real cutoff = upperbound - lpobjval;
   const int halfnedges = graph->edges / 2;

   assert(!SCIPisNegative(scip, cutoff));
   assert(graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || !graph->extended);

   if( initialize )
      for( int k = 0; k < nnodes; k++ )
         replacebounds[k] = -FARAWAY;

   if( extendedsearch )
   {
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &eqmark, halfnedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &eqstack, halfnedges) );
   }

   for( int node = 0; node < nnodes; node++ )
   {
      const int degree = graph->grad[node];

      if( degree >= 3 && !Is_anyTerm(graph->term[node]) )
      {
         SCIP_Real fixbnd;

         /* bound already good enough? */
         if( SCIPisLT(scip, upperbound, replacebounds[node]) )
               continue;

         fixbnd = pathdist[node] + vnoi[node].dist + vnoi[node + nnodes].dist + lpobjval;

         assert(!Is_pseudoTerm(graph->term[node]));

         /* Y-test for small degrees */
         if( degree <= STP_DABD_MAXDEGREE && !SCIPisLT(scip, upperbound, fixbnd) )
         {
            int eqstack_size = 0;
            int edgecount = 0;

            fixbnd = FARAWAY;

            for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
               outedges[edgecount++] = e;

            assert(edgecount == degree);

            /* compute cost for each root and save minimum */
            for( int i = 0; i < degree; i++ )
            {
               const int tmproot = graph->head[outedges[i]];
               const int rootedge = flipedge(outedges[i]);

               /* loop over all combinations of Y with root tmproot */
               for( int j = i + 1; j <= i + degree - 2; j++ )
               {
                  for( int k = j + 1; k <= i + degree - 1; k++ )
                  {
                     const int outedge1 = outedges[j % degree];
                     const int outedge2 = outedges[k % degree];
                     const int leaf1 = graph->head[outedge1];
                     const int leaf2 = graph->head[outedge2];

                     SCIP_Real tmpcostY;

                     assert(leaf1 != leaf2 && tmproot != leaf1 && tmproot != leaf2);
                     assert(vbase[leaf1] >= 0 || vnoi[leaf1].dist >= FARAWAY);
                     assert(vbase[leaf2] >= 0 || vnoi[leaf2].dist >= FARAWAY);

                     if( leaf1 == root || leaf2 == root )
                        continue;

                     tmpcostY = pathdist[tmproot] + cost[rootedge] + cost[outedge1] + cost[outedge2];

                     if( vbase[leaf1] != vbase[leaf2] )
                     {
                        tmpcostY += vnoi[leaf1].dist + vnoi[leaf2].dist;
                     }
                     else
                     {
                        /* also covers the case that leaf is a terminal */
                        const SCIP_Real leaf1far = vnoi[leaf1 + nnodes].dist;
                        const SCIP_Real leaf2far = vnoi[leaf2 + nnodes].dist;

                        assert(vbase[leaf1 + nnodes] >= 0 || leaf1far == FARAWAY);
                        assert(vbase[leaf2 + nnodes] >= 0 || leaf2far == FARAWAY);

                        tmpcostY += MIN(leaf1far + vnoi[leaf2].dist, vnoi[leaf1].dist + leaf2far);
                     }

                     if( tmpcostY < fixbnd )
                     {
                        if( extendedsearch && SCIPisLE(scip, tmpcostY, cutoff) )
                        {
                           int tree3outedges[2];
                           SCIP_Bool ruleout;
#ifndef NDEBUG
                           const SCIP_Real tmpcostYorg = tmpcostY;
#endif
                           tree3outedges[0] = outedge1;
                           tree3outedges[1] = outedge2;

                           SCIP_CALL( reduce_extendedCheck3Tree(scip, graph, root, cost, pathdist, vnoi, vbase, cutoff, tree3outedges, rootedge,
                                       &tmpcostY, &ruleout, eqstack, &eqstack_size, eqmark) );

                           if( ruleout )
                              tmpcostY = FARAWAY;

#ifndef NDEBUG
                           assert(tmpcostY >= tmpcostYorg);
#endif
                        }

                        if( tmpcostY < fixbnd )
                           fixbnd = tmpcostY;
                     }
                  }
               } /* Y loop */
            } /* root loop */

            fixbnd += lpobjval;

            assert(SCIPisGE(scip, fixbnd, pathdist[node] + vnoi[node].dist + vnoi[node + nnodes].dist + lpobjval)
                  || fixbnd >= FARAWAY);

            if( extendedsearch )
            {
               for( int i = 0; i < eqstack_size; i++ )
                  eqmark[eqstack[i]] = FALSE;

               for( int i = 0; i < halfnedges; i++ )
                  assert(eqmark[i] == FALSE);
            }
         }

         if( fixbnd > replacebounds[node] )
            replacebounds[node] = fixbnd;
      }
   }
   if( extendedsearch )
   {
      assert(eqstack != NULL && eqmark != NULL);

      SCIPfreeBufferArray(scip, &eqstack);
      SCIPfreeCleanBufferArray(scip, &eqmark);
   }

   return SCIP_OKAY;
}


/** updates edge fixing bounds for reduced cost fixings */
static
void updateEdgeFixingBounds(
   SCIP_Real*            fixingbounds,       /**< fixing bounds */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   SCIP_Real             lpobjval,           /**< LP objective  */
   int                   extnedges,          /**< number of edges for extended problem */
   SCIP_Bool             initialize,         /**< initialize fixing bounds? */
   SCIP_Bool             undirected          /**< only consider undirected edges */
)
{
   const int nnodes = graph->knots;

   assert(graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || !graph->extended);
   assert(extnedges > 0);

   if( initialize )
      for( int e = 0; e < extnedges; e++ )
         fixingbounds[e] = -FARAWAY;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      if( undirected )
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            const int head = graph->head[e];

            if( graph->mark[head] )
            {
               const int erev = flipedge(e);
               const SCIP_Real fixbnd = pathdist[k] + cost[e] + vnoi[head].dist + lpobjval;
               const SCIP_Real fixbndrev = pathdist[head] + cost[erev] + vnoi[k].dist + lpobjval;
               const SCIP_Real minbnd = MIN(fixbnd, fixbndrev);

               assert(fixingbounds[e] == fixingbounds[erev]);

               if( minbnd > fixingbounds[e] )
               {
                  fixingbounds[e] = minbnd;
                  fixingbounds[erev] = minbnd;
               }
            }
         }
      }
      else
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            const int head = graph->head[e];

            if( graph->mark[head] )
            {
               const SCIP_Real fixbnd = pathdist[k] + cost[e] + vnoi[head].dist + lpobjval;

               if( fixbnd > fixingbounds[e] )
                  fixingbounds[e] = fixbnd;
            }
         }
      }
   }
}

/** eliminate nodes by using fixing-bounds and reduced costs */
static
int reduceWithNodeFixingBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   GRAPH*                transgraph,         /**< graph data structure or NULL */
   const SCIP_Real*      fixingbounds,       /**< fixing bounds */
   SCIP_Real             upperbound          /**< best upperbound */
)
{
   int nfixed = 0;
   const int nnodes = graph->knots;

   assert(graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || !graph->extended);

   graph_mark(graph);

   for( int k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] || Is_term(graph->term[k]) )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      if( SCIPisLT(scip, upperbound, fixingbounds[k]) )
      {
         SCIPdebugMessage("delete knot %d %f < %f %d\n", k, upperbound, fixingbounds[k], graph->grad[k]);
         nfixed += graph->grad[k];

         graph_knot_del(scip, graph, k, TRUE);

         if( transgraph != NULL )
            graph_knot_del(scip, transgraph, k, FALSE);
      }
   }

   return nfixed;
}


/* marks nodes that can be pseudo-eliminated */
static
void markPseudoDeletablesFromBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      replacebounds,      /**< replacement bounds */
   SCIP_Real             upperbound,         /**< best upper bound */
   SCIP_Bool*            pseudoDelNodes      /**< pseudo deletable nodes */
)
{
   const int nnodes = graph->knots;
   const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);

   for( int k = 0; k < nnodes; k++ )
      pseudoDelNodes[k] = FALSE;

   /* main loop */
   for( int degree = 3; degree <= STP_DABD_MAXDEGREE; degree++ )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( rpc && degree == 3 && graph_pc_knotIsNonLeafTerm(graph, k) && graph->grad[k] == 3 )
         {
            SCIPdebugMessage("found rpc deg3 candidate %d \n", k);
         }
         else if( (degree != graph->grad[k] || Is_anyTerm(graph->term[k])) )
         {
            continue;
         }

         if( SCIPisLT(scip, upperbound, replacebounds[k]))
         {
            pseudoDelNodes[k] = TRUE;
         }
      }
   }
}

/** eliminate edges by using fixing-bounds and reduced costs */
static
int reduceWithEdgeFixingBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   GRAPH*                transgraph,         /**< graph data structure or NULL */
   const SCIP_Real*      fixingbounds,       /**< fixing bounds */
   const int*            result,             /**< solution */
   SCIP_Real             upperbound          /**< best upperbound */
)
{
   int nfixed = 0;
   int nnodes = graph->knots;
   const SCIP_Bool solgiven = (result != NULL);

   assert(graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || !graph->extended);
   assert(!solgiven || SCIPisEQ(scip, upperbound, getSolObj(scip, graph, result)));

   for( int k = 0; k < nnodes; k++ )
   {
      int e;
      if( !graph->mark[k] )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      e = graph->outbeg[k];
      while( e != EAT_LAST )
      {
         const int enext = graph->oeat[e];

         if( graph->mark[graph->head[e]] )
         {
            const int erev = flipedge(e);
            SCIP_Bool delete;

            if( !solgiven || result[e] == CONNECT || result[erev] == CONNECT )
               delete = (SCIPisLT(scip, upperbound, fixingbounds[e]) && SCIPisLT(scip, upperbound, fixingbounds[erev]));
            else
               delete = (SCIPisLE(scip, upperbound, fixingbounds[e]) && SCIPisLE(scip, upperbound, fixingbounds[erev]));

            if( delete )
            {
               assert(graph->cost[e] == graph->cost[erev]);

               SCIPdebugMessage("delete edge %d (upperbound=%f fixingbound=%f) \n", e, upperbound, fixingbounds[e]);

               graph_edge_del(scip, graph, e, TRUE);

               if( transgraph != NULL )
                  graph_edge_del(scip, transgraph, e, FALSE);

               nfixed++;
            }
         }

         e = enext;
      }
   }

   return nfixed;
}


/** submethod for daFindRoots */
static
void findRootsMark(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const GRAPH*          transgraph,         /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const SCIP_Real*      cost,               /**< da reduced cost */
   const SCIP_Real*      bestcost,           /**< best incumbent da reduced cost */
   SCIP_Real             lpobjval,           /**< da lower bound */
   SCIP_Real             bestlpobjval,       /**< best da lower bound */
   SCIP_Real             upperbound,         /**< da upper bound */
   SCIP_Bool             rerun,              /**< not the first run? */
   SCIP_Bool             probrooted,         /**< is transgraph a rooted RMW or RPC? */
   int                   pseudoterm,         /**< pseudo terminal */
   int                   pseudoedge,         /**< pseudo terminal edge */
   STP_Bool*             isfixedterm,        /**< bool array to indicate fixed terminals */
   int*                  roots,              /**< root array */
   int*                  rootscount,         /**< number of roots */
   int*                  pathedge,           /**< array */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Real*            dist                /**< distances array, initially set to FARAWAY */
)
{
   int realterm = -1;
   SCIP_Bool mark = FALSE;

   assert(!graph->extended && transgraph->extended);
   assert(graph->grad[pseudoterm] == 2);

   if( probrooted )
   {
      assert(transgraph->tail[transgraph->term2edge[pseudoterm]] == pseudoterm);

      realterm = transgraph->head[transgraph->term2edge[pseudoterm]];
   }
   else
   {
      for( int e2 = transgraph->inpbeg[pseudoterm]; e2 != EAT_LAST; e2 = transgraph->ieat[e2] )
      {
         if( transgraph->cost[e2] == 0.0 )
            realterm = transgraph->tail[e2];
         else
            assert(pseudoedge == e2); /* that holds because of correspondence between graph and transgraph for the pseudo-terminal edge */
      }
   }

   assert(realterm >= 0 && graph->mark[realterm]);
   assert(realterm != graph->source && realterm != transgraph->source);
   assert(Is_pseudoTerm(transgraph->term[realterm]) && Is_term(graph->term[realterm]));

   if( rerun && isfixedterm[realterm] )
      return;

   if( SCIPisGT(scip, cost[pseudoedge], upperbound - lpobjval) || SCIPisGT(scip, bestcost[pseudoedge], upperbound - bestlpobjval) )
   {
      mark = TRUE;
   }
   else
   {
      /* get terminals that imply realterm and add corresponding reduced costs up */
      int nvisits;
      double costsum = cost[pseudoedge];
      double bestcostsum = bestcost[pseudoedge];

      assert(graph->path_heap != NULL);
      mark = graph_sdWalksConnected(scip, graph, termmark, graph->cost, isfixedterm, realterm, 1500, dist, pathedge, &nvisits,
            visited, TRUE);

#ifndef NDEBUG
      for( int k = 0; k < graph->knots; k++ )
         assert(graph->path_state[k] == UNKNOWN && visited[k] == FALSE && dist[k] == FARAWAY);
#endif

      if( !mark )
      {
         for( int k = 0; k < nvisits; k++ )
         {
            const int node = pathedge[k];

            assert((termmark[node] == 2) == (Is_term(graph->term[node]) && !graph_pc_termIsNonLeafTerm(graph, node)));

            if( termmark[node] == 2 && node != realterm )
            {
               const int nodepterm = graph_pc_getTwinTerm(graph, node);
               const int rootedge = graph_pc_getRoot2PtermEdge(graph, nodepterm);

               assert(graph->mark[node]);
               assert(!graph_pc_knotIsFixedTerm(graph, node));
               assert(graph->grad[nodepterm] == 2);
               assert(rootedge >= 0);
               assert(graph->cost[rootedge] == transgraph->cost[rootedge]);
               assert(graph->cost[rootedge] == graph->prize[node]);

               costsum += cost[rootedge];
               bestcostsum += bestcost[rootedge];
            }
         }

         if( SCIPisGT(scip, costsum, upperbound - lpobjval) || SCIPisGT(scip, bestcostsum, upperbound - bestlpobjval) )
            mark = TRUE;
      }
   }

   if( mark )
   {
      assert(realterm >= 0);

      assert(SCIPisPositive(scip, graph->prize[realterm]));
      assert((*rootscount) < graph->terms);

      roots[(*rootscount)++] = realterm;
      isfixedterm[realterm] = TRUE;
   }
}


/** special method for RPC does deletes incident edges of terminal, but not the terminal and the extension itself */
static
void daRpcDeleteTermIncidents(
   SCIP*                 scip,               /**< SCIP data structure */
   const PATH*           vnoi,               /**< Voronoi data structure */
   int                   term,               /**< the terminal */
   GRAPH*                graph,              /**< graph data structure */
   int*                  incidents,          /**< int array */
   int*                  nfixedp             /**< number of fixed edges pointer */
)
{
   const int twinterm = graph_pc_getTwinTerm(graph, term);
   int incidcount = 0;

#ifndef NDEBUG
   const int termedge = graph->term2edge[term];
   assert(termedge >= 0 && Is_pseudoTerm(graph->term[twinterm]) && graph->cost[termedge] == 0.0);
   assert(vnoi[twinterm].dist == 0.0);
#endif

   for( int e = graph->outbeg[term]; e != EAT_LAST; e = graph->oeat[e] )
      incidents[incidcount++] = e;

   assert(incidcount == graph->grad[term]);
   (*nfixedp) += graph->grad[term] - 1;

   for( int e = 0; e < incidcount; e++ )
   {
      const int edge = incidents[e];

      assert(graph->tail[edge] == term);

      if( graph->head[edge] == twinterm )
         continue;

      graph_edge_del(scip, graph, edge, TRUE);
   }

   assert(graph->grad[term] == 1);
   assert(graph->outbeg[term] == graph->term2edge[term] && twinterm == graph_pc_getTwinTerm(graph, term));
}

/** find roots for PC and MW during DA reduction */
static
SCIP_RETCODE daPcFindRoots(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const GRAPH*          transgraph,         /**< graph data structure */
   const SCIP_Real*      cost,               /**< da reduced cost */
   const SCIP_Real*      bestcost,           /**< best incumbent da reduced cost */
   SCIP_Real             lpobjval,           /**< da lower bound */
   SCIP_Real             bestlpobjval,       /**< best da lower bound */
   SCIP_Real             upperbound,         /**< da upper bound */
   SCIP_Bool             rerun,              /**< not the first run? */
   SCIP_Bool             probrooted,         /**< is transgraph a rooted RMW or RPC? */
   PATH*                 vnoi,               /**< SP array */
   int*                  pathedge,           /**< array */
   int*                  vbase,              /**< array */
   STP_Bool*             isfixedterm,        /**< bool array */
   int*                  roots,              /**< roots (fixed terminals) array */
   int*                  rootscount          /**< number of roots */
)
{
   SCIP_Real* dist;
   STP_Bool* visited;
   int* termmark;
   int* const state = graph->path_state;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const SCIP_Bool graphextended = graph->extended;
   int nroots = *rootscount;
   int nvisits;

   assert(state);
   assert(graph_pc_isPcMw(graph));
   assert(!graph_pc_isRootedPcMw(graph));

   SCIP_CALL(SCIPallocBufferArray(scip, &dist, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &visited, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &termmark, nnodes));

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      state[i] = UNKNOWN;
      dist[i] = FARAWAY;
   }

   assert(transgraph->extended);

   if( graphextended )
      graph_pc_2org(scip, graph);

   graph_pc_termMarkProper(graph, termmark);

   assert(rerun || nroots == 0);

   /*
    * get possible roots
    */

   BMSclearMemoryArray(isfixedterm, nnodes);

   if( rerun )
   {
      for( int i = 0; i < nroots; i++ )
         isfixedterm[roots[i]] = TRUE;
   }

   SCIPdebugMessage("before findDaRootsMark: all roots: %d, nodes: %d edges: %d terms: %d \n",
         nroots, nnodes, graph->edges, graph->terms);

   /* has transgraph non-artificial root (and arcs to pseudo-terminals)? */
   if( probrooted )
   {
      const int transroot = transgraph->source;

      assert(transgraph->term2edge != NULL);

      for( int e = transgraph->outbeg[transroot]; e != EAT_LAST; e = transgraph->oeat[e] )
      {
         const int pseudoterm = transgraph->head[e];

         if( Is_term(transgraph->term[pseudoterm]) && transgraph->term2edge[pseudoterm] >= 0 )
         {
            findRootsMark(scip, graph, transgraph, termmark, cost, bestcost, lpobjval, bestlpobjval, upperbound, rerun, probrooted, pseudoterm, e,
                  isfixedterm, roots, &nroots, pathedge, visited, dist);
         }
      }
   }
   /* transgraph has artificial root, so no arcs to pseudo-terminals */
   else
   {
      for( int e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int pseudoterm = graph->head[e];

         if( Is_pseudoTerm(graph->term[pseudoterm]) )
         {
            findRootsMark(scip, graph, transgraph, termmark, cost, bestcost, lpobjval, bestlpobjval, upperbound, rerun, probrooted, pseudoterm, e,
                  isfixedterm, roots, &nroots, pathedge, visited, dist);
         }
      }
   }

   SCIPdebugMessage("...after: new roots in rerun %d all roots: %d, nodes: %d edges: %d terms: %d \n\n", nroots - *rootscount,
         nroots, nnodes, graph->edges, graph->terms);

   /* could more roots be found? */
   if( nroots > *rootscount && graph->terms > 2 )
   {
      /*
       * try to find additional roots by connecting walks
       */
      SCIP_Bool rerunwalks = TRUE;

      for( int rounds = 0; rounds < 3 && rerunwalks; rounds++ )
      {
         rerunwalks = FALSE;

         for( int i = 0; i < nnodes; i++ )
         {
            SCIP_Bool connected;

            if( !Is_term(graph->term[i]) || isfixedterm[i] || graph_pc_knotIsFixedTerm(graph, i) )
               continue;

            if( graph->grad[i] == 0 )
            {
               assert(graph_pc_isPcMw(graph) && graph_pc_knotIsNonLeafTerm(graph, i));
               continue;
            }

            connected = graph_sdWalksConnected(scip, graph, termmark, graph->cost, isfixedterm, i, 1500, dist, pathedge, &nvisits,
                  visited, TRUE);

            if( connected )
            {
               assert(nroots < graph->terms);

               roots[nroots++] = i;
               isfixedterm[i] = TRUE;
               rerunwalks = TRUE;

               SCIPdebugMessage("WALK-CONNECT: added new root %d prize: %f  \n", i, graph->prize[i]);
            }

#ifndef NDEBUG
            for( int k = 0; k < nnodes; k++ )
            {
               assert(state[k] == UNKNOWN);
               assert(visited[k] == FALSE);
               assert(dist[k] == FARAWAY);
            }
#endif
         }
      } /* for rounds < 3 */

      SCIPdebugMessage("number of terminals found by DA: %d \n", nroots);

   } /* look for additional roots */

   SCIPdebugMessage("new roots in rerun %d all roots: %d, nodes: %d  \n", nroots - *rootscount, nroots, nnodes);

   *rootscount = nroots;

   if( graphextended )
      graph_pc_2trans(scip, graph);

   SCIPfreeBufferArray(scip, &termmark);
   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &dist);

   return SCIP_OKAY;
}


/** set prize of marked terminals to blocked */
static
void daPcMarkRoots(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            roots,              /**< root array */
   int                   nrootsold,          /**< old number of roots */
   int                   nrootsnew,          /**< new number of roots */
   SCIP_Real             prizesum,           /**< sum of positive prizes */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Bool*            userec,             /**< recombination? */
   STPSOLPOOL**          solpool             /**< solution pool */
)
{
   const int root = graph->source;
   const SCIP_Bool graphextended = graph->extended;

   if( graphextended )
      graph_pc_2org(scip, graph);

   if( *userec && *solpool != NULL )
   {
      *userec = FALSE;
      solpool_free(scip, solpool);
      *solpool = NULL;
   }

   assert(graph != NULL);
   assert(roots != NULL);
   assert(!graph->extended);
   assert(nrootsnew >= 0 && nrootsold >= 0);

   for( int i = nrootsold; i < nrootsnew; i++ )
   {
      int e;
      const int term = roots[i];
      const int pterm = graph_pc_getTwinTerm(graph, term);

      assert(Is_term(graph->term[term]));
      assert(SCIPisPositive(scip, graph->prize[term]));
      assert(pterm >= 0);
      assert(Is_pseudoTerm(graph->term[pterm]));

      for( e = graph->inpbeg[pterm]; e != EAT_LAST; e = graph->ieat[e] )
         if( root == graph->tail[e] )
            break;

      assert(e != EAT_LAST);
      assert(SCIPisEQ(scip, graph->prize[term], graph->cost[e]));

      graph->prize[term] = prizesum;
      graph->cost[e] = prizesum;
   }

   if( graphextended )
      graph_pc_2trans(scip, graph);
}

/** are the reduced costs still valid, i.e. are there zero cost paths from the root to all terminals? */
static
SCIP_Bool daRedCostIsValid(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                transgraph,         /**< graph data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   int*                  nodearrint,         /**< int array */
   STP_Bool*             nodearrbool         /**< bool array */
)
{
   int* const queue = nodearrint;
   STP_Bool*  visited = nodearrbool;
   int qsize;
   const int root = transgraph->source;
   const int nnodes = transgraph->knots;

   /*
    * construct new graph corresponding to zero cost paths from the root to all terminals
    */

   BMSclearMemoryArray(visited, nnodes);

   qsize = 0;
   visited[root] = TRUE;
   queue[qsize++] = root;

   /* DFS */
   while( qsize )
   {
      const int k = queue[--qsize];

      /* traverse outgoing arcs */
      for( int a = transgraph->outbeg[k]; a != EAT_LAST; a = transgraph->oeat[a] )
      {
         const int head = transgraph->head[a];

         if( !visited[head] && SCIPisZero(scip, cost[a]) )
         {
            visited[head] = TRUE;
            queue[qsize++] = head;
         }
      }
      assert(qsize <= nnodes);
   }

   for( int k = 0; k < nnodes; k++ )
      if( Is_term(transgraph->term[k]) && !visited[k] )
      {
         return FALSE;
      }

   return TRUE;
}


/** pertubate edge costs for PCMW dual-ascent */
static
void daPcPertubateEdgeCosts(
   SCIP* scip,
   const GRAPH* graph,
   GRAPH* transgraph,
   const int* result,
   STP_Bool* nodearrchar,
   int randomize
)
{
   int e;
   const int root = graph->source;
   const int newroot = transgraph->source;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;

   BMSclearMemoryArray(nodearrchar, nnodes);

   /* mark all vertices visited in regular graph */
   for( e = 0; e < nedges; e++ )
      if( result[e] == CONNECT && graph->tail[e] != root )
         nodearrchar[graph->head[e]] = TRUE;
   srand((unsigned)graph->terms);

   if( graph->stp_type != STP_MWCSP )
   {
      SCIP_Real pratio = PERTUBATION_RATIO_PC;

      for( int k = 0; k < nnodes; k++ )
      {
         assert(Is_anyTerm(graph->term[k]) == Is_anyTerm(transgraph->term[k]) || transgraph->grad[k] == 0);

         if( randomize > 8 )
            pratio = ((SCIP_Real)(rand() % 10)) / (100.0) - 1.0 / 100.0;
         else if( randomize > 6 )
            pratio = ((SCIP_Real)(rand() % 10)) / (200.0);
         else if( randomize > 4 )
            pratio = ((SCIP_Real)(rand() % 10)) / (300.0);
         else if( randomize > 0 )
            pratio = ((SCIP_Real)(rand() % 10)) / 1000.0;
         else
            pratio = PERTUBATION_RATIO_PC + ((SCIP_Real)(rand() % 10)) / 1000.0;

         assert(SCIPisPositive(scip, 1.0 - pratio));
         assert(SCIPisPositive(scip, 1.0 + pratio));

         if( !Is_anyTerm(graph->term[k]) )
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
            {
               assert(transgraph->tail[e] != root);

               if( result[e] == CONNECT || result[flipedge(e)] == CONNECT )
                  transgraph->cost[e] *= 1.0 - pratio;
               else
                  transgraph->cost[e] *= 1.0 + pratio;
            }
         }
         else if( Is_term(transgraph->term[k]) && k != root && k != newroot )
         {
            assert(transgraph->grad[k] == 2);

            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               if( SCIPisPositive(scip, transgraph->cost[e]) )
               {
                  assert(!Is_pseudoTerm(transgraph->term[transgraph->tail[e]]));
                  assert(transgraph->tail[e] != root);
                  assert(result[flipedge(e)] != CONNECT);

                  if( result[e] == CONNECT )
                     transgraph->cost[e] *= 1.0 - pratio;
                  else
                     transgraph->cost[e] *= 1.0 + pratio;

                  assert(SCIPisPositive(scip, transgraph->cost[e]));
               }
         }
      }

      return;
   }

   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real pratio = PERTUBATION_RATIO;

      assert(Is_anyTerm(graph->term[k]) == Is_anyTerm(transgraph->term[k]));

      if( randomize > 8 )
         pratio = ((SCIP_Real)(rand() % 10)) / (50.0) - 1.0 / 10.0;
      else if( randomize > 6 )
         pratio = ((SCIP_Real)(rand() % 10)) / (20.0);
      else if( randomize > 4 )
         pratio = ((SCIP_Real)(rand() % 10)) / (30.0);
      else if( randomize > 0 )
         pratio = ((SCIP_Real)(rand() % 10)) / 100.0;
      else
         pratio = PERTUBATION_RATIO + ((SCIP_Real)(rand() % 10)) / 200.0;

      if( !Is_anyTerm(graph->term[k]) )
      {
         if( nodearrchar[k] )
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               transgraph->cost[e] *= 1.0 - pratio;
         }
         else
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               transgraph->cost[e] *= 1.0 + pratio;
         }
      }
      else if( Is_term(transgraph->term[k]) && k != root && k != newroot )
      {
         assert(transgraph->grad[k] == 2);

         if( nodearrchar[k] )
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               if( SCIPisPositive(scip, transgraph->cost[e]) )
               {
                  assert(!Is_pseudoTerm(transgraph->term[transgraph->tail[e]]));

                  transgraph->cost[e] *= 1.0 + pratio;
               }
         }
         else
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               if( SCIPisPositive(scip, transgraph->cost[e]) )
               {
                  assert(!Is_pseudoTerm(transgraph->term[transgraph->tail[e]]));
                  transgraph->cost[e] *= 1.0 - pratio;
               }
         }
      }
   }
}

/** order roots */
static
SCIP_RETCODE daOrderRoots(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   int*                  terms,              /**< sol int array corresponding to upper bound */
   int                   nterms,             /**< number of terminals */
   SCIP_Bool             randomize,          /**< randomize */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
)
{
   int* termdegs;
   int maxdeg = 0;
   const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);

   assert(terms != NULL);
   assert(nterms > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &termdegs, nterms) );

   for( int i = 0; i < nterms; i++ )
   {
      const int grad = graph->grad[terms[i]];
      assert(terms[i] >= 0);

      termdegs[i] = -grad;

      if( grad > maxdeg )
         maxdeg = termdegs[i];

      if( rpc )
      {
         assert(graph_pc_knotIsFixedTerm(graph, terms[i] ));

         /* make sure root is selected for RPC */
         if( terms[i] == graph->source )
            termdegs[i] = -graph->knots;
      }
   }

   if( randomize )
      for( int i = 0; i < nterms; i++ )
         termdegs[i] -= SCIPrandomGetInt(randnumgen, 0, maxdeg);

   SCIPsortIntInt(termdegs, terms, nterms);

   SCIPfreeBufferArray(scip, &termdegs);

   return SCIP_OKAY;
}


/** try to improve both dual and primal solution */
static
SCIP_RETCODE computePertubedSol(
   SCIP* scip,
   GRAPH* graph,
   GRAPH* transgraph,
   STPSOLPOOL* pool,
   PATH* vnoi,
   GNODE** gnodearr,
   SCIP_Real* cost,
   SCIP_Real* bestcost,
   SCIP_Real* pathdist,
   int* state,
   int* vbase,
   int* pathedge,
   int* result,
   int* result2,
   int* transresult,
   STP_Bool* nodearrchar,
   SCIP_Real* upperbound,
   SCIP_Real* lpobjval,
   SCIP_Real* bestlpobjval,
   SCIP_Real* minpathcost,
   SCIP_Bool* apsol,
   SCIP_Real offset,
   int extnedges,
   int pertubation
)
{
   SCIP_Real lb;
   int e;
   const int root = graph->source;
   const int nedges = graph->edges;
   const int transnedges = transgraph->edges;

   assert(graph_pc_isPcMw(graph));

   graph_pc_2transcheck(scip, graph);

   /* pertubate the reduced cost? */
   if( graph->stp_type == STP_MWCSP )
   {
      SCIP_Real* transcost;

      SCIP_CALL( SCIPallocBufferArray(scip, &transcost, transnedges) );
      BMScopyMemoryArray(transcost, transgraph->cost, transnedges);

      /* result contains no valid solution?*/
      if( !(*apsol) )
      {
         /* compute new solution */
         SCIP_Real bound;
         SCIP_Bool success;
         SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, bestcost, result2, vbase, -1, nodearrchar, &success, FALSE) );
         assert(success);

         SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, result2) );

         assert(graph_solIsValid(scip, graph, result2));

         bound = getSolObj(scip, graph, result2);

         if( SCIPisLE(scip, bound, *upperbound) )
         {
            *upperbound = bound;
            *apsol = TRUE;
            BMScopyMemoryArray(result, result2, nedges);
         }
         daPcPertubateEdgeCosts(scip, graph, transgraph, result2, nodearrchar, pertubation);
      }
      else
      {
         daPcPertubateEdgeCosts(scip, graph, transgraph, result, nodearrchar, pertubation);
      }

      /* todo use result as guiding solution? */
      SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, &lb, FALSE, FALSE, gnodearr, NULL, transresult, state, root, TRUE, -1.0) );

      BMScopyMemoryArray(transgraph->cost, transcost, transnedges);

      SCIPfreeBufferArray(scip, &transcost);
   }

   SCIP_CALL( computeSteinerTreeRedCostsPcMw(scip, graph, pool, cost, upperbound, result, result2, vbase, pathedge, nodearrchar, apsol) );

   /* does result not contain a valid solution? */
   if( !(*apsol) )
      BMScopyMemoryArray(result, result2, nedges);

   SCIPdebugMessage("DA: pertubated sol value %f \n", *upperbound);

   /* compute guiding solution */

   BMScopyMemoryArray(transresult, result, nedges);

   for( e = nedges; e < transnedges; e++ )
       transresult[e] = UNKNOWN;

   /* non-trivial solution */
   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      if( Is_term(graph->term[graph->head[e]]) && result[e] != CONNECT )
         break;

   if( e != EAT_LAST)
   {
      int k;
      for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         if( !Is_term(graph->term[graph->head[e]]) && result[e] == CONNECT )
            break;

      assert(e != EAT_LAST);
      k = graph->head[e];

      for( e = transgraph->outbeg[k]; e != EAT_LAST; e = transgraph->oeat[e] )
      {
         if( transgraph->head[e] == graph->knots )
         {
            transresult[e] = CONNECT;
            break;
         }
      }
   }

   assert(!(*apsol) || graph_solIsValid(scip, graph, result));
   assert(graph_valid(scip, transgraph));
   assert(root == transgraph->source);

   SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, &lb, FALSE, FALSE, gnodearr, transresult, NULL, state, root, TRUE, -1.0) );

   lb += offset;
   *lpobjval = lb;

   SCIP_CALL( computeSteinerTreeRedCostsPcMw(scip, graph, pool, cost, upperbound, result, result2, vbase, pathedge, nodearrchar, apsol) );

   assert(!(*apsol) || graph_solIsValid(scip, graph, result));

   if( SCIPisGE(scip, lb, *bestlpobjval) )
   {
      *bestlpobjval = lb;
      BMScopyMemoryArray(bestcost, cost, extnedges);
   }

   /* the required reduced path cost to be surpassed */
   *minpathcost = *upperbound - *lpobjval;

   return SCIP_OKAY;
}

/** compute Voronoi region for dual-ascent elimination for PC/MW */
static
void computeTransVoronoi(
    SCIP* scip,
    GRAPH* transgraph,
    PATH* vnoi,
    const SCIP_Real* cost,
    SCIP_Real* costrev,
    SCIP_Real* pathdist,
    int* vbase,
    int* pathedge
)
{
   const int root = transgraph->source;
   const int transnnodes = transgraph->knots;
   const int transnedges = transgraph->edges;

   for( int k = 0; k < transnnodes; k++ )
      transgraph->mark[k] = (transgraph->grad[k] > 0);

   /* distance from root to all nodes */
   graph_path_execX(scip, transgraph, root, cost, pathdist, pathedge);

   for( int k = 0; k < transnnodes; k++ )
      if( Is_term(transgraph->term[k]) )
         assert(SCIPisEQ(scip, pathdist[k], 0.0));

   for( int e = 0; e < transnedges; e++ )
      costrev[e] = cost[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = transgraph->outbeg[root]; e != EAT_LAST; e = transgraph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build Voronoi diagram wrt incoming arcs */
   graph_voronoiTerms(scip, transgraph, costrev, vnoi, vbase, transgraph->path_heap, transgraph->path_state);
}


/** reduce graph with non-artificial root (SPG, RPC ...) based on information from dual ascent and given upper bound  */
static
SCIP_RETCODE reduceRootedProb(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             marked,             /**< edge array to mark which (directed) edge can be removed */
   STP_Bool*             nodearrchar,        /**< node array storing solution vertices */
   const PATH*           vnoi,               /**< Voronoi data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   const SCIP_Real*      pathdist,           /**< distance array from shortest path calculations */
   const int*            result,             /**< sol int array */
   SCIP_Real             minpathcost,        /**< the required reduced path cost to be surpassed */
   int                   root,               /**< the root */
   SCIP_Bool             solgiven,           /**< is sol given? */
   int*                  nfixedp             /**< number of fixed edges pointer */
)
{
   int* incidents = NULL;
   const int nnodes = graph->knots;
   const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);
   const SCIP_Bool keepsol = (solgiven && SCIPisZero(scip, minpathcost));

   if( rpc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &incidents, nnodes) );

#ifndef NDEBUG
      assert(!graph->extended);
      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(graph->term[k]) )
            assert(!graph->mark[k] && graph->grad[k] == 2);
         else
            assert(graph->mark[k] || graph->grad[k] == 0);
      }
#endif
   }

   if( solgiven )
   {
      graph_solSetVertexFromEdge(graph, result, nodearrchar);
   }

   /* main loop: try to reduce */
   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real redcost;
      int e;

      if( !graph->mark[k] )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      redcost = pathdist[k] + vnoi[k].dist;

      if( rpc && Is_term(graph->term[k]) && !graph_pc_knotIsFixedTerm(graph, k) && !graph_pc_termIsNonLeafTerm(graph, k)
         && SCIPisGT(scip, redcost, minpathcost) )
      {
         daRpcDeleteTermIncidents(scip, vnoi, k, graph, incidents, nfixedp);
         continue;
      }

      /* note: if we want to keep the solution we cannot just delete vertices */
      if( !Is_term(graph->term[k]) && !keepsol &&
         (SCIPisGT(scip, redcost, minpathcost) || (solgiven && SCIPisEQ(scip, redcost, minpathcost) && !nodearrchar[k])) )
      {
         (*nfixedp) += graph->grad[k];
         graph_knot_del(scip, graph, k, TRUE);
         continue;
      }

      e = graph->outbeg[k];
      while( e != EAT_LAST )
      {
         const int head = graph->head[e];
         const int enext = graph->oeat[e];

         /* for rpc no artificial terminal arcs should be deleted */
         if( (rpc && !graph->mark[head])
          || (keepsol && (result[e] == CONNECT || result[flipedge(e)] == CONNECT)) )
         {
            e = enext;
            continue;
         }

         redcost = pathdist[k] + cost[e] + vnoi[head].dist;

         if( SCIPisGT(scip, redcost, minpathcost)
            || (solgiven && SCIPisEQ(scip, redcost, minpathcost) && result[e] != CONNECT && result[flipedge(e)] != CONNECT) )
         {
            if( marked[flipedge(e)] )
            {
               graph_edge_del(scip, graph, e, TRUE);
               (*nfixedp)++;
            }
            else
            {
               marked[e] = TRUE;
            }
         }
         e = enext;
      }
   }

   SCIPfreeBufferArrayNull(scip, &incidents);

   return SCIP_OKAY;
}

/** reduce PCSTP or MWCS graph based on information from dual ascent and given upper bound  */
static
int reducePcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   GRAPH*                transgraph,         /**< graph data structure */
   const PATH*           vnoi,               /**< Voronoi data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   const SCIP_Real*      pathdist,           /**< distance array from shortest path calculations */
   SCIP_Real             minpathcost,        /**< the required reduced path cost to be surpassed */
   const int*            result,             /**< sol int array */
   STP_Bool*             marked,             /**< edge array to mark which (directed) edge can be removed */
   STP_Bool*             nodearrchar,        /**< node array storing solution vertices */
   SCIP_Bool             solgiven            /**< is sol given? */
)
{
   const int nnodes = graph_get_nNodes(graph);
   int nfixed;
   SCIP_Real tmpcost;
   SCIP_Bool keepsol = FALSE;

   assert(SCIPisGE(scip, minpathcost, 0.0));

   if( minpathcost < 0.0 )
      minpathcost = 0.0;

   nfixed = 0;

   graph_pc_2orgcheck(scip, graph);

   if( solgiven )
   {
      assert(graph_solIsValid(scip, graph, result));

      graph_solSetVertexFromEdge(graph, result, nodearrchar);

      if( SCIPisZero(scip, minpathcost) )
         keepsol = TRUE;
   }

   /* try to eliminate vertices and edges */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] )
         continue;

      assert(!Is_pseudoTerm(graph->term[k]));

      if( Is_term(graph->term[k]) )
      {
         int e = graph->outbeg[k];
         while( e != EAT_LAST )
         {
            const int etmp = graph->oeat[e];
            tmpcost = pathdist[k] + cost[e] + vnoi[graph->head[e]].dist;

            if( graph->mark[graph->head[e]] &&
               ((SCIPisGT(scip, tmpcost, minpathcost) && (!keepsol || (result[e] != CONNECT && result[flipedge(e)] != CONNECT))) ||
                  (solgiven && tmpcost >= minpathcost && result[e] != CONNECT && result[flipedge(e)] != CONNECT)) )
            {
               if( marked[flipedge(e)] )
               {
                  assert(!Is_pseudoTerm(graph->term[graph->head[e]]));

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
         continue;
      }

      tmpcost = pathdist[k] + vnoi[k].dist;

      if( (SCIPisGT(scip, tmpcost, minpathcost) && !keepsol) ||
         (solgiven && tmpcost >= minpathcost && !nodearrchar[k]))
      {
         while( transgraph->outbeg[k] != EAT_LAST )
         {
            const int e = transgraph->outbeg[k];

            graph_edge_del(scip, transgraph, e, FALSE);
            graph_edge_del(scip, graph, e, TRUE);
            nfixed++;
         }
         assert(graph->outbeg[k] == EAT_LAST);
      }
      else
      {
         int e = transgraph->outbeg[k];
         while( e != EAT_LAST )
         {
            const int etmp = transgraph->oeat[e];
            tmpcost = pathdist[k] + cost[e] + vnoi[transgraph->head[e]].dist;

            if( (SCIPisGT(scip, tmpcost, minpathcost) && (!keepsol || (result[e] != CONNECT && result[flipedge(e)] != CONNECT))) ||
               (solgiven && tmpcost >= minpathcost && result[e] != CONNECT && result[flipedge(e)] != CONNECT) )
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
   SCIPdebugMessage("DA: eliminations %d \n", nfixed);

   return nfixed;
}

/** try to run reduction method for best known reduced costs (if they are valid) */
static
int reducePcMwTryBest(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   GRAPH*                transgraph,         /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   SCIP_Real*            costrev,            /**< SCIP_Real array */
   SCIP_Real*            bestcost,           /**< best dual ascent costs */
   SCIP_Real*            pathdist,           /**< distance array from shortest path calculations */
   SCIP_Real*            upperbound,         /**< upper bound */
   SCIP_Real*            lpobjval,           /**< reduced cost value */
   SCIP_Real*            bestlpobjval,       /**< best reduced cost value */
   SCIP_Real*            minpathcost,        /**< the required reduced path cost to be surpassed */
   SCIP_Real             oldupperbound,      /**< old upper bound */
   const int*            result,             /**< sol int array */
   int*                  vbase,              /**< array for Voronoi bases */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  pathedge,           /**< array for predecessor edge on a path */
   STP_Bool*             marked,             /**< edge array to mark which (directed) edge can be removed */
   STP_Bool*             nodearrchar,        /**< node array storing solution vertices */
   SCIP_Bool*            solgiven,           /**< is sol given? */
   int                   extnedges           /**< number of edges for transgraph */
)
{
   const SCIP_Bool dualisvalid = daRedCostIsValid(scip, transgraph, bestcost, state, nodearrchar);

   if( !dualisvalid )
   {
      *bestlpobjval = *lpobjval;
      BMScopyMemoryArray(bestcost, cost, extnedges);
   }

   /* has upperbound and changed, and is best reduced cost valid and different from cost? */
   if( SCIPisGT(scip, *bestlpobjval, *lpobjval) && SCIPisLT(scip, *upperbound, oldupperbound) )
   {
      *solgiven = *solgiven && graph_solIsUnreduced(scip, graph, result);

      assert(!(*solgiven) || graph_solIsValid(scip, graph, result));

      *minpathcost = *upperbound - *bestlpobjval;
      assert(SCIPisGE(scip, *minpathcost, 0.0));

      computeTransVoronoi(scip, transgraph, vnoi, bestcost, costrev, pathdist, vbase, pathedge);

      return reducePcMw(scip, graph, transgraph, vnoi, bestcost, pathdist, *minpathcost, result, marked, nodearrchar, *solgiven);
   }
   return 0;
}

/** dual ascent based reductions */
SCIP_RETCODE reduce_da(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const RPDA*           paramsda,           /**< parameters */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   SCIP_Real*            ub,                 /**< pointer to provide upper bound and return upper bound found during ascent and prune (if better) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  edgearrint,         /**< int edges array for internal computations or NULL */
   int*                  vbase,              /**< array for Voronoi bases */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  pathedge,           /**< array for predecessor edge on a path */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
)
{
   STPSOLPOOL* pool = NULL;
   SCIP_Real* edgefixingbounds = NULL;
   SCIP_Real* nodefixingbounds = NULL;
   SCIP_Real* nodereplacebounds = NULL;
   SCIP_Real upperbound;
   SCIP_Real damaxdeviation;
   const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);
   const SCIP_Bool directed = (graph->stp_type == STP_SAP || graph->stp_type == STP_NWSPG);
   int* terms;
   int* result;
   SCIP_Bool* pseudoDelNodes;
   int nruns;
   int nFixedTerms;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   int ndeletions;
   STP_Bool* arcsdeleted;
   const SCIP_Bool extended = paramsda->useExtRed;
   const SCIP_Bool nodereplacing = paramsda->nodereplacing;
   const SCIP_Bool userec = paramsda->useRec;
   const int prevrounds = paramsda->prevrounds;
   REDCOST redcostdata = { .redEdgeCost = NULL, .rootToNodeDist = pathdist, .nodeTo3TermsPaths = vnoi,
      .nodeTo3TermsBases = vbase, .cutoff = -1.0, .dualBound = FARAWAY, .redCostRoot = -1
#ifndef NDEBUG
      , .nnodes = graph->knots, .nedges = graph->edges
#endif
   };

   assert(ub && scip && graph && nelims);
   assert(graph_valid_ancestors(scip, graph) && graph_valid(scip, graph));
   assert(!rpc || !graph->extended);

   if( graph->terms <= 2 )
      return SCIP_OKAY;

#ifdef STP_RPC_FIXEDPROPER
   assert(0 && "check");
   if( rpc )
      return SCIP_OKAY;
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &pseudoDelNodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(redcostdata.redEdgeCost), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &arcsdeleted, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgefixingbounds, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodefixingbounds, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodereplacebounds, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &terms, graph->terms) );

   if( userec )
      SCIP_CALL( solpool_init(scip, &pool, nedges, SOLPOOL_SIZE) );

   if( rpc )
      reduce_identifyNonLeafTerms(scip, graph);

   ndeletions = 0;
   upperbound = SCIPisGE(scip, *ub, 0.0) ? (*ub) : FARAWAY;
   graph_mark(graph);
   collectFixedTerminals(graph, terms, &nFixedTerms);
   assert(nFixedTerms >= 1);

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   if( directed || rpc )
   {
      SCIP_CALL( computeSteinerTreeTM(scip, graph, result, &upperbound) );
   }

   nruns = MIN(nFixedTerms, DEFAULT_DARUNS);
   assert(nruns > 0);

   /* select roots for dual ascent */
   SCIP_CALL( daOrderRoots(scip, graph, terms, nFixedTerms, (prevrounds > 0), randnumgen) );

   // todo maybe inside the loop?
   damaxdeviation = getDaMaxDeviation(paramsda, randnumgen);

   assert(!rpc || graph->extended);

   for( int outerrounds = 0; outerrounds < 2; outerrounds++ )
   {
      SCIP_Real cutoffbound = -1.0;

      /* main reduction loop */
      for( int run = 0; run < nruns; run++ )
      {
         const SCIP_Bool guidedDa = (run > 1) && (SCIPrandomGetInt(randnumgen, 0, 2) < 2) && graph->stp_type != STP_RSMT;
         SCIP_Bool havenewsol = FALSE;
         redcostdata.redCostRoot = terms[run];

         /* graph vanished? */
         if( graph->grad[graph->source] == 0 )
            break;
     //   if( rpc ) {      // int todo; // check for more terminals to be added    }
         if( guidedDa )
         {
            /* run dual-ascent (and possibly re-root solution stored in 'result') */
            SCIP_CALL( computeDualSolutionGuided(scip, graph, damaxdeviation, &redcostdata, state, result) );
         }
         else
         {
            SCIP_CALL( computeDualSolution(scip, graph, damaxdeviation, &redcostdata, state) );
         }

         if( !directed )
         {
            SCIP_CALL( computeSteinerTreeRedCosts(scip, graph, &redcostdata, userec, pool, result, &havenewsol, &upperbound) );
         }

         cutoffbound = upperbound - redcostdata.dualBound;
         redcostdata.cutoff = cutoffbound;

         SCIPdebugMessage("upper=%f lower=%f (round=%d, outerround=%d)\n", upperbound, redcostdata.dualBound, run, outerrounds);
         assert(SCIPisGE(scip, cutoffbound, 0.0));

         if( rpc )
            graph_pc_2org(scip, graph);
         else
            graph_mark(graph);

         for( int e = 0; e < nedges; e++ )
            arcsdeleted[e] = FALSE;

         daInitializeDistances(scip, graph, &redcostdata, pathedge, state);
         updateNodeFixingBounds(nodefixingbounds, graph, redcostdata.rootToNodeDist, redcostdata.nodeTo3TermsPaths, redcostdata.dualBound, (run == 0));
         updateEdgeFixingBounds(edgefixingbounds, graph, redcostdata.redEdgeCost, redcostdata.rootToNodeDist, redcostdata.nodeTo3TermsPaths, redcostdata.dualBound, nedges, (run == 0), TRUE);

         SCIP_CALL( reduceRootedProb(scip, graph, arcsdeleted, nodearrchar, redcostdata.nodeTo3TermsPaths, redcostdata.redEdgeCost, redcostdata.rootToNodeDist, result, cutoffbound, redcostdata.redCostRoot, havenewsol, &ndeletions) );

         if( !SCIPisZero(scip, cutoffbound) )
         {
            ndeletions += reduceWithNodeFixingBounds(scip, graph, NULL, nodefixingbounds, upperbound);
            havenewsol = havenewsol && graph_solIsUnreduced(scip, graph, result);
            ndeletions += reduceWithEdgeFixingBounds(scip, graph, NULL, edgefixingbounds, (havenewsol ? result : NULL), upperbound);
         }

         // todo don't call anymore!
         if( extended && !rpc )
         {
            int extfixed = reduce_extendedEdge(scip, graph, redcostdata.nodeTo3TermsPaths, redcostdata.redEdgeCost, redcostdata.rootToNodeDist, (havenewsol ? result : NULL), cutoffbound, redcostdata.redCostRoot, arcsdeleted, FALSE);
            ndeletions += extfixed;
        //    printf("newly fixedFIRST =%d \n", extfixed);
         }

         // todo call this methods fewer times, at the end maybe
         if( extended && !SCIPisZero(scip, cutoffbound) )
         {
            int extfixed;

            SCIP_CALL( extreduce_deleteEdges(scip, &redcostdata, (havenewsol ? result : NULL), graph, arcsdeleted, &extfixed) );
            ndeletions += extfixed;
//#define EXT_WRITE
         //   graph_printInfo(graph);
          //  printf("newly fixedSECOND =%d \n", extfixed);
          //  exit(1);
#ifdef EXT_WRITE
            {
               FILE *fp;
               fp = fopen("/nfs/optimi/kombadon/bzfrehfe/projects/scip/applications/STP/STATS/stpall_hash.txt", "a+");
               fprintf(fp, "%s %d \n", SCIPgetProbName(scip), extfixed);
               fclose(fp);
               exit(1);
            }
#endif
         }

         if( !directed && !SCIPisZero(scip, cutoffbound) && nodereplacing )
            SCIP_CALL( updateNodeReplaceBounds(scip, nodereplacebounds, graph, redcostdata.redEdgeCost, redcostdata.rootToNodeDist, redcostdata.nodeTo3TermsPaths, redcostdata.nodeTo3TermsBases,
                  redcostdata.dualBound, upperbound, redcostdata.redCostRoot, (run == 0), extended && !rpc));

         if( ndeletions > 0 && !rpc )
            reduceLevel0(scip, graph);

         assert(graph_valid(scip, graph));

         if( !rpc )
            graph_mark(graph);
         else
            graph_pc_2trans(scip, graph);

      } /* root loop */

      /* do pseudo-elimination? */
      if( !directed && !SCIPisZero(scip, cutoffbound) && nodereplacing )
      {
         int nreplacings = 0;
         markPseudoDeletablesFromBounds(scip, graph, nodereplacebounds, upperbound, pseudoDelNodes);

         if( extended )
         {
            const SCIP_Bool havenewsol = graph_solIsUnreduced(scip, graph, result);
            assert(!havenewsol || graph_solIsValid(scip, graph, result));

            SCIP_CALL( extreduce_updatePseudoDeletableNodes(scip, &redcostdata, (havenewsol ? result : NULL), pseudoDelNodes,
                  graph, arcsdeleted) );
         }

         SCIP_CALL( reduce_applyPseudoDeletions(scip, &redcostdata, pseudoDelNodes, graph, offsetp, &nreplacings) );
        // printf("nreplacings=%d \n", nreplacings);
         ndeletions += nreplacings;

         if( nreplacings > 0 && userec )
         {
            /* solutions in pool might not be valid anymore */
            solpool_free(scip, &pool);
            SCIP_CALL(solpool_init(scip, &pool, nedges, SOLPOOL_SIZE));
         }

         assert(graph_valid_ancestors(scip, graph));
         graph_mark(graph);
      }

      if( ndeletions == 0 || !userec )
         break;
      else if( userec )
         damaxdeviation = SCIPrandomGetReal(randnumgen, DAMAXDEVIATION_RANDOM_LOWER, DAMAXDEVIATION_RANDOM_UPPER);
   } /* outerrounds */

   *nelims = ndeletions;

   if( rpc )
      graph_pc_2orgcheck(scip, graph);

   if( SCIPisLT(scip, upperbound, *ub) || SCIPisLT(scip, *ub, 0.0) )
      *ub = upperbound;

   if( userec )
      solpool_free(scip, &pool);

   SCIPfreeBufferArray(scip, &terms);
   SCIPfreeBufferArray(scip, &nodereplacebounds);
   SCIPfreeBufferArray(scip, &nodefixingbounds);
   SCIPfreeBufferArray(scip, &edgefixingbounds);
   SCIPfreeBufferArray(scip, &arcsdeleted);
   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &pseudoDelNodes);

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}


/** dual ascent reduction for slack-and-prune heuristic */
SCIP_RETCODE reduce_daSlackPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< array to store reduced costs */
   SCIP_Real*            costrev,            /**< reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   SCIP_Real*            upperbound,         /**< pointer to store new upper bound */
   int*                  edgearrint,         /**< int edges array to store solution value  */
   int*                  edgearrint2,        /**< int edges array for internal computations */
   int*                  vbase,              /**< array for Voronoi bases */
   int*                  pathedge,           /**< array for predecessor edge on a path */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  solnode,            /**< array of nodes of current solution that is not to be destroyed */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations  */
   STP_Bool*             edgearrchar,        /**< STP_Bool edge array for internal computations  */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   int                   minelims,           /**< minimum number of edges to eliminate */
   SCIP_Bool             solgiven            /**< solution provided? */
   )
{
   IDX** ancestors;
   IDX** revancestors;
   SCIP_Real obj;
   SCIP_Real tmpcost;
   SCIP_Real lpobjval;
   SCIP_Real objprune;
   SCIP_Real minpathcost;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   SCIP_Bool rpc;
   SCIP_Bool success;
   SCIP_Bool eliminate;

   int* grad;
   int* adjvert;
   int* incedge;
   int* reinsert;
   int i;
   int k;
   int e;
   int e2;
   int e3;
   int etmp;
   int root;
   int head;
   int nterms;
   int nedges;
   int nnodes;
   int nfixed;
   int redrounds;
   STP_Bool* marked;

   assert(scip != NULL);
   assert(cost != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);
   assert(solnode != NULL);
   assert(costrev != NULL);
   assert(pathedge != NULL);
   assert(upperbound != NULL);
   assert(edgearrint != NULL);
   assert(edgearrint2 != NULL);
   assert(nodearrchar != NULL);
   assert(edgearrchar != NULL);
   assert(!graph_pc_isPcMw(graph));

   /* 1. step: initialize */

   rpc = (graph->stp_type == STP_RPCSPG);
   grad = graph->grad;
   root = graph->source;
   nfixed = 0;
   nedges = graph->edges;
   nnodes = graph->knots;

   /* graph vanished? */
   if( grad[graph->source] == 0 )
      return SCIP_OKAY;

   marked = edgearrchar;

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

   k = 0;
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( !rpc )
         graph->mark[i] = (grad[i] > 0);
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

   /* 2. step: - if not provided, compute lower bound and reduced costs
    *          - try to eliminate edges and nodes                        */

   for( i = 0; i < nnodes; i++ )
      if( Is_term(graph->term[i]) )
         assert(grad[i] > 0);

   if( rpc )
   {
      graph_pc_2trans(scip, graph);
   }


   if( !solgiven )
   {
      /* try to build MST on solnode nodes */
      for( i = 0; i < nnodes; i++ )
         graph->mark[i] = (solnode[i] == CONNECT);

      for( e = 0; e < nedges; e++ )
         edgearrint[e] = UNKNOWN;

      graph_path_exec(scip, graph, MST_MODE, root, graph->cost, vnoi);

      for( i = 0; i < nnodes; i++ )
      {
         e = vnoi[i].edge;
         if( e >= 0 )
         {
            edgearrint[e] = CONNECT;
         }
         else if( Is_term(graph->term[i]) && i != root )
         {
            break;
         }
      }

      if( i == nnodes )
      {
         int l;
         int count;

         do
         {
            count = 0;

            for( l = nnodes - 1; l >= 0; --l )
            {
               if( (solnode[l] != CONNECT) || Is_term(graph->term[l]) )
                  continue;

               for( e = graph->outbeg[l]; e != EAT_LAST; e = graph->oeat[e] )
                  if( edgearrint[e] == CONNECT )
                     break;

               if( e == EAT_LAST )
               {
                  /* there has to be exactly one incoming edge */
                  for( e = graph->inpbeg[l]; e != EAT_LAST; e = graph->ieat[e] )
                  {
                     if( edgearrint[e] == CONNECT )
                     {
                        edgearrint[e] = UNKNOWN;
                        solnode[l] = UNKNOWN;
                        count++;
                        break;
                     }
                  }
               }
            }
         }
         while( count > 0 );
      }
   }


   if( solgiven || i == nnodes )
   {
      obj = getSolObj(scip, graph, edgearrint);

      SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, &lpobjval, FALSE, FALSE, gnodearr, edgearrint, edgearrint2, state, root, FALSE, -1.0) );
   }
   else
   {
      obj = FARAWAY;
      SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, &lpobjval, FALSE, FALSE, gnodearr, NULL, edgearrint2, state, root, FALSE, -1.0) );
   }

#if 0
      SCIP_QUEUE* queue;
      SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );

      GRAPH* prunegraph = graph;
      int* mark = graph->mark;
      int* pnode;
      int a;
      for( k = 0; k < nnodes; k++ )
         mark[k] = FALSE;


      /* BFS from root along incoming arcs of zero cost */

      mark[prunegraph->source] = TRUE;

      SCIP_CALL( SCIPqueueInsert(queue, &(prunegraph->source)) );

      while( !SCIPqueueIsEmpty(queue) )
      {
         pnode = (SCIPqueueRemove(queue));
         k = *pnode;

         /* traverse outgoing arcs */
         for( a = prunegraph->outbeg[k]; a != EAT_LAST; a = prunegraph->oeat[a] )
         {
            head = prunegraph->head[a];

            if( SCIPisEQ(scip, cost[a], 0.0) )
            {
               /* vertex not labeled yet? */
               if( !mark[head] )
               {
                  mark[head] = TRUE;
                  SCIP_CALL( SCIPqueueInsert(queue, &(prunegraph->head[a])) );
               }
            }
         }
      }
      SCIPqueueFree(&queue);
      for( k = 0; k < nnodes; k++ )
         if( Is_term(prunegraph->term[k]) && !mark[k]  )
            printf("in bnd  FAIL %d not marked, but terminal, \n", k);
#endif

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, edgearrint2, vbase, root, nodearrchar, &success, FALSE) );

   objprune = getSolObj(scip, graph, edgearrint2);

   assert(success);

   if( success && SCIPisLT(scip, objprune, obj ) )
   {
      for( i = 0; i < nnodes; i++ )
         solnode[i] = UNKNOWN;

      for( e = 0; e < nedges; e++ )
      {
         edgearrint[e] = edgearrint2[e];
         if( edgearrint[e] == CONNECT )
         {
            solnode[graph->tail[e]] = CONNECT;
            solnode[graph->head[e]] = CONNECT;
         }
      }
   }

   obj = 0.0;

   for( e = 0; e < nedges; e++ )
   {
      if( edgearrint[e] == CONNECT )
         obj += graph->cost[e];

      marked[e] = FALSE;
      costrev[e] = cost[flipedge(e)];
   }

   *upperbound = obj;

   for( k = 0; k < nnodes; k++ )
      graph->mark[k] = (grad[k] > 0);

   /* distance from root to all nodes */
   graph_path_execX(scip, graph, root, cost, pathdist, pathedge);

   /* no paths should go back to the root */
   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build Voronoi diagram */
   graph_get4nextTerms(scip, graph, costrev, costrev, vnoi, vbase, graph->path_heap, state);

   for( k = 0; k < nnodes; k++ )
      if( !Is_term(graph->term[k]) )
         assert(vbase[k + nnodes] != root );

   /* RPC? If yes, restore original graph */
   if( rpc )
   {
      graph_pc_2org(scip, graph);
      graph->mark[root] = FALSE;
   }

   for( e = 0; e < nedges; e++ )
      costrev[e] = -1.0;

   for( redrounds = 0; redrounds < 3; redrounds++ )
   {
      if( redrounds == 0 )
      {
         eliminate = FALSE;
         minpathcost = FARAWAY;
      }
      else if( redrounds == 1 )
      {
         assert(minelims > 0);
         assert(2 * minelims < nedges);

         eliminate = TRUE;
         SCIPsortReal(costrev, nedges);

         /* the required reduced path cost to be surpassed */
         minpathcost = costrev[nedges - 2 * minelims];

         if( SCIPisLE(scip, minpathcost, 0.0) )
            minpathcost = 2 * SCIPepsilon(scip);

         k = nedges - 2 * minelims;

         /* try to first eliminate edges with higher gap */
         for( e = nedges - 1; e > k && e >= 2; e = e - 2 )
         {
            if( SCIPisLE(scip, costrev[e - 2], minpathcost) )
               break;
         }

         if( SCIPisGT(scip, costrev[e], minpathcost) )
            minpathcost = costrev[e];

         if( SCIPisLE(scip, minpathcost, 0.0) )
            minpathcost = 2 * SCIPepsilon(scip);

         for( e = 0; e < nedges; e++ )
            marked[e] = FALSE;
      }
      else
      {
         eliminate = TRUE;

         /* the required reduced path cost to be surpassed */
         minpathcost = costrev[nedges - 2 * minelims];

         if( SCIPisLE(scip, minpathcost, 0.0) )
            minpathcost = 2 * SCIPepsilon(scip);

         for( e = 0; e < nedges; e++ )
            marked[e] = FALSE;
      }

      for( k = 0; k < nnodes; k++ )
      {
         if( grad[k] <= 0 )
            continue;

         if( nfixed > minelims )
            break;

         if( !Is_term(graph->term[k]) && (!eliminate || pathdist[k] + vnoi[k].dist >= minpathcost) && solnode[k] != CONNECT  )
         {
            if( !eliminate )
            {
               tmpcost = pathdist[k] + vnoi[k].dist;

               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
                  if( SCIPisGT(scip, tmpcost, costrev[e]) )
                     costrev[e] = tmpcost;

                  e2 = flipedge(e);

                  if( SCIPisGT(scip, tmpcost, costrev[e2]) )
                     costrev[e2] = tmpcost;
               }

               continue;
            }
            nfixed += grad[k];

            while( graph->outbeg[k] != EAT_LAST )
               graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
         }
         else
         {
            e = graph->outbeg[k];
            while( e != EAT_LAST )
            {
               etmp = graph->oeat[e];
               head = graph->head[e];

               /* for rpc no artificial terminal arcs should be deleted; in general: delete no solution edges */
               if( (rpc && !graph->mark[head])
                  || (edgearrint[e] == CONNECT) || (edgearrint[flipedge(e)] == CONNECT) )
               {
                  e = etmp;
                  continue;
               }

               tmpcost = pathdist[k] + cost[e] + vnoi[head].dist;

               if( (!eliminate) || tmpcost >= minpathcost )
               {
                  if( marked[flipedge(e)] )
                  {
                     if( eliminate )
                     {
                        graph_edge_del(scip, graph, e, TRUE);
                        nfixed++;
                     }
                     else
                     {
                        if( SCIPisGT(scip, tmpcost, costrev[e]) )
                           costrev[e] = tmpcost;

                        if( SCIPisGT(scip, tmpcost, costrev[flipedge(e)]) )
                           costrev[flipedge(e)] = tmpcost;
                     }
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
         if( nfixed > minelims )
            break;

         if( !graph->mark[k] || Is_term(graph->term[k]) || solnode[k] == CONNECT )
            continue;

         if( grad[k] == 3 )
         {
            tmpcost = pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist;

            if( !eliminate || tmpcost >= minpathcost )
            {
               e = graph->outbeg[k];
               assert(graph->oeat[e] != EAT_LAST);
               e2 = graph->oeat[e];
               assert(graph->oeat[e2] != EAT_LAST);
               e3 = graph->oeat[e2];
               assert(graph->oeat[e3] == EAT_LAST);

               if( SCIPisLE(scip, cost[e], 0.0) || SCIPisLE(scip, cost[e2], 0.0) || SCIPisLE(scip, cost[e3], 0.0) )
                  continue;

               if( !eliminate )
               {
                  for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     if( SCIPisGT(scip, tmpcost, costrev[e]) )
                        costrev[e] = tmpcost;

                     e2 = flipedge(e);

                     if( SCIPisGT(scip, tmpcost, costrev[e2]) )
                        costrev[e2] = tmpcost;
                  }

                  continue;
               }
               nfixed += 3;
               while( graph->outbeg[k] != EAT_LAST )
                  graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
            }
         }
      }
   }
   SCIPdebugMessage("deleted by da: %d \n", nfixed );

   if( rpc )
      graph_pc_2trans(scip, graph);
   assert(graph->mark[root]);

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

   if( edgearrchar == NULL )
      SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}


/** dual ascent based reductions for PCSPG and MWCSP */
SCIP_RETCODE reduce_daPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const RPDA*           paramsda,           /**< parameters */
   PATH*                 vnoi,               /**< Voronoi data structure array */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   int*                  vbase,              /**< Voronoi base array */
   int*                  pathedge,           /**< shortest path incoming edge array for shortest path calculations */
   int*                  state,              /**< int 4 * vertices array  */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             prizesum            /**< sum of positive prizes */
)
{
   STPSOLPOOL* pool = NULL;
   GRAPH* transgraph = NULL;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real* bestcost = NULL;
   SCIP_Real* edgefixingbounds = NULL;
   SCIP_Real* nodefixingbounds = NULL;
   SCIP_Real offset;
   SCIP_Real lpobjval;
   SCIP_Real bestlpobjval;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   int* roots = NULL;
   int* result = NULL;
   int* result2 = NULL;
   int* transresult = NULL;
   STP_Bool* marked = NULL;
   int nroots = 0;
   int nfixed = 0;
   int nusedroots;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   const int extnedges = nedges + 2 * (graph->terms - 1);
   const int root = graph->source;
   SCIP_Bool havenewsol;
   SCIP_Bool success;
   SCIP_Bool userec = paramsda->useRec;
   const SCIP_Bool solbasedda = paramsda->pcmw_solbasedda;
   const SCIP_Bool useDifferentRoots = paramsda->pcmw_useMultRoots;
   const SCIP_Bool markroots = paramsda->pcmw_markroots;
   SCIP_Bool varyroot = useDifferentRoots && userec;
   const SCIP_Real damaxdeviation = paramsda->pcmw_fastDa ? DAMAXDEVIATION_FAST : -1.0;

   assert(scip && nelims && nodearrchar);

   if( graph->terms <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transresult, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result2, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcost, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgefixingbounds, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodefixingbounds, nnodes + 1) );

   if( userec )
      SCIP_CALL( solpool_init(scip, &pool, nedges, SOLPOOL_SIZE) );

   if( graph_pc_isPc(graph) )
      reduce_identifyNonLeafTerms(scip, graph);

   /*
    * 1. step: compute lower bound and reduced costs
    */

   graph_pc_2trans(scip, graph);
   offset = 0.0;

   /* transform the problem to a real SAP */
   SCIP_CALL( graph_transPcGetSap(scip, graph, &transgraph, &offset) );

   /* initialize data structures for shortest paths */
   SCIP_CALL( graph_path_init(scip, transgraph) );

   SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, &lpobjval, FALSE, FALSE,
         gnodearr, NULL, transresult, state, root, TRUE, damaxdeviation) );

   lpobjval += offset;
   bestlpobjval = lpobjval;
   BMScopyMemoryArray(bestcost, cost, extnedges);

   /* compute first primal solution */
   upperbound = FARAWAY;
   havenewsol = FALSE;
   SCIP_CALL( computeSteinerTreeRedCostsPcMw(scip, graph, NULL, cost, &upperbound, result, result2, vbase, pathedge, nodearrchar, &havenewsol) );

   assert(havenewsol && upperbound < FARAWAY);
   assert(graph_solIsValid(scip, graph, result));

   /* the required reduced path cost to be surpassed */
   minpathcost = upperbound - lpobjval;
   if( userec)
      SCIPdebugMessage("DA first minpathcost %f \n", minpathcost);

   /* initialize data structures for transgraph */
   SCIP_CALL( graph_init_history(scip, transgraph) );
   computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);

   /*
    * 2. step: try to eliminate
    */

   /* restore original graph */
   graph_pc_2org(scip, graph);

   for( int e = 0; e < extnedges; e++ )
      marked[e] = FALSE;

   /* try to reduce the graph */
   assert(daRedCostIsValid(scip, transgraph, cost, state, nodearrchar));

   updateEdgeFixingBounds(edgefixingbounds, graph, cost, pathdist, vnoi, lpobjval, extnedges, TRUE, FALSE);
   updateNodeFixingBounds(nodefixingbounds, graph, pathdist, vnoi, lpobjval, TRUE);

   nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, TRUE);

   assert(!graph->extended);

   /* edges from result might have been deleted! */
   havenewsol = havenewsol && graph_solIsUnreduced(scip, graph, result);
   assert(!havenewsol || graph_solIsValid(scip, graph, result));

   if( userec )
      SCIPdebugMessage("DA: 1. NFIXED %d \n", nfixed);

   /* rerun dual ascent? */
   if( solbasedda && graph->terms > 2 && SCIPisGT(scip, minpathcost, 0.0) )
   {
      const SCIP_Real oldupperbound = upperbound;

      /* with recombination? */
      if( userec && graph->stp_type != STP_MWCSP )
      {
         SCIP_Real ub;

         /* compute second solution and add to pool */
         SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fromheurdata,
            graph, NULL, NULL, result2, BND_TMHEUR_NRUNS / 5, root, graph->cost, graph->cost, NULL, NULL, &success) );
         assert(success);

         SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, result2) );
         ub = getSolObj(scip, graph, result2);

         SCIP_CALL( solpool_addSol(scip, ub, result2, pool, &success) );
         SCIPdebugMessage("added initial TM sol to pool? %d , ub %f \n", success, ub);
      }

      /* try to improve both dual and primal bound */
      SCIP_CALL( computePertubedSol(scip, graph, transgraph, pool, vnoi, gnodearr, cost, bestcost, pathdist, state, vbase, pathedge, result, result2,
            transresult, nodearrchar, &upperbound, &lpobjval, &bestlpobjval, &minpathcost, &havenewsol, offset, extnedges, 0) );

      assert(graph_solIsValid(scip, graph, result));
      assert(!havenewsol || SCIPisEQ(scip, getSolObj(scip, graph, result), upperbound));

      graph_pc_2orgcheck(scip, graph);

      assert(daRedCostIsValid(scip, transgraph, cost, state, nodearrchar));
      computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);
      updateEdgeFixingBounds(edgefixingbounds, graph, cost, pathdist, vnoi, lpobjval, extnedges, FALSE, FALSE);
      updateNodeFixingBounds(nodefixingbounds, graph, pathdist, vnoi, lpobjval, FALSE);

      nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, havenewsol);

      nfixed += reducePcMwTryBest(scip, graph, transgraph, vnoi, cost, costrev, bestcost, pathdist, &upperbound,
            &lpobjval, &bestlpobjval, &minpathcost, oldupperbound, result, vbase, state, pathedge, marked, nodearrchar, &havenewsol, extnedges);

      nfixed += reduceWithEdgeFixingBounds(scip, graph, transgraph, edgefixingbounds, NULL, upperbound);
      nfixed += reduceWithNodeFixingBounds(scip, graph, transgraph, nodefixingbounds, upperbound);

      // todo test again once the extended reductions are more mature
#if 0
      if( paramsda->useExtRed && !SCIPisZero(scip, minpathcost) )
      {
         int extfixed;

         REDCOST redcostdata = { .redEdgeCost = cost, .rootToNodeDist = pathdist,
                  .nodeTo3TermsPaths = vnoi, .nodeTo3TermsBases = vbase,
                  .cutoff = minpathcost, dualBound = -1.0, .redCostRoot = graph->source };

         havenewsol = havenewsol && graph_solIsUnreduced(scip, graph, result);

         SCIP_CALL( extreduce_deleteEdges(scip, &redcostdata, (havenewsol ? result : NULL), graph, marked, &extfixed) );
         nfixed += extfixed;

         printf("extfixed=%d \n", extfixed);
      }
#endif

      if( userec )
         SCIPdebugMessage("eliminations after sol based run2 with best dual sol %d bestlb %f newlb %f\n", nfixed, bestlpobjval, lpobjval);
   }

   /* pertubation runs for MWCSP */
   if( varyroot && graph->stp_type == STP_MWCSP )
   {
      for( int run = 0; run < DEFAULT_NMAXROOTS && graph->terms > STP_RED_MINBNDTERMS; run++ )
      {
         SCIP_Real oldupperbound = upperbound;

         graph_pc_2trans(scip, graph);

         havenewsol = havenewsol && graph_solIsUnreduced(scip, graph, result);
         assert(!havenewsol || graph_solIsValid(scip, graph, result));

         assert(SCIPisEQ(scip, upperbound, getSolObj(scip, graph, result)));

         /* try to improve both dual and primal bound */
         SCIP_CALL( computePertubedSol(scip, graph, transgraph, pool, vnoi, gnodearr, cost, bestcost, pathdist, state, vbase, pathedge, result, result2,
               transresult, nodearrchar, &upperbound, &lpobjval, &bestlpobjval, &minpathcost, &havenewsol, offset, extnedges, run) );

         SCIPdebugMessage("DA: pertubated run %d ub: %f \n", run, upperbound);
         SCIPdebugMessage("DA: pertubated run %d minpathcost: %f \n", run, upperbound - lpobjval);

         computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);
         updateEdgeFixingBounds(edgefixingbounds, graph, cost, pathdist, vnoi, lpobjval, extnedges, FALSE, FALSE);
         updateNodeFixingBounds(nodefixingbounds, graph, pathdist, vnoi, lpobjval, FALSE);

         nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, havenewsol);

         nfixed += reducePcMwTryBest(scip, graph, transgraph, vnoi, cost, costrev, bestcost, pathdist, &upperbound,
               &lpobjval, &bestlpobjval, &minpathcost, oldupperbound, result, vbase, state, pathedge, marked, nodearrchar, &havenewsol, extnedges);

         nfixed += reduceWithEdgeFixingBounds(scip, graph, transgraph, edgefixingbounds, NULL, upperbound);
         nfixed += reduceWithNodeFixingBounds(scip, graph, transgraph, nodefixingbounds, upperbound);

         SCIPdebugMessage("DA: pertubated run %d NFIXED %d \n", run, nfixed);
      }
   }

   if( graph->stp_type == STP_MWCSP && graph->terms < STP_RED_MINBNDTERMS  )
      varyroot = FALSE;

   /* change or mark roots? */
   if( varyroot || markroots )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &roots, graph->terms) );

      SCIP_CALL( daPcFindRoots(scip, graph, transgraph, cost, bestcost, lpobjval, bestlpobjval, upperbound, FALSE, FALSE,
            vnoi, pathedge, vbase, nodearrchar, roots, &nroots));

      /* should prize of terminals be changed? */
      if( nroots > 0 && markroots  )
         daPcMarkRoots(scip, roots, 0, nroots, prizesum, graph, &userec, &pool);

      if( nroots > 0 && varyroot )
         SCIP_CALL( daOrderRoots(scip, graph, roots, nroots, TRUE, randnumgen) );
   }

   if( varyroot )
      nusedroots = MIN(DEFAULT_NMAXROOTS, nroots);
   else
      nusedroots = -1;

   graph_path_exit(scip, transgraph);
   graph_free(scip, &transgraph, TRUE);

   /* loop and change root for dual ascent run */
   for( int run = 0; run < nusedroots; run++  )
   {
      const int tmproot = roots[run];
      int transnnodes;
      int transnedges;

      assert(nroots > 0);
      assert(roots != NULL);
      assert(Is_term(graph->term[tmproot]));

      if( graph->terms <= 2 )
         break;

      if( graph_pc_knotIsNonLeafTerm(graph, tmproot) )
      {
         assert(0); // might actually happen!
      }

      SCIP_CALL( graph_transPcGetRsap(scip, graph, &transgraph, roots, nroots, tmproot) );

      assert(graph_valid(scip, transgraph) && STP_SAP == transgraph->stp_type);

      transnnodes = transgraph->knots;
      transnedges = transgraph->edges;

      for( int k = 0; k < transnnodes; k++ )
         transgraph->mark[k] = (transgraph->grad[k] > 0);

      /* init data structures for shortest paths and history */
      SCIP_CALL( graph_path_init(scip, transgraph) );
      SCIP_CALL( graph_init_history(scip, transgraph ) );

      if( havenewsol && run > 1 )
      {
         BMScopyMemoryArray(transresult, result, graph->edges);
         SCIP_CALL(graph_solReroot(scip, transgraph, transresult, tmproot));
         SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, &lpobjval, FALSE, FALSE,
               gnodearr, transresult, result2, state, tmproot, FALSE, -1.0));
      }
      else
      {
         SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, &lpobjval, FALSE, FALSE,
               gnodearr, NULL, transresult, state, tmproot, FALSE, -1.0));
      }

      assert(graph_valid(scip, transgraph));

      for( int e = graph->outbeg[tmproot]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int k = graph->head[e];
         if( Is_term(graph->term[k]) )
         {
            if( k == root )
               cost[flipedge(e)] = 0.0;
            else
               cost[e] = 0.0;
         }
      }

      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(graph->term[k]) && SCIPisGE(scip, graph->prize[k], prizesum) )
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            {
               const int head = graph->head[e];
               if( Is_term(graph->term[head]) && head != root )
                  cost[e] = 0.0;
            }
         }
      }

      havenewsol = FALSE;
      SCIP_CALL( computeSteinerTreeRedCostsPcMw(scip, graph, pool, cost, &upperbound, result, result2, vbase, pathedge, nodearrchar, &havenewsol) );

      SCIPdebugMessage("ROOTRUNS upperbound %f \n", upperbound);
      if( pool )
         SCIPdebugMessage("ROOTRUNS best sol in pool %f \n", pool->sols[0]->obj);

      for( int k = 0; k < transnnodes; k++ )
         transgraph->mark[k] = (transgraph->grad[k] > 0);

      /* the required reduced path cost to be surpassed */
      minpathcost = upperbound - lpobjval;

      if( markroots )
      {
         const int oldnroots = nroots;
         SCIP_CALL( daPcFindRoots(scip, graph, transgraph, cost, cost, lpobjval, lpobjval, upperbound, TRUE, TRUE,
               vnoi, pathedge, vbase, nodearrchar, roots, &nroots) );

         /* should prize of terminals be changed? */
         if( nroots > oldnroots  )
            daPcMarkRoots(scip, roots, oldnroots, nroots, prizesum, graph, &userec, &pool);
      }

      SCIPdebugMessage("ROOTRUNS: minpathcost %f \n", minpathcost);
      SCIPdebugMessage("lb: %f ub: %f \n", lpobjval, upperbound);

      /* distance from root to all nodes */
      graph_path_execX(scip, transgraph, tmproot, cost, pathdist, pathedge);

      for( int e = 0; e < transnedges; e++ )
            costrev[e] = cost[flipedge(e)];

      /* no paths should go back to the root */
      for( int e = transgraph->outbeg[tmproot]; e != EAT_LAST; e = transgraph->oeat[e] )
         costrev[e] = FARAWAY;

      /* build Voronoi diagram */
      graph_voronoiTerms(scip, transgraph, costrev, vnoi, vbase, transgraph->path_heap, state);
      graph_get2next(scip, transgraph, costrev, costrev, vnoi, vbase, transgraph->path_heap, state);

      /* restore original graph */
      graph_pc_2org(scip, graph);

      assert(graph->mark[tmproot]);
      graph->mark[tmproot] = FALSE;

      /* at first run switch to undirected case */
      if( run == 0 )
         for( int e = 0; e < extnedges; e++ )
            edgefixingbounds[e] = MIN(edgefixingbounds[e], edgefixingbounds[flipedge(e)]);

      updateEdgeFixingBounds(edgefixingbounds, graph, cost, pathdist, vnoi, lpobjval, extnedges, FALSE, TRUE);
      updateNodeFixingBounds(nodefixingbounds, graph, pathdist, vnoi, lpobjval, FALSE);

      for( int e = 0; e < transnedges; e++ )
         marked[e] = FALSE;

      /* try to eliminate vertices and edges */
      nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, havenewsol);
      nfixed += reduceWithEdgeFixingBounds(scip, graph, transgraph, edgefixingbounds, NULL, upperbound);
      nfixed += reduceWithNodeFixingBounds(scip, graph, transgraph, nodefixingbounds, upperbound);

      havenewsol = havenewsol && graph_solIsUnreduced(scip, graph, result);
      assert(!havenewsol || graph_solIsValid(scip, graph, result));
      SCIPdebugMessage("FIXED with changed root %d \n\n", nfixed);

      graph->mark[tmproot] = TRUE;

      transgraph->stp_type = STP_RPCSPG;
      graph_path_exit(scip, transgraph);
      graph_free(scip, &transgraph, TRUE);
   }

   *nelims = nfixed;

   if( pool != NULL )
      solpool_free(scip, &pool);

   SCIPfreeBufferArrayNull(scip, &roots);
   SCIPfreeBufferArray(scip, &nodefixingbounds);
   SCIPfreeBufferArray(scip, &edgefixingbounds);
   SCIPfreeBufferArray(scip, &bestcost);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &result2);
   SCIPfreeBufferArray(scip, &marked);
   SCIPfreeBufferArray(scip, &transresult);
   SCIPfreeBufferArray(scip, &result);

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}
