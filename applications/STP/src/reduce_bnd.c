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
#include "heur_ascendprune.h"
#include "cons_stp.h"
#include "heur_local.h"
#include "misc_stp.h"
#include "prop_stp.h"
#include "probdata_stp.h"
#include "heur_rec.h"
#include "heur_slackprune.h"

#define DEFAULT_HEURRUNS 100                  /**< number of runs of constructive heuristic */
#define DEFAULT_DARUNS     10                 /**< number of runs for dual ascent heuristic */
#define DEFAULT_NMAXROOTS  8                  /**< max number of roots to use for new graph in dual ascent heuristic */
#define PERTUBATION_RATIO   0.05              /**< pertubation ratio for dual-ascent primal bound computation */
#define PERTUBATION_RATIO_PC   0.005          /**< pertubation ratio for dual-ascent primal bound computation */
#define SOLPOOL_SIZE 20                       /**< size of presolving solution pool */
#define STP_RED_MINBNDTERMS   500



/** set prize of marked terminals to blocked */
static
void markPcMwRoots(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const int*            roots,              /**< root array */
   int                   nroots,             /**< number of roots */
   SCIP_Real             prizesum            /**< sum of positive prizes */
)
{
   const int root = graph->source;

   assert(graph != NULL);
   assert(roots != NULL);
   assert(!graph->extended);
   assert(nroots >= 0);

   for( int i = 0; i < nroots; i++ )
   {
      int e;
      const int term = roots[i];
      const int pterm = graph->head[graph->term2edge[term]];

      assert(Is_term(graph->term[term]));
      assert(SCIPisPositive(scip, graph->prize[term]));
      assert(pterm >= 0);
      assert(Is_pterm(graph->term[pterm]));

      for( e = graph->inpbeg[pterm]; e != EAT_LAST; e = graph->ieat[e] )
         if( root == graph->tail[e] )
            break;

      assert(e != EAT_LAST);
      assert(SCIPisEQ(scip, graph->prize[term], graph->cost[e]));

      graph->prize[term] = prizesum;
      graph->cost[e] = prizesum;
   }
}

/** are the dual costs still valid, i.e. are there zero cost paths from the root to all terminals? */
static
SCIP_Bool dualCostIsValid(
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
void pertubateEdgeCosts(
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
   srand(graph->terms);

   if( graph->stp_type != STP_MWCSP )
   {
      SCIP_Real pratio = PERTUBATION_RATIO_PC;

      for( int k = 0; k < nnodes; k++ )
      {
         assert(Is_gterm(graph->term[k]) == Is_gterm(transgraph->term[k]) || transgraph->grad[k] == 0);

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

         if( !Is_gterm(graph->term[k]) )
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
                  assert(!Is_pterm(transgraph->term[transgraph->tail[e]]));
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

      assert(Is_gterm(graph->term[k]) == Is_gterm(transgraph->term[k]));

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

      if( !Is_gterm(graph->term[k]) )
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
                  assert(!Is_pterm(transgraph->term[transgraph->tail[e]]));

                  transgraph->cost[e] *= 1.0 + pratio;
               }
         }
         else
         {
            for( e = transgraph->inpbeg[k]; e != EAT_LAST; e = transgraph->ieat[e] )
               if( SCIPisPositive(scip, transgraph->cost[e]) )
               {
                  assert(!Is_pterm(transgraph->term[transgraph->tail[e]]));
                  transgraph->cost[e] *= 1.0 - pratio;
               }
         }
      }
   }
}

/** order roots */
static
void orderDaRoots(
   int*                  terms,              /**< sol int array corresponding to upper bound */
   int*                  termdegs,           /**< sol int array */
   int                   nterms,             /**< number of terminals */
   int                   maxdeg,             /**< max degree */
   SCIP_Bool             randomize,          /**< randomize */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator (or NULL) */
)
{
   assert(terms != NULL);
   assert(termdegs != NULL);

   for( int i = 0; i < nterms; i++ )
   {
      termdegs[i] *= -1;
      if( randomize )
         termdegs[i] -= SCIPrandomGetInt(randnumgen, 0, maxdeg);
   }

   SCIPsortIntInt(termdegs, terms, nterms);
}


/** compute primal solution during dual-ascent routine for PCSTP or MWCSP */
static
SCIP_RETCODE computeDaSolPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   STPSOLPOOL*           pool,               /**< solution pool */
   PATH*                 vnoi,               /**< Voronoi data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   SCIP_Real*            pathdist,           /**< distance array from shortest path calculations */
   SCIP_Real*            upperbound,         /**< upperbound pointer */
   int*                  result1,            /**< sol int array corresponding to upper bound */
   int*                  result2,            /**< sol int array */
   int*                  vbase,              /**< int array */
   int*                  pathedge,           /**< int array */
   int                   root,               /**< da root */
   STP_Bool*             nodearrchar,        /**< node array storing solution vertices */
   SCIP_Bool*            apsol               /**< ascend-prune sol? */
)
{
   SCIP_Real ub;
   const int nedges = graph->edges;
   SCIP_Bool success;

   /* compute new solution and store it in result2 */

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, result2, vbase, -1, nodearrchar, &success, FALSE) );
   assert(success);

   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result2) );

   assert(graph_sol_valid(scip, graph, result2));

   /* compute objective value */
   ub = graph_sol_getObj(graph->cost, result2, 0.0, nedges);
   SCIPdebugMessage("DA: first new sol value in computeDaSolPcMw: %f ... old value: %f \n", ub, *upperbound);

   /* try recombination? */
   if( pool != NULL )
   {
      SCIPdebugMessage("ub %f vs best sol %f\n", ub, pool->sols[0]->obj);
      SCIP_CALL( SCIPStpHeurRecAddToPool(scip, ub, result2, pool, &success) );

#ifdef SCIP_DEBUG
      for( int i = 0; i < pool->size; i++ )
      {
         printf(" %f ", pool->sols[i]->obj);
      }
      printf("\n ");
#endif

      SCIPdebugMessage("DA: added this sol to pool? %d \n", success);

      if( success && pool->size >= 2 )
      {
         /* get index of just added solution */
         int solindex = pool->maxindex;

         SCIP_Bool solfound;

         SCIPdebugMessage("POOLSIZE %d \n", pool->size);

         SCIP_CALL( SCIPStpHeurRecRun(scip, pool, NULL, NULL, graph, NULL, &solindex, 3, pool->size, FALSE, &solfound) );

         if( solfound )
         {
            STPSOL* sol = SCIPStpHeurRecSolfromIdx(pool, solindex);

            assert(pool->size >= 2);
            assert(sol != NULL);

            SCIPdebugMessage("DA: rec found better solution with obj %f vs %f \n", sol->obj, ub);

            if( SCIPisLE(scip, sol->obj, ub) )
            {
               BMScopyMemoryArray(result2, sol->soledges, nedges);

               SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result2) );
               ub = graph_sol_getObj(graph->cost, result2, 0.0, nedges);

               if( SCIPisLT(scip, ub, sol->obj) )
                  SCIP_CALL( SCIPStpHeurRecAddToPool(scip, ub, result2, pool, &success) );
            }
         }
      }
   }

   if( SCIPisLE(scip, ub, *upperbound) )
   {
      SCIPdebugMessage("DA: improved incumbent %f vs %f, return \n", ub, *upperbound);

      *apsol = TRUE;
      *upperbound = ub;
      BMScopyMemoryArray(result1, result2, nedges);
   }

   if( graph->stp_type != STP_MWCSP || !(*apsol) )
     return SCIP_OKAY;

   SCIP_CALL( SCIPStpHeurRecExclude(scip, graph, result1, result2, pathedge, nodearrchar, &success) );

   if( success )
   {
      BMScopyMemoryArray(result1, result2, nedges);
      *upperbound = graph_sol_getObj(graph->cost, result1, 0.0, nedges);
      SCIPdebugMessage("DA: afterLastExclusion %f \n", *upperbound);
   }

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
   SCIP_Real* costrev,
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

   graph_pc_2transcheck(graph);

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

         SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result2) );

         assert(graph_sol_valid(scip, graph, result2));

         bound = graph_sol_getObj(graph->cost, result2, 0.0, nedges);

         if( SCIPisLE(scip, bound, *upperbound) )
         {
            *upperbound = bound;
            *apsol = TRUE;
            BMScopyMemoryArray(result, result2, nedges);
         }
         pertubateEdgeCosts(scip, graph, transgraph, result2, nodearrchar, pertubation);
      }
      else
      {
         pertubateEdgeCosts(scip, graph, transgraph, result, nodearrchar, pertubation);
      }

      /* todo use result as guiding solution? */
      SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lb, FALSE, FALSE, gnodearr, NULL, transresult, state, root, 1, nodearrchar) );

      BMScopyMemoryArray(transgraph->cost, transcost, transnedges);

      SCIPfreeBufferArray(scip, &transcost);
   }

   SCIP_CALL( computeDaSolPcMw(scip, graph, pool, vnoi, cost, pathdist, upperbound, result, result2, vbase, pathedge, root, nodearrchar, apsol) );

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

   assert(!(*apsol) || graph_sol_valid(scip, graph, result));
   assert(graph_valid(transgraph));
   assert(root == transgraph->source);

   SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lb, FALSE, FALSE, gnodearr, transresult, NULL, state, root, 1, nodearrchar) );

   lb += offset;
   *lpobjval = lb;

   SCIP_CALL( computeDaSolPcMw(scip, graph, pool, vnoi, cost, pathdist, upperbound, result, result2, vbase, pathedge, root, nodearrchar, apsol) );

   assert(!(*apsol) || graph_sol_valid(scip, graph, result));

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
   int nnodes;
   int nfixed;
   SCIP_Real tmpcost;
   SCIP_Bool keepsol = FALSE;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(pathdist != NULL);
   assert(result != NULL);
   assert(cost != NULL);
   assert(vnoi != NULL);

   assert(SCIPisGE(scip, minpathcost, 0.0));
   assert(!solgiven || graph_sol_valid(scip, graph, result));

   if( minpathcost < 0.0 )
      minpathcost = 0.0;

   nfixed = 0;
   nnodes = graph->knots;

   graph_pc_2orgcheck(graph);

   if( solgiven )
   {
      const int nedges = graph->edges;

      for( int k = 0; k < nnodes; k++ )
         nodearrchar[k] = FALSE;

      for( int e = 0; e < nedges; e++ )
      {
         if( result[e] == CONNECT )
         {
            assert(graph->oeat[e] != EAT_FREE);

            nodearrchar[graph->head[e]] = TRUE;
            nodearrchar[graph->tail[e]] = TRUE;
         }
      }
      if( SCIPisZero(scip, minpathcost) )
         keepsol = TRUE;
   }

   /* try to eliminate vertices and edges */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !graph->mark[k] )
         continue;

      assert(!Is_pterm(graph->term[k]));

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
   const SCIP_Bool dualisvalid = dualCostIsValid(scip, transgraph, bestcost, state, nodearrchar);

   if( !dualisvalid )
   {
      *bestlpobjval = *lpobjval;
      BMScopyMemoryArray(bestcost, cost, extnedges);
   }

   /* has upperbound and changed, and is best reduced cost valid and different from cost? */
   if( SCIPisGT(scip, *bestlpobjval, *lpobjval) && SCIPisLT(scip, *upperbound, oldupperbound) )
   {
      *solgiven = *solgiven && graph_sol_unreduced(scip, graph, result);

      assert(!(*solgiven) || graph_sol_valid(scip, graph, result));

      *minpathcost = *upperbound - *bestlpobjval;
      assert(SCIPisGE(scip, *minpathcost, 0.0));

      computeTransVoronoi(scip, transgraph, vnoi, bestcost, costrev, pathdist, vbase, pathedge);

      return reducePcMw(scip, graph, transgraph, vnoi, bestcost, pathdist, *minpathcost, result, marked, nodearrchar, *solgiven);
   }
   return 0;
}


/** dual ascent based reductions todo: save best lower bound like in PCMW */
SCIP_RETCODE reduce_da(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   SCIP_Real*            ub,                 /**< pointer to provide upper bound and return upper bound found during ascent and prune (if better) */
   SCIP_Real*            offset,             /**< pointer to store offset */
   int*                  edgearrint,         /**< int edges array for internal computations or NULL */
   int*                  vbase,              /**< array for Voronoi bases */
   int*                  state,              /**< int 4 * nnodes array for internal computations */
   int*                  pathedge,           /**< array for predecessor edge on a path */
   int*                  nodearrint,         /**< int nnodes array for internal computations */
   int*                  heursources,        /**< array to store starting points of TM heuristic */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   int                   prevrounds,         /**< number of reduction rounds that have been performed already */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int*                  edgestate,          /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   SCIP_Bool             userec              /**< use recombination heuristic? */
)
{
   STPSOLPOOL* pool;
   SCIP_Real ecost[4];
   IDX* ancestors[4];
   IDX* revancestors[4];
   int adjvert[4];
   int incedge[4];

   SCIP_Real ubnew;
   SCIP_Real maxcost;
   SCIP_Real lpobjval;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   SCIP_Bool success;
   SCIP_Real hopfactor;
   SCIP_Bool apsol;
   const SCIP_Bool rpc = (graph->stp_type == STP_RPCSPG);
   const SCIP_Bool directed = (graph->stp_type == STP_SAP || graph->stp_type == STP_NWSPG);
   int* const grad = graph->grad;
   int* terms;
   int* result;
   int* starts;
   int i;
   int k;
   int e;
   int etmp;
   int runs;
   int root;
   int nruns;
   int nterms;
   int nedges;
   int nnodes;
   int nfixed;
   int best_start;
   STP_Bool* marked;
   SCIP_Bool externsol;
   SCIP_Bool deletable;
   SCIP_Bool checkstate;

   assert(ub != NULL);
   assert(scip != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);
   assert(nodearrint != NULL);

   checkstate = (edgestate != NULL);
   deletable = TRUE;

   root = graph->source;
   apsol = FALSE;
   nfixed = 0;
   nedges = graph->edges;
   nnodes = graph->knots;
   maxcost = 0.0;
   hopfactor = DEFAULT_HOPFACTOR;

   if( graph->terms <= 1 )
      return SCIP_OKAY;

   if( SCIPisGE(scip, *ub, 0.0) )
      upperbound = *ub;
   else
      upperbound = FARAWAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, nedges) );

   if( !rpc && !directed )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &terms, graph->terms) );
   }
   else
   {
      terms = NULL;
   }

   for( i = 0; i < 4; i++ )
   {
      ancestors[i] = NULL;
      revancestors[i] = NULL;
   }

   /*
    * 1. step: compute upper bound
    */
   if( userec )
      SCIP_CALL( SCIPStpHeurRecInitPool(scip, &pool, nedges, SOLPOOL_SIZE) );
   else
      pool = NULL;

   /* initialize */
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

   /* number of runs should not exceed number of connected vertices */
   if( directed )
      runs = MIN(k, DEFAULT_HEURRUNS);
   else
      runs = MIN(k, DEFAULT_HEURRUNS / 5);

   /* neither PC, MW, RPC nor HC? */
   if( !rpc && heursources != NULL )
   {
      /* choose starting points for TM heuristic */
      starts = heursources;
      SCIPStpHeurTMCompStarts(graph, starts, &runs);
   }
   else
   {
      starts = NULL;
   }

   /* fill result array with current best incumbent solution if possible */
   if( nodereplacing == FALSE )
   {
      SCIP_SOL* bestsol;
      SCIP_Real* xval = NULL;

      externsol = TRUE;
      bestsol = SCIPgetBestSol(scip);

      /* no solution available? */
      if( bestsol == NULL )
         externsol = FALSE;
      else
         xval = SCIPprobdataGetXval(scip, bestsol);

      if( xval == NULL )
      {
         externsol = FALSE;
      }
      else
      {
         for( e = 0; e < nedges; e++ )
         {
            if( SCIPisEQ(scip, xval[e], 1.0) )
               result[e] = CONNECT;
            else
               result[e] = UNKNOWN;

            if( result[e] == CONNECT && graph->oeat[e] == EAT_FREE )
               externsol = FALSE;
         }

         if( externsol )
         {
            assert(graph_sol_valid(scip, graph, result));
            upperbound = graph_sol_getObj(graph->cost, result, 0.0, nedges);
         }
      }
   }
   else
   {
      externsol = FALSE;
   }

   for( e = 0; e < nedges; e++ )
   {
      cost[e] = graph->cost[e];
      if( !externsol )
         result[e] = UNKNOWN;
      costrev[e] = graph->cost[flipedge(e)];

      if( graph->stp_type == STP_DHCSTP && SCIPisLT(scip, graph->cost[e], BLOCKED) && SCIPisGT(scip, graph->cost[e], maxcost) )
         maxcost = graph->cost[e];
   }

   if( (directed || rpc) && !externsol )
   {
      SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, graph, starts, &best_start, result, runs, root, cost, costrev, &hopfactor, NULL, maxcost, &success, FALSE) );

      /* calculate objective value of solution */
      ubnew = graph_sol_getObj(graph->cost, result, 0.0, nedges);

      if( SCIPisLT(scip, ubnew, upperbound) )
         upperbound = ubnew;

      /* no feasible solution found? */
      if( !success )
      {
         goto TERMINATE;
      }
   }

   /*
    * 2. step: repeatedly compute lower bound and reduced costs and try to eliminate
    */

   /* if not RPC, select roots for dual ascent */
   if( !rpc && !directed )
   {
      int* termdegs;
      int maxdeg = 0;

      assert(terms != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &termdegs, graph->terms) );

      k = 0;
      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(graph->term[i]) && (grad[i] > 1) )
         {
            termdegs[k] = grad[i];
            terms[k++] = i;

            if( grad[i] > maxdeg )
               maxdeg = grad[i];
         }
      }

      nruns = MIN(k, DEFAULT_DARUNS);

      orderDaRoots(terms, termdegs, k, maxdeg, (prevrounds > 0), randnumgen);

      SCIPfreeBufferArray(scip, &termdegs);
   }
   else
   {
      nruns = 0;
      if( rpc )
         graph_pc_2transcheck(graph);
   }

   maxcost = -1;

   /*  todo: use second array to store best solution and check after each iteration whether it is valid */
   for( int run = 0, bestroot = graph->source; run <= nruns; run++ )
   {
      SCIP_Bool soladded;
      SCIP_Bool usesol = (run > 1);

      if( rpc )
         usesol = TRUE;

      /* graph vanished? */
      if( grad[graph->source] == 0 )
         break;

      if( run == nruns )
      {
         root = bestroot;
         usesol = FALSE;
      }
      else if( !rpc && !directed )
      {
         assert(terms != NULL);
         root = terms[run];
         usesol = usesol && (SCIPrandomGetInt(randnumgen, 0, 2) < 2) && graph->stp_type != STP_RSMT;
      }

      if( externsol && (run == 0) && graph->stp_type != STP_RSMT )
      {
         SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, result, edgearrint, state, root, 1, nodearrchar) );
      }
      else if( usesol )
      {
         /* solution might not be valid anymore due to pseudo-node elimination */
         if( !graph_sol_valid(scip, graph, result) )
         {
            SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, NULL, edgearrint, state, root, 1, nodearrchar) );
         }
         else
         {
            SCIPdebugMessage("reroot solution \n");
            SCIP_CALL( graph_sol_reroot(scip, graph, result, root) );
#ifndef NDEBUG
            {
            const int realroot = graph->source;
            graph->source = root;
            assert(graph_sol_valid(scip, graph, result));
            graph->source = realroot;
            }
#endif

            SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, result, edgearrint, state, root, 1, nodearrchar) );
         }
      }
      else
      {
         SCIPdebugMessage("no rerooting \n");
         SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, NULL, edgearrint, state, root, 1, nodearrchar) );
      }

      /* perform ascent and prune */
      if( externsol && run == 0 )
      {
         apsol = TRUE;
         assert(graph_sol_valid(scip, graph, result));
      }
      else if( !directed )
      {
         SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, result, nodearrint, root, nodearrchar, &success, FALSE) );
         assert(success);
         assert(graph_sol_valid(scip, graph, result));

         /* calculate objective value of solution */
         ubnew = graph_sol_getObj(graph->cost, result, 0.0, nedges);

         SCIPdebugMessage("new obj: %f \n", ubnew);

         soladded = FALSE;

         if( userec )
         {
            assert(graph_valid(graph));

            SCIPdebugMessage("obj before local %f \n", ubnew);
            assert(graph_sol_valid(scip, graph, result));

            SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result) );
            ubnew = graph_sol_getObj(graph->cost, result, 0.0, nedges);

            SCIPdebugMessage("obj after local  %f \n", ubnew);
            assert(graph_sol_valid(scip, graph, result));

            SCIP_CALL(SCIPStpHeurRecAddToPool(scip, ubnew, result, pool, &soladded));
         }

         if( success )
         {
            if( userec && soladded && pool->size >= 2 && (ubnew < upperbound || maxcost > lpobjval) )
            {
               /* get index of just added solution */
               int solindex = pool->maxindex;
               SCIP_Bool solfound;

               SCIPdebugMessage("POOLSIZE %d \n", pool->size);

               SCIP_CALL(SCIPStpHeurRecRun(scip, pool, NULL, NULL, graph, NULL, &solindex, 1, pool->size, FALSE, &solfound));

               if( solfound )
               {
                  const STPSOL* const sol = SCIPStpHeurRecSolfromIdx(pool, solindex);

                  assert(sol != NULL);

                  SCIPdebugMessage("DA: rec found better solution with obj %f vs %f \n", sol->obj, ubnew);

                  if( sol->obj < ubnew )
                  {
                     BMScopyMemoryArray(result, sol->soledges, nedges);
                     SCIPdebugMessage("obj before local2 %f \n", sol->obj);

                     SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result) );

                     ubnew = graph_sol_getObj(graph->cost, result, 0.0, nedges);
                     SCIPdebugMessage("obj after local2  %f \n", ubnew);

                     assert(SCIPisLE(scip, ubnew, sol->obj));

                     if( ubnew < sol->obj )
                        SCIP_CALL(SCIPStpHeurRecAddToPool(scip, ubnew, result, pool, &solfound));
                  }
               }
            }

            if( SCIPisLE(scip, ubnew, upperbound) )
            {
               apsol = TRUE;
               upperbound = ubnew;
            }
            else
            {
               apsol = FALSE;
            }
         }
         else
         {
            if( rpc )
               graph_pc_2org(graph);

            goto TERMINATE;
         }
      }

      if( maxcost < lpobjval && !usesol )
      {
         maxcost = lpobjval;
         bestroot = root;
      }

      /* the required reduced path cost to be surpassed */
      minpathcost = upperbound - lpobjval;

      SCIPdebugMessage("upper: %f lower: %f \n", upperbound, lpobjval);

      assert(SCIPisGE(scip, minpathcost, 0.0));

      if( minpathcost < 0.0 )
         minpathcost = 0.0;

      for( e = 0; e < nedges; e++ )
      {
         marked[e] = FALSE;
         costrev[e] = cost[flipedge(e)];
      }

      /* todo: all the rest in separate function todo: add check if minpathcost is zero to keep best solution! */
      for( k = 0; k < nnodes; k++ )
         graph->mark[k] = (grad[k] > 0);

      /* distance from root to all nodes */
      graph_path_execX(scip, graph, root, cost, pathdist, pathedge);

      /* no paths should go back to the root */
      for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
         costrev[e] = FARAWAY;

      /* build voronoi diagram */
      if( directed )
      {
         for( k = 0; k < nnodes; k++ )
            graph->mark[k] = (graph->grad[k] > 0);

         graph_voronoiTerms(scip, graph, costrev, vnoi, vbase, graph->path_heap, state);
      }
      else
      {
         graph_get4nextTerms(scip, graph, costrev, costrev, vnoi, vbase, graph->path_heap, state);
      }

      /* RPC? If yes, restore original graph */
      if( rpc )
         graph_pc_2org(graph);

      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( rpc && Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) && k != root )
         {
            (*nelims) += graph_pc_deleteTerm(scip, graph, k);
            (*offset) += graph->prize[k];
            continue;
         }

         if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
         {
            if( checkstate )
            {
               deletable = TRUE;
               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  if( edgestate[e] == EDGE_BLOCKED )
                  {
                     deletable = FALSE;
                     break;
                  }
            }

            if( deletable )
            {
               nfixed += graph->grad[k];
               while( graph->outbeg[k] != EAT_LAST )
                  graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
            }
            else
            {
               e = graph->outbeg[k];
               while( e != EAT_LAST )
               {
                  etmp = graph->oeat[e];

                  /* for rpc no artificial terminal arcs should be deleted */
                  if( (rpc && !graph->mark[graph->head[e]]) || (edgestate[e] == EDGE_BLOCKED) )
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
         else
         {
            e = graph->outbeg[k];
            while( e != EAT_LAST )
            {
               etmp = graph->oeat[e];

               /* for rpc no artificial terminal arcs should be deleted */
               if( (rpc && !graph->mark[graph->head[e]]) || (checkstate && edgestate[e] == EDGE_BLOCKED) )
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

      if( apsol )
      {
         for( k = 0; k < nnodes; k++ )
            nodearrchar[k] = FALSE;
         for( e = 0; e < nedges; e++ )
         {
            if( result[e] == CONNECT )
            {
               nodearrchar[graph->tail[e]] = TRUE;
               nodearrchar[graph->head[e]] = TRUE;
            }
         }

         for( k = 0; k < nnodes; k++ )
         {
            if( !graph->mark[k] )
               continue;

            if( !Is_gterm(graph->term[k]) && pathdist[k] + vnoi[k].dist >= minpathcost && !nodearrchar[k] )
            {
               if( checkstate )
               {
                  deletable = TRUE;
                  for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                     if( edgestate[e] == EDGE_BLOCKED )
                     {
                        deletable = FALSE;
                        break;
                     }
               }

               if( deletable )
               {
                  nfixed += graph->grad[k];
                  while( graph->outbeg[k] != EAT_LAST )
                     graph_edge_del(scip, graph, graph->outbeg[k], TRUE);
               }
               else
               {
                  e = graph->outbeg[k];
                  while( e != EAT_LAST )
                  {
                     etmp = graph->oeat[e];

                     /* for rpc no artificial terminal arcs should be deleted */
                     if( (rpc && !graph->mark[graph->head[e]]) || result[e] == CONNECT || result[flipedge(e)] == CONNECT
                           || (edgestate[e] == EDGE_BLOCKED) )
                     {
                        e = etmp;
                        continue;
                     }

                     if( SCIPisGE(scip, pathdist[k] + cost[e] + vnoi[graph->head[e]].dist, minpathcost) )
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
            else
            {
               e = graph->outbeg[k];
               while( e != EAT_LAST )
               {
                  etmp = graph->oeat[e];

                  /* for rpc no artificial terminal arcs should be deleted */
                  if( (rpc && !graph->mark[graph->head[e]]) || result[e] == CONNECT || result[flipedge(e)] == CONNECT
                        || (checkstate && edgestate[e] == EDGE_BLOCKED) )
                  {
                     e = etmp;
                     continue;
                  }

                  if( pathdist[k] + cost[e] + vnoi[graph->head[e]].dist >= minpathcost )
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
      }

      if( !directed && nodereplacing )
      {
         for( k = 0; k < nnodes; k++ )
         {
            if( !graph->mark[k] || Is_gterm(graph->term[k]) )
               continue;
            if( rpc && k == root )
               continue;

            // todo for degree 4 ... 6  as well!!!
            if( graph->grad[k] == 3 )
            {
               SCIP_Real d = FARAWAY;
               int i1 = graph->tail[vnoi[k].edge];
               int i2;
               int base = vbase[k];

               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
                  i2 = graph->head[e];
                  if( i1 != i2 && (!rpc || i2 != root) )
                  {
                     if( base != vbase[i2] && SCIPisLT(scip, cost[e] + vnoi[i2].dist, d) )
                     {
                        d = cost[e] + vnoi[i2].dist;
                     }
                     else if( base == vbase[i2] && SCIPisLT(scip, cost[e] + vnoi[i2 + nnodes].dist, d) )
                     {
                        d = cost[e] + vnoi[i2 + nnodes].dist;
                     }
                  }
               }

               if( SCIPisGT(scip, vnoi[k + nnodes].dist, d) )
                  d = vnoi[k + nnodes].dist;

               if( SCIPisGT(scip, pathdist[k] + vnoi[k].dist + d, minpathcost) )
               {
                  nfixed++;
                  i = 0;

                  for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     graph->mark[graph->head[e]] = FALSE;
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
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[i]), graph->ancestors[incedge[i]], NULL) );
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[i]), graph->ancestors[Edge_anti(incedge[i])], NULL) );
                  }

                  for( i = 0; i < 3; i++ )
                  {
                     i2 = (i + 1) % 3;
                     SCIP_CALL( graph_edge_reinsert(scip, graph, incedge[i], adjvert[i], adjvert[i2], ecost[i] + ecost[i2], ancestors[i],
                           ancestors[i2], revancestors[i], revancestors[i2]) );
                  }
               }
            }
         }
      }

      assert(graph_valid(graph));

      if( !rpc )
      {
         for( k = 0; k < nnodes; k++ )
            graph->mark[k] = graph->grad[k] > 0;
      }
      else
      {
         graph_pc_2trans(graph);
         graph_pc_2org(graph);
      }
   } /* root loop */

 TERMINATE:
   *nelims = nfixed;

   if( SCIPisLT(scip, upperbound, *ub) || SCIPisLT(scip, *ub, 0.0) )
      *ub = upperbound;

   if( userec )
      SCIPStpHeurRecFreePool(scip, &pool);

   /* free memory */
   for( k = 0; k < 4; k++ )
   {
      SCIPintListNodeFree(scip, &(ancestors[k]));
      SCIPintListNodeFree(scip, &(revancestors[k]));
   }

   SCIPfreeBufferArrayNull(scip, &terms);
   SCIPfreeBufferArray(scip, &marked);
   SCIPfreeBufferArray(scip, &result);

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
      graph_pc_2trans(graph);
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
      obj = graph_sol_getObj(graph->cost, edgearrint, 0.0, nedges);

      SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, edgearrint, edgearrint2, state, root, 1, nodearrchar) );

   }
   else
   {
      obj = FARAWAY;
      SCIP_CALL( SCIPStpDualAscent(scip, graph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, NULL, edgearrint2, state, root, 1, nodearrchar) );
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
      for( k = 0; k < nnodes; k++ ) // TODO
         if( Is_term(prunegraph->term[k]) && !mark[k]  )
            printf("in bnd  FAIL %d not marked, but terminal, \n", k);
#endif

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, edgearrint2, vbase, root, nodearrchar, &success, FALSE) );

   objprune = graph_sol_getObj(graph->cost, edgearrint2, 0.0, nedges);

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
      graph_pc_2org(graph);
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
      graph_pc_2trans(graph);
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
   PATH*                 vnoi,               /**< Voronoi data structure array */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< reduced edge costs */
   SCIP_Real*            costrev,            /**< reduced reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   int*                  vbase,              /**< Voronoi base array */
   int*                  pathedge,           /**< shortest path incoming edge array for shortest path calculations */
   int*                  edgearrint,         /**< int edges array for internal computations or NULL */
   int*                  state,              /**< int 4 * vertices array  */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   SCIP_Bool             solbasedda,         /**< rerun Da based on best primal solution */
   SCIP_Bool             varyroot,           /**< vary root for DA if possible */
   SCIP_Bool             shiftcosts,         /**< should costs be shifted to try to reduce number of terminals? */
   SCIP_Bool             markroots,          /**< should terminals proven to be part of an opt. sol. be marked as such? */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             prizesum            /**< sum of positive prizes */
)
{
   STPSOLPOOL* pool;
   GRAPH* transgraph;
   SCIP_Real* bestcost;
   SCIP_Real ub;
   SCIP_Real offset;
   SCIP_Real lpobjval;
   SCIP_Real bestlpobjval;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   int* roots;
   int* result;
   int* result2;
   int* transresult;
   STP_Bool* marked;
   int nroots;
   int nfixed;
   int nusedroots;
   const int root = graph->source;
   const int nedges = graph->edges;
   const int extnedges = nedges + 2 * (graph->terms - 1);
   const int nnodes = graph->knots;
   SCIP_Bool apsol;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);
   assert(nodearrchar != NULL);

   /* not more than one terminal? */
   if( graph->terms <= 1 )
      return SCIP_OKAY;

   nroots = 0;
   nfixed = 0;
   varyroot = varyroot && userec;

   /* allocate memory */
   if( edgearrint == NULL )
      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   else
      result = edgearrint;

   SCIP_CALL( SCIPallocBufferArray(scip, &transresult, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, extnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result2, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcost, extnedges) );

   if( userec )
      SCIP_CALL( SCIPStpHeurRecInitPool(scip, &pool, nedges, SOLPOOL_SIZE) );
   else
      pool = NULL;

   /*
    * 1. step: compute upper bound
    */

#ifndef NDEBUG
   {
   int nterms = 0;
   for( int i = 0; i < nnodes; i++ )
      if( graph->mark[i] )
         if( Is_term(graph->term[i]) )
            nterms++;

   assert(nterms == (graph->terms - ((graph->stp_type != STP_RPCSPG)? 1 : 0)));
   }
#endif

   /*
    * 2. step: compute lower bound and reduced costs
    */

   graph_pc_2trans(graph);
   offset = 0.0;

   /* transform the problem to a real SAP */
   if( shiftcosts )
      SCIP_CALL( graph_pc_getSapShift(scip, graph, &transgraph, &offset) );
   else
      SCIP_CALL( graph_pc_getSap(scip, graph, &transgraph, &offset) );

   /* initialize data structures for shortest paths */
   SCIP_CALL( graph_path_init(scip, transgraph) );

   SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, NULL, transresult, state, root, 1, nodearrchar) );

   lpobjval += offset;
   bestlpobjval = lpobjval;
   BMScopyMemoryArray(bestcost, cost, extnedges);

   /* compute first primal solution */
   upperbound = FARAWAY;
   apsol = FALSE;
   SCIP_CALL( computeDaSolPcMw(scip, graph, NULL, vnoi, cost, pathdist, &upperbound, result, result2, vbase, pathedge, root, nodearrchar, &apsol) );

   assert(apsol && upperbound < FARAWAY);
   assert(graph_sol_valid(scip, graph, result));

   /* the required reduced path cost to be surpassed */
   minpathcost = upperbound - lpobjval;
   if( userec)
      SCIPdebugMessage("DA first minpathcost %f \n", minpathcost);

   /* initialize data structures for transgraph */
   SCIP_CALL( graph_init_history(scip, transgraph) );
   computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);

   /*
    * 3. step: try to eliminate
    */

   /* restore original graph */
   graph_pc_2org(graph);

   if( shiftcosts )
   {
      for( int i = 0; i < nnodes; i++ )
         if( Is_term(graph->term[i]) && transgraph->term[i] == -1 )
            graph->mark[i] = FALSE;

      varyroot = FALSE;
   }

   for( int e = 0; e < extnedges; e++ )
      marked[e] = FALSE;

   /* try to reduce the graph */
   assert(dualCostIsValid(scip, transgraph, cost, state, nodearrchar));
   nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, TRUE);

   /* edges from result might have been deleted! */
   apsol = apsol && graph_sol_unreduced(scip, graph, result);
   assert(!apsol || graph_sol_valid(scip, graph, result));

   if( userec )
      SCIPdebugMessage("DA: 1. NFIXED %d \n", nfixed);

   /* rerun dual ascent? */
   if( solbasedda && graph->terms > 2 )
   {
      const SCIP_Real oldupperbound = upperbound;

      /* with recombination? */
      if( userec && graph->stp_type != STP_MWCSP )
      {
         int dummy;

         /* compute second solution and add to pool */
         SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, graph, NULL, &dummy, result2, DEFAULT_HEURRUNS / 5, root, graph->cost, graph->cost, NULL, NULL, 0.0, &success, FALSE) );
         assert(success);

         SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, graph->cost, result2) );
         ub = graph_sol_getObj(graph->cost, result2, 0.0, nedges);

         SCIP_CALL( SCIPStpHeurRecAddToPool(scip, ub, result2, pool, &success) );
         SCIPdebugMessage("added initial TM sol to pool? %d , ub %f \n", success, ub);
      }

      /* try to improve both dual and primal bound */
      SCIP_CALL( computePertubedSol(scip, graph, transgraph, pool, vnoi, gnodearr, cost, costrev, bestcost, pathdist, state, vbase, pathedge, result, result2,
            transresult, nodearrchar, &upperbound, &lpobjval, &bestlpobjval, &minpathcost, &apsol, offset, extnedges, 0) );

      assert(graph_sol_valid(scip, graph, result));
      assert(!apsol || SCIPisEQ(scip, graph_sol_getObj(graph->cost, result, 0.0, nedges), upperbound));

      graph_pc_2orgcheck(graph);

      assert(dualCostIsValid(scip, transgraph, cost, state, nodearrchar));
      computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);

      nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, apsol);

      nfixed += reducePcMwTryBest(scip, graph, transgraph, vnoi, cost, costrev, bestcost, pathdist, &upperbound,
            &lpobjval, &bestlpobjval, &minpathcost, oldupperbound, result, vbase, state, pathedge, marked, nodearrchar, &apsol, extnedges);

      if( userec )
         SCIPdebugMessage("eliminations after sol based run2 with best dual sol %d bestlb %f newlb %f\n", nfixed, bestlpobjval, lpobjval);
   }

   /* pertubation runs for MWCSP */
   if( varyroot && graph->stp_type == STP_MWCSP )
   {
      for( int run = 0; run < DEFAULT_NMAXROOTS && graph->terms > STP_RED_MINBNDTERMS; run++ )
      {
         SCIP_Real oldupperbound = upperbound;

         graph_pc_2trans(graph);

         apsol = apsol && graph_sol_unreduced(scip, graph, result);
         assert(!apsol || graph_sol_valid(scip, graph, result));

         assert(SCIPisEQ(scip, upperbound, graph_sol_getObj(graph->cost, result, 0.0, nedges)));

         /* try to improve both dual and primal bound */
         SCIP_CALL( computePertubedSol(scip, graph, transgraph, pool, vnoi, gnodearr, cost, costrev, bestcost, pathdist, state, vbase, pathedge, result, result2,
               transresult, nodearrchar, &upperbound, &lpobjval, &bestlpobjval, &minpathcost, &apsol, offset, extnedges, run) );

         SCIPdebugMessage("DA: pertubated run %d ub: %f \n", run, upperbound);
         SCIPdebugMessage("DA: pertubated run %d minpathcost: %f \n", run, upperbound - lpobjval);

         computeTransVoronoi(scip, transgraph, vnoi, cost, costrev, pathdist, vbase, pathedge);
         nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, apsol);

         nfixed += reducePcMwTryBest(scip, graph, transgraph, vnoi, cost, costrev, bestcost, pathdist, &upperbound,
               &lpobjval, &bestlpobjval, &minpathcost, oldupperbound, result, vbase, state, pathedge, marked, nodearrchar, &apsol, extnedges);

         SCIPdebugMessage("DA: pertubated run %d NFIXED %d \n", run, nfixed);
      }
   }

#if 0
   if( graph->stp_type == STP_MWCSP && graph->terms < STP_RED_MINBNDTERMS  )
      varyroot = FALSE;
#endif

   /* change or mark roots? */
   if( varyroot || markroots )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &roots, graph->terms) );

      /*
       * get possible roots
       */

      BMSclearMemoryArray(nodearrchar, nnodes);

      assert(nroots == 0);

      for( int e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
      {
         int k = graph->head[e];

         if( Is_pterm(graph->term[k]) )
         {
            int l = -1;
            int e3 = -1;

            for( int e2 = transgraph->inpbeg[k]; e2 != EAT_LAST; e2 = transgraph->ieat[e2] )
            {
               if( transgraph->cost[e2] == 0.0 )
                  l = transgraph->tail[e2];
               else
                  e3 = e2;
            }

            assert(l >= 0);
            assert(e3 >= 0);
            assert(l != root);
            assert(SCIPisEQ(scip, graph->cost[e], graph->cost[e3]));
            assert(graph->mark[l]);

            if( SCIPisGT(scip, cost[e3], upperbound - lpobjval) || SCIPisGT(scip, bestcost[e3], upperbound - bestlpobjval) )
            {
               assert(SCIPisPositive(scip, graph->prize[l]));
               assert(nroots < graph->terms);

               roots[nroots++] = l;
               nodearrchar[l] = TRUE;
            }
         }
      }

      if( nroots > 0 && graph->terms > 2 )
      {
         /*
          * try to find additional roots by shortest path computations
          */

         SCIP_Real maxprize = 0.0;
         int qstart = 0;

         assert(!graph->extended);

         for( int i = 0; i < nnodes; i++ )
         {
            state[i] = UNKNOWN;
            vnoi[i].dist = FARAWAY;
            vnoi[i].edge = UNKNOWN;

            if( Is_term(graph->term[i]) && !nodearrchar[i] && graph->mark[i] && maxprize < graph->prize[i] && graph->prize[i] < prizesum )
               maxprize = graph->prize[i];
         }

         for( int rounds = 0; rounds < 10 && qstart != nroots; rounds++ )
         {
            const int oldnroots = nroots;

            assert(!graph->extended);

            SCIPdebugMessage("root-finding round %d \n", rounds);

            for( int i = qstart; i < nroots; i++ )
            {
               int nlabeled;
               const int term = roots[i];
               assert(nodearrchar[term]);
               assert(Is_term(graph->term[term]));
               assert(SCIPisPositive(scip, graph->prize[term]));

               /* look for terminals in neighborhood of term */
               graph_sdPaths(scip, graph, vnoi, graph->cost, maxprize, pathedge, state, vbase, &nlabeled, term, term, 200);

               /* check labelled nodes */
               for( int k = 0; k < nlabeled; k++ )
               {
                  const int l = vbase[k];

                  assert(graph->mark[l]);
                  assert(state[l] != UNKNOWN);

                  if( Is_term(graph->term[l]) && !nodearrchar[l] && vnoi[l].dist <= graph->prize[l] )
                  {
                     roots[nroots++] = l;
                     nodearrchar[l] = TRUE;

                     SCIPdebugMessage("added new root %d dist: %f, prize: %f  \n", l, vnoi[l].dist, graph->prize[l]);
                  }

                  /* restore data structures */
                  state[l] = UNKNOWN;
                  vnoi[l].dist = FARAWAY;
                  vnoi[l].edge = UNKNOWN;
               }
            }
#ifndef NDEBUG
            for( int k = 0; k < nnodes; k++ )
            {
               assert(state[k] == UNKNOWN);
               assert(vnoi[k].dist == FARAWAY);
               assert(vnoi[k].edge == UNKNOWN);
            }
#endif

            qstart = oldnroots;
         }
         SCIPdebugMessage("number of terminals found by DA: %d \n", nroots);
      }

      /* should prize of terminals be changed? */
      if( nroots > 0 && markroots  )
      {
         if( userec )
         {
            userec = FALSE;
            SCIPStpHeurRecFreePool(scip, &pool);
            pool = NULL;
         }

         markPcMwRoots(scip, graph, roots, nroots, prizesum);
      }

      if( nroots > 0 && varyroot )
      {
         int* termdegs;
         int maxdeg = 0;

         SCIP_CALL(SCIPallocBufferArray(scip, &termdegs, nroots));

         for( int i = 0; i < nroots; i++ )
         {
            termdegs[i] = graph->grad[roots[i]];
            if( maxdeg < termdegs[i] )
               maxdeg = termdegs[i];
         }

         orderDaRoots(roots, termdegs, nroots, maxdeg, TRUE, randnumgen);

         SCIPfreeBufferArray(scip, &termdegs);
      }
   }
   else
   {
      roots = NULL;
   }

   if( varyroot )
      nusedroots = MIN(DEFAULT_NMAXROOTS, nroots);
   else
      nusedroots = -1;

   graph_path_exit(scip, transgraph);
   graph_free(scip, transgraph, TRUE);

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

      SCIP_CALL( graph_pc_getRsap(scip, graph, &transgraph, roots, nroots, tmproot) );

      assert(graph_valid(transgraph));

      transnnodes = transgraph->knots;
      transnedges = transgraph->edges;

      for( int k = 0; k < transnnodes; k++ )
         transgraph->mark[k] = (transgraph->grad[k] > 0);

      /* init data structures for shortest paths and history */
      SCIP_CALL( graph_path_init(scip, transgraph) );
      SCIP_CALL( graph_init_history(scip, transgraph ) );

      transgraph->stp_type = STP_SAP;

      if( apsol && run > 1 )
      {
         BMScopyMemoryArray(transresult, result, graph->edges);
         SCIP_CALL(graph_sol_reroot(scip, transgraph, transresult, tmproot));
         SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lpobjval, FALSE, FALSE,
               gnodearr, transresult, result2, state, tmproot, 1, nodearrchar));
      }
      else
      {
         SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lpobjval, FALSE, FALSE,
               gnodearr, NULL, transresult, state, tmproot, 1, nodearrchar));
      }

      assert(graph_valid(transgraph));

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

      apsol = FALSE;
      SCIP_CALL( computeDaSolPcMw(scip, graph, pool, vnoi, cost, pathdist, &upperbound, result, result2, vbase, pathedge, root, nodearrchar, &apsol) );

      SCIPdebugMessage("ROOTRUNS upperbound %f \n", upperbound);
      SCIPdebugMessage("ROOTRUNS best sol in pool %f \n", pool->sols[0]->obj);

      for( int k = 0; k < transnnodes; k++ )
         transgraph->mark[k] = (transgraph->grad[k] > 0);

      /* the required reduced path cost to be surpassed */
      minpathcost = upperbound - lpobjval;

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
      graph_pc_2org(graph);

      assert(graph->mark[tmproot]);
      graph->mark[tmproot] = FALSE;

      for( int e = 0; e < transnedges; e++ )
         marked[e] = FALSE;

      /* try to eliminate vertices and edges */
      nfixed += reducePcMw(scip, graph, transgraph, vnoi, cost, pathdist, minpathcost, result, marked, nodearrchar, apsol);

      apsol = apsol && graph_sol_unreduced(scip, graph, result);
      assert(!apsol || graph_sol_valid(scip, graph, result));
      SCIPdebugMessage("FIXED with changed root %d \n\n", nfixed);

      graph->mark[tmproot] = TRUE;

      transgraph->stp_type = STP_RPCSPG;
      graph_path_exit(scip, transgraph);
      graph_free(scip, transgraph, TRUE);
   }

   *nelims = nfixed;

   /* free memory */
   if( pool != NULL )
      SCIPStpHeurRecFreePool(scip, &pool);
   SCIPfreeBufferArrayNull(scip, &roots);
   SCIPfreeBufferArray(scip, &bestcost);
   SCIPfreeBufferArray(scip, &result2);
   SCIPfreeBufferArray(scip, &marked);
   SCIPfreeBufferArray(scip, &transresult);

   if( edgearrint == NULL )
      SCIPfreeBufferArray(scip, &result);

   assert(graph_valid(graph));

   return SCIP_OKAY;
}


/** dual ascent based heuristic reductions for MWCSP */
SCIP_RETCODE reduce_daSlackPruneMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure array */
   GNODE**               gnodearr,           /**< GNODE* terminals array for internal computations or NULL */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reverse edge costs */
   SCIP_Real*            pathdist,           /**< distance array for shortest path calculations */
   int*                  vbase,              /**< Voronoi base array */
   int*                  pathedge,           /**< shortest path incoming edge array for shortest path calculations */
   int*                  soledge,            /**< edge solution array (CONNECT/UNKNOWN) or NULL; needs to contain solution if solgiven == TRUE */
   int*                  state,              /**< int 4 * vertices array */
   int*                  solnode,            /**< array of nodes of current solution that is not to be destroyed */
   STP_Bool*             nodearrchar,        /**< STP_Bool node array for internal computations */
   int*                  nelims,             /**< pointer to store number of reduced edges */
   int                   minelims,           /**< minimum number of edges to eliminate */
   SCIP_Bool             solgiven            /**< solution provided? */
   )
{
   SCIP_HEURDATA* tmheurdata;
   IDX** ancestors;
   IDX** revancestors;
   GRAPH* transgraph;
   SCIP_Real* sd;
   SCIP_Real* ecost;
   SCIP_Real ub;
   SCIP_Real offset;
   SCIP_Bool success;
   SCIP_Real tmpcost;
   SCIP_Real lpobjval;
   SCIP_Real hopfactor;
   SCIP_Real upperbound;
   SCIP_Real minpathcost;
   SCIP_Bool eliminate;
   int* result;
   int* adjvert;
   int* incedge;
   int* reinsert;
   int* transresult;
   int i;
   int k;
   int e;
   int e2;
   int etmp;
   int runs;
   int root;
   int nterms;
   int nedges;
   int nnodes;
   int nfixed;
   int redrounds;
   int best_start;
   int transnnodes;
   int transnedges;
   STP_Bool* marked;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(nelims != NULL);
   assert(nodearrchar != NULL);

   root = graph->source;
   nfixed = 0;
   nnodes = graph->knots;
   nedges = graph->edges;
   success = FALSE;
   hopfactor = DEFAULT_HOPFACTOR;

   /* not more than two terminals? */
   if( graph->terms <= 3 )
      return SCIP_OKAY;

   /* allocate memory */
   if( soledge == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   }
   else
   {
      result = soledge;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &transresult, nedges + 2 * (graph->terms - 1)) );
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

   runs = 0;
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( graph->grad[i] > 0 )
      {
         runs++;
         if( Is_term(graph->term[i]) )
            nterms++;
      }
   }

   assert(nterms == (graph->terms - ((graph->stp_type != STP_RPCSPG)? 1 : 0)));

   offset = 0.0;

   /* transform the problem to a real SAP */
   SCIP_CALL( graph_pc_getSap(scip, graph, &transgraph, &offset) );

   /* initialize data structures for shortest paths and history */
   SCIP_CALL( graph_path_init(scip, transgraph) );
   SCIP_CALL( graph_init_history(scip, transgraph ) );
   transnnodes = transgraph->knots;
   transnedges = transgraph->edges;

   if( !solgiven )
   {
      for( k = 0; k < nnodes; k++ )
         nodearrchar[k] = (solnode[k] == CONNECT);

      for( e = 0; e < nedges; e++ )
      {
         cost[e] = graph->cost[e];
         costrev[e] = graph->cost[flipedge(e)];
         result[e] = UNKNOWN;
      }

      /* build trivial solution */
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, result, nodearrchar) );
   }
   else
   {
      for( e = 0; e < nedges; e++ )
      {
         cost[e] = graph->cost[e];
         costrev[e] = graph->cost[flipedge(e)];
      }
   }

   upperbound = graph_sol_getObj(graph->cost, result, 0.0, nedges);

   /* number of runs should not exceed number of connected vertices */
   runs = MIN(runs, DEFAULT_HEURRUNS);

   /* get TM heuristic data */
   assert(SCIPfindHeur(scip, "TM") != NULL);
   tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));


   /* compute Steiner tree to obtain upper bound */
   SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, graph, NULL, &best_start, transresult, runs, root, cost, costrev, &hopfactor, NULL, 0.0, &success, FALSE) );

   /* feasible solution found? */
   if( success )
   {
      /* calculate objective value of solution */
      ub = graph_sol_getObj(graph->cost, transresult, 0.0, nedges);

      SCIPdebugMessage("TM upperbound in reduce_daSlackPruneMw  %f \n", ub + SCIPprobdataGetOffset(scip));

      if( SCIPisLE(scip, ub, upperbound) )
      {
         upperbound = ub;
         for( e = 0; e < nedges; e++ )
            result[e] = transresult[e];
      }
   }

   /* compute lower bound and reduced costs todo use SCIPdualAscentStpSol */
   SCIP_CALL( SCIPStpDualAscent(scip, transgraph, cost, pathdist, &lpobjval, FALSE, FALSE, gnodearr, NULL, transresult, state, root, 1, nodearrchar) );

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, cost, transresult, vbase, -1, nodearrchar, &success, FALSE) );

   assert(success);
   assert(graph_sol_valid(scip, graph, transresult));

   if( success )
   {
      ub = graph_sol_getObj(graph->cost, transresult, 0.0, nedges);

      SCIPdebugMessage("AP upperbound in reduce_daSlackPruneMw  %f \n", ub + SCIPprobdataGetOffset(scip));

      if( SCIPisLE(scip, ub, upperbound) )
         for( e = 0; e < nedges; e++ )
            result[e] = transresult[e];
   }

   /*
    * 2. step: try to eliminate
    * */

   for( e = 0; e < transnedges; e++ )
   {
      costrev[e] = cost[flipedge(e)];
      marked[e] = FALSE;
   }

   for( k = 0; k < transnnodes; k++ )
      transgraph->mark[k] = (transgraph->grad[k] > 0);

   /* distance from root to all nodes */
   graph_path_execX(scip, transgraph, root, cost, pathdist, pathedge);

   for( i = 0; i < transnnodes; i++ )
      if( Is_term(transgraph->term[i]) )
         assert(SCIPisEQ(scip, pathdist[i], 0.0));

   /* no paths should go back to the root */
   for( e = transgraph->outbeg[root]; e != EAT_LAST; e = transgraph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build Voronoi diagram */
   graph_voronoiTerms(scip, transgraph, costrev, vnoi, vbase, transgraph->path_heap, transgraph->path_state);

   /* restore original graph */
   graph_pc_2org(graph);

   for( k = 0; k < nnodes; k++ )
      solnode[k] = UNKNOWN;

   for( e = 0; e < nedges; e++ )
   {
      costrev[e] = -1.0;
      if( result[e] == CONNECT )
      {
         solnode[graph->head[e]] = CONNECT;
         solnode[graph->tail[e]] = CONNECT;
      }
   }

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

      /* try to eliminate vertices and edges */
      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k]  )
            continue;

         if( nfixed > minelims )
            break;


         if( Is_term(graph->term[k]) )
            continue;

         tmpcost = pathdist[k] + vnoi[k].dist;

         if( SCIPisGT(scip, tmpcost, minpathcost) ||
            (SCIPisGE(scip, tmpcost, minpathcost) && solnode[k] != CONNECT) || (!eliminate && solnode[k] != CONNECT) )
         {
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
               tmpcost = pathdist[k] + cost[e] + vnoi[transgraph->head[e]].dist;

               if( SCIPisGT(scip, tmpcost, minpathcost) ||
                  ( (SCIPisGE(scip, tmpcost, minpathcost) || !eliminate) && result[e] != CONNECT && result[flipedge(e)] != CONNECT) )
               {
                  if( marked[flipedge(e)] )
                  {
                     if( eliminate )
                     {
                        graph_edge_del(scip, graph, e, TRUE);
                        graph_edge_del(scip, transgraph, e, FALSE);
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
   }
   assert(root == transgraph->source);

   SCIPdebugMessage("fixed by da: %d \n", nfixed );

   graph_path_exit(scip, transgraph);
   graph_free(scip, transgraph, TRUE);

   /* restore transformed graph */
   graph_pc_2trans(graph);

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
   SCIPfreeBufferArray(scip, &transresult);

   if( soledge == NULL )
      SCIPfreeBufferArray(scip, &result);

   return SCIP_OKAY;
}



/** bound-based reductions for the (R)PCSTP, the MWCSP and the STP */
SCIP_RETCODE reduce_bound(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            cost,               /**< edge cost array                    */
   SCIP_Real*            prize,              /**< prize (nodes) array                */
   SCIP_Real*            radius,             /**< radius array                       */
   SCIP_Real*            costrev,            /**< reversed edge cost array           */
   SCIP_Real*            offset,             /**< pointer to the offset              */
   SCIP_Real*            upperbound,         /**< pointer to an upper bound          */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nelims              /**< pointer to store number of eliminated edges */
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
   STP_Bool* stnode;
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
   assert(graph->source >= 0);
   assert(upperbound != NULL);

   mst = NULL;
   obj = DEFAULT_HOPFACTOR;
   perm = NULL;
   root = graph->source;
   nedges = graph->edges;
   nnodes = graph->knots;
   mstobj = 0.0;
   *nelims = 0;
   mstobj2 = 0.0;
   best_start = 0;
   ub = SCIPisGT(scip, *upperbound, 0.0);
   mw = (graph->stp_type == STP_MWCSP);
   pc = (graph->stp_type == STP_RPCSPG) || (graph->stp_type == STP_PCSPG);
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
   if( nterms <= 2 || (pcmw && nterms <= 3) )
   {
      /* free memory and return */
      SCIPfreeBufferArrayNull(scip, &stnode);
      SCIPfreeBufferArrayNull(scip, &result);
      return SCIP_OKAY;
   }

   assert(nterms == (graph->terms - ((graph->stp_type == STP_PCSPG || mw)? 1 : 0)));

   runs = MIN(e, DEFAULT_HEURRUNS);

   /* neither PC, MW, RPC nor HC? */
   if( !pcmw && graph->stp_type != STP_DHCSTP )
   {
      /* choose starting points for TM heuristic */

      SCIP_CALL( SCIPallocBufferArray(scip, &starts, nnodes) );

      SCIPStpHeurTMCompStarts(graph, starts, &runs);
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

      if( graph->stp_type == STP_DHCSTP && SCIPisLT(scip, graph->cost[e], BLOCKED) && SCIPisGT(scip, graph->cost[e], maxcost) )
         maxcost = graph->cost[e];
   }

   /* init auxiliary graph */
   if( !mw )
   {
      SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1) );
   }
   else
   {
      adjgraph = NULL;
   }

   /* build voronoi regions, concomitantly building adjgraph and computing radii values*/
   SCIP_CALL( graph_voronoiWithRadius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

   /* get 2nd next terminals to all non-terminal nodes */
   graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* get 3th next terminals to all non-terminal nodes */
   graph_get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* for (rooted) prize collecting get 4th next terminals to all non-terminal nodes */
   if( pc )
      graph_get4next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

   /* no MWCS problem? */
   if( !mw )
   {
      assert(adjgraph != NULL);
      graph_knot_chg(adjgraph, 0, 0);
      adjgraph->source = 0;

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
   if( graph->stp_type == STP_RPCSPG )
   {
      assert(graph->mark[graph->source]);
      for( k = 0; k < nnodes; k++ )
      {
         if( !Is_term(graph->term[k]) || !graph->mark[k] )
            continue;

         if( SCIPisGT(scip, radius[k], prize[k]) && k != graph->source )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_PCSPG )
   {
      for( k = 0; k < nnodes; k++ )
      {
         if( !graph->mark[k] )
            continue;

         if( Is_term(graph->term[k]) && SCIPisGT(scip, radius[k], prize[k])  )
            radius[k] = prize[k];
      }
   }
   else if( graph->stp_type == STP_MWCSP )
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
         graph_pc_2trans(graph);

      SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, graph, starts, &best_start, result, runs, root, cost, costrev, &obj, NULL, maxcost, &success, FALSE) );

      /* PC or RPC? Then restore oringinal graph */
      if( pcmw )
         graph_pc_2org(graph);

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
      *upperbound = obj;
   }
   else
   {
      obj = *upperbound;
      assert(SCIPisGE(scip, obj, 0.0));
   }

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
               if( graph->stp_type == STP_DHCSTP && SCIPisGT(scip, graph->cost[e], graph->cost[flipedge(e)]) )
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
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[0]), graph->ancestors[edges3[0]], NULL) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[0]), graph->ancestors[Edge_anti(edges3[0])], NULL) );

               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[1]), graph->ancestors[edges3[1]], NULL) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[1]), graph->ancestors[Edge_anti(edges3[1])], NULL) );

               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[2]), graph->ancestors[edges3[2]], NULL) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[2]), graph->ancestors[Edge_anti(edges3[2])], NULL) );

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
      SCIP_CALL( graph_get4nextTTerms(scip, graph, cost, vnoi, vbase, heap, state) );

      for( k = 0; k < nnodes; k++ )
      {
         /* is k a terminal other than the root? */
         if( Is_term(graph->term[k]) && graph->mark[k] && graph->grad[k] > 2 && k != graph->source )
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
               (*nelims) += graph_pc_deleteTerm(scip, graph, k);
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
                     (*nelims) += graph_pc_deleteTerm(scip, graph, k);
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

   assert(graph_valid(graph));

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
   SCIP_Real* cost3;
   SCIP_Real* const prize = graph->prize;
   SCIP_Real obj;
   SCIP_Real max;
   SCIP_Real bound;
   SCIP_Real tmpcost;
   SCIP_Real mstobj;
   SCIP_Real maxcost;
   SCIP_Real radiim2;
   SCIP_Real radiim3;
   int* edges3;
   int* nodes3;
   IDX** ancestors = NULL;
   IDX** revancestors = NULL;
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
   assert(graph_valid(graph));

   mstobj = 0.0;
   *nelims = 0;

   cost3 = NULL;
   edges3 = NULL;
   nodes3 = NULL;

   if( nterms <= 2 )
      return SCIP_OKAY;

   assert( graph->stp_type != STP_RPCSPG || Is_term(graph->term[root]) );

   if( !pcmw )
      for( int k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

   /* initialize cost and costrev array */
   maxcost = 0.0;
   for( int e = 0; e < nedges; e++ )
   {
      cost[e] = graph->cost[e];
      costrev[e] = graph->cost[flipedge(e)];

      assert(SCIPisGE(scip, cost[e], 0.0));

      if( graph->stp_type == STP_DHCSTP && SCIPisLT(scip, graph->cost[e], BLOCKED) && SCIPisGT(scip, graph->cost[e], maxcost) )
         maxcost = graph->cost[e];
   }

   /* no MWCSP? */
   if( !mw )
   {
      GRAPH* adjgraph;
      PATH* mst;

      SCIP_CALL( graph_init(scip, &adjgraph, nterms, MIN(nedges, (nterms - 1) * nterms), 1) );

      /* build Voronoi regions, concomitantly building adjgraph and computing radii values*/
      SCIP_CALL( graph_voronoiWithRadius(scip, graph, adjgraph, vnoi, radius, cost, costrev, vbase, heap, state) );

      /* get 2nd next terminals to all non-terminal nodes */
      graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

      /* get 3th next terminals to all non-terminal nodes */
      graph_get3next(scip, graph, cost, costrev, vnoi, vbase, heap, state);

      assert(adjgraph != NULL);
      graph_knot_chg(adjgraph, 0, 0);
      adjgraph->source = 0;
      assert(graph_valid(adjgraph));

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
      graph_free(scip, adjgraph, TRUE);

   }
   else
   {

      /* build Voronoi regions */
      graph_voronoiMw(scip, graph, costrev, vnoi, vbase, heap, state);

      /* get 2nd next positive node to all non-positive nodes */
      graph_get2next(scip, graph, cost, costrev, vnoi, vbase, heap, state);
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
         bound = mstobj;
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

            if( !graph->mark[k] || graph->grad[k] == 0 || Is_gterm(graph->term[k] ) )
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
                     assert(!Is_pterm(graph->term[head]));
                     assert(!Is_pterm(graph->term[tail]));

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
                     assert(!Is_pterm(graph->term[graph->head[e]]));
                     assert(!Is_pterm(graph->term[graph->tail[e]]));

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
                     assert(!Is_pterm(graph->term[head]));
                     assert(!Is_pterm(graph->term[tail]));

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
         /* traverse all nodes, try to eliminate 3 degree nodes todo: extend for 3 ... 5 */
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
               if( graph->grad[k] == 3 && !Is_term(graph->term[k]) )
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

            if( graph->grad[k] == 3 && !Is_term(graph->term[k]) )
            {
               tmpcost = vnoi[k].dist + vnoi[k + nnodes].dist + vnoi[k + 2 * nnodes].dist + radiim3;
               if( SCIPisGT(scip, tmpcost, obj) )
               {
                  (*nelims)++;

                  /* first 3-node elimination? */
                  if( ancestors == NULL )
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &cost3, 3) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &edges3, 3) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &nodes3, 3) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, 3) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, 3) );
                     for( int l = 0; l < 3; l++ )
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
                  for( int e = graph->outbeg[k], l = 0; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     assert(l < 3);
                     edges3[l] = e;
                     nodes3[l] = graph->head[e];
                     cost3[l++] = graph->cost[e];
                  }

                  /* clear */
                  for( int l = 0; l < 3; l++ )
                  {
                     SCIPintListNodeFree(scip, &(ancestors[l]));
                     SCIPintListNodeFree(scip, &(revancestors[l]));
                  }

                  /* store ancestors of incident edges */
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[0]), graph->ancestors[edges3[0]], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[0]), graph->ancestors[Edge_anti(edges3[0])], NULL) );

                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[1]), graph->ancestors[edges3[1]], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[1]), graph->ancestors[Edge_anti(edges3[1])], NULL) );

                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[2]), graph->ancestors[edges3[2]], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[2]), graph->ancestors[Edge_anti(edges3[2])], NULL) );

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
            for( int k = 0; k < 3; k++ )
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

   assert(graph_valid(graph));

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
   root = graph->source;
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
   graph_voronoiTerms(scip, graph, costrev, vnoi, vbase, heap, state);

   if( SCIPisLT(scip, objval, 0.0) )
   {
      /* get TM heuristic data */
      assert(SCIPfindHeur(scip, "TM") != NULL);
      tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

      /* compute UB */
      SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, graph, NULL, &best_start, result, 50, root, cost, costrev, &hopfactor, NULL, maxcost, &success, FALSE) );

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
