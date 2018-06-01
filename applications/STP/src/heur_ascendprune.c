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

/**@file   heur_ascendprune.c
 * @brief  reduction-based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reduction and dual-cost based heuristic for Steiner problems. See
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" (2016) by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_ascendprune.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_ascendprune.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "grph.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"

#define HEUR_NAME             "ascendprune"
#define HEUR_DESC             "Dual-cost reduction heuristic for Steiner problems"
#define HEUR_DISPCHAR         'A'
#define HEUR_PRIORITY         2
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_MAXFREQPRUNE     FALSE         /**< executions of the heuristic at maximum frequency?                               */
#define ASCENPRUNE_MINLPIMPROVE     0.05          /**< minimum percentual improvement of dual bound (wrt to gap) mandatory to execute heuristic */

#ifdef WITH_UG
extern
int getUgRank(void);
#endif

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             lastdualbound;      /**< dual bound after the previous run                                 */
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   nfailures;          /**< number of failures since last successful call                     */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAscendPrune)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurAscendPrune(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */

static
SCIP_DECL_HEURINIT(heurInitAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   /* initialize data */
   heurdata->nfailures = 0;
   heurdata->bestsolindex = -1;
   heurdata->lastdualbound = 0.0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA*    heurdata;
   SCIP_PROBDATA*    probdata;
   SCIP_VAR**        vars;
   SCIP_SOL*         bestsol;                        /* best known solution */
   GRAPH*            graph;
   SCIP_Real         dualbound;
   SCIP_Real         gap;
   SCIP_Real*        redcosts;
   SCIP_Bool         success;
   int       e;
   int       nnodes;
   int       nedges;
   int*      edgearrint;
   int*      nodearrint;
   STP_Bool*     nodearrchar;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);

   *result = SCIP_DIDNOTRUN;

   /* todo: delete this file and move to slack-prune */

   nedges = graph->edges;
   nnodes = graph->knots;
   success = FALSE;

   /* get best current solution */
   bestsol = SCIPgetBestSol(scip);

   /* no solution available? */
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* get dual bound */
   dualbound = SCIPgetDualbound(scip);

   /* no new best solution available? */
   if( heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) && !(heurdata->maxfreq) )
   {
      /* current optimality gap */
      gap = SCIPgetSolOrigObj(scip, bestsol) - dualbound;

      if( SCIPisLT(scip, dualbound - heurdata->lastdualbound, gap * ASCENPRUNE_MINLPIMPROVE ) )
         return SCIP_OKAY;
   }

   heurdata->lastdualbound = dualbound;

   /* allocate memory for ascent and prune */
   SCIP_CALL( SCIPallocBufferArray(scip, &redcosts, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );

   for( e = 0; e < nedges; e++ )
   {
      assert(SCIPvarIsBinary(vars[e]));

      /* variable is already fixed, we must not trust the reduced cost */
      if( SCIPvarGetLbLocal(vars[e]) + 0.5 > SCIPvarGetUbLocal(vars[e]) )
      {
         if( SCIPvarGetLbLocal(vars[e]) > 0.5 )
            redcosts[e] = 0.0;
         else
         {
            assert(SCIPvarGetUbLocal(vars[e]) < 0.5);
            redcosts[e] = FARAWAY;
         }
      }
      else
      {
         if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[e])) )
         {
            assert(!SCIPisDualfeasNegative(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = SCIPgetVarRedcost(scip, vars[e]);
         }
         else
         {
            assert(!SCIPisDualfeasPositive(scip, SCIPgetVarRedcost(scip, vars[e])));
            assert(SCIPisFeasEQ(scip, SCIPgetSolVal(scip, NULL, vars[e]), 1.0) || SCIPisDualfeasZero(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = 0.0;
         }
      }

      if( SCIPisLT(scip, redcosts[e], 0.0) )
         redcosts[e] = 0.0;

      assert(SCIPisGE(scip, redcosts[e], 0.0));
      assert(!SCIPisEQ(scip, redcosts[e], SCIP_INVALID));
   }

   /* perform ascent and prune */
   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, heur, graph, redcosts, edgearrint, nodearrint, graph->source, nodearrchar, &success, TRUE) );

   if( success )
   {
      heurdata->nfailures = 0;
      *result = SCIP_FOUNDSOL;
   }
   else
   {
      heurdata->nfailures++;
   }

   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free memory */
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArray(scip, &redcosts);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */


/** ascent and prune */
SCIP_RETCODE SCIPStpHeurAscendPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure or NULL */
   const GRAPH*          g,                  /**< the graph */
   const SCIP_Real*      redcosts,           /**< the reduced costs */
   int*                  edgearrint,         /**< int edges array to store solution */
   int*                  nodearrint,         /**< int vertices array for internal computations */
   int                   root,               /**< the root (used for dual ascent) */
   STP_Bool*             nodearrchar,        /**< STP_Bool vertices array for internal computations */
   SCIP_Bool*            solfound,           /**< has a solution been found? */
   SCIP_Bool             addsol              /**< should the solution be added to SCIP by this method? */
   )
{
   GRAPH* newgraph;
   SCIP_Real* nval;
   int* const mark = g->mark;
   int* const newedges = edgearrint;
   int* const nodechild = nodearrint;
   int* edgeancestor;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const int probtype = g->stp_type;
   int nnewnodes = 0;
   int nnewedges = 0;
   const SCIP_Bool pcmw = (probtype == STP_PCSPG || probtype == STP_MWCSP || probtype == STP_RPCSPG || probtype == STP_RMWCSP);
   SCIP_Bool success;

   assert(g != NULL);
   assert(scip != NULL);
   assert(redcosts != NULL);
   assert(edgearrint != NULL);
   assert(nodearrint != NULL);
   assert(nodearrchar != NULL);

   if( root < 0 )
      root = g->source;

   assert(Is_term(g->term[root]));
   assert(graph_valid(g));

   if( addsol )
   {
      const int nvars = SCIPprobdataGetNVars(scip);
      SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );
   }
   else
   {
      nval = NULL;
   }

   /* DFS to identify 0-redcost subgraph */
   {
      int* const queue = nodearrint;
      STP_Bool* const scanned = nodearrchar;
      int qsize;

      /*
       * construct new graph corresponding to zero cost paths from the root to all terminals
       */

      BMSclearMemoryArray(mark, nnodes);
      BMSclearMemoryArray(scanned, nnodes);

      qsize = 0;
      mark[root] = TRUE;
      queue[qsize++] = root;
      nnewnodes++;

      /* DFS */
      while( qsize )
      {
         const int k = queue[--qsize];
         scanned[k] = TRUE;

         /* traverse outgoing arcs */
         for( int a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];

            if( SCIPisZero(scip, redcosts[a]) )
            {
               if( pcmw && k == root && Is_term(g->term[head]) )
                  continue;

               /* vertex not labeled yet? */
               if( !mark[head] )
               {
                  mark[head] = TRUE;
                  nnewnodes++;
                  queue[qsize++] = head;
               }

               if( (!scanned[head] || !SCIPisZero(scip, redcosts[flipedge(a)])) && k != root )
               {
                  assert(g->tail[a] != root);
                  assert(g->head[a] != root);

                  newedges[nnewedges++] = a;
               }
            }
         }
      }

#ifndef NDEBUG
      for( int k = 0; k < nnewedges && pcmw; k++ )
      {
         const int e = newedges[k];
         assert(!(g->tail[e] == root && Is_pterm(g->term[g->head[e]])));
         assert(!(g->head[e] == root && Is_pterm(g->term[g->tail[e]])));
      }
#endif

      for( int a = g->outbeg[root]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];
         if( mark[head] )
            newedges[nnewedges++] = a;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &edgeancestor, 2 * nnewedges) );

   /* initialize new graph */
   SCIP_CALL( graph_init(scip, &newgraph, nnewnodes, 2 * nnewedges, 1) );

   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      newgraph->stp_type = STP_SPG;
   else
      newgraph->stp_type = probtype;

   if( pcmw )
      SCIP_CALL( graph_pc_init(scip, newgraph, nnewnodes, nnewnodes) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( mark[k] )
      {
         if( pcmw )
         {
            if( (!Is_term(g->term[k])) )
               newgraph->prize[newgraph->knots] = g->prize[k];
            else
               newgraph->prize[newgraph->knots] = 0.0;
         }

         nodechild[k] = newgraph->knots;
         graph_knot_add(newgraph, g->term[k]);
      }
   }

   if( pcmw )
   {
      newgraph->norgmodelknots = nnewnodes;
      newgraph->extended = TRUE;
   }

   assert(nnewnodes == newgraph->knots);

   /* set root of new graph */
   newgraph->source = nodechild[root];
   assert(newgraph->source >= 0);

   if( g->stp_type == STP_RPCSPG || g->stp_type == STP_RMWCSP )
      newgraph->prize[newgraph->source] = FARAWAY;

   /* add edges to new graph */
   for( int a = 0; a < nnewedges; a++ )
   {
      int i;
      const int e = newedges[a];
      const int tail = nodechild[g->tail[e]];
      const int head = nodechild[g->head[e]];

      assert(tail >= 0);
      assert(head >= 0);

      for( i = newgraph->outbeg[tail]; i != EAT_LAST; i = newgraph->oeat[i] )
         if( newgraph->head[i] == head )
            break;

      if( i == EAT_LAST )
      {
         edgeancestor[newgraph->edges] = e;
         edgeancestor[newgraph->edges + 1] = flipedge(e);

         if( pcmw )
            graph_pc_updateTerm2edge(newgraph, g, tail, head, g->tail[e], g->head[e]);

         graph_edge_add(scip, newgraph, tail, head, g->cost[e], g->cost[flipedge(e)]);
      }
   }

   nnewedges = newgraph->edges;
   newgraph->norgmodeledges = nnewedges;

   assert(!pcmw || -1 == newgraph->term2edge[newgraph->source]);

   /* initialize ancestors of new graph edges */
   SCIP_CALL( graph_init_history(scip, newgraph) );

   /* initialize shortest path algorithm */
   SCIP_CALL( graph_path_init(scip, newgraph) );

   SCIP_CALL( level0(scip, newgraph) );

#ifdef DEBUG_ASCENDPRUNE
   for( int k = 0; k < nnodes && !pcmw; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         const int i = nodechild[k];
         if( i < 0 )
         {
            printf("k %d root %d \n", k, root);
            printf("FAIL in AP \n\n\n");
            return SCIP_ERROR;
         }

         if( newgraph->grad[i] == 0 && newgraph->knots > 1 )
         {
            printf("FAIL GRAD \n\n\n");
            return SCIP_ERROR;

         }
      }
   }
#endif
   assert(graph_valid(newgraph));

   /* get solution on new graph by PRUNE heuristic */
   SCIP_CALL( SCIPStpHeurPruneRun(scip, NULL, newgraph, newedges, &success, FALSE, TRUE) );

#ifdef DEBUG_ASCENDPRUNE
   for( int k = 0; k < newgraph->knots; ++k )
   {
      if( Is_term(newgraph->term[k]) && newgraph->grad[k] == 0 && k != newgraph->source )
      {
         printf("after i %d r %d \n", k, root);
         return SCIP_ERROR;
      }
   }

   if( !graph_sol_valid(scip, newgraph, newedges) )
   {
      printf("not valid %d \n", 0);
      return SCIP_ERROR;
   }
#endif
   if( !success )
   {
      SCIPdebugMessage("failed to build tree in ascend-prune (by prune) \n");
      goto TERMINATE;
   }

   assert(success && graph_sol_valid(scip, newgraph, newedges));

   SCIPdebugMessage("obj after prune %f \n", graph_sol_getObj(newgraph->cost, newedges, 0.0, newgraph->edges));

   SCIP_CALL( SCIPStpHeurLocalRun(scip, newgraph, newgraph->cost, newedges) );

   SCIPdebugMessage("obj after local %f \n", graph_sol_getObj(newgraph->cost, newedges, 0.0, newgraph->edges));

   assert(graph_sol_valid(scip, newgraph, newedges));
   graph_path_exit(scip, newgraph);


    /*
    * prune solution (in the original graph)
    */

   BMSclearMemoryArray(nodearrchar, nnodes);

   for( int e = 0; e < nnewedges; e++ )
      if( newedges[e] == CONNECT )
      {
         const int eorg = edgeancestor[e];
         nodearrchar[g->tail[eorg]] = TRUE;
         nodearrchar[g->head[eorg]] = TRUE;
      }

   for( int e = 0; e < nedges; e++ )
      newedges[e] = UNKNOWN;

   if( pcmw )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, g, g->cost, newedges, nodearrchar) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, g, g->cost, 0, newedges, nodearrchar) );

   assert(graph_sol_valid(scip, g, newedges));

   if( addsol )
   {
      assert(nval != NULL);
      for( int e = 0; e < nedges; e++ )
      {
         if( newedges[e] == CONNECT )
            nval[e] = 1.0;
         else
            nval[e] = 0.0;
      }
   }

   success = graph_sol_valid(scip, g, newedges);

   if( success && addsol )
   {
      /* add solution */
      SCIP_SOL* sol = NULL;
      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );
      SCIPdebugMessage("Ascend-and-prune added solution \n");
   }

   *solfound = success;

 TERMINATE:

   for( int k = 0; k < nnodes; k++ )
      mark[k] = (g->grad[k] > 0);

   /* free memory */
   graph_free(scip, &newgraph, TRUE);
   SCIPfreeBufferArray(scip, &edgeancestor);
   SCIPfreeBufferArrayNull(scip, &nval);

   return SCIP_OKAY;
}


/** creates the prune primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurAscendPrune(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create prune primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAscendPrune, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAscendPrune) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAscendPrune) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAscendPrune) );

   /* add ascend and prune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_MAXFREQPRUNE, NULL, NULL) );

   return SCIP_OKAY;
}
