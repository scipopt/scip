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

/**@file   heur_prune.c
 * @brief  reduction-based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reduction based heuristic for Steiner problems. See
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" (2016) by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_prune.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_prune.h"
#include "heur_local.h"
#include "grph.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"

#define HEUR_NAME             "prune"
#define HEUR_DESC             "Reduction based heuristic for Steiner problems"
#define HEUR_DISPCHAR         'P'
#define HEUR_PRIORITY         2
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE            /**< does the heuristic use a secondary SCIP instance?                  */

#define DEFAULT_PRUNE_MAXFRQ      FALSE       /**< executions of the heuristic at maximum frequency?                  */
#define DEFAULT_PRUNE_TMRUNS     100          /**< number of runs in TM heuristic when called by prune heuristic      */
#define PRUNE_MINREDELIMS        2            /**< maximum number of eliminations for reduction package when called by prune heuristic */
#define PRUNE_MAXREDROUNDS       6            /**< maximum number of reduction rounds in prune heuristic */
#define BREAKONERROR  FALSE
#define MAXNTERMINALS 500
#define MAXNEDGES     10000

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
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   ntmruns;            /**< number of runs in TM heuristic                                    */
   int                   nfailures;          /**< number of failures since last successful call                     */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};

/*
 * Local methods
 */



/* set node solution array and get solution value */
static
void setNodeSolArray(
   const GRAPH* g,
   SCIP_Real* uborg,
   int* solnode,
   const int* soledge
   )
{
   SCIP_Real ub;
   int e;
   int k;
   int nedges;
   int nnodes;

   assert(g != NULL);
   assert(solnode != NULL);
   assert(soledge != NULL);

   ub = 0.0;
   nedges = g->edges;
   nnodes = g->knots;

   for( k = 0; k < nnodes; k++ )
      solnode[k] = UNKNOWN;

   for( e = 0; e < nedges; e++ )
   {
      if( soledge[e] == CONNECT )
      {
         ub += g->cost[e];
         solnode[g->tail[e]] = CONNECT;
         solnode[g->head[e]] = CONNECT;
      }
   }

   if( uborg != NULL)
      *uborg = ub;
}

/** computes new solution during execution of prune and updates best global one if possible */
static
SCIP_RETCODE computeNewSols(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   GRAPH*                prunegraph,         /**< pruned graph data structure */
   PATH*                 path,               /**< shortest path struct */
   int*                  nodearrint,         /**< array */
   int*                  edgearrint,         /**< array */
   int*                  solnode,            /**< array for best solution nodes wrt prunegraph */
   int*                  soledge,            /**< array for best solution edges wrt prunegraph */
   int*                  globalsoledge,      /**< array storing best solution wrt g */
   STP_Bool*             nodearrchar,        /**< array */
   SCIP_Real*            globalobj,          /**< pointer to objective value of best solution wrt g */
   SCIP_Bool             incumbentgiven,     /**< incumbent solution for pruned graph given? */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
   int* const pmark = prunegraph->mark;
   SCIP_Real dummyhop = 0.0;
   int nruns;
   int best_start;
   const int nnodes = g->knots;

   assert(graph_valid(prunegraph));
   assert(g->edges == prunegraph->edges);

   /*
    * compute new solution on pruned graph
    */

   nruns = 0;
   for( int k = 0; k < nnodes; k++ )
   {
      pmark[k] = (prunegraph->grad[k] > 0);
      if( pmark[k] )
         nruns++;
   }

   nruns = MIN(nruns, DEFAULT_PRUNE_TMRUNS);
   SCIPStpHeurTMCompStarts(prunegraph, nodearrint, &nruns);

   /* run shortest path heuristic */
   SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, prunegraph, nodearrint, &best_start, soledge, nruns,
         prunegraph->source, prunegraph->cost, prunegraph->cost, &dummyhop, NULL, 0.0, success, FALSE));

   SCIP_CALL( SCIPStpHeurLocalRun(scip, prunegraph, prunegraph->cost, soledge) );

   SCIP_CALL( SCIPStpHeurPruneUpdateSols(scip, g, prunegraph, path, nodearrint, edgearrint, solnode, soledge,
         globalsoledge, nodearrchar, globalobj, incumbentgiven, success) );

   return SCIP_OKAY;
}

/* get reduction bound */
static
int getRedBound(
   int nround,
   int nedges
   )
{
   if( nround == 0 )
      return(MAX(nedges / 1000, PRUNE_MINREDELIMS));
   if( nround == 1 )
      return(MAX(nedges / 500, PRUNE_MINREDELIMS));
   return(MAX(nedges / 200, PRUNE_MINREDELIMS));
}


static
void setMinMaxElims(
   SCIP*                 scip,
   int*                  minnelims,
   int*                  lminnelims,
   int                   annodes,
   int                   anedges,
   int                   anterms,
   int                   nround,
   int                   maxrounds
   )
{
   int min;
   int totminnelims;
   SCIP_Real factor;

   assert(scip != NULL);
   assert(minnelims != NULL);
   assert(lminnelims != NULL);
   assert(annodes > 0);
   assert(nround >= 0);
   assert(maxrounds >= 1);

   anedges = anedges / 2;

   totminnelims = MAX(PRUNE_MINREDELIMS, (anedges / 20));

   min = (anedges / 10);

   min -= (int) ( ((SCIP_Real) min * anterms) / ( (SCIP_Real) annodes) );
   min = MAX(min, 1);

   factor = (SCIP_Real) anedges / min;
   factor = ((SCIP_Real) nround / (3.0 * maxrounds)) * factor;

   if( SCIPisGT(scip, factor, 1.0) )
      min = (int) ((SCIP_Real) min * factor);

   min = MAX(totminnelims, min);
   min = MIN(min, (anedges - 1));
   min = MAX(min, 1);

   *lminnelims = min / 2;
   *lminnelims = MAX(*lminnelims, 1);

   *minnelims = min;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyPrune)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurPrune(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreePrune)
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
SCIP_DECL_HEURINIT(heurInitPrune)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   /* initialize data */
   heurdata->ntmruns = DEFAULT_PRUNE_TMRUNS;
   heurdata->nfailures = 0;
   heurdata->bestsolindex = -1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolPrune)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;                        /* best known solution */
   GRAPH* graph;
   SCIP_Real* xval;
   SCIP_Bool success;
   int e;
   int nvars;
   int nedges;
   int*  soledge;

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

   nedges = graph->edges;
   *result = SCIP_DIDNOTRUN;

   /* if not STP like variant, return */
   if( graph->stp_type != STP_SPG && graph->stp_type != STP_RSMT && graph->stp_type != STP_OARSMT &&
      graph->stp_type != STP_GSTP )
      return SCIP_OKAY;

   if( (graph->edges > MAXNEDGES) && (graph->terms > MAXNTERMINALS) )
      return SCIP_OKAY;

   /* get best current solution */
   bestsol = SCIPgetBestSol(scip);

   /* no solution available? */
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* has the new solution been found by this very heuristic? */
   if( SCIPsolGetHeur(bestsol) == heur || heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   xval = SCIPprobdataGetXval(scip, bestsol);

   if( xval == NULL )
      return SCIP_OKAY;

   vars = SCIPprobdataGetVars(scip);
   nvars = SCIPprobdataGetNVars(scip);

   SCIPdebugMessage("START prune heuristic \n");

   assert(vars != NULL);

   /* allocate array to store primal solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &soledge, nedges) );

   for( e = 0; e < nedges; e++ )
   {
      if( SCIPisEQ(scip, xval[e], 1.0) )
         soledge[e] = CONNECT;
      else
         soledge[e] = UNKNOWN;
   }

   /* execute prune heuristic */
   SCIP_CALL( SCIPStpHeurPruneRun(scip, vars, graph, soledge, &success, TRUE, FALSE) );

   /* solution found by prune heuristic? */
   if( success )
   {
      SCIP_SOL* sol;
      SCIP_Real pobj;
      SCIP_Real* nval;

      SCIPdebugMessage("ADDED valid solution in prune \n");

      assert(graph_sol_valid(scip, graph, soledge));
      pobj = 0.0;

      /* allocate memory to store solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

      for( e = 0; e < nedges; e++ )
      {
         if( soledge[e] == CONNECT )
         {
            nval[e] = 1.0;
            pobj += graph->cost[e];
         }
         else
         {
            nval[e] = 0.0;
         }
      }

      SCIPdebugMessage("prune, best: old %f, new %f \n  \n",  SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj);

      /* try to add new solution to pool */
      sol = NULL;
      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

      /* has solution been added? */
      if( success )
      {
         SCIPdebugMessage("solution added by PRUNE \n  \n");

         *result = SCIP_FOUNDSOL;

         assert(graph_sol_valid(scip, graph, soledge));

         /* is the solution the new incumbent? */
         if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj) )
         {
            heurdata->nfailures = 0;
         }
         else
         {
            success = FALSE;
         }
      }

      /* solution could not be added or is not best? */
      if( !success )
         heurdata->nfailures++;

      /* free memory */
      SCIPfreeBufferArray(scip, &nval);
   }

   /* store index of incumbent solution */
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free memory */
   SCIPfreeBufferArray(scip, &soledge);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** updates solutions for pruned graph */
SCIP_RETCODE SCIPStpHeurPruneUpdateSols(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   GRAPH*                prunegraph,         /**< pruned graph data structure */
   PATH*                 path,               /**< shortest path struct */
   int*                  nodearrint,         /**< array */
   int*                  edgearrint,         /**< array */
   int*                  solnode,            /**< array for best solution nodes wrt prunegraph */
   int*                  soledge,            /**< array for best solution edges wrt prunegraph */
   int*                  globalsoledge,      /**< array storing best solution wrt g */
   STP_Bool*             nodearrchar,        /**< array */
   SCIP_Real*            globalobj,          /**< pointer to objective value of best solution wrt g */
   SCIP_Bool             incumbentgiven,     /**< incumbent solution for pruned graph given? */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   )
{
   SCIP_Real objnew;
   int* const pmark = prunegraph->mark;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const int probtype = g->stp_type;
   const SCIP_Bool pcmw = (probtype == STP_MWCSP || probtype == STP_RMWCSP || probtype == STP_RPCSPG || probtype == STP_PCSPG);

   assert(g != NULL);
   assert(scip != NULL);
   assert(soledge != NULL);
   assert(path != NULL);
   assert(solnode != NULL);
   assert(edgearrint != NULL);
   assert(nodearrint != NULL);
   assert(nodearrchar != NULL);

   /*
    * compare new solution on pruned graph with (reconstructed) incumbent
    */

   if( incumbentgiven )
   {
      SCIP_Real objold;

      objnew = graph_sol_getObj(prunegraph->cost, soledge, 0.0, nedges);

      if( pcmw )
         SCIP_CALL( SCIPStpHeurTMBuildTreePcMw(scip, prunegraph, path, prunegraph->cost, &objold, solnode) );
      else
         SCIP_CALL( SCIPStpHeurTMBuildTree(scip, prunegraph, path, prunegraph->cost, &objold, solnode) );

      SCIPdebugMessage("objold %f objnew %f \n", objold, objnew);

      /* keep (reconstructed) old solution? */
      if( objold < objnew )
      {
         for( int e = 0; e < nedges; e++ )
            soledge[e] = UNKNOWN;
         for( int k = 0; k < nnodes; k++ )
         {
            const int e = path[k].edge;

            if( e >= 0 )
               soledge[e] = CONNECT;
         }
      }
   }

   assert(graph_sol_valid(scip, prunegraph, soledge));

   setNodeSolArray(prunegraph, NULL, solnode, soledge);

   /*
    * retransform new solution and compare with best global one
    */

   SCIP_CALL( graph_sol_getOrg(scip, prunegraph, g, soledge, edgearrint) );

   assert(graph_sol_valid(scip, g, edgearrint));

#if BREAKONERROR
   if( !graph_sol_valid(scip, g, edgearrint) )
   {
      printf("TM orig sol in prune not valid %d \n", 0);
      exit(1);
   }
#endif

   objnew = graph_sol_getObj(g->cost, edgearrint, 0.0, nedges);

   SCIPdebugMessage("old global obj: %f ... new global obj: %f \n", *globalobj, objnew);

   if( objnew < (*globalobj) )
   {
      SCIPdebugMessage("new global solution is better \n");
      *globalobj = objnew;

      SCIP_CALL( SCIPStpHeurLocalRun(scip, g, g->cost, edgearrint) );

      objnew = graph_sol_getObj(g->cost, edgearrint, 0.0, nedges);

      assert(SCIPisLE(scip, objnew, *globalobj));

      if( objnew < (*globalobj) )
         *globalobj = objnew;

      BMScopyMemoryArray(globalsoledge, edgearrint, nedges);
   }

   assert(*globalobj < FARAWAY);

   for( int k = 0; k < nnodes; k++ )
      pmark[k] = (prunegraph->grad[k] > 0);

   *success = TRUE;

   return SCIP_OKAY;
}

/** execute prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL (need to be provided whenever available) */
   GRAPH*                g,                  /**< the graph */
   int*                  soledge,            /**< array to store primal solution (if no solution is provided,
                                                withinitialsol must be set to FALSE) */
   SCIP_Bool*            success,            /**< feasible solution found? */
   const SCIP_Bool       withinitialsol,     /**< solution given? */
   const SCIP_Bool       reducegraph         /**< try to reduce graph initially? */
   )
{
   GRAPH* prunegraph;
   PATH* vnoi;
   PATH* path;
   SCIP_Real globalobj;
   SCIP_Real uborg = 0.0;
   SCIP_Real offset;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real* nodearrreal;
   const int probtype = g->stp_type;
   const SCIP_Bool pc = (probtype == STP_RPCSPG || probtype == STP_PCSPG);
   const SCIP_Bool mw = (probtype == STP_MWCSP || probtype == STP_RMWCSP);
   const SCIP_Bool pcmw = (pc || mw);
   const int nnodes = g->knots;
   const int nedges = g->edges;
   int annodes;
   int anedges;
   int anterms;
   int reductbound;
   int* heap;
   int* state;
   int* vbase;
   int* solnode;
   int* nodearrint;
   int* edgearrint;
   int* nodearrint2;
   int* globalsoledge;
   STP_Bool* nodearrchar;
   SCIP_Bool solgiven;

   assert(g != NULL);
   assert(scip != NULL);
   assert(soledge != NULL);
   assert(probtype != STP_DHCSTP);
   assert(!(withinitialsol && reducegraph));

   *success = TRUE;

   if( reducegraph )
      solgiven = FALSE;
   else
      solgiven = withinitialsol;

   if( g->terms <= 1)
   {
      for( int e = 0; e < nedges; e++ )
         soledge[e] = UNKNOWN;
      return SCIP_OKAY;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &globalsoledge, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes + 1) );

   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      g->stp_type = STP_SPG;

   /* copy the graph */
   SCIP_CALL( graph_copy(scip, g, &prunegraph) );

   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      g->stp_type = probtype;

   if( vars != NULL )
   {
      prunegraph->norgmodeledges = prunegraph->edges;
      prunegraph->norgmodelknots = prunegraph->knots;
   }

   /* set ancestors of the new graph */
   SCIP_CALL( graph_init_history(scip, prunegraph) );

   /* set offset (within new graph) to 0.0 */
   offset = 0.0;

   reductbound = getRedBound(0, nedges);

   SCIPdebugMessage("Starting prune... \n");

   assert(!pcmw || g->extended);

   /* problem variables given? */
   if( vars != NULL )
   {
      int nfixedges = 0;

      /* delete fixed edges from the new graph */
      for( int e = 0; e < nedges; e += 2 )
      {
         /* both e and its anti-parallel edge fixed to zero? */
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5 )
         {
            if( pcmw )
            {
               if( (!solgiven || (soledge[e] != CONNECT && soledge[e + 1] != CONNECT))
                     && !Is_term(prunegraph->head[e]) && !Is_term(prunegraph->tail[e]) )
               {
                  graph_edge_del(scip, prunegraph, e, TRUE);
                  nfixedges++;
               }
            }
            else
            {
               if( !solgiven || (soledge[e] != CONNECT && soledge[e + 1] != CONNECT) )
               {
                  graph_edge_del(scip, prunegraph, e, TRUE);
                  nfixedges++;
               }
            }
         }
      }
      SCIPdebugMessage("fixed edges in prune: %d \n", nfixedges);

      if( nfixedges >= reductbound )
      {
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
         reductbound = getRedBound(0, anedges);
      }
   }

   SCIP_CALL(graph_path_init(scip, prunegraph));

   /* perform initial reductions? */
   if( reducegraph )
   {
      if( pc )
         SCIP_CALL( redLoopPc(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
               vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, &offset, FALSE, FALSE, reductbound, FALSE) );
      else if( mw )
         SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, state,
               vbase, nodearrint, edgearrint, nodearrint2, heap, NULL, nodearrchar, &offset, FALSE, FALSE, FALSE, reductbound, FALSE) );
      else
         SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
               vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, &offset, -1.0, FALSE, FALSE, TRUE, reductbound, FALSE, FALSE) );
   }

   /* get number of remaining nodes, edges and terminals */
   graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

   if( solgiven )
   {
      BMScopyMemoryArray(globalsoledge, soledge, nedges);
      globalobj = graph_sol_getObj(g->cost, soledge, 0.0, nedges);
      setNodeSolArray(g, &uborg, solnode, soledge);
   }
   else
   {
      globalobj = FARAWAY;
      SCIP_CALL( computeNewSols(scip, g, prunegraph, path, nodearrint, edgearrint, solnode, soledge, globalsoledge, nodearrchar, &globalobj, FALSE, success) );

      assert(success);
   }

   SCIPdebugMessage("initially reduced graph: |V| = %d, |E| = %d, |T| = %d  \n", annodes, anedges, anterms);
   SCIPdebugMessage("entering prune \n\n");

   if( annodes > 0 )
   {
      /* main prune reduction loop */
      for( int i = 0; i < PRUNE_MAXREDROUNDS; i++ )
      {
         int minnelims = 0;
         int brednelims = 0;
         int lminnelims = 0;

         /* get number of remaining edges */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         if( anterms <= 2 )
         {
            SCIPdebugMessage("less than two terminals, break !! \n");
            break;
         }

         setMinMaxElims(scip, &minnelims, &lminnelims, annodes, anedges, anterms, i + 1, PRUNE_MAXREDROUNDS);

         if( i > 0 )
         {
            SCIP_CALL( computeNewSols(scip, g, prunegraph, path, nodearrint, edgearrint, solnode, soledge, globalsoledge,
                  nodearrchar, &globalobj, TRUE, success) );
         }

         if( pcmw )
            graph_pc_2org(prunegraph);

         SCIP_CALL( reduce_boundPrune(scip, prunegraph, vnoi, cost, nodearrreal, costrev,
               &offset, heap, state, vbase, solnode, soledge, &brednelims, minnelims));

         if( pcmw )
            graph_pc_2trans(prunegraph);

#ifdef SCIP_DEBUG
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
         printf("PRUNE round: %d edges %d  nodes: %d \n", i, anedges / 2, annodes);
         printf("PRUNE round: %d minelims %d  really reduced: %d \n", i, minnelims, brednelims);
#endif

         /* not enough reductions possible? */
         if( brednelims < lminnelims  )
         {
            SCIPdebugMessage("not enough elims in PRUNE; exit %d \n\n", brednelims);
            i = PRUNE_MAXREDROUNDS;
         }

         /* delete all vertices not reachable from the root */
         SCIP_CALL( level0(scip, prunegraph) );

         assert(graph_valid(prunegraph));

         /* is the reduced graph still feasible? */
         if( !graph_valid(prunegraph) )
            break;

         /* get number of remaining edges */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         reductbound = getRedBound(i + 1, anedges);

         /* reduce graph */
         if( pc )
            SCIP_CALL( redLoopPc(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
                  vbase, nodearrint, edgearrint, nodearrint2, solnode, nodearrchar, &offset, FALSE, FALSE, reductbound, FALSE) );
         else if( mw )
            SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, state,
                  vbase, nodearrint, edgearrint, nodearrint2, heap, solnode, nodearrchar, &offset, FALSE, FALSE, FALSE, reductbound, FALSE) );
         else
            SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state, vbase, nodearrint, edgearrint,
                  nodearrint2, solnode, nodearrchar, &offset, -1.0, FALSE, FALSE, TRUE, reductbound, FALSE, FALSE));

         /* delete all vertices not reachable from the root */
         SCIP_CALL( level0(scip, prunegraph) );

         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         if( anterms <= 2 )
         {
            SCIPdebugMessage("less than two terminals, break !! \n");
            SCIP_CALL( computeNewSols(scip, g, prunegraph, path, nodearrint, edgearrint, solnode, soledge, globalsoledge,
                              nodearrchar, &globalobj, TRUE, success) );
            break;
         }
      } /* reduction loop */
   } /* main prune loop */

   for( int k = 0; k < nnodes; k++ )
      nodearrchar[k] = FALSE;

   graph_path_exit(scip, prunegraph);

   *success = graph_sol_valid(scip, g, globalsoledge);
   assert(*success);

   SCIPdebugMessage("obj of best prune sol: %f \n", graph_sol_getObj(g->cost, globalsoledge, 0.0, nedges));

   BMScopyMemoryArray(soledge, globalsoledge, nedges);

#if BREAKONERROR
   if( !(*success) )
   {
      printf("final sol not valid %d \n", 0);
      return SCIP_ERROR;
   }
#endif

   /* free memory */
   graph_free(scip, &prunegraph, TRUE);

   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &solnode);
   SCIPfreeBufferArray(scip, &globalsoledge);

   return SCIP_OKAY;
}


/** creates the prune primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurPrune(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPrune, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPrune) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePrune) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitPrune) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolPrune) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolPrune) );

   /* add prune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_PRUNE_MAXFRQ, NULL, NULL) );

   return SCIP_OKAY;
}
