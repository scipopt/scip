/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE            /**< does the heuristic use a secondary SCIP instance?                  */

#define DEFAULT_PRUNE_MAXFRQ      FALSE       /**< executions of the heuristic at maximum frequency?                  */
#define DEFAULT_PRUNE_TMRUNS     100          /**< number of runs in TM heuristic when called by prune heuristic      */
#define PRUNE_MINREDELIMS        2            /**< maximum number of eliminations for reduction package when called by prune heuristic */
#define PRUNE_MAXREDROUNDS       6            /**< maximum number of reduction rounds in prune heuristic */
#define BREAKONERROR FALSE
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

/* get reduction bound */
static
inline int getRedBound(
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

/* set solution array and get solution value */
static
void setSolArray(
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
   assert(uborg != NULL);
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

   *uborg += ub;
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
   STP_Bool* stnodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, graph->knots ) );

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->knots; e++ )
   {
      stnodes[e] = FALSE;
   }
   for( e = 0; e < graph->edges; e++ )
   {
      if( result[e] == CONNECT )
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
         if( n == graph->source )
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
      if( result[e] == CONNECT )
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
#ifdef SCIP_DEBUG
      printf("prune, best: old %f, new %f \n  \n",  SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj);
#endif
      /* try to add new solution to pool */
      sol = NULL;
      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

      /* has solution been added? */
      if( success )
      {
#ifdef SCIP_DEBUG
         printf("solution added by PRUNE \n  \n");
#endif
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

/** execute prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL (need to be provided whenever available) */
   GRAPH*                g,                  /**< the graph */
   int*                  soledge,            /**< array to store primal solution (if no solution is provided,
                                                solgiven must be set to FALSE) */
   SCIP_Bool*            success,            /**< feasible solution found? */
   const SCIP_Bool       solgiven,           /**< solution given? */
   const SCIP_Bool       reducegraph         /**< try to reduce graph initially? */
   )
{
   GRAPH* prunegraph;
   PATH* vnoi;
   PATH* path;
   SCIP_Real bestobj;
   SCIP_Real uborg;
   SCIP_Real hopfactor;
   SCIP_Real offsetnew;
   SCIP_Real offsetold;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real* nodearrreal;
   IDX** ancestors;
   const int probtype = g->stp_type;
   const SCIP_Bool pc = (probtype == STP_RPCSPG || probtype == STP_PCSPG);
   const SCIP_Bool mw = (probtype == STP_MWCSP || probtype == STP_RMWCSP);
   const SCIP_Bool pcmw = (pc || mw);
   SCIP_Bool solvalid;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   int annodes;
   int anedges;
   int anterms;
   int minnelims;
   int lminnelims;
   int best_start;
   int reductbound;
   int* heap;
   int* state;
   int* vbase;
   int* solnode;
   int* nodearrint;
   int* edgearrint;
   int* nodearrint2;
   int* bestsoledges;
   STP_Bool* nodearrchar;

   solvalid = TRUE;
   *success = TRUE;
   hopfactor = DEFAULT_HOPFACTOR;

   assert(g != NULL);
   assert(scip != NULL);
   assert(soledge != NULL);
   assert(probtype != STP_DHCSTP);

   if( g->terms <= 1)
   {
      for( int e = 0; e < nedges; e++ )
         soledge[e] = UNKNOWN;
      return SCIP_OKAY;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestsoledges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

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
   offsetold = 0.0;

   reductbound = getRedBound(0, nedges);

   /* problem variables given? */
   if( vars != NULL )
   {
      int nfixedges = 0;

      /* delete fixed edges from the new graph */

      if( pcmw )
      {
         for( int e = 0; e < nedges; e += 2 )
            /* both e and its anti-parallel edge fixed to zero? */
            if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5 )
            {
               if( (!solgiven || (soledge[e] != CONNECT && soledge[e + 1] != CONNECT)) && !Is_term(prunegraph->head[e]) && !Is_term(prunegraph->tail[e]) )
                  graph_edge_del(scip, prunegraph, e, TRUE);
               nfixedges++;
            }
      }
      else
      {
         for( int e = 0; e < nedges; e += 2 )
            /* both e and its anti-parallel edge fixed to zero? */
            if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5 )
            {
               if( !solgiven || (soledge[e] != CONNECT && soledge[e + 1] != CONNECT) )
                  graph_edge_del(scip, prunegraph, e, TRUE);
               nfixedges++;
            }
      }

      SCIPdebugMessage("fixed edges in prune: %d \n", nfixedges);

      if( nfixedges >= reductbound )
      {
         /* get number of remaining edges */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
         reductbound = getRedBound(0, anedges);
      }
   }

   SCIP_CALL(graph_path_init(scip, prunegraph));

   /* perform initial reductions? */
   if( reducegraph )
   {
      if( solgiven )
         /* set solution node array and get solution value */
         setSolArray(g, &uborg, solnode, soledge);

      if( pc )
      {
         SCIP_CALL( redLoopPc(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
               vbase, nodearrint, edgearrint, nodearrint2, (solgiven) ? solnode : NULL, nodearrchar, &offsetold, FALSE, FALSE, reductbound, FALSE) );
      }
      else if( mw )
      {
         SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, state,
               vbase, nodearrint, edgearrint, nodearrint2, heap, (solgiven) ? solnode : NULL, nodearrchar, &offsetold, FALSE, FALSE, FALSE, reductbound, FALSE) );
      }
      else
      {
         SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
               vbase, nodearrint, edgearrint, nodearrint2, (solgiven) ? solnode : NULL, nodearrchar, &offsetold, -1.0, FALSE, FALSE, TRUE, reductbound, NULL, FALSE) );
      }
   }

   /* get number of remaining nodes, edges and terminals */
   graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

   SCIPdebugMessage("entering prune \n\n");

   /* no solution provided and graph not completely reduced? */
   if( !solgiven && annodes > 0 )
   {
      /*
       * compute new solution on given graph
       */

      if( reducegraph )
      {
         int nruns = 0;
         int* tmstarts;

         for( int k = 0; k < nnodes; k++ )
         {
            prunegraph->mark[k] = (prunegraph->grad[k] > 0);
            if( prunegraph->mark[k] )
               nruns++;
         }

         nruns = MIN(nruns, DEFAULT_PRUNE_TMRUNS);

         if( !pcmw )
         {
            SCIPStpHeurTMCompStarts(prunegraph, nodearrint, &nruns);

            for( int k = 0; k < nruns; k++ )
               assert(prunegraph->grad[nodearrint[k]] > 0 );
            tmstarts = nodearrint;
         }
         else
         {
            if( nruns > 1 )
               nruns--;
            tmstarts = NULL;
         }

         assert(graph_valid(prunegraph));

         /* run shortest path heuristic */
         SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, prunegraph, tmstarts, &best_start, edgearrint, nruns,
               prunegraph->source, prunegraph->cost, prunegraph->cost, &hopfactor, NULL, 0.0, success, FALSE) );

         if( pcmw )
         {
            SCIP_Bool dummy;
            SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, prunegraph, prunegraph->cost, path, edgearrint, nodearrint, nodearrchar, &dummy) );
         }

         setSolArray(prunegraph, &uborg, solnode, edgearrint);

         /* transform solution to original graph */

         ancestors = prunegraph->ancestors;

         for( int k = 0; k < nnodes; k++ )
            nodearrchar[k] = FALSE;

         for( int e = 0; e < nedges; e++ )
         {
            soledge[e] = UNKNOWN;

            if( edgearrint[e] == CONNECT )
               graph_sol_setNodeList(g, nodearrchar, ancestors[e]);
         }

         /* retransform edges fixed during graph reduction */
         graph_sol_setNodeList(g, nodearrchar, prunegraph->fixedges);

         if( pcmw )
         {
            for( int k = 0; k < nnodes; k++ )
               if( nodearrchar[k] )
               {
                  IDX* curr = prunegraph->pcancestors[k];
                  while( curr != NULL )
                  {
                     const int idx = curr->index;
                     nodearrchar[g->tail[idx]] = TRUE;
                     nodearrchar[g->head[idx]] = TRUE;

                     curr = curr->parent;
                  }
               }
         }

         /* prune solution (in the original graph) */
         if( pcmw )
            SCIP_CALL( SCIPStpHeurTMPrunePc(scip, g, g->cost, soledge, nodearrchar) );
         else
            SCIP_CALL( SCIPStpHeurTMPrune(scip, g, g->cost, 0, soledge, nodearrchar) );

#ifdef SCIP_DEBUG
         SCIP_Real tmpub = graph_sol_getObj(g->cost, soledge, 0.0, nedges);
         printf("first obj %f \n", tmpub);
#endif
         assert(graph_sol_valid(scip, g, soledge));
      }
      else
      {
         int nruns;

         nruns = MIN(annodes, DEFAULT_PRUNE_TMRUNS);

         assert(graph_valid(prunegraph));

         /* run shortest path heuristic */
         SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, prunegraph, NULL, &best_start, soledge, nruns,
               prunegraph->source, prunegraph->cost, prunegraph->cost, &hopfactor, NULL, 0.0, success, FALSE) );

         if( pcmw )
         {
            SCIP_Bool dummy;
            SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, prunegraph, prunegraph->cost, path, soledge, nodearrint, nodearrchar, &dummy) );
         }

         setSolArray(prunegraph, &uborg, solnode, soledge);

         for( int e = 0; e < nedges; e++ )
            edgearrint[e] = soledge[e];
      }

      assert(*success);
      assert(graph_sol_valid(scip, g, soledge));

      // todo
      if( !(*success) )
      {
         printf("failed to build initial tree in prune \n");
         return SCIP_ERROR;
      }
   }
   else if( annodes > 0 && !reducegraph )
   {
      BMScopyMemoryArray(edgearrint, soledge, nedges);
   }

   offsetnew = offsetold;

   if( annodes > 0 )
   {
      int brednelims;
      int* const pmark = prunegraph->mark;

      if( solgiven && !reducegraph )
         setSolArray(g, &uborg, solnode, soledge);

      /* main reduction loop */
      for( int i = 0; i < PRUNE_MAXREDROUNDS; i++ )
      {
         brednelims = 0;

         offsetold = offsetnew;

         /* get number of remaining edges */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         if( anterms <= 2 )
         {
            SCIPdebugMessage("less than two terminals, break !! \n");
            break;
         }

         setMinMaxElims(scip, &minnelims, &lminnelims, annodes, anedges, anterms, i + 1, PRUNE_MAXREDROUNDS);

         /* no initial solution available? */
         if( i > 0 || (reducegraph && solgiven) )
         {
            SCIP_Real objold;
            SCIP_Real objnew;
            int nruns = 0;

            SCIPdebugMessage("compute new temporary solution in prune  \n");

            for( int k = 0; k < nnodes; k++ )
            {
               pmark[k] = (prunegraph->grad[k] > 0);
               if( pmark[k] )
                  nruns++;
            }

            nruns = MIN(nruns, DEFAULT_PRUNE_TMRUNS);
            SCIPStpHeurTMCompStarts(prunegraph, nodearrint, &nruns);

            assert(graph_valid(prunegraph));

            /* run shortest path heuristic */
            SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, prunegraph, nodearrint, &best_start, edgearrint, nruns,
                  prunegraph->source, prunegraph->cost, prunegraph->cost, &hopfactor, NULL, 0.0, success, FALSE) );

            if( pcmw )
            {
               SCIP_Bool dummy;
               SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, prunegraph, prunegraph->cost, path, edgearrint, nodearrint, nodearrchar, &dummy) );
            }

            objnew = graph_sol_getObj(prunegraph->cost, edgearrint, 0.0, nedges);

            assert(graph_sol_valid(scip, prunegraph, edgearrint));

#if BREAKONERROR
            if( !graph_sol_valid(scip, prunegraph, edgearrint) )
            {
               printf("TM sol in prune not valid %d \n", 0);
               return SCIP_ERROR;
            }
#endif
            if( pcmw )
            {
               SCIP_CALL( SCIPStpHeurTMBuildTreePcMw(scip, prunegraph, path, cost, &objold, solnode) );
            }
            else
            {
               SCIP_CALL( SCIPStpHeurTMBuildTree(scip, prunegraph, path, cost, &objold, solnode) );
            }

#ifdef SCIP_DEBUG
            printf("objold %f objnew %f \n", objold, objnew);
#endif
            if( SCIPisLT(scip, objold, objnew ) )
            {
               int e;
               for( e = 0; e < nedges; e++ )
                  edgearrint[e] = UNKNOWN;
               for( int k = 0; k < nnodes; k++ )
               {
                  e = path[k].edge;

                  if( e >= 0 )
                     edgearrint[e] = CONNECT;
               }
            }

            assert(graph_sol_valid(scip, prunegraph, edgearrint));

            for( int k = 0; k < nnodes; k++ )
               pmark[k] = (prunegraph->grad[k] > 0);

            setSolArray(prunegraph, &hopfactor, solnode, edgearrint);
         }
         else
         {
            int count;

#if BREAKONERROR
            if( !graph_sol_valid(scip, prunegraph, edgearrint) )
            {
               printf("sol not valid %d \n", 0);
               return SCIP_ERROR;
            }
#endif
#ifdef SCIP_DEBUG
               printf("continue with temporary solution in prune  \n");
#endif
            /* prune solution */
            for( int k = 0; k < nnodes; k++ )
               solnode[k] = UNKNOWN;

            for( int e = 0; e < nedges; e++ )
            {
               if( edgearrint[e] == CONNECT)
               {
                  solnode[prunegraph->head[e]] = CONNECT;
                  solnode[prunegraph->tail[e]] = CONNECT;
               }
            }

            do
            {
               count = 0;
               for( int k = nnodes - 1; k >= 0; --k )
               {
                  int e;
                  if( (solnode[k] != CONNECT) || Is_term(prunegraph->term[k]) )
                     continue;

                  for( e = prunegraph->outbeg[k]; e != EAT_LAST; e = prunegraph->oeat[e] )
                     if( edgearrint[e] == CONNECT )
                        break;

                  if( e == EAT_LAST )
                  {
                     /* there has to be exactly one incoming edge */
                     for( e = prunegraph->inpbeg[k]; e != EAT_LAST; e = prunegraph->ieat[e] )
                     {
                        if( edgearrint[e] == CONNECT )
                        {
                           edgearrint[e] = UNKNOWN;
                           solnode[k] = UNKNOWN;
                           count++;
                           break;
                        }
                     }
                     assert(e != EAT_LAST);
                  }
               }
            }
            while( count > 0 );

            for( int k = 0; k < nnodes; k++ )
               pmark[k] = (prunegraph->grad[k] > 0);
         }

#if BREAKONERROR
         if( !graph_sol_valid(scip, prunegraph, edgearrint) )
         {
            printf("sol 2 not valid %d \n", 0);
            return SCIP_ERROR;
         }
#endif
#ifdef SCIP_DEBUG
         if( !graph_valid(prunegraph) )
            printf("reduced graph in prune not valid \n");
#endif

         if( pcmw )
            SCIP_CALL( graph_pc_2org(scip, prunegraph) );

         SCIP_CALL( bound_reducePrune(scip, prunegraph, vnoi, cost, (pcmw) ? prunegraph->prize : NULL, nodearrreal, costrev,
               &offsetnew, heap, state, vbase, solnode, edgearrint, &brednelims, minnelims));

         if( pcmw )
            SCIP_CALL( graph_pc_2trans(scip, prunegraph) );

#ifdef SCIP_DEBUG
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
         printf("PRUNE round: %d edges %d  nodes: %d \n", i, anedges / 2, annodes);
         printf("PRUNE round: %d minelims %d  really reduced: %d \n", i, minnelims, brednelims);
#endif

         /* not enough reductions possible? */
         if( brednelims < lminnelims  )
         {
#ifdef SCIP_DEBUG
            printf("too little elims in PRUNE; out !! %d \n\n", brednelims);
#endif
            i = PRUNE_MAXREDROUNDS;
         }

         /* delete all vertices not reachable from the root */
         SCIP_CALL( level0(scip, prunegraph) );

         /* is the reduced graph still feasible? */
         if( !graph_valid(prunegraph) )
         {
            solvalid = FALSE;
#if BREAKONERROR
            printf("reduced graph in prune not valid X \n \n \n");
            return SCIP_ERROR;
#endif
            break;
         }

         /* get number of remaining edges */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         reductbound = getRedBound(i + 1, anedges);

         /* reduce graph */
         if( pc )
         {
            SCIP_CALL( redLoopPc(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
                  vbase, nodearrint, edgearrint, nodearrint2, solnode, nodearrchar, &offsetnew, FALSE, FALSE, reductbound, FALSE) );
         }
         else if( mw )
         {
            SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, state,
                  vbase, nodearrint, edgearrint, nodearrint2, heap, solnode, nodearrchar, &offsetnew, FALSE, FALSE, FALSE, reductbound, FALSE) );
         }
         else
         {
            SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state, vbase, nodearrint, edgearrint,
                  nodearrint2, solnode, nodearrchar, &offsetnew, -1.0, FALSE, FALSE, TRUE, reductbound, NULL, FALSE));
         }


         /* delete all vertices not reachable from the root */
         SCIP_CALL( level0(scip, prunegraph) );

         /* get number of remaining terminals */
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

         if( anterms <= 2 )
         {
            SCIPdebugMessage("less than two terminals, break !! \n");
            break;
         }

      } /* reduction loop */

      if( solvalid )
      {
         SCIP_CALL( level0(scip, prunegraph) );
      }
   }

   assert(solvalid);


   for( int k = 0; k < nnodes; k++ )
      nodearrchar[k] = FALSE;

   graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
   SCIPdebugMessage("Xin prune grad: %d , nedges: %d nodes: %d \n", prunegraph->grad[prunegraph->source], anedges, annodes);

   /* if graph not vanished, compute solution */
   if( prunegraph->grad[prunegraph->source] > 0 )
   {
      GRAPH* pprunegraph;
      SCIP_Real objsph;
      SCIP_Real objprune = 0.0;
      int npruneedges;

      ancestors = prunegraph->ancestors;

      assert(soledge != NULL);

      /* try to build MST on solnode nodes */

      if( pcmw )
      {
         SCIP_CALL( SCIPStpHeurTMBuildTreePcMw(scip, prunegraph, path, prunegraph->cost, &objprune, solnode) );
      }
      else
      {
         SCIP_CALL( SCIPStpHeurTMBuildTree(scip, prunegraph, path, prunegraph->cost, &objprune, solnode) );
      }

      /* solution valid? */
      if( SCIPisLT(scip, objprune, FARAWAY) )
      {
         SCIPdebugMessage("solution in prune finally still valid \n");
         for( int k = 0; k < nnodes; k++ )
            if( path[k].edge >= 0 )
               graph_sol_setNodeList(g, nodearrchar, ancestors[path[k].edge]);
      }

      graph_path_exit(scip, prunegraph);

      /* pack the graph */
      SCIP_CALL( graph_pack(scip, prunegraph, &pprunegraph, FALSE) );

      prunegraph = pprunegraph;

      assert(graph_valid(prunegraph));

      ancestors = prunegraph->ancestors;
      npruneedges = prunegraph->edges;

      for( int e = 0; e < npruneedges; e++ )
         soledge[e] = UNKNOWN;

      /* initialize shortest path algorithm */
      SCIP_CALL( graph_path_init(scip, prunegraph) );

      /* compute new solution on heuristically reduced graph */

      /* run TM heuristic */
      SCIP_CALL( SCIPStpHeurTMRun(scip, NULL, prunegraph, NULL, &best_start, soledge, DEFAULT_PRUNE_TMRUNS,
            prunegraph->source, prunegraph->cost, prunegraph->cost, &hopfactor, NULL, 0.0, success, FALSE) );

      if( pcmw )
      {
         SCIP_CALL( SCIPStpHeurLocalRun(scip, prunegraph, cost, soledge) );
      }

#if BREAKONERROR
      if( !graph_sol_valid(scip, prunegraph, soledge))
      {
         printf("sol 3 not valid %d \n", 0);
         return SCIP_ERROR;
      }
#endif

      /* free shortest path memory */
      graph_path_exit(scip, prunegraph);

      if( !(*success) || !graph_sol_valid(scip, prunegraph, soledge) )
      {
         SCIPdebugMessage("failed to build reduced tree in prune \n");
         goto TERMINATE;
      }

      /* retransform solution found by TM heuristic */

      objsph = graph_sol_getObj(prunegraph->cost, soledge, 0.0, npruneedges);
#ifdef SCIP_DEBUG
      printf("prune: sph weight %f \n", objsph);
      printf("prune: recovered weight %f \n", objprune);
#endif

      if( SCIPisLT(scip, objsph, objprune) )
      {
         for( int k = 0; k < nnodes; k++ )
            nodearrchar[k] = FALSE;

         for( int e = 0; e < npruneedges; e++ )
            if( soledge[e] == CONNECT )
               graph_sol_setNodeList(g, nodearrchar, ancestors[e]);
      }
   }
   else
   {
      graph_path_exit(scip, prunegraph);
   }

   /* retransform edges fixed during graph reduction */
   graph_sol_setNodeList(g, nodearrchar, prunegraph->fixedges);

   if( pcmw )
   {
      IDX* curr;
      nodearrchar[prunegraph->source] = TRUE;

      for( int k = 0; k < nnodes; k++ )
      {
         if( nodearrchar[k] )
         {
            curr = prunegraph->pcancestors[k];
            while( curr != NULL )
            {
               assert(prunegraph->orgtail[curr->index] < nnodes);
               assert(prunegraph->orghead[curr->index] < nnodes);

               nodearrchar[prunegraph->orgtail[curr->index]] = TRUE;
               nodearrchar[prunegraph->orghead[curr->index]] = TRUE;

               curr = curr->parent;
            }
         }
      }
   }

   /* prune solution (in the original graph) */

   for( int e = 0; e < nedges; e++ )
      soledge[e] = UNKNOWN;

   if( pcmw )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, g, g->cost, soledge, nodearrchar) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, g, g->cost, 0, soledge, nodearrchar) );

   *success = graph_sol_valid(scip, g, soledge);
   assert(*success);

#ifdef SCIP_DEBUG
   SCIP_Real x = 0.0;
   int s = 0;
   for( e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
      {
         s++;
         x += g->cost[e];
      }
   printf("x %f s: %d \n", x, s);
#endif


#if BREAKONERROR
   if( !(*success) )
   {
      printf("final sol not valid %d \n", 0);
      return SCIP_ERROR;
   }
#endif
 TERMINATE:

   /* free memory */
   graph_free(scip, prunegraph, TRUE);

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
   SCIPfreeBufferArray(scip, &bestsoledges);

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


   /* add prune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_PRUNE_MAXFRQ, NULL, NULL) );

   return SCIP_OKAY;
}
