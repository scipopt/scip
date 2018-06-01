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

/**@file   heur_slackprune.c
 * @brief  dual-ascent and reduction based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a dual-ascent and reduction based heuristic for Steiner problems. It is based on an approach
 * described in T. Polzin's "Algorithms for the Steiner problem in networks".
 *
 * A list of all interface methods can be found in heur_slackprune.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_slackprune.h"
#include "heur_ascendprune.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "grph.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"
#include "prop_stp.h"

#define HEUR_NAME             "slackprune"
#define HEUR_DESC             "Reduction based heuristic for Steiner problems"
#define HEUR_DISPCHAR         'S'
#define HEUR_PRIORITY         1
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_SLACKPRUNE_MAXFREQ   FALSE       /**< executions of the heuristic at maximum frequency?                             */
#define SLACKPRUNE_MINREDELIMS       2           /**< minimum number of eliminations for reduction package when called by slack-and-prune heuristic */
#define SLACKPRUNE_MAXREDROUNDS      10          /**< maximum number of reduction rounds in slack-prune heuristic */
#define SLACKPRUNE_MINSTALLPROPORTION   0.25      /**< minimum proportion of arcs to be fixed before restarting slack-prune heuristic */
#define SLACKPRUNE_MAXSTALLPROPORTION   0.5       /**< maximum proportion of arcs to be fixed before restarting slack-prune heuristic */
#define BREAKONERROR FALSE
#define MAXNTERMINALS 500
#define MAXNEDGES     10000
#define SLACK_MAXTOTNEDGES 5000

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
   int                   lastnfixededges;    /**< number of fixed edges before the previous run                     */
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   nfailures;          /**< number of failures since last successful call                     */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};

/*
 * Local methods
 */

/* get reduction bound */
static
int getRedBound(
   int nrounds,
   int nedges
   )
{
   if( nrounds == 0)
      return (MAX(nedges / 2000, SLACKPRUNE_MINREDELIMS));
   if( nrounds == 1)
      return (MAX(nedges / 1000, SLACKPRUNE_MINREDELIMS));
   return(MAX(nedges / 500, SLACKPRUNE_MINREDELIMS));
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
   int                   maxnrounds
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
   assert(maxnrounds >= 1);

   anedges = anedges / 2;

   totminnelims = MAX(SLACKPRUNE_MINREDELIMS, (anedges / 25));

   min = (int) (anedges * 0.15);
   min -= (int) (((double) min * anterms) / (annodes));
   min = MAX(min, 1);

   factor = (double) anedges / min;
   factor = ((double) nround / (2.5 * maxnrounds)) * factor;
#if 1
   if( SCIPisGT(scip, factor, 1.0) )
   {
      SCIP_Real tmp = min * factor;
      min = (int) tmp;
   }
#endif
   min = MAX(totminnelims, min);

   min = MIN(min, (anedges - 1));
   min = MAX(min, 1);

   *lminnelims = min / 2;
   *lminnelims = MAX(*lminnelims, 1);

   *minnelims = min;

}

static
void updateSolNodeArray(
   GRAPH*                graph,
   int*                  soledge,
   int*                  solnode,
   int                   nnodes,
   int                   nnedges
   )
{
   int j;
   int e;

   for( j = 0; j < nnodes; j++ )
      solnode[j] = UNKNOWN;

   /* reset vertices that are to be kept */
   for( e = 0; e < nnedges; e++ )
   {
      if( soledge[e] == CONNECT )
      {
         solnode[graph->tail[e]] = CONNECT;
         solnode[graph->head[e]] = CONNECT;
      }
   }
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
   const char           filename,           /**< Name of the output file */
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
SCIP_DECL_HEURCOPY(heurCopySlackPrune)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurSlackPrune(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSlackPrune)
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
SCIP_DECL_HEURINIT(heurInitSlackPrune)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolSlackPrune)
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
   heurdata->lastnfixededges = -1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSlackPrune)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSlackPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;
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
   if( graph->stp_type != STP_SPG && graph->stp_type != STP_RSMT && graph->stp_type != STP_OARSMT && graph->stp_type != STP_GSTP )
      return SCIP_OKAY;

   if( (graph->edges > MAXNEDGES) && (graph->terms > MAXNTERMINALS) )
      return SCIP_OKAY;

   if( !(heurdata->maxfreq) && heurdata->nfailures > 0 )
      return SCIP_OKAY;

   /* get best current solution */
   bestsol = SCIPgetBestSol(scip);

   /* no solution available? */
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* heuristic not at maximum or ...*/
   if( !(heurdata->maxfreq)
      /* has the new solution been found by this very heuristic or is new best solution available? */
      || (SCIPsolGetHeur(bestsol) == heur || heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip))) )
   {
      if( heurdata->lastnfixededges >= 0 )
      {
         SCIP_Real stallproportion;

         stallproportion = (1.0 + heurdata->nfailures) * SLACKPRUNE_MINSTALLPROPORTION;

         if( SCIPisGT(scip, stallproportion, SLACKPRUNE_MAXSTALLPROPORTION) )
            stallproportion = SLACKPRUNE_MAXSTALLPROPORTION;

         if( (SCIPStpNfixedEdges(scip) - heurdata->lastnfixededges) < (int) (stallproportion * nedges) )
            return SCIP_OKAY;
      }
   }

   /* deactivate as expensive is too expensive to be reiterated todo: do this properly */
   heurdata->nfailures = 1;

   xval = SCIPprobdataGetXval(scip, bestsol);

   if( xval == NULL )
      return SCIP_OKAY;

   heurdata->lastnfixededges = SCIPStpNfixedEdges(scip);

   vars = SCIPprobdataGetVars(scip);
   nvars = SCIPprobdataGetNVars(scip);

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

   /* execute slackprune heuristic */
   SCIP_CALL( SCIPStpHeurSlackPruneRun(scip, vars, graph, soledge, &success, FALSE, (graph->edges < SLACK_MAXTOTNEDGES)) );

   /* solution found by slackprune heuristic? */
   if( success )
   {
      SCIP_SOL* sol;
      SCIP_Real pobj;
      SCIP_Real* nval;

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

      SCIPdebugMessage("SP final solution: best: old %f, new %f \n",  SCIPgetSolOrigObj(scip, bestsol), pobj + SCIPprobdataGetOffset(scip));

      /* try to add new solution to pool */
      sol = NULL;
      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

      /* has solution been added? */
      if( success )
      {
         SCIPdebugMessage("better solution added by SLACKPRUNE %f \n", pobj + SCIPprobdataGetOffset(scip));
         *result = SCIP_FOUNDSOL;

         assert(graph_sol_valid(scip, graph, soledge));

         /* is the solution the new incumbent? */
         if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj) )
         {
            /* deactivated subsequent runs since heuristics seems to be to expensive */
            heurdata->nfailures = 1;
         }
         else
         {
            success = FALSE;
         }
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &nval);
   }
   else
   {
      SCIPdebugMessage("slack-prune: solution not valid \n");
   }

   /* solution could not be found, added or is not best? */
   if( !success )
      heurdata->nfailures++;

   /* store index of incumbent solution */
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* free memory */
   SCIPfreeBufferArray(scip, &soledge);
   SCIPdebugMessage("leaving slack and prune heuristic \n");

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** execute slack-and-prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurSlackPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   GRAPH*                g,                  /**< graph data structure */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   SCIP_Bool*            success,            /**< feasible solution found? */
   SCIP_Bool             reducegraph,        /**< try to reduce graph initially? */
   SCIP_Bool             fullreduce          /**< use full reduction techniques? */
   )
{
   GRAPH*    prunegraph;
   PATH*    vnoi;
   PATH*    path;
   GNODE** gnodearr;
   SCIP_Real    ubnew;
   SCIP_Real    ubbest;
   SCIP_Real    offsetnew;
   SCIP_Real    globalobj;
   SCIP_Real*    cost;
   SCIP_Real*    costrev;
   SCIP_Real*    nodearrreal;
   int     i;
   int     k;
   int     e;
   int     nterms;
   int     nnodes;
   int     nedges;
   int     anterms;
   int     anedges;
   int     annodes;
   int     probtype;
   int     minnelims;
   int     lminnelims;
   int     reductbound;
   int*     heap;
   int*     state;
   int*     vbase;
   int*     solnode;
   int*     edgearrint2;
   int*     nodearrint;
   int*     nodearrint2;
   int*     globalsoledge;
   STP_Bool*     nodearrchar;
   STP_Bool*     edgearrchar;

   assert(g != NULL);
   assert(scip != NULL);
   assert(soledge != NULL);
   assert(graph_sol_valid(scip, g, soledge));

   nterms = g->terms;
   nedges = g->edges;
   nnodes = g->knots;
   probtype = g->stp_type;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   /* mark solution vertices */
   for( k = 0; k < nnodes; k++ )
      solnode[k] = UNKNOWN;

   ubbest = 0.0;

   /* set solution array and get solution value */
   for( e = 0; e < nedges; e++ )
   {
      if( soledge[e] == CONNECT )
      {
         ubbest += g->cost[e];
         solnode[g->tail[e]] = CONNECT;
         solnode[g->head[e]] = CONNECT;
      }
   }

   globalobj = ubbest;

   /* set offset (within new graph) to 0.0 */
   offsetnew = 0.0;

   /* allocate memory for reduction methods */
   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &(gnodearr[i])) ); /*lint !e866*/
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &globalsoledge, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrchar, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint2, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   BMScopyMemoryArray(globalsoledge, soledge, nedges);


   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      g->stp_type = STP_SPG;

   /* copy the graph */
   SCIP_CALL( graph_copy(scip, g, &prunegraph) );

   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      g->stp_type = probtype;

   /* set ancestors of the new graph */
   SCIP_CALL( graph_init_history(scip, prunegraph) );

   reductbound = getRedBound(0, nedges);

   /* variables given? */
   if( vars != NULL )
   {
      int nfixedges = 0;

      /* delete fixed edges from the new graph */
      for( e = 0; e < nedges; e += 2 )
      {
         /* both e and its anti-parallel edge fixed to zero? */
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5
            && soledge[e] != CONNECT && soledge[e + 1] != CONNECT )
         {
            graph_edge_del(scip, prunegraph, e, TRUE);
            nfixedges++;
         }
      }
      SCIPdebugMessage("fixed edges in slack and prune: %d \n", nfixedges);

      if( nfixedges > reductbound && reducegraph )
      {
         graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
         reductbound = getRedBound(0, anedges);
      }
   }

   SCIP_CALL( graph_path_init(scip, prunegraph) );

   /* perform initial reductions? */
   if( reducegraph )
   {
      SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
            vbase, nodearrint, soledge, nodearrint2, solnode, nodearrchar, &offsetnew, -1.0, TRUE, FALSE, TRUE, reductbound, TRUE, fullreduce) );
   }

   /* get number of remaining vertices, edges and terminals */
   graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

   /* main reduction loop */
   for( i = 0; i < SLACKPRUNE_MAXREDROUNDS && anterms > 2; i++ )
   {
      SCIP_Real obj;
      SCIP_Bool apsuccess;
      int danelims;

      /* update reduction bounds */
      setMinMaxElims(scip, &minnelims, &lminnelims, annodes, anedges, anterms, i + 1, SLACKPRUNE_MAXREDROUNDS);
#if BREAKONERROR
      if( minnelims > (anedges / 2) )
         return SCIP_ERROR;
#endif

      /*  perform reductions */
      SCIP_CALL( reduce_daSlackPrune(scip, vars, prunegraph, vnoi, gnodearr, cost, costrev, nodearrreal, &ubnew,
            soledge, edgearrint2, vbase, nodearrint, state, solnode, nodearrchar, edgearrchar, &danelims, minnelims, ((i == 0) && !reducegraph)) );

      /* delete all vertices not reachable from the root */
      SCIP_CALL( level0(scip, prunegraph) );

      assert(graph_valid(prunegraph));

      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

      if( anterms <= 2 )
         break;

      SCIPdebugMessage("minelims %d really reduced: %d \n", minnelims, danelims);

      /* not enough reductions possible? */
      if( danelims < lminnelims )
      {
         SCIPdebugMessage("too little elims in DA, break %d < %d\n\n", danelims, lminnelims);
         i = SLACKPRUNE_MAXREDROUNDS;
      }

      /* compute potential new guiding solution */
      SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, prunegraph, cost, soledge, nodearrint, prunegraph->source, nodearrchar, &apsuccess, FALSE) );

      /* solution found by ascend and prune? */
      if( apsuccess )
      {
         SCIP_CALL( SCIPStpHeurLocalRun(scip, prunegraph, prunegraph->cost, soledge) );

         assert(graph_sol_valid(scip, prunegraph, soledge));

         SCIP_CALL( SCIPStpHeurPruneUpdateSols(scip, g, prunegraph, path, nodearrint, edgearrint2, solnode, soledge,
               globalsoledge, nodearrchar, &globalobj, TRUE, success) );

         /* calculate objective value of solution */
         obj = graph_sol_getObj(prunegraph->cost, soledge, offsetnew, nedges);

         /* obj <= incumbent objective value? */
         if( SCIPisLE(scip, obj, ubbest) )
            ubbest = obj;

         SCIPdebugMessage("old solution: %f new solution %f \n\n", ubbest + SCIPprobdataGetOffset(scip), obj + SCIPprobdataGetOffset(scip));
      }

      /* get number of remaining edges */
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

      reductbound = getRedBound(i, anedges);

      /* reduce graph, using the new upper bound and not letting BND eliminate solution edges */
      SCIP_CALL( redLoopStp(scip, prunegraph, vnoi, path, NULL, nodearrreal, cost, costrev, heap, state,
            vbase, nodearrint, soledge, nodearrint2, solnode, nodearrchar, &offsetnew, -1.0, TRUE, FALSE, TRUE, reductbound, TRUE, fullreduce) );

      /* graph vanished? */
      if( prunegraph->grad[prunegraph->source] == 0 )
         break;

      /* get number of remaining edges */
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

   } /* reduction loop */

   SCIP_CALL(SCIPStpHeurLocalRun(scip, g, g->cost, globalsoledge));

   graph_path_exit(scip, prunegraph);

   *success = graph_sol_valid(scip, g, globalsoledge);
   BMScopyMemoryArray(soledge, globalsoledge, nedges);

#if BREAKONERROR
   if( !(*success) )
   {
      printf("slack prune solution not valid %d \n", 0);
      return SCIP_ERROR;
   }
#endif

   /* free memory */
   graph_free(scip, &prunegraph, TRUE);

   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &edgearrint2);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &edgearrchar);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &globalsoledge);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBlockMemory(scip, &(gnodearr[i]));
   SCIPfreeBufferArray(scip, &gnodearr);

   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &solnode);

   return SCIP_OKAY;
}


/** execute slack-and-prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurSlackPruneRunPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   GRAPH*                g,                  /**< graph data structure */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   SCIP_Bool*            success             /**< feasible solution found? */
   )
{
   GRAPH*    prunegraph;
   PATH*    vnoi;
   PATH*    path;
   GNODE** gnodearr;
   SCIP_Real    ubbest;
   SCIP_Real    offsetnew;
   SCIP_Real*    cost;
   SCIP_Real*    costrev;
   SCIP_Real*    nodearrreal;
   int     i;
   int     k;
   int     e;
   int     nterms;
   int     nnodes;
   int     nedges;
   int     anterms;
   int     anedges;
   int     annodes;
   int     minnelims;
   int     lminnelims;
   int     reductbound;
   int     nprunenodes;
   int     npruneedges;
   int*     state;
   int*     vbase;
   int*     solnode;
   int*     edgearrint;
   int*     nodearrint;
   int*     nodearrint2;
   int*     nodearrint3;
   STP_Bool*     nodearrchar;

   assert(g != NULL);
   assert(scip != NULL);
   assert(soledge != NULL);
   assert(graph_sol_valid(scip, g, soledge));

   nterms = g->terms;
   nedges = g->edges;
   nnodes = g->knots;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   /* mark solution vertices */
   for( k = 0; k < nnodes; k++ )
      solnode[k] = UNKNOWN;

   ubbest = 0.0;

   /* set solution array and get solution value (with respect to current, possibly reduced, graph) */
   for( e = 0; e < nedges; e++ )
   {
      if( soledge[e] == CONNECT )
      {
         ubbest += g->cost[e];
         solnode[g->tail[e]] = CONNECT;
         solnode[g->head[e]] = CONNECT;
      }
   }

   /* set offset (within new graph) to 0.0 */
   offsetnew = 0.0;

   /* allocate memory for reduction methods */
   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint3, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   /* copy the graph */
   SCIP_CALL( graph_copy(scip, g, &prunegraph) );

   /* set ancestors of the new graph */
   SCIP_CALL( graph_init_history(scip, prunegraph) );

   edgearrint = soledge;

   /* variables given? */
   if( vars != NULL )
   {
      /* delete fixed edges from the new graph */
      for( e = 0; e < nedges; e += 2 )
      {
         /* both e and its anti-parallel edge fixed to zero? */
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5
            && soledge[e] != CONNECT && soledge[e + 1] != CONNECT && !Is_term(g->term[g->tail[e]]) && !Is_term(g->term[g->head[e]]) )
         {
            graph_edge_del(scip, prunegraph, e, TRUE);
         }
      }
   }

   SCIP_CALL( graph_path_init(scip, prunegraph) );

   npruneedges = prunegraph->edges;
   nprunenodes = prunegraph->knots;
   prunegraph->norgmodeledges = npruneedges;
   prunegraph->norgmodelknots = nprunenodes;

   /* get number of remaining vertices, edges and terminals */
   graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

   /* main reduction loop */
   for( i = 0; i < SLACKPRUNE_MAXREDROUNDS && anterms > 3; i++ )
   {
      SCIP_Real obj;
      SCIP_Bool apsuccess;
      int danelims;

      /* update reduction bounds */
      setMinMaxElims(scip, &minnelims, &lminnelims, annodes, anedges, anterms, i + 1, SLACKPRUNE_MAXREDROUNDS);

#if BREAKONERROR
      if( minnelims > (anedges / 2) )
      {
         printf("too many elim %d \n", minnelims);
         return SCIP_ERROR;
      }
#endif

      /*  perform heuristic reductions */
      SCIP_CALL( reduce_daSlackPruneMw(scip, prunegraph, vnoi, gnodearr, cost, costrev, nodearrreal, vbase, nodearrint,
            edgearrint, nodearrint2, solnode, nodearrchar, &danelims, minnelims, ((i == 0))) );

      updateSolNodeArray(prunegraph, edgearrint, solnode, nprunenodes, npruneedges);

      /* calculate objective value of solution */
      obj = graph_sol_getObj(prunegraph->cost, edgearrint, offsetnew, npruneedges);

      ubbest = obj;

      /* delete all vertices not reachable from the root */
      SCIP_CALL( level0(scip, prunegraph) );

      assert(graph_valid(prunegraph));

      /* get number of remaining vertices, edges and terminals */
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

      /* update reduction bounds */
      setMinMaxElims(scip, &minnelims, &lminnelims, annodes, anedges, anterms, i + 1, SLACKPRUNE_MAXREDROUNDS);

      if( anterms <= 3 )
      {
         SCIPdebugMessage("graph vanished in SLACKPRUNE \n\n");
         break;
      }

      SCIPdebugMessage("SP: minelimsx %d really reduced: %d \n", minnelims, danelims);

      /* not enough reductions possible? */
      if( danelims < lminnelims )
      {
         SCIPdebugMessage("too little elims in DA OUT !! %d < %d\n\n", danelims, lminnelims);
         i = SLACKPRUNE_MAXREDROUNDS;
      }

      /* compute new guiding solution */
      SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, prunegraph, cost, edgearrint, vbase, -1, nodearrchar, &apsuccess, FALSE) );

      /* solution found by ascend and prune? */
      if( apsuccess )
      {
         /* calculate objective value of solution */
         obj = graph_sol_getObj(prunegraph->cost, edgearrint, offsetnew, npruneedges);

         SCIPdebugMessage(" old solution: %f AP solution %f \n", ubbest + SCIPprobdataGetOffset(scip), obj + SCIPprobdataGetOffset(scip));
         SCIPdebugMessage("offsetnew %f \n", offsetnew);

         /* obj <= incumbent objective value? */
         if( SCIPisLT(scip, obj, ubbest) )
            updateSolNodeArray(prunegraph, edgearrint, solnode, nprunenodes, npruneedges);
      }

      /* get number of remaining edges */
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);

      reductbound = getRedBound(i, anedges);

      /* reduction loop */
      SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, path, NULL, NULL, NULL, NULL, state,
            vbase, nodearrint, NULL, nodearrint2, nodearrint3, solnode, nodearrchar, &offsetnew, FALSE, FALSE, FALSE, reductbound, TRUE) );

      assert(graph_valid(prunegraph));

#if BREAKONERROR
      if( !graph_valid(prunegraph) )
      {
         printf("SP: not valid2 %d \n", 0);
         return SCIP_ERROR;
      }
#endif

      /* get number of remaining edges */
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
   } /* reduction loop */

   /* if graph not vanished, compute solution */
   if( prunegraph->grad[prunegraph->source] > 0 )
   {
      IDX** ancestors = prunegraph->ancestors;
      SCIP_Real objorg;
      SCIP_Real objprune;

      /* build solution on solnode nodes */

      for( e = 0; e < nedges; e++ )
         soledge[e] = UNKNOWN;

      for( k = 0; k < nnodes; k++ )
         nodearrchar[k] = (solnode[k] == CONNECT);

      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, prunegraph, prunegraph->cost, soledge, nodearrchar) );

      for( k = 0; k < nnodes; k++ )
         nodearrchar[k] = FALSE;

      for( e = 0; e < nedges; e++ )
         if( soledge[e] == CONNECT )
            graph_sol_setNodeList(g, nodearrchar, ancestors[e]);

         /* calculate objective value of solution */
         objorg = graph_sol_getObj(prunegraph->cost, soledge, offsetnew, nedges);

      /* compute new solution on heuristically reduced graph */

#ifdef SCIP_DEBUG
      graph_get_NVET(prunegraph, &annodes, &anedges, &anterms);
      printf("final SP anedges: %d anterms: %d \n", anedges, anterms);
#endif

      SCIP_CALL( SCIPStpHeurPruneRun(scip, NULL, prunegraph, soledge, success, FALSE, TRUE) );

      if( !(*success) )
      {
#if BREAKONERROR
         printf("Xfailed to build tree\n");
         return SCIP_ERROR;
#endif
         goto TERMINATE;
      }

      /* calculate objective value of solution */
      objprune = graph_sol_getObj(prunegraph->cost, soledge, offsetnew, nedges);

      if( SCIPisLT(scip, objprune, objorg) )
      {
         /* mark vertices of solution found by prune heuristic */

         for( k = 0; k < nnodes; k++ )
            nodearrchar[k] = FALSE;

         for( e = 0; e < nedges; e++ )
            if( soledge[e] == CONNECT )
               graph_sol_setNodeList(g, nodearrchar, ancestors[e]);

      }
   }
   else
   {
      for( k = 0; k < nnodes; k++ )
         nodearrchar[k] = FALSE;
   }

   /* retransform edges fixed during graph reduction */
   graph_sol_setNodeList(g, nodearrchar, prunegraph->fixedges);

   SCIP_CALL( graph_sol_markPcancestors(scip, prunegraph->pcancestors, prunegraph->orgtail, prunegraph->orghead, nnodes,
         nodearrchar, NULL, NULL, NULL, NULL) );

   for( e = 0; e < nedges; e++ )
      soledge[e] = UNKNOWN;

   SCIP_CALL( SCIPStpHeurTMPrunePc(scip, g, g->cost, soledge, nodearrchar) );
   *success = graph_sol_valid(scip, g, soledge);

 TERMINATE:

   /* free memory */
   graph_path_exit(scip, prunegraph);
   graph_free(scip, &prunegraph, TRUE);

   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint3);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &nodearrreal);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBuffer(scip, &gnodearr[i]);
   SCIPfreeBufferArray(scip, &gnodearr);

   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &solnode);

   return SCIP_OKAY;
}


/** creates the slackprune primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurSlackPrune(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create slackprune primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSlackPrune, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySlackPrune) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSlackPrune) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSlackPrune) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolSlackPrune) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSlackPrune) );

   /* add slackprune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_SLACKPRUNE_MAXFREQ, NULL, NULL) );


   return SCIP_OKAY;
}
