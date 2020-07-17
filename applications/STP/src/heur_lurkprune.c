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

/**@file   heur_lurkprune.c
 * @brief  lurking-bounds reduction based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reduction based heuristic for Steiner problems that makes use of lurking bounds.
 *
 * A list of all interface methods can be found in heur_lurkprune.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_lurkprune.h"
#include "heur_ascendprune.h"
#include "dualascent.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "solstp.h"
#include "cons_stp.h"
#include "probdata_stp.h"
#include "prop_stp.h"

#define HEUR_NAME             "lurkprune"
#define HEUR_DESC             "Reduction based heuristic for Steiner problems"
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         1
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_LURKPRUNE_MAXFREQ   FALSE       /**< executions of the heuristic at maximum frequency?                             */
#define LURKPRUNE_NTMRUNS           30          /**< TM heur runs */
#define LURKPRUNE_MINREDELIMS       10           /**< minimum number of eliminations for reduction package when called by lurk-and-prune heuristic */
#define LURKPRUNE_MAXREDROUNDS      10           /**< maximum number of reduction rounds in lurk-prune heuristic */
#define LURKPRUNE_MINLURKEDGE_RATIO 0.05
#define LURKPRUNE_MINSTALLPROPORTION   0.25      /**< minimum proportion of arcs to be fixed before restarting lurk-prune heuristic */
#define LURKPRUNE_MAXSTALLPROPORTION   0.5       /**< maximum proportion of arcs to be fixed before restarting lurk-prune heuristic */
#define LURK_MAXTOTNEDGES 5000

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastnfixededges;    /**< number of fixed edges before the previous run                     */
   int                   nfailures;          /**< number of failures since last successful call                     */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};


/** data */
typedef struct lurking_prune
{
   int*                  solnode;
   int*                  soledge;
   int*                  globalsoledge;
   SCIP_Real*            lurkingbounds_half;
   SCIP_Real             ubnew;
   SCIP_Real             ubbest;
   SCIP_Real             offsetnew;
   SCIP_Real             globalobj;
   SCIP_Real             obj_old;
   int                   minlurkelims;
} LURKPRUNE;



/*
 * Local methods
 */


/** initializes prune graph */
static
SCIP_RETCODE prunegraphInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph */
   GRAPH**               prunegraph          /**< graph */
)
{
   SCIP_CALL( graph_copy(scip, g, prunegraph) );
   SCIP_CALL( graph_init_history(scip, *prunegraph) );
   if( graph_pc_isPcMw(g) )
      (*prunegraph)->norgmodelknots = (*prunegraph)->knots;
   SCIP_CALL( graph_path_init(scip, *prunegraph) );
   graph_mark(*prunegraph);

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE lurkpruneInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      lurkingbounds,      /**< lurking edge bounds */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   LURKPRUNE*            lurkprune,          /**< lurking-prune data */
   GRAPH**               prunegraph          /**< graph */
   )
{
   SCIP_Real* lurkingbounds_half;
   int* solnode;
   int* globalsoledge;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   int nedges_red;

   SCIP_CALL( prunegraphInit(scip, g, prunegraph) );

   SCIP_CALL( SCIPallocBufferArray(scip, &globalsoledge, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lurkingbounds_half, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );

   BMScopyMemoryArray(globalsoledge, soledge, nedges);

   for( int k = 0; k < nnodes; k++ )
      solnode[k] = UNKNOWN;

   for( int e = 0; e < nedges; e++ )
   {
      if( soledge[e] == UNKNOWN )
         continue;

      assert(CONNECT == soledge[e]);

      solnode[g->tail[e]] = CONNECT;
      solnode[g->head[e]] = CONNECT;
   }

   for( int e = 0; e < nedges; e += 2 )
   {
      const SCIP_Real min = MIN(lurkingbounds[e], lurkingbounds[e + 1]);
      assert(GE(min, 0.0) || EQ(min, -FARAWAY));

      lurkingbounds_half[e / 2] = min;
   }

   graph_get_nVET(g, NULL, &nedges_red, NULL);

   lurkprune->minlurkelims = (nedges_red / 2) * LURKPRUNE_MINLURKEDGE_RATIO;
   lurkprune->ubbest = solstp_getObjBounded(g, soledge, 0.0, nedges);
   lurkprune->globalobj = lurkprune->ubbest;
   lurkprune->offsetnew = 0.0;

   lurkprune->globalsoledge = globalsoledge;
   lurkprune->lurkingbounds_half = lurkingbounds_half;
   lurkprune->solnode = solnode;
   lurkprune->soledge = soledge;
   lurkprune->obj_old = solstp_getObj(g, soledge, 0.0);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("lurking-prune starts \n");
   graph_printInfoReduced(g);
#endif

   return SCIP_OKAY;
}


/** finalizes */
static
void lurkpruneFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph */
   GRAPH*                prunegraph,         /**< graph */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   LURKPRUNE*            lurkprune,          /**< lurking-prune data */
   SCIP_Bool*            solimproved         /**< could a better solution be found? */
   )
{
   const int nedges = graph_get_nEdges(prunegraph);
   const SCIP_Real obj_new = solstp_getObj(g, lurkprune->globalsoledge, 0.0);
   const SCIP_Real obj_old = lurkprune->obj_old;

   SCIPdebugMessage("lurk-prune: obj_old=%f obj_new=%f \n", obj_old, obj_new);
   *solimproved = (LT(obj_new, obj_old));

   graph_path_exit(scip, prunegraph);
   BMScopyMemoryArray(soledge, lurkprune->globalsoledge, nedges);

   graph_free(scip, &prunegraph, TRUE);
   SCIPfreeBufferArray(scip, &(lurkprune->solnode));
   SCIPfreeBufferArray(scip, &(lurkprune->lurkingbounds_half));
   SCIPfreeBufferArray(scip, &(lurkprune->globalsoledge));
}


/** does exact reductions */
static
SCIP_RETCODE reduceExact(
   SCIP*                 scip,               /**< SCIP data structure */
   LURKPRUNE*            lurkprune,          /**< lurking-prune data */
   GRAPH*                prunegraph          /**< graph data structure */
   )
{
   PATH* vnoi;
   PATH* path;
   SCIP_Real* nodearrreal;
   int* heap;
   int* state;
   int* vbase;
   int* nodearrint;
   int* edgearrint;
   int* nodearrint2;
   STP_Bool* nodearrchar;
   const int nnodes = graph_get_nNodes(prunegraph);
   const int nedges = graph_get_nEdges(prunegraph);

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

   if( graph_pc_isPc(prunegraph) )
   {
      SCIP_CALL( redLoopPc(scip, NULL, prunegraph, vnoi, path, nodearrreal, heap, state,
            vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, &(lurkprune->offsetnew),
            FALSE, FALSE, FALSE, LURKPRUNE_MINREDELIMS, FALSE, TRUE) );
   }
   else if( graph_pc_isMw(prunegraph) )
   {
      SCIP_CALL( redLoopMw(scip, prunegraph, vnoi, nodearrreal, state,
            vbase, nodearrint, NULL, nodearrchar, &(lurkprune->offsetnew),
            FALSE, FALSE, FALSE, LURKPRUNE_MINREDELIMS, FALSE) );
   }
   else
   {
      const RPARAMS parameters = { .dualascent = FALSE, .boundreduce = FALSE, .nodereplacing = TRUE,
                                                    .reductbound = LURKPRUNE_MINREDELIMS, .userec = FALSE, .fullreduce = FALSE };

      SCIP_CALL( redLoopStp(scip, &parameters, prunegraph, vnoi, path, nodearrreal, heap, state,
            vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, &(lurkprune->offsetnew)) );
   }

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

   return SCIP_OKAY;
}


/** does heuristic reductions */
static
SCIP_RETCODE reduceLurk(
   SCIP*                 scip,               /**< SCIP data structure */
   LURKPRUNE*            lurkprune,          /**< lurking-prune data */
   GRAPH*                prunegraph,         /**< graph data structure */
   SCIP_Bool*            success             /**< could enough edges be removed? */
   )
{
   IDX** ancestors = prunegraph->ancestors;
   int* lurkhalfedges;
   SCIP_Real* lurkingbounds_local;
   const SCIP_Real* lurkingbounds_half = lurkprune->lurkingbounds_half;
   const int* soledge = lurkprune->soledge;
   const int nedges = graph_get_nEdges(prunegraph);
   int i;
   int nelims = 0;
   const int nelims_min = lurkprune->minlurkelims;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(prunegraph);

   assert(ancestors && lurkingbounds_half && soledge);
   assert(solstp_isValid(scip, prunegraph, soledge));

   *success = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &lurkhalfedges, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lurkingbounds_local, nedges / 2) );

   /* get averaged lurking ancestor bounds */
   for( i = 0; i < nedges / 2; i++ )
   {
      SCIP_Real bound = 0;
      int nancestors = 0;
      const int edge = i * 2;

      for( IDX* curr = ancestors[edge]; curr != NULL; curr = curr->parent )
      {
         const int edge_ancestor = curr->index;
         assert(graph_edge_isInRange(prunegraph, edge_ancestor));

         bound += lurkingbounds_half[edge_ancestor / 2];
         nancestors++;
      }

      if( nancestors != 0 )
      {
         bound /= (SCIP_Real) nancestors;
      }
      else
      {
         assert(graph_edge_isDeleted(prunegraph, edge));
      }

      lurkingbounds_local[i] = bound;
      lurkhalfedges[i] = i;
   }

   SCIPsortDownRealInt(lurkingbounds_local, lurkhalfedges, nedges / 2);

   SCIPdebugMessage("first/last %f %f \n", lurkingbounds_local[0], lurkingbounds_local[nedges / 2 - 1]);

   for( i = 0; i < nedges / 2; i++ )
   {
      const int edge = lurkhalfedges[i] * 2;
      const int edge_rev = edge + 1;

      assert(flipedge(edge) == edge_rev);
      assert(graph_edge_isDeleted(prunegraph, edge) == graph_edge_isDeleted(prunegraph, edge_rev));

      if( graph_edge_isDeleted(prunegraph, edge) )
         continue;

      if( soledge[edge] == CONNECT || soledge[edge_rev] == CONNECT )
         continue;

      if( isPcMw && graph_pc_edgeIsExtended(prunegraph, edge) )
         continue;

      graph_edge_del(scip, prunegraph, edge, TRUE);
      nelims++;

      if( nelims >= nelims_min )
      {
         *success = TRUE;
         break;
      }
   }

   SCIP_CALL( reduceLevel0(scip, prunegraph) );
   assert(graph_valid(scip, prunegraph));

   SCIPfreeBufferArray(scip, &lurkingbounds_local);
   SCIPfreeBufferArray(scip, &lurkhalfedges);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("reduce-lurk: nelims=%d nelims_min=%d \n", nelims, nelims_min);
   graph_printInfoReduced(prunegraph);
#endif

   return SCIP_OKAY;
}


/** delete fixed edges */
static
void reduceFixedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   const int*            soledge,            /**< primal solution */
   GRAPH*                prunegraph          /**< graph data structure */
   )
{
   const int nedges = graph_get_nEdges(prunegraph);
   int nfixededges = 0;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(prunegraph);

   /* delete fixed edges from the new graph */
   for( int e = 0; e < nedges; e += 2 )
   {
      /* both e and its anti-parallel edge fixed to zero? */
      if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[e + 1]) < 0.5
         && soledge[e] != CONNECT && soledge[e + 1] != CONNECT )
      {
         if( isPcMw && graph_pc_edgeIsExtended(prunegraph, e) )
            continue;

         graph_edge_del(scip, prunegraph, e, TRUE);
         nfixededges++;
      }
   }
   SCIPdebugMessage("fixed edges in lurk and prune: %d \n", nfixededges);
}


/** updates */
static
SCIP_RETCODE updateSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   LURKPRUNE*            lurkprune,          /**< lurking-prune data */
   GRAPH*                prunegraph          /**< graph data structure */
)
{
   /* compute potential new guiding solution */
      //
   SCIP_Bool success;
   SCIP_Real obj;
   int* soledge = lurkprune->soledge;
   int* solnode = lurkprune->solnode;

   if( graph_pc_isRootedPcMw(g) || graph_typeIsSpgLike(g) )
   {
      SCIP_Real* redcost;
      DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = FALSE, .root = prunegraph->source,
                      .is_pseudoroot = FALSE, .damaxdeviation = -1.0 };
      obj = 0.0;


      SCIP_CALL( SCIPallocBufferArray(scip, &redcost, g->edges) );

      SCIP_CALL( dualascent_exec(scip, g, NULL, &daparams, redcost, &obj) );

      SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, prunegraph, redcost, soledge, -1, &success, FALSE) );
      {
         int todo; // do reducctions based on redcost
      }
      SCIPfreeBufferArray(scip, &redcost);
   }
   else
   {
      int* starts;
      int nruns = LURKPRUNE_NTMRUNS;

      SCIP_CALL( SCIPallocBufferArray(scip, &starts, prunegraph->knots) );
      SCIPStpHeurTMCompStarts(prunegraph, starts, &nruns);

      SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fromheurdata, prunegraph, starts, NULL, soledge, nruns,
            prunegraph->source, prunegraph->cost, prunegraph->cost, NULL, NULL, &success));

      SCIPfreeBufferArray(scip, &starts);
   }

   SCIP_CALL( SCIPStpHeurLocalRun(scip, prunegraph, soledge) );

   assert(solstp_isValid(scip, prunegraph, soledge));

   SCIP_CALL( SCIPStpHeurPruneUpdateSols(scip, g, prunegraph, solnode, soledge,
         lurkprune->globalsoledge, &(lurkprune->globalobj), TRUE, &success) );

   obj = solstp_getObj(prunegraph, soledge, lurkprune->offsetnew);

   /* obj < incumbent objective value? */
   if( LT(obj, lurkprune->ubbest) )
      lurkprune->ubbest = obj;

   SCIPdebugMessage("old solution: %f new solution %f \n\n", lurkprune->ubbest, obj);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLurkPrune)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurLurkPrune(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLurkPrune)
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
SCIP_DECL_HEURINIT(heurInitLurkPrune)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolLurkPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   heurdata->nfailures = 0;
   heurdata->lastnfixededges = -1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolLurkPrune)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLurkPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;
   GRAPH* graph;
   SCIP_Real* xval;
   SCIP_Bool success;
   int nedges;
   int* soledge;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   probdata = SCIPgetProbData(scip);
   graph = SCIPprobdataGetGraph(probdata);

   assert(heurdata && probdata && graph);

   nedges = graph->edges;
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;

   if( !graph_pc_isPcMw(graph) && !graph_typeIsSpgLike(graph) )
      return SCIP_OKAY;

   if( !(heurdata->maxfreq) && heurdata->nfailures > 0 )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);

   /* no solution available? */
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* heuristic not at maximum or ...*/
   if( !heurdata->maxfreq )
   {
      if( heurdata->lastnfixededges >= 0 )
      {
         SCIP_Real stallproportion;

         stallproportion = (1.0 + heurdata->nfailures) * LURKPRUNE_MINSTALLPROPORTION;

         if( SCIPisGT(scip, stallproportion, LURKPRUNE_MAXSTALLPROPORTION) )
            stallproportion = LURKPRUNE_MAXSTALLPROPORTION;

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
   assert(vars != NULL);

   /* allocate array to store primal solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &soledge, nedges) );

   for( int e = 0; e < nedges; e++ )
   {
      if( EQ(xval[e], 1.0) )
         soledge[e] = CONNECT;
      else
         soledge[e] = UNKNOWN;
   }

   /* execute lurkprune heuristic */
   SCIP_CALL( SCIPStpHeurLurkPruneRun(scip, vars, NULL, graph, FALSE, (graph->edges < LURK_MAXTOTNEDGES), soledge, &success) );

   /* solution found by lurkprune heuristic? */
   if( success )
   {
      SCIP_Real pobj = 0.0;
      SCIP_Real* nval;

      assert(SCIPprobdataGetNVars(scip) == nedges);

      /* allocate memory to store solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );

      for( int e = 0; e < nedges; e++ )
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
      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, heur, &success) );

      /* has solution been added? */
      if( success )
      {
         SCIPdebugMessage("better solution added by LURKPRUNE %f \n", pobj + SCIPprobdataGetOffset(scip));
         *result = SCIP_FOUNDSOL;
         assert(solstp_isValid(scip, graph, soledge));

         /* is the solution the new incumbent? */
         if( SCIPisLT(scip, SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj) )
            success = FALSE;
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &nval);
   }
   else
   {
      SCIPdebugMessage("lurk-prune: solution not valid \n");
   }

   /* solution could not be found, added or is not best? */
   if( !success )
      heurdata->nfailures++;

   /* free memory */
   SCIPfreeBufferArray(scip, &soledge);
   SCIPdebugMessage("leaving lurk and prune heuristic \n");

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** execute lurk-and-prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurLurkPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   const SCIP_Real*      lurkingbounds,      /**< lurking edge bounds */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool             initialreduce,      /**< try to reduce graph initially? */
   SCIP_Bool             ascendprune,        /**< use ascend-prune? */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   SCIP_Bool*            solimproved         /**< could a better solution be found? */
   )
{
   LURKPRUNE lurkprune;
   GRAPH *prunegraph;

   assert(scip && soledge && solimproved && lurkingbounds);
   assert(!graph_pc_isPcMw(g) || g->extended);
   assert(solstp_isValid(scip, g, soledge));

   *solimproved = FALSE;
   SCIP_CALL( lurkpruneInit(scip, g, lurkingbounds, soledge, &lurkprune, &prunegraph) );

   if( vars != NULL )
   {
      reduceFixedVars(scip, vars, soledge, prunegraph);
   }

   /* perform initial reductions? */
   if( initialreduce )
   {
      SCIP_CALL( reduceExact(scip, &lurkprune, prunegraph) );
      SCIP_CALL( updateSolution(scip, g, &lurkprune, prunegraph) );
   }

   /* main reduction loop */
   for( int i = 0; i < LURKPRUNE_MAXREDROUNDS && g->terms > 2; i++ )
   {
      SCIP_Bool lurksuccess;

      SCIPdebugMessage("starting round %d \n", i);

      SCIP_CALL( reduceLurk(scip, &lurkprune, prunegraph, &lurksuccess) );
      SCIP_CALL( reduceExact(scip, &lurkprune, prunegraph) );
      SCIP_CALL( updateSolution(scip, g, &lurkprune, prunegraph) );

      if( !lurksuccess )
      {
         SCIPdebugMessage("breaking early \n");
         break;
      }
   }

   lurkpruneFinalize(scip, g, prunegraph, soledge, &lurkprune, solimproved);
   assert(solstp_isValid(scip, g, soledge));

   //exit(1);

   return SCIP_OKAY;
}


/** creates the lurkprune primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurLurkPrune(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create lurkprune primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLurkPrune, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLurkPrune) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLurkPrune) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLurkPrune) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLurkPrune) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLurkPrune) );

   /* add lurkprune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_LURKPRUNE_MAXFREQ, NULL, NULL) );


   return SCIP_OKAY;
}
