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

/**@file   cons_stpcomponents.c
 * @brief  Component constraint handler for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file checks for biconnected components and solves them separately.
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

#include "cons_stpcomponents.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "bidecomposition.h"
#include "cons_stp.h"
#include "reduce.h"
#include "heur_tm.h"
#include "reader_stp.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "heur_ascendprune.h"
#include "heur_slackprune.h"
#include "heur_lurkprune.h"
#include "heur_rec.h"
#include "pricer_stp.h"
#include "probdata_stp.h"
#include "dialog_stp.h"
#include "prop_stp.h"
#include "branch_stp.h"
#include "solhistory.h"
#include "graph.h"
#include "solstp.h"

/**@name Constraint handler properties
 *
 * @{
 */

#define CONSHDLR_NAME          "stpcomponents"
#define CONSHDLR_DESC          "steiner tree components constraint handler"
#define CONSHDLR_SEPAPRIORITY   9999999 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             0 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_BEFORELP // SCIP_PROPTIMING_DURINGLPLOOP //  SCIP_PROPTIMING_BEFORELP
#define DAQ_MINCOMPRATIO 0.9


/**@} */

/*
 * Data structures
 */

/** @brief Constraint data for  \ref cons_stp.c "Stp" constraints */
struct SCIP_ConsData
{
   GRAPH*                graph;              /**< graph data structure */
};


/** @brief Constraint handler data for \ref cons_stp.c "Stp" constraint handler */
struct SCIP_ConshdlrData
{
   SCIP_Bool             probWasDecomposed;  /**< already decomposed? */
};


/** representation of a biconnected component */
typedef struct sub_component
{
   int*                  subedgesSol;        /**< solution: UNKNOWN/CONNECT */
   int                   subroot;            /**< root for orientation of solution */
   int                   nsubedges;          /**< number of edges */
   int                   nsubnodes;          /**< number of nodes */
} SUBCOMP;


/*
 * Local methods
 */



/** gets optimal sub-problem solution */
static
SCIP_RETCODE subsolGet(
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBCOMP*              subcomp             /**< component */
   )
{

   SOLHISTORY* solhistory;
   SCIP_SOL* subsol = SCIPgetBestSol(subscip);
   GRAPH* subgraph = SCIPprobdataGetGraph2(subscip);
   const STP_Bool* edges_isInSol;
   int* subedgesSol = subcomp->subedgesSol;
   const int nsubedges = subcomp->nsubedges;

   assert(subgraph && subsol);
   assert(nsubedges > 0);

   SCIP_CALL( solhistory_init(subscip, subgraph, &solhistory) );
   SCIP_CALL( solhistory_computeHistory(subscip, subsol, subgraph, solhistory) );

   assert(nsubedges == solhistory->norgedges);
   edges_isInSol = solhistory->orgedges_isInSol;

   for( int i = 0; i < nsubedges; i++ )
   {
      if( edges_isInSol[i] )
      {
         subedgesSol[i] = CONNECT;
      }
      else
      {
         subedgesSol[i] = UNKNOWN;
      }
   }

   solhistory_free(subscip, &solhistory);

   return SCIP_OKAY;
}


/** fixes original edges */
static
SCIP_RETCODE subsolFixOrgEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBINOUT*       subinout,
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBCOMP*              subcomp             /**< component */
   )
{
   SCIP_VAR** orgvars = SCIPprobdataGetVars(scip);
   const int* const subedgesSol = subcomp->subedgesSol;
   const int* const subedgesToOrg = graph_subinoutGetSubToOrgEdgeMap(subinout);
   const int nsubedges = subcomp->nsubedges;
   SCIP_Bool success;
#ifndef NDEBUG
   GRAPH* orggraph = SCIPprobdataGetGraph2(scip);
#endif
   assert(orggraph);
   assert(orgvars);
   assert(subedgesSol && subedgesToOrg);
   assert(nsubedges > 0);

   for( int i = 0; i < nsubedges; i++ )
   {
      const int orgedge = subedgesToOrg[i];
      assert(graph_edge_isInRange(orggraph, orgedge));

      if( subedgesSol[i] == CONNECT )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("fix to 1: ");
         graph_edge_printInfo(SCIPprobdataGetGraph2(scip), orgedge);
#endif
         SCIP_CALL( SCIPStpFixEdgeVarTo1(scip, orgvars[orgedge], &success) );
      }
      else
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("fix to 0: ");
         graph_edge_printInfo(SCIPprobdataGetGraph2(scip), orgedge);
#endif
         assert(subedgesSol[i] == UNKNOWN);
         SCIP_CALL( SCIPStpFixEdgeVarTo0(scip, orgvars[orgedge], &success) );
      }
   }

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE subcompInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBCOMP**             subcomponent        /**< to initialize */
   )
{
   SUBCOMP* subcomp;
   const GRAPH* const subgraph = SCIPprobdataGetGraph2(subscip);
   const int nsubnodes = graph_get_nNodes(subgraph);
   const int nsubedges = graph_get_nEdges(subgraph);

   assert(scip && subscip);
   assert(nsubnodes > 0);
   assert(nsubedges > 0);

   SCIP_CALL( SCIPallocMemory(subscip, subcomponent) );
   subcomp = *subcomponent;

   SCIP_CALL( SCIPallocMemoryArray(subscip, &(subcomp->subedgesSol), nsubedges) );

   subcomp->nsubnodes = nsubnodes;
   subcomp->nsubedges = nsubedges;
   subcomp->subroot = subgraph->source;

   return SCIP_OKAY;
}


/** fixes original edges with sub-solution */
static
SCIP_RETCODE subcompFixOrgEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBINOUT*       subinout,
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBCOMP*              subcomp             /**< component */
   )
{
   assert(scip && subscip && subcomp);

   SCIP_CALL( subsolGet(subscip, subcomp) );
   SCIP_CALL( subsolFixOrgEdges(scip, subinout, subscip, subcomp) );

   return SCIP_OKAY;
}


/** exits */
static
void subcompFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBCOMP**             subcomponent        /**< to initialize */
   )
{
   SUBCOMP* subcomp;

   assert(scip && subcomponent);
   subcomp = *subcomponent;
   assert(subcomp);

   SCIPfreeMemoryArray(subscip, &(subcomp->subedgesSol));
   SCIPfreeMemory(subscip, subcomponent);
}


/** helper */
static
SCIP_RETCODE subScipSetupCallbacks(
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludePricerStp(subscip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   SCIP_CALL( SCIPincludeReaderStp(subscip) );

   SCIP_CALL( SCIPincludeDialogStp(subscip) );

   SCIP_CALL( SCIPincludeConshdlrStp(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurTM(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurLocal(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurRec(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurPrune(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurAscendPrune(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurSlackPrune(subscip) );

   SCIP_CALL( SCIPStpIncludeHeurLurkPrune(subscip) );

   SCIP_CALL( SCIPincludeBranchruleStp(subscip) );

   SCIP_CALL( SCIPincludePropStp(subscip) );

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE subScipSetupParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_Real timelimit;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   timelimit -= SCIPgetSolvingTime(scip);

   assert(GT(timelimit, 0.0));


   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   //SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   /* disable expensive resolving */
   // SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* disable STP presolving */
   SCIP_CALL( SCIPsetIntParam(subscip, "stp/reduction", 0) );

   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );

   /* set hard-coded default parameters */
   SCIP_CALL( SCIPprobdataSetDefaultParams(subscip) );

   return SCIP_OKAY;
}


/** sets up the sub-SCIP */
static
SCIP_RETCODE subScipSetup(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_CALL( subScipSetupCallbacks(subscip) );

   SCIP_CALL( subScipSetupParameters(scip, subscip) );

   return SCIP_OKAY;
}


/** is promising? */
static
SCIP_Bool decomposeIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const BIDECOMP*       bidecomp
   )
{
   const SCIP_Real mincompratio = DAQ_MINCOMPRATIO;
   const SCIP_Real maxratio = bidecomposition_getMaxcompNodeRatio(bidecomp);

   assert(GT(maxratio, 0.0));

   return (maxratio < mincompratio);
}


/** gets subgraph */
static
SCIP_RETCODE decomposeGetSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                orggraph,           /**< graph data structure */
   GRAPH**               subgraph
   )
{
/*
   SCIP_CALL( graph_copy(scip, orggraph, subgraph) );
   (*subgraph)->is_packed = FALSE;
   return SCIP_OKAY;
*/
   int subroot;
   GRAPH* subg;
   SUBINOUT* subinout = bidecomp->subinout;
   assert(graph_valid(scip, orggraph));

   bidecomposition_markSub(bidecomp, compindex, orggraph);
   SCIP_CALL( graph_subgraphExtract(scip, orggraph, subinout, subgraph) );
   subg = *subgraph;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("original nodes of connected components: \n");

   for( int i = 0; i < orggraph->knots; i++ )
   {
      if( orggraph->mark[i] )
         graph_knot_printInfo(orggraph, i);
   }
#endif

   SCIP_CALL( bidecomposition_getMarkedSubRoot(scip, bidecomp, orggraph, subg, &subroot) );
   assert(graph_knot_isInRange(subg, subroot));
   assert(Is_term(subg->term[subroot]));
   subg->source = subroot;

   return SCIP_OKAY;
}


/** solves subproblem */
static
SCIP_RETCODE decomposeSolveSub(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                orggraph,           /**< graph data structure */
   SCIP_Bool*            success             /**< solved? */
   )
{
   SCIP* subscip;
   SUBCOMP* subcomp;
   GRAPH* subgraph;

   assert(scip && orggraph && bidecomp && success);

   *success = TRUE;

   if( bidecomposition_componentIsTrivial(bidecomp, compindex) )
   {
      SCIPdebugMessage("trivial component, skip! \n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("building subproblem... \n");

   SCIP_CALL( decomposeGetSubgraph(scip, bidecomp, compindex, orggraph, &subgraph) );

   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( subScipSetup(scip, subscip) );
   SCIP_CALL( graph_subinoutCompleteNewHistory(subscip, orggraph, bidecomp->subinout, subgraph) );
   /* NOTE: subgraph will be moved into subscip probdata! */
   SCIP_CALL( SCIPprobdataCreateFromGraph(subscip, 0.0, "subproblem", TRUE, subgraph) );
   SCIP_CALL( subcompInit(scip, subscip, &subcomp) );
   SCIP_CALL( SCIPsolve(subscip) );

   SCIPdebugMessage("subproblem has been solved \n");

   if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
   {
      /* here we basically solve the entire problem */
      SCIP_CALL( subcompFixOrgEdges(scip, bidecomp->subinout, subscip, subcomp) );
   }
   else
   {
      *success = FALSE;

      printf("solving of sub-problem interrupted (status=%d, time=%.2f)\n",
          SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));
   }

   graph_subinoutClean(scip, bidecomp->subinout);
   subcompFree(scip, subscip, &subcomp);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** tries to decompose and solve */
static
SCIP_RETCODE decomposeExec(
   SCIP*                 scip,               /**< SCIP data structure */
   CUTNODES*             cutnodes,
   GRAPH*                orggraph,           /**< graph to decompose */
   SCIP_Bool*            success             /**< decomposed? */
   )
{
   BIDECOMP* bidecomp;

   assert(graph_valid(scip, orggraph));
   assert(*success == FALSE);

   SCIP_CALL( bidecomposition_init(scip, cutnodes, orggraph, &bidecomp) );

   if( decomposeIsPromising(orggraph, bidecomp) )
   {
      SCIP_CALL( bidecomposition_initSubInOut(scip, orggraph, bidecomp) );
      SCIP_CALL( graph_subinoutActivateEdgeMap(orggraph, bidecomp->subinout) );
      graph_subinoutActivateNewHistory(bidecomp->subinout);

      printf("solving problem by decomposition (%d components) \n",  bidecomp->nbicomps);

      *success = TRUE;

      /* solve each biconnected component individually */
      for( int i = 0; i < bidecomp->nbicomps; i++ )
      {
         SCIP_CALL( decomposeSolveSub(scip, bidecomp, i, orggraph, success) );

         if( *success == FALSE )
         {
            printf("could not solve component %d; aborting decomposition now \n", i);
            break;
         }
      }
   }

   bidecomposition_free(scip, &bidecomp);

   return SCIP_OKAY;
}


/** decomposes and solves */
static
SCIP_RETCODE divideAndConquer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            success             /**< decomposed? */
   )
{
   GRAPH* orggraph = SCIPprobdataGetGraph2(scip);
   CUTNODES* cutnodes;

   assert(success);
   assert(orggraph->terms > 1);
   assert(graph_typeIsSpgLike(orggraph) && "only SPG decomposition supported yet");

   *success = FALSE;
   graph_mark(orggraph);

   if( !bidecomposition_isPossible(orggraph) )
   {
      SCIPdebugMessage("graph is too large...don't decompose \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( bidecomposition_cutnodesInit(scip, orggraph, &cutnodes) );
   bidecomposition_cutnodesCompute(orggraph, cutnodes);

   if( cutnodes->biconn_ncomps > 0 )
   {
      SCIP_Real fixed = 0.0;
      int nelims = 0;

      // todo this method changes the graph! probably only good for propgraph, otherwise gets into troubles!
      SCIP_CALL( reduce_nonTerminalComponents(scip, cutnodes, orggraph, &fixed, &nelims) );

      SCIPdebugMessage("simple reductions: %d \n", nelims);

      /* try to decompose and reduce recursively */
      SCIP_CALL( decomposeExec(scip, cutnodes, orggraph, success) );
   }

   bidecomposition_cutnodesFree(scip, &cutnodes);

   return SCIP_OKAY;
}


/*
 * Callbacks
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyStpcomponents)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrStpcomponents(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeStpcomponents)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolStpcomponents)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXITSOL(consExitsolStpcomponents)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteStpcomponents)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransStpcomponents)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   targetdata->graph = sourcedata->graph;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpStpcomponents)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropStpcomponents)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   GRAPH* graph = SCIPprobdataGetGraph2(scip);
   SCIP_Bool success;

   assert(conshdlrdata);
   assert(graph);

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   if( graph->terms == 1 )
      return SCIP_OKAY;

   if( !graph_typeIsSpgLike(graph) )
      return SCIP_OKAY;

   assert(!conshdlrdata->probWasDecomposed);

   *result = SCIP_DIDNOTFIND;

   // todo call update propgraph from prop_stp
   // todo ... get the graph from prop_stp

   SCIP_CALL( divideAndConquer(scip, &success) );

   if( success )
   {
      printf("problem solved by decomposition \n");
      conshdlrdata->probWasDecomposed = TRUE;
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockStpcomponents)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


#define consEnfolpStpcomponents NULL
#define consEnfopsStpcomponents NULL
#define consCheckStpcomponents NULL

/*
 * Interface methods
 */

/** creates the handler for stp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrStpcomponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create stp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->probWasDecomposed = FALSE;

   conshdlr = NULL;
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpStpcomponents, consEnfopsStpcomponents, consCheckStpcomponents, consLockStpcomponents,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyStpcomponents, NULL) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteStpcomponents) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransStpcomponents) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolStpcomponents) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolStpcomponents) );

   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpStpcomponents) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropStpcomponents, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeStpcomponents) );

   return SCIP_OKAY;
}
