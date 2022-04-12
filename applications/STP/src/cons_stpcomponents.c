/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
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
#include "bidecomposition.h"
#include "prop_stp.h"
#include "substpsolver.h"
#include "graph.h"
#include "solstp.h"

/**@name Constraint handler properties
 *
 * @{
 */

#define CONSHDLR_NAME          "stpcomponents"
#define CONSHDLR_DESC          "steiner tree components constraint handler"
#define CONSHDLR_SEPAPRIORITY        -1 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY        -1 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY       -1 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             0 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_BEFORELP // SCIP_PROPTIMING_DURINGLPLOOP //  SCIP_PROPTIMING_BEFORELP
#define DECOMP_MAXCOMPRATIO 0.95
#define DECOMP_MINNNODES     10


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
   CUTNODES*             cutnodes;
   BIDECOMP*             bidecomposition;
   SCIP_Bool             probDecompIsReady;  /**< decomposition available?*/
   SCIP_Bool             probWasDecomposed;  /**< already decomposed? */
};


/** sub-solution */
typedef struct sub_solution
{
   int*                  subedgesSol;        /**< solution: UNKNOWN/CONNECT */
   int                   nsubedges;          /**< number of edges of subproblem */
} SUBSOL;


/*
 * Local methods
 */



/** gets optimal sub-problem solution */
static
SCIP_RETCODE subsolGet(
   SUBSTP*               substp,             /**< sub-problem data structure */
   SUBSOL*               subcomp             /**< component */
   )
{
   SCIP_CALL( substpsolver_getSolution(substp, subcomp->subedgesSol) );

   return SCIP_OKAY;
}


/** fixes original edges */
static
SCIP_RETCODE subsolFixOrgEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBINOUT*       subinout,           /**< helper for problem mapping */
   SUBSOL*               subsol              /**< solution */
   )
{
   SCIP_VAR** orgvars = SCIPprobdataGetVars(scip);
   const int* const subedgesSol = subsol->subedgesSol;
   const int* const subedgesToOrg = graph_subinoutGetSubToOrgEdgeMap(subinout);
   const int nsubedges = subsol->nsubedges;
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
SCIP_RETCODE subsolInit(
   SCIP*                 scip,               /**< sub-SCIP data structure */
   const SUBSTP*         substp,             /**< sub problem */
   SUBSOL**              subsolution        /**< to initialize */
   )
{
   SUBSOL* subsol;

   assert(scip && substp);

   SCIP_CALL( SCIPallocMemory(scip, subsolution) );
   subsol = *subsolution;

   subsol->nsubedges = substpsolver_getNsubedges(substp);
   assert(subsol->nsubedges >= 2);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(subsol->subedgesSol), subsol->nsubedges) );

   return SCIP_OKAY;
}


/** exits */
static
void subsolFree(
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SUBSOL**              subsolution         /**< to initialize */
   )
{
   SUBSOL* subsol;

   assert(subsolution);
   subsol = *subsolution;
   assert(subsol);

   SCIPfreeMemoryArray(subscip, &(subsol->subedgesSol));
   SCIPfreeMemory(subscip, subsolution);
}


/** fixes original edges with sub-solution */
static
SCIP_RETCODE subcompFixOrgEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBINOUT*       subinout,           /**< sub-problem insertion/extraction data structure */
   SUBSTP*               substp              /**< sub-problem data structure */
   )
{
   SUBSOL* subcomp;

   assert(scip && subinout && substp);

   SCIP_CALL( subsolInit(scip, substp, &subcomp) );
   SCIP_CALL( subsolGet(substp, subcomp) );
   SCIP_CALL( subsolFixOrgEdges(scip, subinout, subcomp) );

   subsolFree(scip, &subcomp);

   return SCIP_OKAY;
}

/** is promising? */
static
SCIP_Bool decomposeIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const BIDECOMP*       bidecomp            /**< bidecomposition data structure */
   )
{
   if( g->knots < DECOMP_MINNNODES )
   {
      return FALSE;
   }
   else
   {
      const SCIP_Real maxcompratio = DECOMP_MAXCOMPRATIO;
      const SCIP_Real maxratio = bidecomposition_getMaxcompNodeRatio(bidecomp);

      assert(GT(maxratio, 0.0));

      return (maxratio < maxcompratio);
   }
}


/** gets subgraph */
static
SCIP_RETCODE decomposeGetSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                orggraph,           /**< graph data structure */
   GRAPH**               subgraph            /**< subgraph */
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
   SUBINOUT* subinout;
   SUBSTP* substp;
   GRAPH* subgraph;

   assert(scip && orggraph && bidecomp && success);
   assert(graph_subinoutUsesNewHistory(bidecomp->subinout));

   *success = TRUE;
   subinout = bidecomp->subinout;

   if( bidecomposition_componentIsTrivial(bidecomp, compindex) )
   {
      SCIPdebugMessage("trivial component, skip! \n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("building subproblem... \n");

   SCIP_CALL( decomposeGetSubgraph(scip, bidecomp, compindex, orggraph, &subgraph) );

   /* NOTE: subgraph will be moved into substp! */
   SCIP_CALL( substpsolver_init(scip, subgraph, &substp) );
   SCIP_CALL( substpsolver_transferHistory(graph_subinoutGetSubToOrgEdgeMap(subinout),
      orggraph, substp) );

   {
      int verblevel;
      SCIP_CALL( SCIPgetIntParam(scip, "display/verblevel", &verblevel) );
      if( verblevel == 0 )
         SCIP_CALL( substpsolver_setMute(substp) );
   }

   SCIP_CALL( substpsolver_solve(scip, substp, success) );

   if( *success )
   {
      /* here we basically fix the entire sub-problem */
      SCIP_CALL( subcompFixOrgEdges(scip, subinout, substp) );
   }

   graph_subinoutClean(scip, subinout);
   substpsolver_free(scip, &substp);

   return SCIP_OKAY;
}


/** tries to decompose and solve */
static
SCIP_RETCODE decomposeExec(
   SCIP*                 scip,               /**< SCIP data structure */
   BIDECOMP*             bidecomp,           /**< bidecomposition data structure */
   CUTNODES*             cutnodes,           /**< cut nodes data structure */
   GRAPH*                orggraph,           /**< graph to decompose */
   SCIP_Bool*            success             /**< decomposed? */
   )
{
   assert(bidecomp && cutnodes && success);
   assert(graph_valid(scip, orggraph));
   assert(*success == FALSE);

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

   return SCIP_OKAY;
}


/** initializes for conshdlrdata */
static
SCIP_RETCODE initDecompose(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< cosntraint handler data */
   GRAPH*                orggraph,           /**< graph to decompose */
   SCIP_Bool*            isPromsing          /**< promising decomposition? */
   )
{
   BIDECOMP* bidecomp;
   CUTNODES* cutnodes;

   assert(scip && conshdlrdata && orggraph && isPromsing);
   assert(!conshdlrdata->cutnodes);
   assert(!conshdlrdata->bidecomposition);
   assert(graph_typeIsSpgLike(orggraph) && "only SPG decomposition supported yet");

   *isPromsing = FALSE;

   if( orggraph->terms == 1 )
   {
      SCIPdebugMessage("only one terminal...don't decompose \n");
      return SCIP_OKAY;
   }

   graph_mark(orggraph);

   if( !bidecomposition_isPossible(orggraph) )
   {
      SCIPdebugMessage("graph is too large...don't decompose \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( bidecomposition_cutnodesInit(scip, orggraph, &cutnodes) );
   conshdlrdata->cutnodes = cutnodes;
   bidecomposition_cutnodesCompute(orggraph, cutnodes);

   if( cutnodes->biconn_ncomps > 0 )
   {
      SCIP_CALL( bidecomposition_init(scip, cutnodes, orggraph, &bidecomp) );
      conshdlrdata->bidecomposition = bidecomp;

      *isPromsing = decomposeIsPromising(orggraph, bidecomp);
   }

   return SCIP_OKAY;
}



/** tries to decompose and solve */
static
void freeDecompose(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   assert(scip && conshdlrdata);

   if( conshdlrdata->cutnodes )
   {
      bidecomposition_cutnodesFree(scip, &(conshdlrdata->cutnodes));
   }

   if( conshdlrdata->bidecomposition )
   {
      bidecomposition_free(scip, &(conshdlrdata->bidecomposition));
   }
}

/** decomposes and solves */
static
SCIP_RETCODE divideAndConquer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraints handler data */
   SCIP_Bool*            success             /**< decomposed? */
   )
{
   GRAPH* orggraph = SCIPprobdataGetGraph2(scip);

   assert(conshdlrdata && success);
   assert(conshdlrdata->bidecomposition && conshdlrdata->cutnodes);
   assert(orggraph->terms > 1);
   assert(decomposeIsPromising(orggraph, conshdlrdata->bidecomposition));
   assert(graph_typeIsSpgLike(orggraph) && "only SPG bidecomposition supported yet");
   assert(conshdlrdata->cutnodes->biconn_ncomps > 0);

   *success = FALSE;
   graph_mark(orggraph);

   SCIP_CALL( decomposeExec(scip, conshdlrdata->bidecomposition, conshdlrdata->cutnodes, orggraph, success) );

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
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata);

   freeDecompose(scip, conshdlrdata);

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

   if( !conshdlrdata->probDecompIsReady )
      return SCIP_OKAY;

   assert(!conshdlrdata->probWasDecomposed);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( divideAndConquer(scip, conshdlrdata, &success) );

   if( success )
   {
      printf("problem solved by bidecomposition \n");
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


/** sets the data for bidecomposition up  */
SCIP_RETCODE SCIPStpcomponentsSetUp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< graph data */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool isPromising = FALSE;

   assert(scip && graph);

   conshdlr = SCIPfindConshdlr(scip, "stpcomponents");
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata);

   SCIP_CALL( initDecompose(scip, conshdlrdata, graph, &isPromising) );

   conshdlrdata->probDecompIsReady = isPromising;

   if( !isPromising )
   {
      freeDecompose(scip, conshdlrdata);
   }

   return SCIP_OKAY;
}


/** is a promising bidecomposition available? */
SCIP_Bool SCIPStpcomponentsAllowsDecomposition(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip);

   conshdlr = SCIPfindConshdlr(scip, "stpcomponents");
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata);

   return conshdlrdata->probDecompIsReady;
}

/** creates the handler for stp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrStpcomponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create stp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->cutnodes = NULL;
   conshdlrdata->bidecomposition = NULL;
   conshdlrdata->probWasDecomposed = FALSE;
   conshdlrdata->probDecompIsReady = FALSE;

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
