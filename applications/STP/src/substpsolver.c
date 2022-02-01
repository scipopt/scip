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

/**@file   substpsolver.c
 * @brief  Solver for Steiner tree (sub-) problems
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "substpsolver.h"
#include "solhistory.h"

#include "pricer_stp.h"
#include "probdata_stp.h"
#include "cons_stp.h"
#include "reduce.h"
#include "cons_stpcomponents.h"
#include "heur_tm.h"
#include "reader_stp.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "heur_ascendprune.h"
#include "heur_slackprune.h"
#include "heur_lurkprune.h"
#include "heur_rec.h"
#include "dialog_stp.h"
#include "prop_stp.h"
#include "branch_stp.h"
#include "dpterms.h"
#include "solstp.h"
#include "relax_stp.h"
#include "relax_stpdp.h"


//#define CHECK_BCOPT


/*
 * Data structures
 */


/** sub-problem */
struct sub_steiner_tree_problem
{
   SCIP*                 subscip;            /**< SCIP or NULL (if dynamic programming is used) */
   GRAPH*                subgraph;           /**< subgraph; OWNED! */
   int*                  dpsubsol;           /**< solution storage for dynamic programming case */
   int                   nsubedges;          /**< number of edges */
   SCIP_Bool             useDP;
   SCIP_Bool             useOutput;
   SCIP_Bool             subprobIsIndependent;
};



/*
 * Local methods
 */

#ifdef CHECK_BCOPT
#define CHECK_BCOPT_MAXNTERMS 50

/** checks solution with dynamic programming */
static inline
SCIP_RETCODE checkSolWithDP(
   GRAPH*                graph,              /**< sub-problem */
   const int*            edgesol             /**< array to store edge solution; CONNECT/UNKNOWN */
)
{
   return SCIP_ERROR;
}

#endif




/** helper */
static
SCIP_RETCODE subscipSetupCallbacks(
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludePricerStp(subscip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   SCIP_CALL( SCIPincludeReaderStp(subscip) );

   SCIP_CALL( SCIPincludeDialogStp(subscip) );

   SCIP_CALL( SCIPincludeConshdlrStp(subscip) );

   SCIP_CALL( SCIPincludeConshdlrStpcomponents(subscip) );

   SCIP_CALL( SCIPincludeRelaxStp(subscip) );

   SCIP_CALL( SCIPincludeRelaxStpdp(subscip) );

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
SCIP_RETCODE subscipSetupParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBSTP*         substp,             /**< sub-problem */
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_Real timelimit;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   timelimit -= SCIPgetSolvingTime(scip);

   assert(GT(timelimit, 0.0));

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console? */
   if( !substp->useOutput )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   }

   if( substp->subprobIsIndependent )
   {
      /* disable statistic timing inside sub SCIP */
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
   }
   else
   {
      /* disable STP presolving */
      SCIP_CALL( SCIPsetIntParam(subscip, "stp/reduction", 0) );
   }

   /* disable expensive resolving */
   // SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );

   /* set hard-coded default parameters */
   SCIP_CALL( SCIPprobdataSetDefaultParams(subscip) );

   return SCIP_OKAY;
}


/** sets up the sub-SCIP */
static
SCIP_RETCODE subscipSetup(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBSTP*         substp,             /**< sub-problem */
   SCIP*                 subscip             /**< sub-SCIP data structure */
   )
{
   SCIP_CALL( subscipSetupCallbacks(subscip) );

   SCIP_CALL( subscipSetupParameters(scip, substp, subscip) );

   return SCIP_OKAY;
}


/** solves sub-problem by using branch-and-cut */
static
SCIP_RETCODE subscipSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBSTP*               substp,             /**< sub-problem */
   SCIP_Bool*            success             /**< was sub-problem solved to optimality? */
)
{
   SCIP* subscip = substp->subscip;
   GRAPH* subgraph = substp->subgraph;
   const SCIP_Bool isSubProb = !substp->subprobIsIndependent;

   assert(subscip);
   assert(subgraph);

   *success = TRUE;

   SCIP_CALL( SCIPprobdataCreateFromGraph(subscip, 0.0, "subproblem", isSubProb, subgraph) );
   SCIP_CALL( SCIPsolve(subscip) );

   if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
   {
      SCIPdebugMessage("subproblem has been solved to optimality \n");
   }
   else
   {
      *success = FALSE;

      printf("solving of sub-problem interrupted (status=%d, time=%.2f)\n",
          SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));
   }

   /* NOTE: already moved into subscip */
   substp->subgraph = NULL;

   return SCIP_OKAY;
}


/** gets solution */
static
SCIP_RETCODE subscipGetSol(
   SUBSTP*               substp,             /**< sub-problem */
   int*                  subedgesSol         /**< solution */
   )
{
   SCIP* subscip = substp->subscip;
   SOLHISTORY* solhistory;
   SCIP_SOL* subsol = SCIPgetBestSol(subscip);
   GRAPH* subgraph = SCIPprobdataGetGraph2(subscip);
   const STP_Bool* edges_isInSol;
   const int nsubedges = substp->nsubedges;

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

/** Initializes default stuff */
static
SCIP_RETCODE initDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< subproblem to initialize from; NOTE: will be moved */
   SUBSTP**              substp              /**< initialize */
)
{
   SUBSTP* sub;

   assert(scip && subgraph);
   assert(subgraph->edges >= 2);

   SCIP_CALL( SCIPallocMemory(scip, substp) );
   sub = *substp;

   sub->subscip = NULL;
   sub->subgraph = subgraph;
   sub->nsubedges = subgraph->edges;
   sub->useOutput = TRUE;
   sub->subprobIsIndependent = FALSE;
   sub->dpsubsol = NULL;

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** Initializes from given sub-graph. Decides automatically whether to use DP or B&C.
 *  NOTE: Sub-graph will be moved into internal data structure and also released
 *  once substp is freed! */
SCIP_RETCODE substpsolver_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< subproblem to initialize from; NOTE: will be moved */
   SUBSTP**              substp              /**< initialize */
)
{
   /* decide whether to do dynamic programming or branch-and-cut */
   if( dpterms_isPromisingFully(subgraph) )
   {
      SCIP_CALL( substpsolver_initDP(scip, subgraph, substp) );
   }
   else
   {
      SCIP_CALL( substpsolver_initBC(scip, subgraph, substp) );

   }

   return SCIP_OKAY;
}


/** Initializes from given sub-graph, and always uses dynamic programming for solving
 * NOTE: Sub-graph will be moved into internal data structure and also released
 *  once substp is freed! */
SCIP_RETCODE substpsolver_initDP(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< subproblem to initialize from; NOTE: will be moved */
   SUBSTP**              substp              /**< initialize */
)
{
   SUBSTP* sub;

   SCIP_CALL( initDefault(scip, subgraph, substp) );

   sub = *substp;
   sub->useDP = TRUE;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sub->dpsubsol), subgraph->edges) );

   return SCIP_OKAY;
}



/** Initializes from given sub-graph, and always uses B&C for solving
 * NOTE: Sub-graph will be moved into internal data structure and also released
 *  once substp is freed! */
SCIP_RETCODE substpsolver_initBC(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< subproblem to initialize from; NOTE: will be moved */
   SUBSTP**              substp              /**< initialize */
)
{
   SUBSTP* sub;

   SCIP_CALL( initDefault(scip, subgraph, substp) );

   sub = *substp;
   sub->useDP = FALSE;
   SCIP_CALL( SCIPcreate(&(sub->subscip)) );
   SCIP_CALL( subscipSetup(scip, sub, sub->subscip) );

   return SCIP_OKAY;
}



/** frees */
void substpsolver_free(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBSTP**              substp              /**< to be freed */
   )
{
   SUBSTP* sub;

   assert(scip && substp);

   sub = *substp;

   SCIPfreeMemoryArrayNull(scip, &(sub->dpsubsol));

   if( sub->subscip )
   {
      SCIPfree(&(sub->subscip));
   }

   SCIPfreeMemory(scip, substp);
}


/** Transfers history. Useful in order to have (fixed and edge) pseudo-ancestors up-to-date.
 *  NOTE: fixed edges are not transferred! Only pseudo-ancestors */
SCIP_RETCODE substpsolver_transferHistory(
   const int*            edgeMapToOrg,       /**< maps edges of subgraph to original graph */
   GRAPH*                orggraph,           /**< original graph */
   SUBSTP*               substp              /**< sub-problem*/
)
{
   assert(edgeMapToOrg && orggraph && substp);
   assert(substp->subgraph);

   /* NOTE: no history needed in dynamic programming algorithm
    * todo might want to change that to exploit ancestor conflicts */
   if( substp->useDP )
   {
      assert(!substp->subscip);
      return SCIP_OKAY;
   }

   assert(substp->subscip);

   SCIP_CALL( graph_subgraphCompleteNewHistory(substp->subscip, edgeMapToOrg, orggraph, substp->subgraph) );

   return SCIP_OKAY;
}


/** initializes history, but does not transfer any information! */
SCIP_RETCODE substpsolver_initHistory(
   SUBSTP*               substp              /**< sub-problem*/
)
{
   assert(substp);
   assert(substp->subgraph);

   if( substp->useDP )
   {
      assert(!substp->subscip);
      return SCIP_OKAY;
   }

   assert(substp->subscip);

   SCIP_CALL( graph_initHistory(substp->subscip, substp->subgraph) );

   return SCIP_OKAY;
}


/** solves sub-problem */
SCIP_RETCODE substpsolver_solve(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBSTP*               substp,             /**< sub-problem */
   SCIP_Bool*            success             /**< was sub-problem solved to optimality? */
)
{
   assert(scip && substp && success);
   assert(substp->subgraph);

   if( substp->useDP )
   {
      assert(substp->dpsubsol);

      if( !graph_path_exists(substp->subgraph) )
      {
         SCIP_CALL( graph_path_init(scip, substp->subgraph) );
      }

      SCIP_CALL( dpterms_solve(scip, substp->subgraph, substp->dpsubsol, success) );

      // todo add an extra flag to avoid deletion of graph in some cases?
      graph_free(scip, &(substp->subgraph), TRUE);
   }
   else
   {
      SCIP_CALL( subscipSolve(scip, substp, success) );

   }

   assert(substp->subgraph == NULL);

   return SCIP_OKAY;
}


/** gets number of edges of subproblem */
int substpsolver_getNsubedges(
   const SUBSTP*         substp              /**< sub-problem */
   )
{
   assert(substp);
   assert(substp->nsubedges >= 2);

   return substp->nsubedges;
}


/** gets solution to sub-problem */
SCIP_RETCODE substpsolver_getSolution(
   SUBSTP*               substp,             /**< sub-problem */
   int*                  edgesol             /**< array to store edge solution; CONNECT/UNKNOWN */
)
{
   assert(substp && edgesol);
   assert(!substp->subgraph);

   if( !substp->useDP  )
   {
      assert(substp->subscip);

      SCIP_CALL( subscipGetSol(substp, edgesol) );
   }
   else
   {
      assert(substp->dpsubsol);
      BMScopyMemoryArray(edgesol, substp->dpsubsol, substp->nsubedges);
   }

   return SCIP_OKAY;
}


/** sets to no output */
SCIP_RETCODE substpsolver_setMute(
   SUBSTP*               substp              /**< sub-problem */
)
{
   assert(substp);
   substp->useOutput = FALSE;

   if( substp->subscip )
   {
      SCIP_CALL( SCIPsetIntParam(substp->subscip, "display/verblevel", 0) );
   }

   return SCIP_OKAY;
}


/** sets to independent problem */
SCIP_RETCODE substpsolver_setProbIsIndependent(
   SUBSTP*               substp              /**< sub-problem */
)
{
   assert(substp);
   substp->subprobIsIndependent = TRUE;
   return SCIP_OKAY;
}


/** disables DP for solving sub-problems */
SCIP_RETCODE substpsolver_setProbNoSubDP(
   SUBSTP*               substp              /**< sub-problem */
)
{
   assert(substp);
   if( substp->subscip )
   {
      SCIP_CALL( SCIPsetIntParam(substp->subscip, "stp/usedp", 0) );
   }

   return SCIP_OKAY;
}


/** sets full presolving */
SCIP_RETCODE substpsolver_setProbFullPresolve(
   SUBSTP*               substp              /**< sub-problem */
)
{
   assert(substp);
   if( substp->subscip )
   {
      SCIP_CALL( SCIPsetIntParam(substp->subscip, "stp/reduction", 2) );
   }

   return SCIP_OKAY;
}


/** Gets optimal objective by B&C.
 * NOTE: Meant for debugging only!
 * NOTE: Need to be careful with recursion! Might run into infinite loop. */
SCIP_RETCODE substpsolver_getObjFromGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   SCIP_Real*            obj                 /**< to store solution */
   )
{
   int* soledges;
   GRAPH* graph_copied;
   SUBSTP* substp;
   SCIP_Bool success;

   assert(graph);

   SCIP_CALL( SCIPallocMemoryArray(scip, &soledges, graph->edges) );

   SCIP_CALL( graph_copy(scip, graph, &graph_copied) );
   graph_copied->is_packed = FALSE;

   SCIP_CALL( substpsolver_initBC(scip, graph_copied, &substp) );
   SCIP_CALL( substpsolver_initHistory(substp) );

   SCIP_CALL( substpsolver_setMute(substp) );
   SCIP_CALL( substpsolver_setProbNoSubDP(substp) );
   SCIP_CALL( substpsolver_solve(scip, substp, &success) );

   assert(success);

   SCIP_CALL( substpsolver_getSolution(substp, soledges) );

   substpsolver_free(scip, &substp);

   *obj = solstp_getObj(graph, soledges, 0.0);

   SCIPfreeMemoryArray(scip, &soledges);

   return SCIP_OKAY;
}
