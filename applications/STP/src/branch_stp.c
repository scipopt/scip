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
/**@file   branch_stp.c
 * @brief  Steiner vertex branching rule
 * @author Daniel Rehfeldt
 *
 * The Steiner branching rule implemented in this file is described in
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 * It includes and excludes Steiner vertices during branching.
 *
*/
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/branch_fullstrong.h"
#include "scip/cons_linear.h"
#include "scip/var.h"
#include "scip/set.h"
#include "scip/pub_tree.h"
#include "scip/struct_scip.h"
#include "scip/clock.h"
#include "grph.h"
#include "branch_stp.h"
#include "probdata_stp.h"

#define BRANCHRULE_NAME            "stp"
#define BRANCHRULE_DESC            "stp branching on vertices"
#define BRANCHRULE_PRIORITY        1000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
};


/*
 * Local methods
 */

static
SCIP_RETCODE selectBranchingVertex(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  vertex              /**< the vertex to branch on */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VAR** edgevars;
   GRAPH* g;
   SCIP_Real maxflow;
   SCIP_Real* inflow;
   int a;
   int k;
   int nnodes;
   int branchvert;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   /* LP has not been solved */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      *vertex = UNKNOWN;
      return SCIP_OKAY;
   }

   edgevars = SCIPprobdataGetEdgeVars(scip);
   assert(edgevars != NULL);

   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &inflow, nnodes) );

   branchvert = UNKNOWN;
   maxflow = 1.0;
   for( k = 0; k < nnodes; k++ )
   {
      inflow[k] = 0.0;
      for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
         inflow[k] += SCIPvarGetLPSol(edgevars[a]);

      if( !Is_term(g->term[k]) && SCIPisLT(scip, inflow[k], 1.0) && SCIPisLT(scip, fabs(inflow[k] - 0.5), maxflow) )
      {
         branchvert = k;
         maxflow = fabs(inflow[k] - 0.5);
         SCIPdebugMessage("new maxflow %f on vertex %d \n", inflow[k], branchvert );
      }
   }
   SCIPdebugMessage("maxflow %f on vertex %d, term? %d \n", maxflow, branchvert, Is_term(g->term[branchvert])  );
   (*vertex) = branchvert;

   SCIPfreeBufferArray(scip, &inflow);

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleStp(scip) ) ;

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeStp)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitStp)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->lastcand = 0;

   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitStp)
{  /*lint --e{715}*/
#if 0
   SCIP_BRANCHRULEDATA* branchruledata;
#endif
   SCIPstatistic(int j = 0);

#if 0
   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
#endif
   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpStp)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_CONS* consin;
   SCIP_CONS* consout;
   SCIP_NODE* vertexin;
   SCIP_NODE* vertexout;
   SCIP_VAR** edgevars;
   SCIP_Real estimatein;
   SCIP_Real estimateout;
   GRAPH* g;
   int e;
   int probtype;
   int branchvertex;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of Stp branching\n ");
   estimatein = SCIPgetUpperbound(scip);
   estimateout = SCIPgetUpperbound(scip);
   *result = SCIP_DIDNOTRUN;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   probtype = g->stp_type;

   if( probtype != STP_SPG && probtype != STP_RSMT && probtype != STP_OARSMT )
      return SCIP_OKAY;

   /* get vertex to branch on */
   SCIP_CALL( selectBranchingVertex(scip, &branchvertex) );

   if( branchvertex == UNKNOWN )
   {
      SCIPdebugMessage("Branching did not run \n");
      return SCIP_OKAY;
   }

   edgevars = SCIPprobdataGetEdgeVars(scip);

   /* create constraints */
   SCIP_CALL( SCIPcreateConsLinear(scip, &consin, "consin", 0,
         NULL, NULL, 1.0, 1.0,
         TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcreateConsLinear(scip, &consout, "consout", 0,
         NULL, NULL, 0.0, 0.0,
         TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   for( e = g->inpbeg[branchvertex]; e != EAT_LAST; e = g->ieat[e] )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, consin,  edgevars[e], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consout, edgevars[e], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consout, edgevars[flipedge(e)], 1.0) );
   }

   /* create the child nodes */
   SCIP_CALL( SCIPcreateChild(scip, &vertexin, 1.0, estimatein) );

   SCIP_CALL( SCIPcreateChild(scip, &vertexout, 1.0, estimateout) );

   assert(vertexin != NULL);
   assert(vertexout != NULL);

   SCIP_CALL( SCIPaddConsNode(scip, vertexin, consin, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, vertexout, consout, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consin) );
   SCIP_CALL( SCIPreleaseCons(scip, &consout) );

   SCIPdebugMessage("Branched on stp vertex %d \n", branchvertex);

   *result = SCIP_BRANCHED;


   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the multi-aggregated branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create stp branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->lastcand = 0;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyStp) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeStp) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitStp) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitStp) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpStp) );

   return SCIP_OKAY;
}
