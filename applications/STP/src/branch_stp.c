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
#include "scip/cons_setppc.h"
#include "scip/var.h"
#include "scip/set.h"
#include "scip/pub_tree.h"
#include "scip/struct_scip.h"
#include "scip/clock.h"
#include "grph.h"
#include "heur_tm.h"
#include "heur_local.h"
#include "branch_stp.h"
#include "prop_stp.h"
#include "probdata_stp.h"
#include "scip/tree.h"

#define BRANCHRULE_NAME            "stp"
#define BRANCHRULE_DESC            "stp branching on vertices"
#define BRANCHRULE_PRIORITY        10000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0
#define BRANCHRULE_TMRUNS          20

#define DEFAULT_USELOCALSOL      FALSE

#define BRANCH_STP_VERTEX_KILLED      -1
#define BRANCH_STP_VERTEX_NONTERM      0
#define BRANCH_STP_VERTEX_TERM         1

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   SCIP_Bool             uselocalsol;        /**< branch on vertex of highest degree in best local solution? */
};


/*
 * Local methods
 */


/** check whether branching-rule is compatible with given problem type */
static
SCIP_Bool isProbCompatible(
   int                   probtype            /**< the problem type */
)
{
   return (probtype == STP_SPG || probtype == STP_RSMT || probtype == STP_OARSMT);
}

/** select vertex to branch on by chosing vertex of highest degree */
static
SCIP_RETCODE selectBranchingVertexByDegree(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  vertex,             /**< the vertex to branch on */
   const GRAPH*          g                   /**< graph */
   )
{
   int* nodestate;
   int maxdegree = 0;
   const int nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );

   *vertex = UNKNOWN;
   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         nodestate[k] = BRANCH_STP_VERTEX_KILLED;
      else if( Is_term(g->term[k]) )
         nodestate[k] = BRANCH_STP_VERTEX_TERM;
      else
         nodestate[k] = BRANCH_STP_VERTEX_NONTERM;
   }

   SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, nodestate, NULL) );

   /* */
   for( int k = 0; k < nnodes; k++ )
   {
      if( nodestate[k] == BRANCH_STP_VERTEX_NONTERM && g->grad[k] > maxdegree )
      {
         assert(!Is_term(g->term[k]));

         maxdegree = g->grad[k];
         *vertex = k;
      }
   }

   printf("vertex selected by degree branching rule: %d \n\n", *vertex);

   SCIPfreeBufferArray(scip, &nodestate);

   return SCIP_OKAY;
}

/** select vertex to branch on by using primal solution */
static
SCIP_RETCODE selectBranchingVertexBySol(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  vertex,             /**< the vertex to branch on */
   SCIP_Bool             addsol              /**< add new solution to pool? */
   )
{
   GRAPH* propgraph;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   int* soledges;
   SCIP_Longint graphnodenumber;
   SCIP_Bool success;
   int nedges;
   int nnodes;
   int maxdeg;

   *vertex = UNKNOWN;

   /* check whether LP solution is available */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get local graph from propagator */
   SCIPStpPropGetGraph(scip, &propgraph, &graphnodenumber);

   if( propgraph == NULL || graphnodenumber != SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      return SCIP_OKAY;

   nedges = propgraph->edges;
   nnodes = propgraph->knots;

   /*
    * compute locally feasible solution (SPH + local)
    */

   SCIP_CALL(SCIPallocBufferArray(scip, &cost, nedges));
   SCIP_CALL(SCIPallocBufferArray(scip, &costrev, nedges));
   SCIP_CALL(SCIPallocBufferArray(scip, &soledges, nedges));

   SCIP_CALL( SCIPStpHeurTMRunLP(scip, propgraph, NULL, soledges, BRANCHRULE_TMRUNS, cost, costrev, &success) );
   assert(success);

   SCIP_CALL( SCIPStpHeurLocalRun(scip, propgraph, propgraph->cost, soledges) );

   if( addsol )
   {
      SCIP_SOL* sol = NULL;

      /* use cost array to store solution */
      for( int e = 0; e < nedges; e++ )
         if( soledges[e] == CONNECT )
            cost[e] = 1.0;
         else
            cost[e] = 0.0;

      SCIP_CALL( SCIPprobdataAddNewSol(scip, cost, sol, NULL, &success) );
   }

   /* get vertex with highest degree in solution */
   maxdeg = -1;
   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(propgraph->term[i]) && propgraph->grad[i] != 0 )
      {
         int soldeg = 0;
         for( int e = propgraph->outbeg[i]; e != EAT_LAST; e = propgraph->oeat[e] )
            if( soledges[e] == CONNECT || soledges[flipedge(e)] == CONNECT )
               soldeg++;
         if( soldeg > maxdeg )
         {
            maxdeg = soldeg;
            *vertex = i;
         }
      }
   }

   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}

/** select vertex to branch on by using LP */
static
SCIP_RETCODE selectBranchingVertexByLp(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  vertex,             /**< the vertex to branch on */
   const GRAPH*          g                   /**< graph */
   )
{
   SCIP_VAR** edgevars;
   SCIP_Real maxflow;
   SCIP_Real* inflow;
   int a;
   int k;
   int branchvert;
   const int nnodes = g->knots;

   /* LP has not been solved */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      *vertex = UNKNOWN;
      return SCIP_OKAY;
   }

   edgevars = SCIPprobdataGetEdgeVars(scip);
   assert(edgevars != NULL);

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

/** branch on specified (Steiner) vertex */
static
SCIP_RETCODE branchOnVertex(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   branchvertex        /**< the vertex to branch on */
   )
{
   char consnamein[SCIP_MAXSTRLEN];
   char consnameout[SCIP_MAXSTRLEN];
   SCIP_CONS* consin = NULL;
   SCIP_CONS* consout = NULL;
   SCIP_NODE* vertexin = NULL;
   SCIP_NODE* vertexout = NULL;
   SCIP_VAR** const edgevars = SCIPprobdataGetEdgeVars(scip);
   const SCIP_Real estimatein = SCIPgetUpperbound(scip);
   const SCIP_Real estimateout = SCIPgetUpperbound(scip);

   assert(scip != NULL);
   assert(g != NULL);
   assert(branchvertex >= 0 && branchvertex < g->knots);
   assert(!Is_term(g->term[branchvertex]));

   (void) SCIPsnprintf(consnamein, SCIP_MAXSTRLEN, "consin%d", branchvertex);
   (void) SCIPsnprintf(consnameout, SCIP_MAXSTRLEN, "consout%d", branchvertex);

   /* create constraints */
   SCIP_CALL( SCIPcreateConsSetpart(scip, &consin,
         consnamein, 0, NULL, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcreateConsLinear(scip, &consout,
         consnameout, 0, NULL, NULL, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   for( int e = g->inpbeg[branchvertex]; e != EAT_LAST; e = g->ieat[e] )
   {
      SCIP_CALL(SCIPaddCoefSetppc(scip, consin, edgevars[e]));
      SCIP_CALL(SCIPaddCoefLinear(scip, consout, edgevars[e], 1.0));
      SCIP_CALL(SCIPaddCoefLinear(scip, consout, edgevars[flipedge(e)], 1.0));
   }

   /* create the child nodes */
   SCIP_CALL(SCIPcreateChild(scip, &vertexin, 1.0, estimatein));
   SCIP_CALL(SCIPcreateChild(scip, &vertexout, 1.0, estimateout));

   assert(vertexin != NULL);
   assert(vertexout != NULL);

   SCIP_CALL(SCIPaddConsNode(scip, vertexin, consin, NULL));
   SCIP_CALL(SCIPaddConsNode(scip, vertexout, consout, NULL));

   /* release constraints */
   SCIP_CALL(SCIPreleaseCons(scip, &consin));
   SCIP_CALL(SCIPreleaseCons(scip, &consout));

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
   SCIP_BRANCHRULEDATA* branchruledata;
   GRAPH* g;
   int branchvertex;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of STP branching\n ");
   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   if( !isProbCompatible(g->stp_type) )
      return SCIP_OKAY;

   /* get vertex to branch on */
   if( branchruledata->uselocalsol )
      SCIP_CALL( selectBranchingVertexBySol(scip, &branchvertex, TRUE) );
   else
      SCIP_CALL( selectBranchingVertexByLp(scip, &branchvertex, g) );

   if( branchvertex == UNKNOWN )
      SCIP_CALL( selectBranchingVertexByDegree(scip, &branchvertex, g) );

   if( branchvertex == UNKNOWN )
   {
      SCIPdebugMessage("Stp branching did not run \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( branchOnVertex(scip, g, branchvertex) );

   SCIPdebugMessage("Branched on stp vertex %d \n", branchvertex);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsStp)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   GRAPH* g;
   int branchvertex;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execps method of STP branching\n");
   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   if( !isProbCompatible(g->stp_type) )
      return SCIP_OKAY;

   /* select vertex to branch on */
   SCIP_CALL( selectBranchingVertexByDegree(scip, &branchvertex, g) );

   if( branchvertex == UNKNOWN )
   {
      SCIPdebugMessage("Stp branching did not run \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( branchOnVertex(scip, g, branchvertex) );

   *result = SCIP_BRANCHED;

   printf("branchExecpsStp %d \n", 0);

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** applies vertex changes caused by this branching rule, either on a graph or on an array */
SCIP_RETCODE SCIPStpBranchruleApplyVertexChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vertexchgs,         /**< array to store changes or NULL */
   GRAPH*                graph               /**< graph to apply changes on or NULL */
   )
{
   SCIP_CONS* parentcons;
   int naddedconss;

   assert(scip != NULL);
   assert(graph != NULL || vertexchgs != NULL);

   printf("added %d \n", SCIPnodeGetNAddedConss(SCIPgetCurrentNode(scip)));
   if( SCIPnodeGetNAddedConss(SCIPgetCurrentNode(scip)) != 1 )
      exit(1);

   assert(SCIPnodeGetNAddedConss(SCIPgetCurrentNode(scip)) == 1);

   /* move up branch-and-bound path and check constraints */
   for( SCIP_NODE* node = SCIPgetCurrentNode(scip); SCIPnodeGetDepth(node) > 0; node = SCIPnodeGetParent(node) )
   {
      char* consname;

      if( SCIPnodeGetNAddedConss(node) != 1 )
         continue;

      /* get constraints */
      SCIPnodeGetAddedConss(node, &parentcons, &naddedconss, 1);
      consname = (char*) SCIPconsGetName(parentcons);

      /* terminal inclusion constraint? */
      if( strncmp(consname, "consin", 6) == 0 )
      {
         char* tailptr;
         const int term = (int) strtol(consname + 6, &tailptr, 10);

         SCIPdebugMessage("make terminal %d \n", term);

         if( graph != NULL)
            graph_knot_chg(graph, term, 0);

         if( vertexchgs != NULL )
            vertexchgs[term] = BRANCH_STP_VERTEX_TERM;
      }
      /* node removal constraint? */
      else if( strncmp(consname, "consout", 7) == 0 )
      {
         char* tailptr;
         const int vert = (int) strtol(consname + 7, &tailptr, 10);

         SCIPdebugMessage("delete vertex %d \n", vert);

         if( graph != NULL)
            graph_knot_del(scip, graph, vert, TRUE);

         if( vertexchgs != NULL )
            vertexchgs[vert] = BRANCH_STP_VERTEX_KILLED;
      }
      else
      {
         printf("found unknown constraint at b&b node \n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

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
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsStp) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/stp/uselocalsol",
         "branch on vertex of highest degree in best local solution?",
         &branchruledata->uselocalsol, FALSE, DEFAULT_USELOCALSOL, NULL, NULL) );

   return SCIP_OKAY;
}
