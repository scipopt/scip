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
#include "graph.h"
#include "heur_tm.h"
#include "heur_local.h"
#include "branch_stp.h"
#include "prop_stp.h"
#include "probdata_stp.h"
#include "scip/pub_tree.h"

#define BRANCHRULE_NAME            "stp"
#define BRANCHRULE_DESC            "stp branching on vertices"
#define BRANCHRULE_PRIORITY        10000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0
#define BRANCHRULE_TMRUNS          20

#define BRANCH_STP_ON_LP   0
#define BRANCH_STP_ON_LP2  1
#define BRANCH_STP_ON_SOL  2
#define BRANCH_STP_ON_DEGREE  3


/*
 * Data structures
 */



/** branching rule data */
struct SCIP_BranchruleData
{
   int                   lastcand;           /**< last evaluated candidate of last branching rule execution */
   int                   branchtype;         /**< type of branching */
   SCIP_Bool             active;             /**< is branch-rule being used? */
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
   return (probtype == STP_SPG || probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_RPCSPG || probtype == STP_PCSPG);
}

/** select vertex to branch on by choosing vertex of highest degree */
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
   SCIP_Bool ptermselected;

   assert(g != NULL && scip != NULL && vertex != NULL);
   assert(!graph_pc_isPcMw(g) || g->extended);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );

   *vertex = UNKNOWN;

   SCIPStpBranchruleInitNodeState(g, nodestate);

   /* get vertex changes from branch-and-bound */
   if( SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) != 0 )
      SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, nodestate, NULL) );

   ptermselected = FALSE;
   for( int k = 0; k < nnodes; k++ )
   {
      if( nodestate[k] == BRANCH_STP_VERTEX_NONTERM )
      {
         assert(!Is_term(g->term[k]));
         assert(!ptermselected || graph_pc_isPcMw(g));

         /* first pterm? Then update */
         if( !ptermselected && Is_pterm(g->term[k]) )
         {
            assert(graph_pc_isPcMw(g));
            maxdegree = g->grad[k];
            *vertex = k;
            ptermselected = TRUE;
         }
         else if( g->grad[k] > maxdegree && (!ptermselected || Is_pterm(g->term[k])) )
         {
            maxdegree = g->grad[k];
            *vertex = k;
         }
      }
   }

   assert(maxdegree > 0);

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
   GRAPH* graph;
   SCIP_Real* cost;
   SCIP_Real* costorg;
   SCIP_Real* costrev;
   SCIP_Real* prizeorg;
   int* soledges;
   int* nodestatenew;
   int* termorg;
   SCIP_Bool pcmw;
   SCIP_Bool ptermselected;
   SCIP_Bool success;
   int nedges;
   int nnodes;
   int maxdeg;

   assert(vertex != NULL);
   *vertex = UNKNOWN;

   /* check whether LP solution is available */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   graph = SCIPprobdataGetGraph2(scip);
   assert(graph != NULL);
   assert(!graph_pc_isPcMw(graph) || graph->extended);

   pcmw = graph_pc_isPcMw(graph);
   nedges = graph->edges;
   nnodes = graph->knots;

   /*
    * compute locally feasible solution (SPH + local)
    */

   SCIP_CALL( SCIPallocBufferArray(scip, &costorg, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodestatenew, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termorg, nnodes) );

   if( pcmw )
      SCIP_CALL( SCIPallocBufferArray(scip, &prizeorg, nnodes) );
   else
      prizeorg = NULL;

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) )
         nodestatenew[k] = BRANCH_STP_VERTEX_TERM;
      else
         nodestatenew[k] = BRANCH_STP_VERTEX_NONTERM;
   }

   BMScopyMemoryArray(costorg, graph->cost, nedges);
   BMScopyMemoryArray(termorg, graph->term, nnodes);

   if( prizeorg != NULL )
   {
      assert(graph->prize != NULL);
      BMScopyMemoryArray(prizeorg, graph->prize, nnodes);
   }

   /* get vertex changes from branch-and-bound */
   if( SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) != 0 )
      SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, nodestatenew, NULL) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) )
         continue;

      if( nodestatenew[k] == BRANCH_STP_VERTEX_TERM )
      {
         assert(graph->grad[k] > 0);

         if( pcmw && Is_pterm(graph->term[k]) )
         {
            graph_pc_enforcePterm(graph, k);
         }
         else
         {
            if( graph_pc_isRootedPcMw(graph) )
            {
               assert(graph->term2edge[k] == -1);
               graph->prize[k] = FARAWAY;
            }
            else if( pcmw )
            {
               continue;
            }

            graph_knot_chg(graph, k, 0);
         }
      }
      else if( nodestatenew[k] == BRANCH_STP_VERTEX_KILLED )
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            if( pcmw && graph->term2edge[k] == e ) /* do not change edge going to pseudo-terminal */
            {
               assert(Is_pterm(graph->term[k]));
               assert(Is_term(graph->term[graph->head[e]]));
               continue;
            }

            graph->cost[e] += BLOCKED;
            graph->cost[flipedge(e)] += BLOCKED;
         }
      }
   }

   SCIP_CALL( SCIPStpHeurTMRunLP(scip, graph, NULL, soledges, BRANCHRULE_TMRUNS, cost, costrev, &success) );
   assert(success);
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, soledges) );

   assert(graph_sol_valid(scip, graph, soledges));

   for( int k = 0; k < nnodes; k++ )
      if( graph->term[k] != termorg[k] )
         graph_knot_chg(graph, k, termorg[k]);

   BMScopyMemoryArray(graph->cost, costorg, nedges);

   if( prizeorg != NULL )
   {
      assert(graph->prize != NULL);
      BMScopyMemoryArray(graph->prize, prizeorg, nnodes);
   }

   assert(graph_sol_valid(scip, graph, soledges));

   if( addsol )
   {
      SCIP_SOL* sol = NULL;

      const int nvars = SCIPprobdataGetNVars(scip);

      assert(graph->edges == nvars);

      /* use cost array to store solution */
      for( int e = 0; e < nvars; e++ )
         if( soledges[e] == CONNECT )
            cost[e] = 1.0;
         else
            cost[e] = 0.0;

      SCIP_CALL( SCIPprobdataAddNewSol(scip, cost, sol, NULL, &success) );
      SCIPdebugMessage("solution added? %d \n", success);
   }

   /* get vertex with highest degree in solution */
   maxdeg = -1;
   ptermselected = FALSE;
   for( int i = 0; i < nnodes; i++ )
   {
      if( nodestatenew[i] == BRANCH_STP_VERTEX_NONTERM && graph->grad[i] != 0 )
      {
         int soldeg = 0;
         assert(!Is_term(graph->term[i]));
         assert(!ptermselected || pcmw);

         for( int e = graph->outbeg[i]; e != EAT_LAST; e = graph->oeat[e] )
            if( soledges[e] == CONNECT || soledges[flipedge(e)] == CONNECT )
               soldeg++;

         /* first pterm? Then update */
         if( !ptermselected && Is_pterm(graph->term[i]) )
         {
            assert(pcmw);
            maxdeg = soldeg;
            *vertex = i;
            ptermselected = TRUE;
         }
         else if( soldeg > maxdeg && (!ptermselected || Is_pterm(graph->term[i])) )
         {
            maxdeg = soldeg;
            *vertex = i;
         }
      }
   }

   assert(maxdeg >= 0);

   SCIPfreeBufferArrayNull(scip, &prizeorg);
   SCIPfreeBufferArray(scip, &termorg);
   SCIPfreeBufferArray(scip, &nodestatenew);
   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &costorg);

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
   SCIP_Real bestflow;
   const int nnodes = g->knots;

   assert(g != NULL && vertex != NULL);

   *vertex = UNKNOWN;

   /* LP has not been solved */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   edgevars = SCIPprobdataGetEdgeVars(scip);
   assert(edgevars != NULL);

   bestflow = 1.0;
   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_Real inflow = 0.0;
      for( int a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
         inflow += SCIPvarGetLPSol(edgevars[a]);

      if( !Is_term(g->term[k]) && SCIPisLT(scip, inflow, 1.0) && fabs(inflow - 0.5) < bestflow && SCIPisPositive(scip, inflow) )
      {
         *vertex = k;
         bestflow = fabs(inflow - 0.5);
         SCIPdebugMessage("new maxflow %f on vertex %d \n", inflow, k );
      }
   }

   SCIPdebugMessage("maxflow %f on vertex %d, term? %d \n", bestflow, *vertex, Is_term(g->term[*vertex])  );

   return SCIP_OKAY;
}

/** select vertices with flow on at least two incoming arcs */
static
SCIP_RETCODE selectBranchingVertexByLp2Flow(
   SCIP*                 scip,               /**< original SCIP data structure */
   int*                  vertex,             /**< the vertex to branch on */
   const GRAPH*          g                   /**< graph */
   )
{
   SCIP_VAR** edgevars;
   int* nodestate;
   SCIP_Real bestflow;
   const int nnodes = g->knots;

   assert(g != NULL && vertex != NULL);

   *vertex = UNKNOWN;

   /* has LP not been solved? */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
         nodestate[k] = BRANCH_STP_VERTEX_TERM;
      else
         nodestate[k] = BRANCH_STP_VERTEX_NONTERM;
   }

   /* get vertex changes from branch-and-bound */
   if( SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) != 0 )
      SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, nodestate, NULL) );

   edgevars = SCIPprobdataGetEdgeVars(scip);
   assert(edgevars != NULL);

   bestflow = 1.0;
   for( int k = 0; k < nnodes; k++ )
   {
      if( nodestate[k] == BRANCH_STP_VERTEX_NONTERM )
      {
         int ninarcs = 0;
         assert(g->grad[k] > 0);
         assert(!Is_term(g->term[k]));

         for( int a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
            if( SCIPisPositive(scip, SCIPvarGetLPSol(edgevars[a])) )
               ninarcs++;

         if( ninarcs > 1 )
         {
            SCIP_Real inflow = 0.0;
            for( int a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
               inflow += SCIPvarGetLPSol(edgevars[a]);

            if( fabs(inflow - 0.5) < bestflow )
            {
               *vertex = k;
               bestflow = fabs(inflow - 0.5);
               SCIPdebugMessage("new maxflow %f on vertex %d \n", inflow, k );
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &nodestate);

   SCIPdebugMessage("LP2: branch on vertex %d with abs(flow - 0.5): %f \n", *vertex, bestflow);

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
   int branchruletype;
   int branchvertex = UNKNOWN;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of STP branching\n ");
   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( !branchruledata->active )
      return SCIP_OKAY;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   if( !isProbCompatible(g->stp_type) )
      return SCIP_OKAY;

   branchruletype = branchruledata->branchtype;

   /* get vertex to branch on */
   if( graph_pc_isPcMw(g) && (branchruletype == BRANCH_STP_ON_LP || branchruletype == BRANCH_STP_ON_LP2) )
   {
      SCIP_CALL( selectBranchingVertexBySol(scip, &branchvertex, TRUE) );
   }

   if( branchruletype == BRANCH_STP_ON_LP )
      SCIP_CALL( selectBranchingVertexByLp(scip, &branchvertex, g) );
   else if( branchruletype == BRANCH_STP_ON_LP2 )
      SCIP_CALL( selectBranchingVertexByLp2Flow(scip, &branchvertex, g) );
   else if( branchruletype == BRANCH_STP_ON_SOL )
      SCIP_CALL( selectBranchingVertexBySol(scip, &branchvertex, TRUE) );

   /* fall-back strategy */
   if( branchvertex == UNKNOWN )
      SCIP_CALL( selectBranchingVertexByDegree(scip, &branchvertex, g) );

   /* we should only have terminals in this case */
   if( branchvertex == UNKNOWN )
   {
      printf("Stp branching could not run \n");
      return SCIP_ERROR;
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
   SCIP_BRANCHRULEDATA* branchruledata;
   GRAPH* g;
   int branchvertex;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPdebugMessage("Execps method of STP branching\n");
   *result = SCIP_DIDNOTRUN;

   if( !branchruledata->active )
      return SCIP_OKAY;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   g = SCIPprobdataGetGraph(probdata);
   assert(g != NULL);

   if( !isProbCompatible(g->stp_type) )
      return SCIP_OKAY;

   branchvertex = UNKNOWN;

   if( branchruledata->branchtype != BRANCH_STP_ON_DEGREE )
      SCIP_CALL( selectBranchingVertexBySol(scip, &branchvertex, TRUE) );

   /* fall-back strategy */
   if( branchvertex == UNKNOWN )
      SCIP_CALL( selectBranchingVertexByDegree(scip, &branchvertex, g) );

   if( branchvertex == UNKNOWN )
   {
      SCIPdebugMessage("Stp branching did not run \n");
      return SCIP_OKAY;
   }

   branchruledata->active = TRUE;
   SCIP_CALL( branchOnVertex(scip, g, branchvertex) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** parse constraint name and apply changes to graph or array */
SCIP_RETCODE STPStpBranchruleParseConsname(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vertexchgs,         /**< array to store changes or NULL */
   GRAPH*                graph,              /**< graph to modify or NULL */
   const char*           consname,           /**< constraint name */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
)
{

   /* terminal inclusion constraint? */
   if( strncmp(consname, "consin", 6) == 0 )
   {
      char* tailptr;
      const int term = (int) strtol(consname + 6, &tailptr, 10);

      SCIPdebugMessage("make terminal %d \n", term);

      if( graph != NULL && !Is_term(graph->term[term]) )
      {
         if( graph_pc_isPcMw(graph) )
         {
            if( Is_pterm(graph->term[term]) )
            {
               graph_pc_enforcePterm(graph, term);
            }
            else if( graph_pc_isRootedPcMw(graph) )
            {
               graph->prize[term] = FARAWAY;
               graph_knot_chg(graph, term, 0);
            }
         }
         else
         {
            graph_knot_chg(graph, term, 0);
         }
      }

      if( vertexchgs != NULL )
         vertexchgs[term] = BRANCH_STP_VERTEX_TERM;
   }
   /* node removal constraint? */
   else if( strncmp(consname, "consout", 7) == 0 )
   {
      char* tailptr;
      const int vert = (int) strtol(consname + 7, &tailptr, 10);

      SCIPdebugMessage("delete vertex %d \n", vert);
      assert(graph == NULL || (vert >= 0 && vert < graph->knots));

      if( graph != NULL && graph->grad[vert] > 0 )
      {
         assert(!Is_term(graph->term[vert]));

         if( Is_pterm(graph->term[vert]) )
         {
            assert(graph_pc_isPcMw(graph));
            graph_pc_deleteTerm(scip, graph, graph_pc_getTwinTerm(graph, vert));
         }
         else
         {
            graph_knot_del(scip, graph, vert, freeancestors);
         }
      }
      if( vertexchgs != NULL )
         vertexchgs[vert] = BRANCH_STP_VERTEX_KILLED;
   }
   else
   {
      printf("found unknown constraint at b&b node \n");
      return SCIP_ERROR;
   }
   return SCIP_OKAY;
}

/** applies vertex changes caused by this branching rule, either on a graph or on an array */
void SCIPStpBranchruleInitNodeState(
   const GRAPH*          g,                  /**< graph data structure */
   int*                  nodestate           /**< node state array */
   )
{
   const int nnodes = g->knots;

   assert(g != NULL);
   assert(nodestate != NULL);

   assert(nnodes > 0);

   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         nodestate[k] = BRANCH_STP_VERTEX_KILLED;
      else if( Is_term(g->term[k]) )
         nodestate[k] = BRANCH_STP_VERTEX_TERM;
      else
         nodestate[k] = BRANCH_STP_VERTEX_NONTERM;
   }
}


/** applies vertex changes caused by this branching rule, either on a graph or on an array */
SCIP_RETCODE SCIPStpBranchruleApplyVertexChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vertexchgs,         /**< array to store changes or NULL */
   GRAPH*                graph               /**< graph to apply changes on or NULL */
   )
{
   SCIP_CONS* parentcons;
   SCIP_BRANCHRULE* branchrule = NULL;
   SCIP_BRANCHRULEDATA* branchruledata;
   int naddedconss;

   assert(scip != NULL);
   assert(graph != NULL || vertexchgs != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( !(branchruledata->active) )
   {
      printf("STP branchrule not active! \n");
      return SCIP_OKAY;
   }

   assert(SCIPnodeGetNAddedConss(SCIPgetCurrentNode(scip)) == 1);

   /* not very clean, but this should only happen when the whole b&b tree is explored */
   if( SCIPnodeGetNAddedConss(SCIPgetCurrentNode(scip)) != 1 )
      return SCIP_OKAY;

   /* move up branch-and-bound path and check constraints */
   for( SCIP_NODE* node = SCIPgetCurrentNode(scip); SCIPnodeGetDepth(node) > 0; node = SCIPnodeGetParent(node) )
   {
      char* consname;

      assert(SCIPnodeGetNAddedConss(node) == 1);

      if( SCIPnodeGetNAddedConss(node) != 1 )
         break;

      /* get constraints */
      SCIPnodeGetAddedConss(node, &parentcons, &naddedconss, 1);
      consname = (char*) SCIPconsGetName(parentcons);

      SCIP_CALL( STPStpBranchruleParseConsname(scip, vertexchgs, graph, consname, TRUE) );
   }

   return SCIP_OKAY;
}

/** creates the stp branching rule and includes it in SCIP */
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

   SCIP_CALL( SCIPaddIntParam(scip, "branching/stp/type",
         "Branching: 0 based on LP, 1 based on LP and with indegree > 1, 2 based on best solution,"
         "3 based on degree",
         &(branchruledata->branchtype), FALSE, BRANCH_STP_ON_LP, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/stp/vertexbranching",
         "use branching on vertices?",
         &(branchruledata->active), FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
