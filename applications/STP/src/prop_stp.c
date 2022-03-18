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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_stp.c
 * @brief  propagator for Steiner tree problems, using the LP reduced costs
 * @author Daniel Rehfeldt
 *
 * This propagator makes use of the reduced cost of an optimally solved LP relaxation to propagate the variables
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "prop_stp.h"
#include "graph.h"
#include "reduce.h"
#include "solstp.h"
#include "redcosts.h"
#include "extreduce.h"
#include "cons_stp.h"
#include "branch_stp.h"
#include "portab.h"
#include "scip/tree.h"


/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "stp"
#define PROP_DESC              "stp propagator"
#define PROP_TIMING            SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           1000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_STP_EDGE_KILLED -1
#define PROP_STP_EDGE_UNSET   0
#define PROP_STP_EDGE_SET     1
#define PROP_STP_EDGE_FIXED   2

#define PROP_STP_REDCOST_LEVELS 5

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_MAXNWAITINGROUNDS        2     /**< maximum number of rounds to wait until propagating again */
#define REDUCTION_WAIT_RATIO_INITIAL             0.01  /**< ratio of edges to be newly fixed before performing reductions for additional fixing */
#define REDUCTION_WAIT_FACTOR                     2.0
/**@} */

/*
 * Data structures
 */


/** propagator data */
struct SCIP_PropData
{
   REDCOST*              redcostdata;        /**< reduced cost data */
   GRAPH*                propgraph;          /**< graph data */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                                           */
   SCIP_Real*            fixingbounds;       /**< saves largest upper bound to each variable that would allow to fix it */
   SCIP_Real*            deg2bounds;         /**< saves largest upper bound to each variable that would allow to set degree 2 constraint */
   SCIP_Bool*            deg2bounded;        /**< maximum degree of vertex is 2? */
   SCIP_Longint          nfails;             /**< number of failures since last successful call */
   SCIP_Longint          ncalls;             /**< number of calls */
   SCIP_Longint          nlastcall;          /**< number of last call */
   SCIP_Longint          nlastnlpiter;       /**< number of last LP iterations */
   SCIP_Longint          lastnodenumber_red; /**< number for last reduction call */
   SCIP_Longint          lastnodenumber;     /**< number nodes for last propagator call */
   SCIP_Longint          propgraphnodenumber;/**< B&B node number at which propgraph was updated */
   SCIP_Real             lpobjval;           /**< value of current LP */
   SCIP_Real             lpobjval_last;      /**< value of last LP */
   SCIP_Real             redwaitratio;       /**< ratio */
   int                   redcostnupdates;    /**< number of reduced costs updates */
   int                   nfixededges_bicurr; /**< number of arcs fixed by 'fixedgevar' method of this propagator in current run */
   int                   nfixededges_curr;   /**< number of arcs fixed by 'fixedgevar' method of this propagator in current run */
   int                   nfixededges_all;    /**< total number of arcs fixed by 'fixedgevar' method of this propagator */
   int                   nfixededges_bipost; /**< total number of arcs fixed by 'fixedgevar' method of this propagator after the last reductions */
   int                   maxnwaitrounds;     /**< maximum number of rounds to wait until propagating again */
   SCIP_Bool             isInitialized;      /**< already initialized? */
   SCIP_Bool             aggressive;         /**< be aggressive? */
};

/**@name Local methods
 *
 * @{
 */

#if 0
static
SCIP_RETCODE fixedgevarTo1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   int*                  nfixed              /**< counter that is incriminated if variable could be fixed */
   )
{
   assert(scip != NULL);

   if( SCIPvarGetLbLocal(edgevar) < 0.5 && SCIPvarGetUbLocal(edgevar) > 0.5 )
   {
      SCIP_PROPDATA* propdata;

      printf("FIXED TO ONE %d \n", 0);

      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(edgevar), 1.0));

      /* get propagator data */
      assert(SCIPfindProp(scip, "stp") != NULL);
      propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));
      assert(propdata != NULL);

      SCIP_CALL( SCIPchgVarLb(scip, edgevar, 1.0) );
      (*nfixed)++;
      propdata->nfixededges_all++;
      propdata->nfixededges_bipost++;
   }
   return SCIP_OKAY;
}
#endif



/** fix a variable (corresponding to an edge) to zero */
static
SCIP_RETCODE fixEdgeVar(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge,               /**< the edge */
   SCIP_VAR**            vars,               /**< variables  */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip && vars && propdata);

   if( SCIPvarGetLbLocal(vars[edge]) < 0.5 && SCIPvarGetUbLocal(vars[edge]) > 0.5 )
   {
      SCIP_CALL( SCIPchgVarUb(scip, vars[edge], 0.0) );
      propdata->nfixededges_curr++;

      /* opposite edge already fixed? */
      if( SCIPvarGetUbLocal(vars[flipedge(edge)]) < 0.5 )
         propdata->nfixededges_bicurr++;
   }

   return SCIP_OKAY;
}


/** gets cutoff bound */
static
SCIP_Real getCutoffbound(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_Real             lpbound            /**< LP bound */
   )
{
   SCIP_Real cutoffbound;
   const SCIP_Real cutoffbound_scip = SCIPgetCutoffbound(scip);
   SCIP_Real cutoffbound_presol = SCIPprobdataGetPresolUpperBound(scip);

   // todo why doesnt that work???
   /*
   const SCIP_Bool objIsIntegral = SCIPprobdataObjIsIntegral(scip);

   if( objIsIntegral && LT(lpbound, cutoffbound_presol - 1.0) )
   {
      cutoffbound_presol = SCIPprobdataGetPresolCutoffbound(scip);
   }
   */

   SCIPdebugMessage("Cutoffbound SCIP vs Presol: %f vs %f \n", cutoffbound_scip, cutoffbound_presol);

   cutoffbound = MIN(cutoffbound_scip, cutoffbound_presol);

   return cutoffbound;
}


/** gets bound changes specific for PC/MW */
static
void getBoundchangesPcMW(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          propgraph,          /**< graph data structure */
   int*                  nodestate,          /**< state */
   SCIP_Bool*            conflictFound       /**< conflict found? */
)
{
   const int nnodes = propgraph->knots;

   assert(graph_pc_isPcMw(propgraph));
   assert(propgraph->extended);

   *conflictFound = FALSE;

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph_pc_knotIsPropPotTerm(propgraph, k) )
      {
         const int pterm2term = propgraph->term2edge[k];
         const int twinterm = graph_pc_getTwinTerm(propgraph, k);
         const int root2term = graph_pc_getRoot2PtermEdge(propgraph, twinterm);

         assert(pterm2term >= 0);
         assert(root2term >= 0);
         assert(graph_pc_knotIsDummyTerm(propgraph, twinterm));
         assert(SCIPisEQ(scip, propgraph->prize[k], propgraph->cost[root2term]));

         /* root->dummyterm edge fixed to 0? Then take terminal */
         if( SCIPvarGetUbLocal(vars[root2term]) < 0.5 )
         {
            if( BRANCH_STP_VERTEX_KILLED == nodestate[k] )
            {
               *conflictFound = TRUE;
               return;
            }

            nodestate[k] = BRANCH_STP_VERTEX_TERM;
         }
         /* root->dummyterm edge fixed to 1? Then delete terminal */
         else if( SCIPvarGetLbLocal(vars[root2term]) > 0.5 )
         {
            if( BRANCH_STP_VERTEX_TERM == nodestate[k] )
            {
               *conflictFound = TRUE;
               return;
            }

            nodestate[k] = BRANCH_STP_VERTEX_KILLED;
         }

         /* term->dummyterm edge fixed to 0? Then delete terminal */
         if( SCIPvarGetUbLocal(vars[pterm2term]) < 0.5 )
         {
            if( BRANCH_STP_VERTEX_TERM == nodestate[k] )
            {
               *conflictFound = TRUE;
               return;
            }

            nodestate[k] = BRANCH_STP_VERTEX_KILLED;
         }
         /* term->dummyterm edge fixed to 1? Then take terminal */
         else if( SCIPvarGetLbLocal(vars[pterm2term]) > 0.5 )
         {
            if( BRANCH_STP_VERTEX_KILLED == nodestate[k] )
            {
               *conflictFound = TRUE;
               return;
            }

            nodestate[k] = BRANCH_STP_VERTEX_TERM;
         }
      }
   }
}


/** gets bound changes */
static
SCIP_RETCODE getBoundchanges(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          propgraph,          /**< graph data structure */
   int*                  nodestate,          /**< node state (uninitialized) */
   int*                  edgestate,          /**< edge state (uninitialized) */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
)
{
   const int nnodes = graph_get_nNodes(propgraph);
   const int nedges = graph_get_nEdges(propgraph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(propgraph);

   assert(!pcmw || propgraph->extended);
   assert(!(*probisinfeas));

   for( int k = 0; k < nnodes; k++ )
      nodestate[k] = BRANCH_STP_VERTEX_UNSET;

   for( int e = 0; e < nedges; e++ )
      edgestate[e] = PROP_STP_EDGE_UNSET;

   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      if( pcmw && graph_pc_edgeIsExtended(propgraph, e) )
         continue;

      /* e OR its anti-parallel edge fixed to 1? */
      if( SCIPvarGetLbLocal(vars[e]) > 0.5 || SCIPvarGetLbLocal(vars[erev]) > 0.5 )
      {
         const int tail = propgraph->tail[e];
         const int head = propgraph->head[e];

         edgestate[e] = PROP_STP_EDGE_FIXED;
         edgestate[erev] = PROP_STP_EDGE_FIXED;

         nodestate[tail] = BRANCH_STP_VERTEX_TERM;
         nodestate[head] = BRANCH_STP_VERTEX_TERM;

         assert(!( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[erev]) < 0.5 ));
      }
      /* both e AND its anti-parallel edge fixed to 0? */
      else if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[erev]) < 0.5 )
      {
         assert(SCIPvarGetLbLocal(vars[e]) < 0.5 && SCIPvarGetLbLocal(vars[erev]) < 0.5);

         edgestate[e] = PROP_STP_EDGE_KILLED;
         edgestate[erev] = PROP_STP_EDGE_KILLED;
      }
   }

   if( pcmw )
   {
      SCIP_Bool conflict = FALSE;

      getBoundchangesPcMW(scip, vars, propgraph, nodestate, &conflict);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      assert(!(BRANCH_STP_VERTEX_TERM == nodestate[k] && graph_pc_knotIsDummyTerm(propgraph, k)));
#endif

      if( conflict )
         *probisinfeas = TRUE;
   }

   /* not at root? */
   if( SCIPgetDepth(scip) > 0 && !(*probisinfeas) )
   {
      SCIP_Bool conflict = FALSE;

      SCIP_CALL( SCIPStpBranchruleGetVertexChgs(scip, nodestate, &conflict) );

      if( conflict )
         *probisinfeas = TRUE;
   }

   return SCIP_OKAY;
}



/** gets state of graph at current B&B node, including bound changes;
 *  takes direction of edges into account! */
static
SCIP_RETCODE getGraphStatesDirected(
   SCIP*                 scip,               /**< SCIP structure */
   const GRAPH*          graph,              /**< graph data structure */
   int*                  nodestate,          /**< node state (uninitialized) */
   int*                  edgestate,          /**< edge state (uninitialized) */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
)
{
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Bool isPcMw = graph_pc_isPcMw(graph);

   assert(vars);
   assert(!isPcMw || graph->extended);
   assert(!(*probisinfeas));

   SCIPStpBranchruleInitNodeState(graph, nodestate);

   for( int e = 0; e < nedges; e++ )
      edgestate[e] = PROP_STP_EDGE_UNSET;

   for( int e = 0; e < nedges; e++ )
   {
      if( SCIPvarGetLbLocal(vars[e]) > 0.5 )
      {
         assert(SCIPvarGetUbLocal(vars[e]) > 0.5);
         edgestate[e] = PROP_STP_EDGE_FIXED;
         nodestate[graph->tail[e]] = BRANCH_STP_VERTEX_TERM;
         nodestate[graph->head[e]] = BRANCH_STP_VERTEX_TERM;
      }
      else if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
      {
         assert(SCIPvarGetLbLocal(vars[e]) < 0.5);
         edgestate[e] = PROP_STP_EDGE_KILLED;
      }
   }

   if( isPcMw )
   {
      SCIP_Bool conflict = FALSE;

      getBoundchangesPcMW(scip, vars, graph, nodestate, &conflict);

      if( conflict )
         *probisinfeas = TRUE;
   }

   /* not at root? */
   if( SCIPgetDepth(scip) > 0 && !(*probisinfeas) && SCIPStpBranchruleIsActive(scip) )
   {
      SCIP_Bool conflict = FALSE;

      SCIP_CALL( SCIPStpBranchruleGetVertexChgs(scip, nodestate, &conflict) );

      if( conflict )
         *probisinfeas = TRUE;
   }

   return SCIP_OKAY;
}



/** trails graph and checks for infeasibility */
static
SCIP_RETCODE trailGraphWithStates(
   SCIP*                 scip,               /**< SCIP structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            nodestate,          /**< node state  */
   const int*            edgestate,          /**< edge state  */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
)
{
   STP_Vectype(int) queue = NULL;
   SCIP_Bool* nodes_isVisited;
   int termcount = 0;
   const int nnodes = graph_get_nNodes(graph);
   const int root = graph->source;

   assert(!(*probisinfeas));
   assert(nodestate[graph->source] == BRANCH_STP_VERTEX_TERM);

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodestate[i] == BRANCH_STP_VERTEX_TERM )
         termcount++;
   }

   assert(termcount >= graph->terms);

   /* now do a BFS from the root */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodes_isVisited, nnodes) );

   assert(!nodes_isVisited[root]);
   nodes_isVisited[root] = TRUE;
   StpVecReserve(scip, queue, nnodes);
   StpVecPushBack(scip, queue, root);

   /* BFS loop stopping at roots */
   for( int i = 0; i < StpVecGetSize(queue); i++ )
   {
      const int node = queue[i];
      assert(nodes_isVisited[node]);

      if( nodestate[node] == BRANCH_STP_VERTEX_TERM )
         termcount--;

      for( int e = graph->outbeg[node]; e >= 0; e = graph->oeat[e] )
      {
         const int head = graph->head[e];
         if( !nodes_isVisited[head] && edgestate[e] != PROP_STP_EDGE_KILLED )
         {
            nodes_isVisited[head] = TRUE;
            StpVecPushBack(scip, queue, head);
            assert(nodestate[head] != BRANCH_STP_VERTEX_KILLED);
         }
      }
   }

   assert(termcount >= 0);

   if( termcount != 0 )
   {
      *probisinfeas = TRUE;
   }

   StpVecFree(scip, queue);
   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}


/** initialize reduced cost distances */
static
SCIP_RETCODE getRedCost2ndNextDistances(
   SCIP*                 scip,               /**< SCIP structure */
   const SCIP_Real*      redcost,            /**< reduced costs */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi paths  */
   SCIP_Real*            pathdist,           /**< path distance */
   int*                  vbase,              /**< Voronoi base */
   int*                  state               /**< state  */
)
{
   SCIP_Real* redcostrev = NULL;
   int* pathedge = NULL;
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);

   assert(graph_isMarked(g));

   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcostrev, nedges) );

   /* distance from root to all nodes */
   graph_path_execX(scip, g, g->source, redcost, pathdist, pathedge);

   for( unsigned int e = 0; e < (unsigned) nedges; e++ )
      redcostrev[e] = redcost[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = g->outbeg[g->source]; e != EAT_LAST; e = g->oeat[e] )
      redcostrev[e] = FARAWAY;

   graph_get2nextTermPaths(g, redcostrev, redcostrev, vnoi, vbase, state);

   SCIPfreeBufferArray(scip, &redcostrev);
   SCIPfreeBufferArray(scip, &pathedge);

   return SCIP_OKAY;
}


/** marks arcs fixed to 0 */
static inline
void mark0FixedArcs(
   const GRAPH*          graph,              /**< graph structure */
   SCIP_VAR**            vars,               /**< variables (in) */
   STP_Bool*             arcIs0Fixed         /**< array (out) */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(vars && arcIs0Fixed);

#ifdef USE_EXTRED_FULL
   for( int e = 0; e < nedges; e++ )
   {
      arcIs0Fixed[e] = (SCIPvarGetUbGlobal(vars[e]) < 0.5);
   }
#else
   for( int e = 0; e < nedges; e++ )
   {
      arcIs0Fixed[e] = (SCIPvarGetUbLocal(vars[e]) < 0.5);
   }
#endif
}


/** some checks */
static
void validateEdgestate(
   const GRAPH*          graph,              /**< graph structure */
   const GRAPH*          propgraph,          /**< propagator graph */
   /*const*/ SCIP_VAR**  vars,               /**< variables */
   const int*            edgestate,          /**< edge state array */
   SCIP_Bool*            error               /**< error during update? */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   *error = FALSE;

   /* 1-fixed edge to be deleted? */
   for( int e = 0; e < nedges; e++ )
   {
      if( (edgestate[e] == PROP_STP_EDGE_UNSET || edgestate[e] == PROP_STP_EDGE_KILLED)
            && (SCIPvarGetLbLocal(vars[e]) > 0.5) )
      {
         if( pcmw && graph_pc_edgeIsExtended(graph, e) )
            continue;

         graph_edge_printInfo(propgraph, e);
         printf(
               "1-fixed arc deleted by reduction methods ... can't propagate  \n \n \n");
         *error = TRUE;
         return;
      }
   }

   /* 0-fixed edge been contracted? */
   for( int e = 0; e < nedges; e++ )
   {
      if( edgestate[e] == PROP_STP_EDGE_FIXED
            && SCIPvarGetUbLocal(vars[e]) < 0.5
            && SCIPvarGetUbLocal(vars[flipedge(e)]) < 0.5 )
      {
         graph_edge_printInfo(propgraph, e);
         printf(
               "0-fixed arc contracted by reduction methods ... can't propagate  \n \n \n");
         *error = TRUE;
         return;
      }
   }
}


/** helper method for reduction based variable fixings */
static inline
void setEdgestate(
   const GRAPH*          graph,              /**< graph structure */
   IDX*                  curr,               /**< current ancestor */
   int* RESTRICT         edgestate           /**< edge state array */
)
{
   while( curr != NULL )
   {
      const int i = curr->index;

      assert(i < graph->edges);
      assert(edgestate[i] != PROP_STP_EDGE_KILLED && edgestate[flipedge(i)] != PROP_STP_EDGE_KILLED);

      if( edgestate[i] == PROP_STP_EDGE_UNSET )
      {

#ifdef SCIP_DEBUG
         printf("set remain edge ");
         graph_edge_printInfo(graph, i);
#endif

         edgestate[i] = PROP_STP_EDGE_SET;
         edgestate[flipedge(i)] = PROP_STP_EDGE_SET;
      }
      curr = curr->parent;
   }
}


/** helper method for reduction based variable fixings */
static inline
void fixEdgestate(
   const GRAPH*          graph,              /**< graph structure */
   IDX*                  curr,               /**< current ancestor */
   int* RESTRICT         edgestate           /**< edge state array */
)
{
   while( curr != NULL )
   {
      const int e = curr->index;
      assert(e < graph->edges);

#ifdef SCIP_DEBUG
         printf("fix remain edge ");
         graph_edge_printInfo(graph, e);
#endif

      edgestate[e] = PROP_STP_EDGE_FIXED;
      edgestate[flipedge(e)] = PROP_STP_EDGE_FIXED;

      curr = curr->parent;
   }
}


/** method for reduction based variable fixings */
static
SCIP_RETCODE applyEdgestateToProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   SCIP_VAR**            vars,               /**< variables */
   const int*            edgestate,          /**< edge state array */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   /* fix edge variables to 0 */
   for( int e = 0; e < nedges; e += 2 )
   {
      /* edge not set yet? */
      if( edgestate[e] == PROP_STP_EDGE_UNSET )
      {
         const int erev = e + 1;

         assert(edgestate[erev] == PROP_STP_EDGE_UNSET);

         if( pcmw  )
         {
            assert(graph->extended);

            if( SCIPisGE(scip, graph->cost[e], FARAWAY) || SCIPisGE(scip, graph->cost[erev], FARAWAY) )
               continue;
         }

#ifdef SCIP_DEBUG
         printf("red-based fixing of ");
         graph_edge_printInfo(graph, e);
#endif

         SCIP_CALL( fixEdgeVar(scip, e, vars, propdata) );
         SCIP_CALL( fixEdgeVar(scip, erev, vars, propdata) );
      }
   }

   return SCIP_OKAY;
}


/** update method for reduction based variable fixings */
static
void updateEdgestateFromRed(
   const GRAPH*          graph,              /**< graph structure */
   const GRAPH*          propgraph,          /**< propagator graph */
   SCIP_VAR**            vars,               /**< variables */
   const int*            nodestate,          /**< node state array */
   int*                  edgestate,          /**< edge state array */
   SCIP_Bool*            error               /**< error during update? */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(graph_pc_isPcMw(graph) == graph_pc_isPcMw(propgraph));

   if( pcmw )
   {
      assert(graph->extended);
      setEdgestate(propgraph, propgraph->pcancestors[propgraph->source], edgestate);
   }

   for( int e = 0; e < nedges; e++ )
   {
      if( propgraph->ieat[e] != EAT_FREE )
      {
         assert(propgraph->ieat[flipedge(e)] != EAT_FREE);
         setEdgestate(propgraph, propgraph->ancestors[e], edgestate);
         if( pcmw )
         {
            setEdgestate(propgraph, propgraph->pcancestors[propgraph->head[e]], edgestate);
            setEdgestate(propgraph, propgraph->pcancestors[propgraph->tail[e]], edgestate);
         }
      }
   }

   fixEdgestate(propgraph, graph_getFixedges(propgraph), edgestate);

   validateEdgestate(graph, propgraph, vars, edgestate, error);
}


/** update method for reduction based variable fixings */
static
void updateEdgestateFromRedPcmw(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph structure */
   const GRAPH*          propgraph,          /**< propagator graph */
   SCIP_VAR**            vars,               /**< variables */
   const int*            nodestate,          /**< node state array */
   int*                  edgestate,          /**< edge state array */
   SCIP_Bool*            error               /**< error during update? */
)
{
   const int nnodes = graph_get_nNodes(propgraph);
   const int nedges = graph_get_nEdges(propgraph);
   IDX** ancestors = propgraph->ancestors;
   int nsetnodes = 0;
   int nsetedges = 0;
   STP_Bool* node_isSet;
   STP_Bool* edge_isSet;
   int* setnodesqueue;

   assert(graph_pc_isPcMw(propgraph));
   assert(nnodes == graph->knots);
   assert(nedges == graph->edges);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &node_isSet, nnodes) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &edge_isSet, nedges) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &setnodesqueue, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      node_isSet[k] = FALSE;

   for( int e = 0; e < nedges; e++ )
      edge_isSet[e] = FALSE;

   if( graph_pc_isRootedPcMw(propgraph) )
   {
      setnodesqueue[nsetnodes++] = propgraph->source;
      node_isSet[propgraph->source] = TRUE;
   }

   for( int e = 0; e <= nedges; e++ )
   {
      if( e == nedges || propgraph->ieat[e] != EAT_FREE )
      {
         // ugly as hell, use extra method!
         IDX* curr = (e < nedges) ? ancestors[e] : graph_getFixedges(propgraph);

         while( curr != NULL )
         {
            const int ancestoredge = curr->index;

            if( !edge_isSet[ancestoredge] )
            {
               edge_isSet[ancestoredge] = TRUE;
               nsetedges++;
            }
            if( !node_isSet[propgraph->tail[ancestoredge]] )
            {
               node_isSet[propgraph->tail[ancestoredge]] = TRUE;
               setnodesqueue[nsetnodes++] = propgraph->tail[ancestoredge];
            }
            if( !node_isSet[propgraph->head[ancestoredge]] )
            {
               node_isSet[propgraph->head[ancestoredge]] = TRUE;
               setnodesqueue[nsetnodes++] = propgraph->head[ancestoredge];
            }

            curr = curr->parent;
         }
      }
   }

   SCIP_CALL_ABORT( solstp_markPcancestors(scip, propgraph->pcancestors, propgraph->tail, propgraph->head,
         nnodes, node_isSet, edge_isSet, setnodesqueue, &nsetnodes, &nsetedges ) );

   for( int e = 0; e < nedges; e++ )
   {
      if( edge_isSet[e] )
      {
         if( edgestate[e] == PROP_STP_EDGE_UNSET )
         {
            edgestate[e] = PROP_STP_EDGE_SET;
            edgestate[flipedge(e)] = PROP_STP_EDGE_SET;
         }
      }
   }

   SCIPfreeBufferArray(scip, &setnodesqueue);
   SCIPfreeBufferArray(scip, &edge_isSet);
   SCIPfreeBufferArray(scip, &node_isSet);

   fixEdgestate(propgraph, graph_getFixedges(propgraph), edgestate);

   validateEdgestate(graph, propgraph, vars, edgestate, error);
}



/** updates fixing bounds for reduced cost fixings */
static
void updateEdgeLurkingBounds(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   SCIP_Real             lpobjal,            /**< LP objective  */
   SCIP_Real*            fixingbounds        /**< fixing bounds */
)
{
   const int nnodes = graph->knots;

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) && (graph->stp_type == STP_MWCSP || graph->stp_type == STP_PCSPG) )
         continue;

      for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const SCIP_Real fixbnd = pathdist[k] + cost[e] + vnoi[graph->head[e]].dist + lpobjal;
         if( fixbnd > fixingbounds[e] )
            fixingbounds[e] = fixbnd;
      }
   }
}


/* updates fixing bounds for reduced cost based constraints */
static
void updateDeg2LurkingBounds(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   SCIP_Real             lpobjal,            /**< LP objective  */
   SCIP_Real*            deg2bounds          /**< bounds */
)
{
   const int nnodes = graph->knots;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) )
      {
         const SCIP_Real fixbnd = pathdist[k] + vnoi[k].dist + vnoi[k + nnodes].dist + lpobjal;
         if( fixbnd > deg2bounds[k] )
            deg2bounds[k] = fixbnd;
      }
   }
}


/** fixes node of propgraph (ignores dummy terminals!) */
static inline
void propgraphFixNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node,               /**< the node */
   GRAPH*                propgraph,          /**< propagator graph */
   SCIP_Real*            offset              /**< pointer to the offset (in/out) */
   )
{
   assert(node >= 0 && node < propgraph->knots);
   assert(!propgraph->extended);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("try to fix ");
   graph_knot_printInfo(propgraph, node);
#endif


   if( graph_pc_isPcMw(propgraph) )
   {
      assert(!graph_pc_knotIsDummyTerm(propgraph, node));

      if( graph_pc_isRootedPcMw(propgraph) )
      {
         if( !graph_pc_knotIsFixedTerm(propgraph, node) )
         {
            graph_pc_knotToFixedTerm(scip, propgraph, node, offset);

            SCIPdebugMessage("...fixed \n");
         }
      }
      else
      {
         graph_pc_enforceNode(scip, propgraph, node, offset);

         SCIPdebugMessage("...enforced \n");
      }
   }
   else
   {
      graph_knot_chg(propgraph, node, STP_TERM);

      SCIPdebugMessage("...made standard terminal \n");
   }
}


/** deletes node of propgraph (ignores dummy terminals!) */
static inline
void propgraphDeleteNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node,               /**< the node */
   GRAPH*                propgraph,          /**< propagator graph */
   SCIP_Real*            offset              /**< offset pointer */
   )
{
   assert(scip && propgraph && offset);
   assert(node >= 0 && node < propgraph->knots);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("try to delete ");
   graph_knot_printInfo(propgraph, node);
#endif

   if( graph_pc_isPcMw(propgraph) )
   {
      assert(!graph_pc_knotIsDummyTerm(propgraph, node));
      assert(!propgraph->extended);
      assert(!graph_pc_knotIsFixedTerm(propgraph, node));

      if( Is_term(propgraph->term[node]) )
      {
         assert(graph_pc_knotIsPropPotTerm(propgraph, node) || graph_pc_knotIsNonLeafTerm(propgraph, node));

         graph_pc_deleteTerm(scip, propgraph, node, offset);

         SCIPdebugMessage("deleted term \n");
      }
      else
      {
         assert(!Is_anyTerm(propgraph->term[node]));

         graph_knot_del(scip, propgraph, node, TRUE);

         SCIPdebugMessage("deleted node \n");
      }
   }
   else
   {
      assert(!Is_term(propgraph->term[node]));

      graph_knot_del(scip, propgraph, node, TRUE);

      SCIPdebugMessage("deleted node \n");
   }
}


/** fixes edge of propgraph */
static
void propgraphFixEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   e,                  /**< the edge */
   GRAPH*                propgraph           /**< propagator graph */
   )
{
   const SCIP_Bool pcmw = graph_pc_isPcMw(propgraph);
   const int erev = flipedge(e);
#ifndef NDEBUG
   const int tail = propgraph->tail[e];
   const int head = propgraph->head[e];
#endif

   assert(!pcmw || !propgraph->extended);
   assert(e >= 0 && e < propgraph->edges);

   if( pcmw && (SCIPisGE(scip, propgraph->cost[e], FARAWAY) || SCIPisGE(scip, propgraph->cost[erev], FARAWAY)) )
   {
      assert(!SCIPisEQ(scip, propgraph->cost[e], propgraph->cost[erev]));
      return;
   }

   /* nothing has to be done here, because tail and head will be fixed
    * ...and later merged during presolving */
   if( graph_pc_isMw(propgraph) )
   {
      return;
   }

   assert(!pcmw || !graph_pc_knotIsDummyTerm(propgraph, tail));
   assert(!pcmw || !graph_pc_knotIsDummyTerm(propgraph, head));

   assert(graph_pc_isMw(propgraph) || SCIPisEQ(scip, propgraph->cost[e], propgraph->cost[erev]));

   propgraph->cost[e] = 0.0;
   propgraph->cost[erev] = 0.0;
}


/** deletes edge of propgraph */
static
void propgraphDeleteEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   e,                  /**< the edge */
   GRAPH*                propgraph           /**< propagator graph */
   )
{
   const SCIP_Bool pcmw = graph_pc_isPcMw(propgraph);
   const int erev = flipedge(e);
#ifndef NDEBUG
   const int tail = propgraph->tail[e];
   const int head = propgraph->head[e];
#endif

   assert(!pcmw || !propgraph->extended);
   assert(e >= 0 && e < propgraph->edges);

   if( pcmw )
   {
      assert(!propgraph->extended);

      if( SCIPisGE(scip, propgraph->cost[e], FARAWAY) || SCIPisGE(scip, propgraph->cost[erev], FARAWAY) )
      {
         assert(Is_term(propgraph->term[tail]) || Is_term(propgraph->term[head]));

         return;
      }

      assert(!graph_pc_knotIsDummyTerm(propgraph, tail));
      assert(!graph_pc_knotIsDummyTerm(propgraph, head));
   }

   graph_edge_del(scip, propgraph, e, TRUE);

}


/** apply current bound changes to propgraph */
static
void propgraphMarkFixedTermsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            nodestate,          /**< node state*/
   GRAPH*                propgraph,          /**< propagator graph */
   SCIP_Bool*            hasFixedTerms       /**< have terminals been fixed? */
)
{
   assert(hasFixedTerms && nodestate && scip && propgraph);
   assert(graph_pc_isPcMw(propgraph));
   assert(!propgraph->extended);

   *hasFixedTerms = FALSE;

   if( !graph_pc_isRootedPcMw(propgraph) )
   {
      const int nnodes = graph_get_nNodes(propgraph);

      for( int k = 0; k < nnodes; k++ )
      {
         if( nodestate[k] == BRANCH_STP_VERTEX_TERM && graph_pc_knotIsPropPotTerm(propgraph, k) )
         {
            graph_pc_knotChgPrize(propgraph, BLOCKED_MINOR, k);
            *hasFixedTerms = TRUE;
         }
      }
   }
}

/** helper method for propgraphApplyBoundchanges */
static inline
void propgraphApplyStates(
   SCIP*                 scip,               /**< SCIP data structure (in) */
   const int*            nodestate,          /**< node state (in) */
   const int*            edgestate,          /**< edge state (in) */
   GRAPH*                propgraph,          /**< propagator graph (in/out) */
   SCIP_Real*            offset              /**< pointer to the offset (in/out) */
   )
{
   const int nnodes = graph_get_nNodes(propgraph);
   const int nedges = graph_get_nEdges(propgraph);

   for( int e = 0; e < nedges; e += 2 )
   {
      assert(edgestate[e] == edgestate[e + 1]);
      assert(e + 1 == flipedge(e));

      if( PROP_STP_EDGE_FIXED == edgestate[e] )
      {
#ifndef NDEBUG
         const int tail = propgraph->tail[e];
         const int head = propgraph->tail[e];
         assert(BRANCH_STP_VERTEX_TERM == nodestate[tail]);
         assert(BRANCH_STP_VERTEX_TERM == nodestate[head]);
#endif

         propgraphFixEdge(scip, e, propgraph);
      }
      else if( PROP_STP_EDGE_KILLED == edgestate[e] )
      {
         propgraphDeleteEdge(scip, e, propgraph);
      }
   }

   for( int k = 0; k < nnodes; k++ )
   {
      if( BRANCH_STP_VERTEX_TERM == nodestate[k] )
      {
         propgraphFixNode(scip, k, propgraph, offset);
      }
      else if( BRANCH_STP_VERTEX_KILLED == nodestate[k] )
      {
         propgraphDeleteNode(scip, k, propgraph, offset);
      }
   }
}


/** Apply current implications resulting from bound changes to propgraph. */
static
SCIP_RETCODE propgraphApplyImplicationsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the original graph */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   SCIP_Real*            offset              /**< pointer to the offset */
   )
{
   const int* verts = SCIPStpGetPcImplVerts(scip);
   const int* start = SCIPStpGetPcImplStarts(scip);
   GRAPH* const propgraph = propdata->propgraph;

   assert(g && vars && propgraph && probisinfeas && offset);
   assert(graph_pc_isPcMw(propgraph));
   assert(g->extended);
   assert(!propgraph->extended);

   if( verts != NULL )
   {
      const int nnodes = graph_get_nNodes(propgraph);
      int ptermcount = 0;

      assert(start != NULL);
      assert(propgraph->knots == g->knots);

      for( int i = 0; i < nnodes && !(*probisinfeas); i++ )
      {
         if( !Is_pseudoTerm(g->term[i]) )
            continue;

         assert(graph_pc_knotIsPropPotTerm(g, i));

         ptermcount++;

         assert(ptermcount <= SCIPStpGetPcImplNstarts(scip));

         /* has the vertex 'i' been killed in propgraph? */
         if( propgraph->grad[i] == 0 )
         {
            assert(g->grad[i] > 0);

            for( int j = start[ptermcount - 1]; j < start[ptermcount]; j++ )
            {
               const int vert = verts[j];

               if( propgraph->grad[vert] == 0 )
                  continue;

               if( graph_pc_knotIsFixedTerm(propgraph, vert) )
               {
                  *probisinfeas = TRUE;
                  SCIPdebugMessage("implication tries to kill fixed terminal -> problem is infeasible \n");
                  break;
               }

               assert(!graph_pc_knotIsDummyTerm(propgraph, vert));

               if( Is_anyTerm((propgraph->term[vert])) )
               {
                  assert(graph_pc_knotIsPropPotTerm(propgraph, vert) || graph_pc_knotIsNonLeafTerm(propgraph, vert));

                  graph_pc_deleteTerm(scip, propgraph, vert, offset);
               }
               else
               {
                  assert(!graph_pc_knotIsPropPotTerm(propgraph, vert) && !graph_pc_knotIsNonLeafTerm(propgraph, vert));

                  graph_knot_del(scip, propgraph, vert, TRUE);
               }
            }
         }
         else if( graph_pc_knotIsPropPotTerm(propgraph, i) )
         {
            assert(SCIPisLT(scip, propgraph->prize[i], BLOCKED_MINOR)); // should have been fixed already!

            for( int j = start[ptermcount - 1]; j < start[ptermcount]; j++ )
            {
               const int vert = verts[j];

               assert(SCIPisLT(scip, propgraph->prize[vert], BLOCKED_MINOR) || graph_pc_knotIsFixedTerm(propgraph, vert) );

               /* is 'vert' fixed? */
               if( graph_pc_knotIsFixedTerm(propgraph, vert) )
               {
                  const int twin = graph_pc_getTwinTerm(propgraph, i);
                  const int rootedge = graph_pc_getRoot2PtermEdge(g, twin);
                  SCIP_CALL( fixEdgeVar(scip, rootedge, vars, propdata ) );

                  graph_pc_knotToFixedTerm(scip, propgraph, i, offset);
                  break;
               }
            }
         }
      }

      assert(*probisinfeas || SCIPStpGetPcImplNstarts(scip) == ptermcount);
      assert(*probisinfeas || graph_pc_nNonFixedTerms(g) == ptermcount);
   } /* verts != NULL */

   return SCIP_OKAY;
}

/** Prunes unconnected vertices and edges.
 *  Also checks for resulting infeasibility */
static
SCIP_RETCODE propgraphPruneUnconnected(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                propgraph,          /**< propagator graph (in) */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   SCIP_Real*            offset              /**< pointer to the offset */
   )
{
   assert(probisinfeas);
   assert(!(*probisinfeas));

   if( graph_pc_isRootedPcMw(propgraph) )
      SCIP_CALL( reduce_unconnectedRpcRmwInfeas(scip, propgraph, offset, probisinfeas) );
   else
      SCIP_CALL( reduce_unconnectedInfeas(scip, FALSE, propgraph, probisinfeas) );

   return SCIP_OKAY;
}

/** Apply current bound changes to propgraph.
 *  Note that the graph type of propgraph might be changed! */
static
SCIP_RETCODE propgraphApplyBoundchanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          g,                  /**< graph data structure */
   int*                  nodestate,          /**< node state (uninitialized) */
   int*                  edgestate,          /**< edge state (uninitialized) */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   SCIP_Real*            offset              /**< pointer to the offset */
   )
{
   GRAPH* const propgraph = propdata->propgraph;
   const SCIP_Bool pcmw = graph_pc_isPcMw(propgraph);

   assert(scip && vars && g && edgestate && nodestate && probisinfeas && offset);
   assert(g->stp_type == propgraph->stp_type);
   assert(propgraph->knots == g->knots);
   assert(propgraph->edges == g->edges);
   assert(!(*probisinfeas));

   SCIP_CALL( getBoundchanges(scip, vars, propgraph, nodestate, edgestate, probisinfeas) );

   if( *probisinfeas )
   {
      SCIP_CALL( graph_path_init(scip, propgraph) );
      return SCIP_OKAY;
   }

   if( pcmw )
   {
      SCIP_Bool rootable = FALSE;
      graph_pc_2org(scip, propgraph);

      propgraphMarkFixedTermsPcMw(scip, nodestate, propgraph, &rootable);

      if( rootable )
      {
         assert(!graph_pc_isRootedPcMw(propgraph));

         graph_transPcmw2rooted(scip, propgraph, BLOCKED_MINOR, FALSE);

         assert(graph_pc_isRootedPcMw(propgraph));
      }
   }

   SCIP_CALL( graph_path_init(scip, propgraph) );

   propgraphApplyStates(scip, nodestate, edgestate, propgraph, offset);

   *probisinfeas = FALSE;

   if( pcmw )
   {
      SCIP_CALL( propgraphApplyImplicationsPcMw(scip, g, vars, propdata, probisinfeas, offset) );

      graph_pc_2trans(scip, propgraph);
   }

   return SCIP_OKAY;
}


/** initializes */
static inline
SCIP_RETCODE initPropgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   assert(scip && graph && propdata);

   assert(propdata->propgraph == NULL);

   propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   SCIP_CALL( graph_copy(scip, graph, &(propdata->propgraph)) );

   propdata->propgraph->norgmodeledges = propdata->propgraph->edges;
   propdata->propgraph->norgmodelknots = propdata->propgraph->knots;

   SCIP_CALL( graph_initHistory(scip, propdata->propgraph) );
#ifndef WITH_UG
   SCIP_CALL( graph_copyPseudoAncestors(scip, graph, propdata->propgraph) );
#endif

   assert(propdata->nfixededges_all == 0);
   assert(propdata->nfixededges_curr == 0);
   assert(propdata->nfixededges_bicurr == 0);
   assert(propdata->propgraph != NULL);

   return SCIP_OKAY;
}


/** initializes */
static inline
SCIP_Bool useRedcostdata(
   const GRAPH*          graph               /**< graph structure to use for the update */
)
{
   return (graph_typeIsSpgLike(graph) || graph_pc_isPc(graph));
}


/** writes reduced costs into given level */
static
void writeRedcostdata(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   level,              /**< the level */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   REDCOST* redcostdata = propdata->redcostdata;
   SCIP_Real* const redcost = redcosts_getEdgeCosts(redcostdata, level);
   const SCIP_Real lpobjval = propdata->lpobjval;
   const SCIP_Real cutoffbound = getCutoffbound(scip, lpobjval);
   const SCIP_Real cutoffgap = cutoffbound - lpobjval;

   assert(GE(lpobjval, 0.0) && LE(lpobjval, cutoffbound));
   assert(GE(cutoffgap, 0.0));

   redcosts_forLPget(scip, vars, graph, redcost);
   redcosts_setDualBound(lpobjval, level, redcostdata);
   redcosts_setCutoff(cutoffgap, level, redcostdata);
   redcosts_setRoot(graph->source, level, redcostdata);
}


/** initializes */
static inline
SCIP_RETCODE initRedcostdata(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
#ifdef USE_EXTRED_FULL
   int nlevels = PROP_STP_REDCOST_LEVELS;
#else
   int nlevels = 1;
#endif
   RCPARAMS rcparams = { .cutoff = -1.0, .nLevels = nlevels, .nCloseTerms = 3,
                              .nnodes = nnodes, .nedges = nedges, .redCostRoot = -1 };

   assert(scip && graph && propdata);
   assert(propdata->redcostdata == NULL);
   assert(propdata->redcostnupdates == 0);

   SCIP_CALL( redcosts_initFromParams(scip, &rcparams, &(propdata->redcostdata)) );
   writeRedcostdata(scip, 0, graph, vars, propdata);

   return SCIP_OKAY;
}


/** updates */
static inline
void updateRedcostdata(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
#ifdef USE_EXTRED_FULL
   REDCOST* redcostdata = propdata->redcostdata;;
#endif

   assert(propdata->redcostdata);

   /* NOTE: ugly workaround */
   if( propdata->redcostnupdates == 0 )
   {
      assert(redcosts_getNlevels(propdata->redcostdata) == 1);
      propdata->redcostnupdates++;
      return;
   }

#ifdef USE_EXTRED_FULL
   assert(redcosts_getNlevels(redcostdata) >= 1);

   if( redcosts_getNlevels(redcostdata) < PROP_STP_REDCOST_LEVELS )
   {
      redcosts_addLevel(redcostdata);
      writeRedcostdata(scip, redcosts_getLevel(redcostdata), graph, vars, propdata);

      propdata->redcostnupdates++;
   }
   else
   {
      const SCIP_Bool overwrite = SCIPisGT(scip, propdata->lpobjval, propdata->lpobjval_last); //(SCIPrandomGetInt(propdata->randnumgen, 0, 1) == 0);

      assert(SCIPisGE(scip, propdata->lpobjval, propdata->lpobjval_last));

      if( overwrite )
      {
         const int level = SCIPrandomGetInt(propdata->randnumgen, 1, PROP_STP_REDCOST_LEVELS - 1);

         SCIPdebugMessage("overwrite level %d \n", level);

         writeRedcostdata(scip, level, graph, vars, propdata);

         propdata->redcostnupdates++;
      }
   }
#endif

   propdata->lpobjval_last = propdata->lpobjval;
}


/** update the graph from propdata from given graph */
static inline
SCIP_RETCODE updatePropgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   GRAPH* propgraph = propdata->propgraph;

   assert(propgraph != NULL);

   graph_freeHistory(scip, propgraph);
   graph_freeHistoryDeep(scip, propgraph);

   SCIP_CALL( graph_copyData(scip, graph, propgraph) );

   propgraph->norgmodeledges = propgraph->edges;
   propgraph->norgmodelknots = propgraph->knots;
   propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   assert(!graph_pc_isRootedPcMw(graph) || graph_pc_nFixedTerms(graph) == graph_pc_nFixedTerms(propgraph));

   SCIP_CALL( graph_initHistory(scip, propdata->propgraph) );

   return SCIP_OKAY;
}


/** try to make global fixings based on lurking bounds */
static
SCIP_RETCODE fixVarsDualcostLurking(
   SCIP*                 scip,               /**< SCIP structure */
   const GRAPH*          graph,              /**< graph structure */
   SCIP_Real             cutoffbound,        /**< cutoff bound  */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR**            vars                /**< variables */
)
{
   const SCIP_Real* fixingbounds = propdata->fixingbounds;
   const SCIP_Real* deg2bounds = propdata->deg2bounds;
   SCIP_Bool* deg2bounded = propdata->deg2bounded;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(fixingbounds != NULL && deg2bounds != NULL && deg2bounded != NULL);

   for( int e = 0; e < nedges; e++ )
   {
      if( SCIPisLT(scip, cutoffbound, fixingbounds[e]) )
      {
         SCIP_VAR* const edgevar = vars[e];

         if( SCIPvarGetLbGlobal(edgevar) < 0.5 && SCIPvarGetUbGlobal(edgevar) > 0.5 )
         {
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(edgevar), 1.0));

            SCIPchgVarUbGlobal(scip, edgevar, 0.0);
            propdata->nfixededges_curr++;
         }
      }
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, cutoffbound, deg2bounds[i]) && !Is_term(graph->term[i]) )
         deg2bounded[i] = TRUE;
   }

#ifdef SCIP_DEBUG
   {
      int d2bounded = 0;
      for( int i = 0; i < nnodes; i++ )
         if( deg2bounded[i] )
            d2bounded++;

      printf("deg2bounded=%d \n", d2bounded);

   }
#endif

   return SCIP_OKAY;
}


/** dual cost based fixing of variables */
static
SCIP_RETCODE fixVarsDualcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
)
{
   PATH* vnoi;
   SCIP_Real* redcost;
   SCIP_Real* pathdist;
   const SCIP_Real lpobjval = propdata->lpobjval;
   const SCIP_Real cutoffbound = getCutoffbound(scip, lpobjval);
   const SCIP_Real cutoffgap = cutoffbound - lpobjval;
   int* vbase;
   int* state;
   int* termorg = NULL;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(FALSE == *probisinfeas);
 //  printf("minpathcost=%f \n", minpathcost);

   if( SCIPisLT(scip, cutoffgap, 0.0) )
   {
      SCIPdebugMessage("infeasible sub-problem! \n");
      *probisinfeas = TRUE;
      return SCIP_OKAY;
   }

   if( propdata->fixingbounds == NULL )
   {
      assert(propdata->deg2bounds == NULL && propdata->deg2bounded == NULL);

      SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->deg2bounds), nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->deg2bounded), nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->fixingbounds), nedges) );

      for( int i = 0; i < nedges; i++ )
         propdata->fixingbounds[i] = -FARAWAY;

      for( int i = 0; i < nnodes; i++ )
      {
         propdata->deg2bounds[i] = -FARAWAY;
         propdata->deg2bounded[i] = FALSE;
      }

#ifndef WITH_UG
      /* first call, so we can also fix incoming arcs of root to zero */
      for( int e = graph->inpbeg[graph->source]; e != EAT_LAST; e = graph->ieat[e] )
         SCIP_CALL( fixEdgeVar(scip, e, vars, propdata) );
#endif
   }

#ifdef WITH_UG
   if( !redcosts_forLPareReliable(scip, vars, graph) )
   {
      graph_mark(graph);
      return SCIP_OKAY;
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &state, 2 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 2 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 2 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );

   if( SCIPgetDepth(scip) > 0 && SCIPStpBranchruleIsActive(scip) )
   {
      int* nodestatenew;
      SCIP_Bool conflict = FALSE;

      SCIP_CALL( SCIPallocBufferArray(scip, &termorg, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nodestatenew, nnodes) );

      BMScopyMemoryArray(termorg, graph->term, nnodes);
      SCIPStpBranchruleInitNodeState(graph, nodestatenew);
      SCIP_CALL( SCIPStpBranchruleGetVertexChgs(scip, nodestatenew, &conflict) );
      assert(!conflict);

      for( int k = 0; k < nnodes; k++ )
      {
         if( nodestatenew[k] == BRANCH_STP_VERTEX_TERM && !Is_term(graph->term[k]) )
            graph_knot_chg(graph, k, STP_TERM);

      }

      SCIPfreeBufferArray(scip, &nodestatenew);
   }

   graph_mark(graph);

   redcosts_forLPget(scip, vars, graph, redcost);
   SCIP_CALL( getRedCost2ndNextDistances(scip, redcost, graph, vnoi, pathdist, vbase, state) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, cutoffgap) )
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            /* try to fix edge and reversed one */
            SCIP_CALL( fixEdgeVar(scip, e, vars, propdata) );
            SCIP_CALL( fixEdgeVar(scip, flipedge(e), vars, propdata) );
         }
      }
      else
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            if( SCIPisGT(scip, pathdist[k] + redcost[e] + vnoi[graph->head[e]].dist, cutoffgap) )
               SCIP_CALL( fixEdgeVar(scip, e, vars, propdata) );
         }
      }
   }

   if( termorg )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( graph->term[k] != termorg[k] )
            graph_knot_chg(graph, k, termorg[k]);
      }
      SCIPfreeBufferArray(scip, &termorg);
   }

   /* at root? */
   if( SCIPgetDepth(scip) == 0 )
   {
      updateEdgeLurkingBounds(graph, redcost, pathdist, vnoi, lpobjval, propdata->fixingbounds);
      updateDeg2LurkingBounds(graph, pathdist, vnoi, lpobjval, propdata->deg2bounds);
   }

   SCIP_CALL( fixVarsDualcostLurking(scip, graph, cutoffbound, propdata, vars) );

   SCIPfreeBufferArray(scip, &pathdist);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &vbase);

   return SCIP_OKAY;
}


/** extended reduced cost reduction */
static
SCIP_RETCODE fixVarsExtendedRed(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   REDCOST* const redcostdata = propdata->redcostdata;
   GRAPH* const propgraph = propdata->propgraph;
   STP_Bool* arcdeleted = NULL;
   int nfixededges = 0;
   const int nedges = propgraph->edges;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(propgraph);
#ifdef USE_EXTRED_FULL
   const SCIP_Real cutoffbound = getCutoffbound(scip, propdata->lpobjval);
#endif

   assert(redcostdata);

   /* in this case the reduced cost reductions are not valid anymore! (because the root has changed) */
   if( graph->stp_type != propgraph->stp_type )
   {
      assert(graph_pc_isRootedPcMw(propgraph));
      return SCIP_OKAY;
   }
#ifdef WITH_UG
   return SCIP_OKAY;
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &arcdeleted, nedges) );
   mark0FixedArcs(graph, vars, arcdeleted);

   graph_mark(propgraph);

   if( isPcMw )
      graph_pc_2org(scip, propgraph);
#ifdef USE_EXTRED_FULL
   SCIPdebugMessage("starting extended reductions with %d red. cost levels \n", redcosts_getNlevels(redcostdata));

   for( int i = 0; i < redcosts_getNlevels(redcostdata); i++ )
   {
      redcosts_setCutoffFromBound(cutoffbound, i, redcostdata);

      redcosts_increaseOnDeletedArcs(propgraph, arcdeleted, i, redcostdata);

      redcosts_unifyBlockedEdgeCosts(propgraph, i, redcostdata);

      SCIP_CALL( redcosts_initializeDistances(scip, i, propgraph, redcostdata) );

      SCIPdebugMessage("cutoff for level %d: %f \n", i, redcosts_getCutoff(redcostdata, i));
   }

   /* Note that all in-root arcs will be set to 0! (w.r.t redCostRoot) */
   /* reduce graph and mark deletable arcs */
   SCIP_CALL( extreduce_deleteArcs(scip, redcostdata, NULL, propgraph,
         arcdeleted, &nfixededges) );
#else
   writeRedcostdata(scip, 0, graph, vars, propdata);
   redcosts_increaseOnDeletedArcs(propgraph, arcdeleted, 0, redcostdata);
   SCIP_CALL( redcosts_initializeDistances(scip, 0, propgraph, redcostdata) );

   nfixededges = reduce_extendedEdge(scip, propgraph, redcosts_getNodeToTermsPathsTop(redcostdata),
        redcosts_getEdgeCostsTop(redcostdata), redcosts_getRootToNodeDistTop(redcostdata),
        NULL, redcosts_getCutoffTop(redcostdata), propgraph->source, arcdeleted, TRUE);
#endif

   if( isPcMw )
      graph_pc_2trans(scip, propgraph);

   SCIPdebugMessage("extended-reduction graph deletions: %d \n", nfixededges);

   for( int e = 0; e < nedges; e++ )
   {
      if( SCIPvarGetUbLocal(vars[e]) > 0.5 && arcdeleted[e] )
      {
         if( isPcMw && graph_pc_edgeIsExtended(propgraph, e) )
            continue;

         SCIPdebugMessage("fix edge %d to 0 \n", e);

         SCIP_CALL(fixEdgeVar(scip, e, vars, propdata));
      }
   }

   SCIPdebugMessage("after extended-reduction number of fixed variables: %d \n", propdata->nfixededges_curr);

   SCIPfreeBufferArray(scip, &arcdeleted);

   assert(nfixededges >= 0);

   return SCIP_OKAY;
}


/** This methods tries to fix edges by performing reductions on the current graph.
 *  To this end, the already 0-fixed edges are (temporarily) removed from the underlying graph to strengthen the reduction methods. */
static
SCIP_RETCODE fixVarsRedbased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
   )
{
   REDSOL* redsol;
   GRAPH* propgraph = NULL;
   int* nodestate = NULL;
   int* edgestate = NULL;
   SCIP_Real offset = -1.0;  /* don't use it for pcmw! */
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   SCIP_Bool error;

   assert(propdata && scip && vars && probisinfeas);
   assert(!graph_pc_isPcMw(graph) || graph->extended);

   SCIPdebugMessage("start redbasedVarfixing, fixings: %d \n", propdata->nfixededges_curr);

   *probisinfeas = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgestate, nedges) );

   SCIP_CALL( updatePropgraph(scip, graph, propdata) );
   propgraph = propdata->propgraph;

   /* note: graph type might change (from PC/MW to RPC/RMW) */
   SCIP_CALL( propgraphApplyBoundchanges(scip, vars, graph, nodestate, edgestate, propdata, probisinfeas,
        &offset ) );

   if( *probisinfeas )
   {
      SCIPdebugMessage("problem has become infeasible after applying bound changes: terminating! \n");
      goto TERMINATE;
   }

   /* remove unconnected vertices and edges */
   SCIP_CALL( propgraphPruneUnconnected(scip, propgraph, probisinfeas, &offset ) );

   if( *probisinfeas )
   {
      SCIPdebugMessage("problem has become infeasible after pruning (terminals not connected): terminating! \n");
      goto TERMINATE;
   }

   assert(graph_valid(scip, propgraph));

   /* start with extended reductions based on reduced costs of LP */
   if( !graph_pc_isMw(propgraph) )
   {
      SCIP_CALL( fixVarsExtendedRed(scip, graph, vars, propdata) );
   }

   SCIP_CALL( reduce_solInit(scip, propgraph, FALSE, &redsol) );

   /* now reduce the graph by standard reductions */
   if( graph_pc_isPc(propgraph) )
   {
      SCIP_CALL( reduce_pc(scip, redsol, propgraph, 2, FALSE, FALSE, FALSE, FALSE) );
   }
   else if( graph_pc_isMw(propgraph) )
   {
      SCIPdebugMessage("starting MW reductions \n");
      SCIP_CALL( reduce_mw(scip, redsol, propgraph, 2, FALSE, FALSE, FALSE) );
   }
   else
   {

#ifdef SCIP_DISABLED
      {
      DISTDATA* distdata;
      EXTPERMA* extpermanent;
      REDCOST* redcostdata = propdata->redcostdata;
      const SCIP_Real cutoffbound = getCutoffbound(scip, propdata->lpobjval);

      int ndeleted;


      for( int i = 0; i < redcosts_getNlevels(redcostdata); i++ )
      {
         redcosts_setCutoffFromBound(cutoffbound, i, redcostdata);

//         redcosts_increaseOnDeletedArcs(propgraph, arcdeleted, i, redcostdata);

         redcosts_unifyBlockedEdgeCosts(propgraph, i, redcostdata);

         SCIP_CALL( redcosts_initializeDistances(scip, i, propgraph, redcostdata) );

         SCIPdebugMessage("cutoff for level %d: %f \n", i, redcosts_getCutoff(redcostdata, i));
      }

      SCIP_CALL( graph_init_dcsr(scip, propgraph) );
      SCIP_CALL( extreduce_distDataInit(scip, propgraph, 40, TRUE, FALSE, &distdata) );
      SCIP_CALL( extreduce_extPermaInit(scip, extred_fast, propgraph, NULL, &extpermanent) );

      assert(!extpermanent->distdata_default);

      extpermanent->distdata_default = distdata;
      extpermanent->redcostdata = propdata->redcostdata;
      extpermanent->redcostEqualAllow = FALSE;

      SCIP_CALL( extreduce_pseudoDeleteNodes(scip, propdata->deg2bounded, extpermanent, propgraph, NULL, &ndeleted) );

      /* clean up */
      extreduce_distDataFree(scip, propgraph, &distdata);

      if( extpermanent->distdata_biased )
         extreduce_distDataFree(scip, propgraph, &(extpermanent->distdata_biased));

      extreduce_extPermaFree(scip, &extpermanent);

      printf("ndeleted=%d \n", ndeleted);

     graph_free_dcsr(scip, propgraph);
      }
#endif

      // todo Call two times, and with node-replacing!
      // todo: before make all the node replacements from lurking bounds!
      assert(graph_typeIsSpgLike(propgraph));
      SCIP_CALL( reduce_unconnected(scip, propgraph) );
      SCIP_CALL( reduce_stp(scip, propgraph, redsol, 2, FALSE, FALSE, FALSE, FALSE) );

    //  SCIP_CALL( reduceStp(scip, propgraph, redsol, 2, FALSE, TRUE, FALSE) );

   }
   offset += reduce_solGetOffset(redsol);
   reduce_solFree(scip, &redsol);

   assert(graph_valid(scip, propgraph));

   /* mark surviving original edges of propgraph reductions in array 'edgestate' */
   if( graph_pc_isPcMw(propgraph) )
   {
      updateEdgestateFromRedPcmw(scip, graph, propgraph, vars, nodestate, edgestate, &error);
   }
   else
   {
      updateEdgestateFromRed(graph, propgraph, vars, nodestate, edgestate, &error);
   }

   if( error )
   {
      return SCIP_ERROR; // todo check! is terminating enough?
   }

   SCIP_CALL( applyEdgestateToProb(scip, graph, vars, edgestate, propdata) );

   SCIPdebugMessage("number of reduction-based variable fixings: %d \n", propdata->nfixededges_curr);

TERMINATE:
   graph_path_exit(scip, propgraph);
   SCIPfreeBufferArray(scip, &edgestate);
   SCIPfreeBufferArray(scip, &nodestate);

   return SCIP_OKAY;
}



/** promising? */
static
SCIP_Bool fixVarsRedbasedIsPromising(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_Bool callreduce = FALSE;

   assert(scip && graph && vars && propdata);

   if( graph_typeIsSpgLike(graph) || (graph_pc_isPcMw(graph) && graph->stp_type != STP_BRMWCSP) )
   {
      /* in the tree? */
      if( SCIPgetDepth(scip) > 0 )
      {
         SCIP_NODE* const currnode = SCIPgetCurrentNode(scip);
         const SCIP_Longint nodenumber = SCIPnodeGetNumber(currnode);

         if( nodenumber != propdata->lastnodenumber_red || propdata->aggressive )
         {
            propdata->lastnodenumber_red = nodenumber;
            callreduce = TRUE;
         }
      }
      /* at root (== 0 is necessary, could be -1) */
      else if( SCIPgetDepth(scip) == 0 )
      {
         const SCIP_Real redratio = (2.0 * (SCIP_Real) propdata->nfixededges_bipost) / ((SCIP_Real) graph->edges);

         assert(GE(propdata->redwaitratio, REDUCTION_WAIT_RATIO_INITIAL));

         if( useRedcostdata(graph) )
            updateRedcostdata(scip, graph, vars, propdata);

         if( SCIPisGT(scip, redratio, propdata->redwaitratio) || propdata->redcostnupdates == PROP_STP_REDCOST_LEVELS )
         {
            callreduce = TRUE;
            propdata->redwaitratio *= REDUCTION_WAIT_FACTOR;
         }
      }

#ifdef WITH_UG
      if( propdata->ncalls == 1 && SCIPgetDepth(scip) == 0 )
      {
         SCIPdebugMessage("trigger UG root reductions \n");
         callreduce = TRUE;
      }
#endif
   }

   return callreduce;
}

/** block edges of the underlying graphs by using global fixings */
static inline
void blockEdgesWithGlobalFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   GRAPH*                graph               /**< graph data structure */
)
{
   const int nedges = graph_get_nEdges(graph);

   assert(vars && scip);

   if( !graph_typeIsSpgLike(graph) && !graph_pc_isPc(graph) && graph->stp_type != STP_DCSTP )
   {
      return;
   }

   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      /* both e and its anti-parallel edge fixed to zero? */
      if( SCIPvarGetUbGlobal(vars[e]) < 0.5 && SCIPvarGetUbGlobal(vars[erev]) < 0.5 && LT(graph->cost[e], BLOCKED) )
      {
         assert(SCIPvarGetLbLocal(vars[e]) < 0.5 && SCIPvarGetLbLocal(vars[erev]) < 0.5);

         // todo for pc/rpc can we maybe also block edges with unequal costs? Seems to lead to bugs so far...
         if( EQ(graph->cost[e], graph->cost[erev]) )
         {
            graph->cost[e] = BLOCKED;
            graph->cost[erev] = BLOCKED;
         }
      }
   }
}


/** initializes for first call */
static
SCIP_RETCODE initPropAtFirstCall(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   assert(scip && graph && vars && propdata);

   if( graph_typeIsSpgLike(graph) || (graph_pc_isPcMw(graph) && graph->stp_type != STP_BRMWCSP) )
   {
      if( !propdata->propgraph )
      {
         assert(!propdata->isInitialized);

         SCIP_CALL( initPropgraph(scip, graph, propdata) );

         assert(propdata->propgraph);
      }

      if( !propdata->isInitialized && useRedcostdata(graph) )
      {
         SCIP_CALL( initRedcostdata(scip, graph, vars, propdata) );
      }
   }

   if( !propdata->isInitialized && graph->stp_type == STP_DCSTP )
   {
      const int nnodes = graph_get_nNodes(graph);

      assert(graph->maxdeg);

      for( int k = 0; k < nnodes; k++ )
      {
         if( graph->maxdeg[k] != 1 )
            continue;

         if( !Is_term(graph->term[k]) || k == graph->source )
            continue;

         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            if( SCIPvarGetUbLocal(vars[e]) > 0.5 )
            {
               assert(SCIPvarGetLbLocal(vars[e]) < 0.5);
               SCIP_CALL( SCIPchgVarUb(scip, vars[e], 0.0) );
            }
         }
      }
   }

   propdata->isInitialized = TRUE;

   return SCIP_OKAY;
}



/** applies vertex branching changes.
 *  NOTE: this is necessary, because SCIP is crazy slow in propagating constraints,
 *  so we do it ourselves */
static
SCIP_RETCODE applyLastVertexBranch(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   assert(SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0);

   if( SCIPStpBranchruleIsActive(scip) )
   {
      int branchvertex;
      SCIP_Bool isDeleted = FALSE;
      SCIP_CALL( SCIPStpBranchruleGetVertexChgLast(scip, &branchvertex, &isDeleted) );

      if( isDeleted )
      {
         for( int e = graph->outbeg[branchvertex]; e != EAT_LAST; e = graph->oeat[e] )
         {
            SCIP_CALL( fixEdgeVar(scip, e, vars, propdata) );
            SCIP_CALL( fixEdgeVar(scip, flipedge(e), vars, propdata) );
         }
      }
   }

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of propagator
 *
 * @{
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropStp(scip) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeStp)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeRandom(scip, &propdata->randnumgen);

   SCIPfreeMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** reduced cost propagation method for an LP solution */
static
SCIP_DECL_PROPEXEC(propExecStp)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   GRAPH* graph;
   SCIP_Bool probisinfeas = FALSE;
   SCIP_NODE* currnode;

   *result = SCIP_DIDNOTRUN;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   vars = SCIPprobdataGetVars(scip);
   assert(vars);

   propdata = SCIPpropGetData(prop);
   graph = SCIPprobdataGetGraph(probdata);
   assert(graph);

   currnode = SCIPgetCurrentNode(scip);

   if( SCIPnodeGetDepth(currnode) > 0 && propdata->lastnodenumber != SCIPnodeGetNumber(currnode) )
   {
      /* propagate based on last vertex branching decision */
      SCIP_CALL( applyLastVertexBranch(scip, graph, vars, propdata) );

      SCIP_CALL( SCIPStpPropCheckForInfeas(scip, &probisinfeas) );

      propdata->lastnodenumber = SCIPnodeGetNumber(currnode);
   }

   if( probisinfeas )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* check if all integral variables are fixed */
   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if( !redcosts_forLPareAvailable(scip) )
      return SCIP_OKAY;

   propdata->ncalls++;

   if( propdata->nfails > 0 && (propdata->nlastcall + propdata->maxnwaitrounds >= propdata->ncalls)
     && (propdata->nlastcall + propdata->nfails > propdata->ncalls) )
      return SCIP_OKAY;

   /* new LP needs to have been solved to avoid conflicts between reduction based and red. cost propagation */
   if( propdata->nlastnlpiter == SCIPgetNLPIterations(scip) )
      return SCIP_OKAY;

   propdata->nlastnlpiter = SCIPgetNLPIterations(scip);
   propdata->nlastcall = propdata->ncalls;

   propdata->nfixededges_curr = 0;
   propdata->nfixededges_bicurr = 0;
   *result = SCIP_DIDNOTFIND;

   propdata->lpobjval = SCIPgetLPObjval(scip);

   SCIP_CALL( initPropAtFirstCall(scip, graph, vars, propdata) );

   /* call dual cost based variable fixing */
   SCIP_CALL( fixVarsDualcost(scip, vars, propdata, graph, &probisinfeas) );

   if( !probisinfeas && fixVarsRedbasedIsPromising(scip, graph, vars, propdata) )
   {
      assert(!probisinfeas);
      SCIPdebugMessage("use reduction techniques \n");

      /* call reduced cost based based variable fixing */
      SCIP_CALL( fixVarsRedbased(scip, graph, vars, propdata, &probisinfeas) );

      propdata->nfixededges_bipost = 0;
   }

   if( probisinfeas )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( propdata->nfixededges_curr > 0 )
   {
      SCIPdebugMessage("newly fixed by STP propagator: %d of (%d) \n", propdata->nfixededges_curr, propdata->nfixededges_all);

      propdata->nfails = 0;
      propdata->nfixededges_all += propdata->nfixededges_curr;
      propdata->nfixededges_bipost += propdata->nfixededges_bicurr;

      *result = SCIP_REDUCEDDOM;

      /* translate the global fixings of variables into blocking of graph edges */
      blockEdgesWithGlobalFixings(scip, vars, graph);
   }
   else
   {
      propdata->nfails++;
   }

   return SCIP_OKAY;
}

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolStp)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->nfails = 0;
   propdata->ncalls = 0;
   propdata->nlastcall = 0;
   propdata->nlastnlpiter = 0;
   propdata->lastnodenumber_red = -1;
   propdata->lastnodenumber = -1;
   propdata->nfixededges_all = 0;
   propdata->redwaitratio = REDUCTION_WAIT_RATIO_INITIAL;
   propdata->nfixededges_curr = 0;
   propdata->nfixededges_bicurr = 0;
   propdata->nfixededges_bipost = 0;
   propdata->fixingbounds = NULL;
   propdata->deg2bounded = NULL;
   propdata->deg2bounds = NULL;
   propdata->propgraph = NULL;
   propdata->redcostdata = NULL;
   propdata->redcostnupdates = 0;
   propdata->propgraphnodenumber = -1;
   propdata->lpobjval = -1.0;
   propdata->lpobjval_last = -1.0;
   propdata->isInitialized = FALSE;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolStp)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &(propdata->fixingbounds));
   SCIPfreeMemoryArrayNull(scip, &(propdata->deg2bounded));
   SCIPfreeMemoryArrayNull(scip, &(propdata->deg2bounds));

   if( propdata->propgraph )
      graph_free(scip, &(propdata->propgraph), TRUE);

   if( propdata->redcostdata )
      redcosts_free(scip, &(propdata->redcostdata));

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */


/** fix a variable (corresponding to an edge) to 0 */
SCIP_RETCODE SCIPStpFixEdgeVarTo0(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   SCIP_Bool*            success             /**< could variable be fixed? */
   )
{
   assert(scip && edgevar && success);

   /* not fixed yet? */
   if( SCIPvarGetLbLocal(edgevar) < 0.5 && SCIPvarGetUbLocal(edgevar) > 0.5 )
   {
      SCIP_CALL( SCIPchgVarUb(scip, edgevar, 0.0) );
      *success = TRUE;
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}


/** fix a variable (corresponding to an edge) to 1 */
SCIP_RETCODE SCIPStpFixEdgeVarTo1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   SCIP_Bool*            success             /**< could variable be fixed? */
   )
{
   assert(scip && edgevar && success);

   /* not fixed yet? */
   if( SCIPvarGetLbLocal(edgevar) < 0.5 && SCIPvarGetUbLocal(edgevar) > 0.5 )
   {
      SCIP_CALL( SCIPchgVarLb(scip, edgevar, 1.0) );
      *success = TRUE;
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}


/** return total number of arcs fixed by 'fixedgevar' method of this propagator */
int SCIPStpNfixedEdges(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);

   /* get propagator data */
   assert(SCIPfindProp(scip, "stp") != NULL);
   propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));

   assert(propdata != NULL);

   return (propdata->nfixededges_all);
}


/** checks whether problem has become infeasible at current node
 *  NOTE: we basically check whether all terminals (at given B&B node) are reachable from root,
 *  taking bound changes into account */
SCIP_RETCODE SCIPStpPropCheckForInfeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            probisinfeas        /**< is infeasible? */
)
{
   const GRAPH* orggraph = SCIPprobdataGetGraph2(scip);
   int* nodestate;
   int* edgestate;
   const int nnodes = graph_get_nNodes(orggraph);
   const int nedges = graph_get_nEdges(orggraph);

   assert(probisinfeas && orggraph);

   *probisinfeas = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgestate, nedges) );

   SCIP_CALL( getGraphStatesDirected(scip, orggraph, nodestate, edgestate, probisinfeas) );

   if( *probisinfeas == FALSE )
   {
      SCIP_CALL( trailGraphWithStates(scip, orggraph, nodestate, edgestate, probisinfeas) );
   }

   SCIPfreeBufferArray(scip, &edgestate);
   SCIPfreeBufferArray(scip, &nodestate);

   return SCIP_OKAY;
}


/** gives propagator graph  */
SCIP_RETCODE SCIPStpPropGetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data */
   SCIP_Longint*         graphnodenumber,    /**< point to b&b node for which graph is valid */
   SCIP_Bool*            probisinfeas,       /**< infeasible problem? */
   SCIP_Real*            offset              /**< needed for PC/MW */
   )
{
   SCIP_PROPDATA* propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));
   SCIP_VAR** vars = SCIPprobdataGetVars(scip);
   const GRAPH* orggraph = SCIPprobdataGetGraph2(scip);
   int* nodestate ;
   int* edgestate;
   const int nnodes = graph_get_nNodes(orggraph);
   const int nedges = graph_get_nEdges(orggraph);

   assert(probisinfeas && offset);
   assert(vars && orggraph && propdata);

   *probisinfeas = FALSE;
   *graphnodenumber = -1;
   *graph = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgestate, nedges) );

   if( !propdata->propgraph )
   {
      SCIP_CALL( initPropgraph(scip, orggraph, propdata) );
   }
   else
   {
      SCIP_CALL( updatePropgraph(scip, orggraph, propdata) );
   }

   assert(propdata->propgraph);

   /* note: graph type might change (from PC/MW to RPC/RMW) */
   SCIP_CALL( propgraphApplyBoundchanges(scip, vars, orggraph, nodestate, edgestate, propdata, probisinfeas,
        offset ) );

   if( *probisinfeas )
   {
      SCIPdebugMessage("problem has become infeasible after applying bound changes: terminating! \n");
      *probisinfeas = TRUE;
   }
   else
   {
      SCIP_CALL( propgraphPruneUnconnected(scip, propdata->propgraph, probisinfeas, offset ) );

      if( *probisinfeas )
      {
         SCIPdebugMessage("problem has become infeasible after pruning (terminals not connected): terminating! \n");
         *probisinfeas = TRUE;
      }
   }

   if( FALSE == *probisinfeas )
   {
      assert(graph_valid(scip, propdata->propgraph));
      *graph = propdata->propgraph;
      *graphnodenumber = propdata->propgraphnodenumber;
   }

   graph_path_exit(scip, propdata->propgraph);
   SCIPfreeBufferArray(scip, &edgestate);
   SCIPfreeBufferArray(scip, &nodestate);

   return SCIP_OKAY;
}

/** gives array indicating which nodes are degree-2 bounded */
const SCIP_Bool* SCIPStpPropGet2BoundedArr(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PROPDATA* propdata;
   assert(scip != NULL);

   /* get propagator data */
   assert(SCIPfindProp(scip, "stp") != NULL);
   propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));
   assert(propdata != NULL);

   return propdata->deg2bounded;
}


/** creates the stp propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   /* create redcost propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecStp, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyStp) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeStp) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolStp) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolStp) );

   /* add parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/nwaitingrounds",
         "maximum number of rounds before propagating again",
         &propdata->maxnwaitrounds, FALSE, DEFAULT_MAXNWAITINGROUNDS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/aggressive",
         "should the heuristic be executed at maximum frequeny?",
         &propdata->aggressive, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPcreateRandom(scip, &propdata->randnumgen, 0, TRUE) );

   return SCIP_OKAY;
}

/**@} */
