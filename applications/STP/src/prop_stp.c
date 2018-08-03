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

/**@file   prop_stp.c
 * @brief  propagator for Steiner tree problems, using the LP reduced costs
 * @author Daniel Rehfeldt
 *
 * This propagator makes use of the reduced cost of an optimally solved LP relaxation to propagate the variables
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "prop_stp.h"
#include "grph.h"
#include "branch_stp.h"
#include "scip/tree.h"


/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "stp"
#define PROP_DESC              "stp propagator"
#define PROP_TIMING            SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY          1000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_STP_EDGE_KILLED -1
#define PROP_STP_EDGE_UNSET   0
#define PROP_STP_EDGE_SET     1
#define PROP_STP_EDGE_FIXED   2

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_MAXNWAITINGROUNDS        2     /**< maximum number of rounds to wait until propagating again */
#define REDUCTION_WAIT_RATIO             0.02  /**< ratio of edges to be newly fixed before performing reductions for additional fixing */

/**@} */

/*
 * Data structures
 */


/** propagator data */
struct SCIP_PropData
{
   GRAPH*                propgraph;          /**< graph data */
   SCIP_Real*            fixingbounds;       /**< saves largest upper bound to each variable that would allow to fix it */
   SCIP_Longint          nfails;             /**< number of failures since last successful call */
   SCIP_Longint          ncalls;             /**< number of calls */
   SCIP_Longint          nlastcall;          /**< number of last call */
   SCIP_Longint          lastnodenumber;     /**< number of last call */
   SCIP_Longint          propgraphnodenumber;/**< b&b node number at which propgraph was updated */
   int                   nfixededges;        /**< total number of arcs fixed by 'fixedgevar' method of this propagator */
   int                   postrednfixededges; /**< total number of arcs fixed by 'fixedgevar' method of this propagator after the last reductions */
   int                   maxnwaitrounds;     /**< maximum number of rounds to wait until propagating again */
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

      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(edgevar), 1.0));

      /* get propagator data */
      assert(SCIPfindProp(scip, "stp") != NULL);
      propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));
      assert(propdata != NULL);

      SCIP_CALL( SCIPchgVarLb(scip, edgevar, 1.0) );
      (*nfixed)++;
      propdata->nfixededges++;
      propdata->postrednfixededges++;
   }
   return SCIP_OKAY;
}
#endif

/** try to make global fixings */
static
SCIP_RETCODE globalfixing(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   int*                  nfixededges,        /**< points to number of fixed edges */
   const SCIP_Real*      fixingbounds,       /**< fixing bounds */
   const GRAPH*          graph,              /**< graph structure */
   SCIP_Real             cutoffbound,        /**> cutoff bound  */
   int                   nedges              /**< number of edges */
)
{
   for( int e = 0; e < nedges; e++ )
   {
      if( SCIPisLT(scip, cutoffbound, fixingbounds[e]) )
      {
         SCIP_VAR* const edgevar = vars[e];

         if( SCIPvarGetLbGlobal(edgevar) < 0.5 && SCIPvarGetUbGlobal(edgevar) > 0.5 )
         {
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(edgevar), 1.0));

            SCIPchgVarUbGlobal(scip, edgevar, 0.0);
            (*nfixededges)++;
         }
      }
   }

   return SCIP_OKAY;
}

static
SCIP_Bool redcostAvailable(
   SCIP*                 scip                /**< SCIP structure */
)
{
   /* only execute dualcostVarfixing if current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
      return FALSE;

   /* only execute dualcostVarfixing if optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return FALSE;

   /* only execute dualcostVarfixing if current LP is valid relaxation */
   if( !SCIPisLPRelax(scip) )
      return FALSE;

   /* we cannot apply reduced cost strengthening if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
      return FALSE;

   /* reduced cost strengthening can only be applied if cutoff is finite */
   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
      return FALSE;

   return TRUE;
}

/* initialize reduced costs*/
static
void setRedcosts(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   int                   nedges,             /**< nedges */
   SCIP_Real*            cost                /**< reduced costs */
)
{
   assert(nedges >= 0);

   for( unsigned int e = 0; e < (unsigned) nedges; e++ )
   {
      assert(SCIPvarIsBinary(vars[e]));

      /* variable is already fixed, we must not trust the reduced cost */
      if( SCIPvarGetLbLocal(vars[e]) + 0.5 > SCIPvarGetUbLocal(vars[e]) )
      {
         if( SCIPvarGetLbLocal(vars[e]) > 0.5 )
            cost[e] = 0.0;
         else
         {
            assert(SCIPvarGetUbLocal(vars[e]) < 0.5);
            cost[e] = FARAWAY;
         }
      }
      else
      {
         if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[e])) )
         {
            assert(!SCIPisDualfeasNegative(scip, SCIPgetVarRedcost(scip, vars[e])));
            cost[e] = SCIPgetVarRedcost(scip, vars[e]);
         }
         else
         {
            assert(!SCIPisDualfeasPositive(scip, SCIPgetVarRedcost(scip, vars[e])));
            assert(SCIPisFeasEQ(scip, SCIPgetSolVal(scip, NULL, vars[e]), 1.0) || SCIPisDualfeasZero(scip, SCIPgetVarRedcost(scip, vars[e])));
            cost[e] = 0.0;
         }
      }

      if( cost[e] < 0.0 )
         cost[e] = 0.0;
   }
}


/* initialize (Voronoi) distances */
static
void setVnoiDistances(
   SCIP*                 scip,               /**< SCIP structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const GRAPH*          graph,              /**< graph data structure */
   PATH*                 vnoi,               /**> Voronoi paths  */
   SCIP_Real*            costrev,            /**< reversed reduced costs */
   SCIP_Real*            pathdist,           /**< path distance */
   int*                  pathedge,           /**< path edge */
   int*                  vbase               /**< Voronoi base */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;

   for( int k = 0; k < nnodes; k++ )
      graph->mark[k] = (graph->grad[k] > 0);

   /* distance from root to all nodes */
   graph_path_execX(scip, graph, graph->source, cost, pathdist, pathedge);

   for( unsigned int e = 0; e < (unsigned) nedges; e++ )
      costrev[e] = cost[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = graph->outbeg[graph->source]; e != EAT_LAST; e = graph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build voronoi diagram */
   graph_voronoiTerms(scip, graph, costrev, vnoi, vbase, graph->path_heap, graph->path_state);
}

/* updates fixing bounds for reduced cost fixings */
static
void updateFixingBounds(
   SCIP_Real*            fixingbounds,       /**< fixing bounds */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**> shortest path distances  */
   const PATH*           vnoi,               /**> Voronoi paths  */
   SCIP_Real             lpobjal             /**> LP objective  */
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

/** dual cost based fixing of variables */
static
SCIP_RETCODE dualcostVarfixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   const GRAPH*          graph               /**< graph data structure */
      )
{
   PATH* vnoi;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real* pathdist;
   const SCIP_Real cutoffbound = SCIPgetCutoffbound(scip);
   const SCIP_Real minpathcost = cutoffbound - lpobjval;
   int* vbase;
   int* pathedge;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(SCIPisGE(scip, minpathcost, 0.0));

   if( propdata->fixingbounds == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->fixingbounds), nedges) );
      for( int i = 0; i < nedges; i++ )
         propdata->fixingbounds[i] = -FARAWAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   /* initialize reduced costs*/
   setRedcosts(scip, vars, nedges, cost);

   /* initialize Voronoi structures */
   setVnoiDistances(scip, cost, graph, vnoi, costrev, pathdist, pathedge, vbase);

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            /* try to fix edge and reversed one */
            SCIP_CALL( fixedgevar(scip, vars[e], nfixed) );
            SCIP_CALL( fixedgevar(scip, vars[flipedge(e)], nfixed) );
         }
      }
      else
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            if( SCIPisGT(scip, pathdist[k] + cost[e] + vnoi[graph->head[e]].dist, minpathcost) )
               SCIP_CALL( fixedgevar(scip, vars[e], nfixed) );
      }
   }

   /* at root? */
   if( SCIPgetDepth(scip) == 0 )
      updateFixingBounds(propdata->fixingbounds, graph, cost, pathdist, vnoi, lpobjval);

   SCIP_CALL( globalfixing(scip, vars, nfixed, propdata->fixingbounds, graph, cutoffbound, nedges) );

   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &pathdist);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/** extended reduced cost reduction */
static
SCIP_RETCODE reduceRedcostExtended(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   GRAPH*                propgraph           /**< graph data structure */
)
{
   PATH* vnoi;
   SCIP_Real* redcost;
   SCIP_Real* redcostrev;
   SCIP_Real* pathdist;
   int* nodearr;
   int* pathedge;
   int* vbase;
   const SCIP_Real cutoffbound = SCIPgetCutoffbound(scip);
   SCIP_Real minpathcost;
   int extnfixed;
   const int nnodes = propgraph->knots;
   const int nedges = propgraph->edges;

   minpathcost = cutoffbound - lpobjval;

   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcostrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr, nnodes) );

   /* initialize reduced costs*/
   setRedcosts(scip, vars, nedges, redcost);

   /* initialize Voronoi structures */
   setVnoiDistances(scip, redcost, propgraph, vnoi, redcostrev, pathdist, pathedge, vbase);

   /* reduce graph */
   extnfixed = reduce_extendedEdge(scip, propgraph, vnoi, redcost, pathdist, NULL, minpathcost, propgraph->source, nodearr, NULL);
   SCIPdebugMessage("extended fixes: %d \n", extnfixed);

   SCIPfreeBufferArray(scip, &nodearr);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathdist);
   SCIPfreeBufferArray(scip, &redcostrev);
   SCIPfreeBufferArray(scip, &redcost);

   if( extnfixed < 0 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/** This methods tries to fix edges by performing reductions on the current graph.
 *  To this end, the already 0-fixed edges are (temporarily) removed from the underlying graph to strengthen the reduction methods. */
static
SCIP_RETCODE redbasedVarfixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   const GRAPH*          g                   /**< graph data structure */
)
{
   GRAPH* propgraph;
   IDX* curr;
   SCIP_Real offset;
   int* remain;
   const int nedges = g->edges;
   const SCIP_Bool pc = (g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG);

   assert(propdata != NULL);
   assert(scip != NULL);
   assert(g != NULL);
   assert(vars != NULL);

   *probisinfeas = FALSE;
   offset = 0.0;

   SCIP_CALL( SCIPallocBufferArray(scip, &remain, nedges) );

   propgraph = propdata->propgraph;
   assert(propgraph != NULL);

   graph_free_history(scip, propgraph);
   graph_free_historyDeep(scip, propgraph);

   /* copy the graph */
   SCIP_CALL( graph_copy_data(scip, g, propgraph) );
   propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   SCIP_CALL( graph_init_history(scip, propdata->propgraph) );

   for( int e = 0; e < nedges; e++ )
      remain[e] = PROP_STP_EDGE_UNSET;

   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      /* e or its anti-parallel edge fixed to one? */
      if( SCIPvarGetLbLocal(vars[e]) > 0.5 || SCIPvarGetLbLocal(vars[erev]) > 0.5 )
      {
         const int tail = propgraph->tail[e];
         const int head = propgraph->head[e];

         graph_knot_chg(propgraph, tail, 0);
         graph_knot_chg(propgraph, head, 0);

         propgraph->cost[e] = 0.0;
         propgraph->cost[erev] = 0.0;

         remain[e] = PROP_STP_EDGE_FIXED;
         remain[erev] = PROP_STP_EDGE_FIXED;
      }

      /* both e and its anti-parallel edge fixed to zero? */
      if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[erev]) < 0.5 )
      {
         assert(SCIPvarGetLbLocal(vars[e]) < 0.5 && SCIPvarGetLbLocal(vars[erev]) < 0.5);

         graph_edge_del(scip, propgraph, e, TRUE);
         remain[e] = PROP_STP_EDGE_KILLED;
         remain[erev] = PROP_STP_EDGE_KILLED;
      }
   }

   if( pc )
      SCIP_CALL( reduce_contractZeroEdges(scip, propgraph, TRUE) );

   /* not at root? */
   if( SCIPgetDepth(scip) > 0 )
   {
      /* then modify the graph according to vertex kills/fixes during branch-and-bound */
      SCIP_CALL( SCIPStpBranchruleApplyVertexChgs(scip, NULL, propgraph) );
   }

   if( !graph_valid(propgraph) )
   {
      *probisinfeas = TRUE;
      goto TERMINATE;
   }

   SCIP_CALL( graph_path_init(scip, propgraph) );


   if( propdata->fixingbounds == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->fixingbounds), nedges) );
      for( int i = 0; i < nedges; i++ )
         propdata->fixingbounds[i] = -FARAWAY;
   }

   /* reduce graph */
   SCIP_CALL( reduceRedcostExtended(scip, lpobjval, vars, propgraph) );

   SCIP_CALL( level0(scip, propgraph) );
   SCIP_CALL( reduceStp(scip, &propgraph, &offset, 2, FALSE, FALSE, FALSE) );

   assert(graph_valid(propgraph));

   graph_path_exit(scip, propgraph);

   /* try to fix edges ... */

   for( int e = 0; e < nedges; e++ )
   {
      if( propgraph->ieat[e] != EAT_FREE )
      {
         assert(propgraph->ieat[flipedge(e)] != EAT_FREE);

         curr = propgraph->ancestors[e];

         while( curr != NULL )
         {
            const int i = curr->index;
            assert(i < nedges);
            assert(remain[i] != PROP_STP_EDGE_KILLED);

            if( remain[i] == PROP_STP_EDGE_UNSET )
            {
               remain[i] = PROP_STP_EDGE_SET;
               remain[flipedge(i)] = PROP_STP_EDGE_SET;
            }

            curr = curr->parent;
         }
      }
   }

   curr = propgraph->fixedges;

   while( curr != NULL )
   {
      const int e = curr->index;
      assert(e < nedges);

      remain[e] = PROP_STP_EDGE_FIXED;
      remain[flipedge(e)] = PROP_STP_EDGE_FIXED;

      curr = curr->parent;
   }

   /* 1-fixed edge to be deleted? */
   for( int e = 0; e < nedges; e++ )
      if( (remain[e] == PROP_STP_EDGE_UNSET || remain[e] == PROP_STP_EDGE_KILLED) && (SCIPvarGetLbLocal(vars[e]) > 0.5) )
      {
         SCIPdebugMessage("1-fixed arc deleted by reduction methods ... can't propagate  \n \n \n");
         goto TERMINATE;
      }

   /* fix to zero and one */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;
      /* edge not set yet? */
      if( remain[e] == PROP_STP_EDGE_UNSET )
      {
         assert(remain[erev] == PROP_STP_EDGE_UNSET);

         SCIP_CALL( fixedgevar(scip, vars[e], nfixed) );
         SCIP_CALL( fixedgevar(scip, vars[erev], nfixed) );
      }
   }

   SCIPdebugMessage("reduction based fixings: %d \n", *nfixed);

TERMINATE:

   SCIPfreeBufferArray(scip, &remain);

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
   SCIP_Real lpobjval;
   int nfixed;
   SCIP_Bool callreduce;

   *result = SCIP_DIDNOTRUN;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get all variables (corresponding to the edges) */
   vars = SCIPprobdataGetVars(scip);

   if( vars == NULL )
      return SCIP_OKAY;

   /* check if all integral variables are fixed */
   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if( !redcostAvailable(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->ncalls++;

   if( propdata->nfails > 0 && (propdata->nlastcall + propdata->maxnwaitrounds >= propdata->ncalls)
     && (propdata->nlastcall + propdata->nfails > propdata->ncalls) )
      return SCIP_OKAY;

   propdata->nlastcall = propdata->ncalls;

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   nfixed = 0;
   *result = SCIP_DIDNOTFIND;

   lpobjval = SCIPgetLPObjval(scip);

   /* call dual cost based variable fixing */
   SCIP_CALL( dualcostVarfixing(scip, lpobjval, vars, propdata, &nfixed, graph) );

   callreduce = FALSE;

   if( graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT )
   {
      const SCIP_Real redratio = ((SCIP_Real) propdata->postrednfixededges ) / (graph->edges);

      if( propdata->propgraph == NULL )
      {
         propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
         SCIP_CALL( graph_copy(scip, graph, &(propdata->propgraph)) );
         SCIP_CALL( graph_init_history(scip, propdata->propgraph) );
         assert(propdata->nfixededges == 0);
      }

      /* in the tree? */
      if( SCIPgetDepth(scip) > 0 )
      {
         SCIP_NODE* const currnode = SCIPgetCurrentNode(scip);
         const SCIP_Longint nodenumber = SCIPnodeGetNumber(currnode);

         if( nodenumber != propdata->lastnodenumber || propdata->aggressive )
         {
            propdata->lastnodenumber = nodenumber;
            callreduce = TRUE;
         }
      }
      /* and root and is ratio of newly fixed edges higher than bound? */
      else if( SCIPisGT(scip, redratio, REDUCTION_WAIT_RATIO) && SCIPgetDepth(scip) == 0 )
      {
         callreduce = TRUE;
      }
#ifdef WITH_UG
   if( propdata->ncalls == 1 && SCIPgetDepth(scip) == 0 )
   {
      SCIPdebugMessage("trigger UG root reductions \n");
      callreduce = TRUE;
   }
#endif
   }

   if( callreduce )
   {
      SCIP_Bool probisinfeas;

      SCIPdebugMessage("use reduction techniques \n");

      /* call reduced cost based based variable fixing */
      SCIP_CALL( redbasedVarfixing(scip, lpobjval, vars, propdata, &nfixed, &probisinfeas, graph) );

      propdata->postrednfixededges = 0;

      if( probisinfeas )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   if( nfixed > 0 )
   {
      SCIPdebugMessage("newly fixed by STP propagator: %d \n", nfixed );

      propdata->nfails = 0;
      propdata->nfixededges += nfixed;
      propdata->postrednfixededges += nfixed;

      *result = SCIP_REDUCEDDOM;

      if( graph->stp_type == STP_SPG || graph->stp_type == STP_RSMT || graph->stp_type == STP_RPCSPG ||
          graph->stp_type == STP_PCSPG || graph->stp_type == STP_DCSTP )
      {
         for( int e = 0; e < graph->edges; e += 2 )
         {
            const int erev = e + 1;

            /* both e and its anti-parallel edge fixed to zero? */
            if( SCIPvarGetUbGlobal(vars[e]) < 0.5 && SCIPvarGetUbGlobal(vars[erev]) < 0.5 && graph->cost[e] < BLOCKED )
            {
               assert(SCIPvarGetLbLocal(vars[e]) < 0.5 && SCIPvarGetLbLocal(vars[erev]) < 0.5);
               if( graph->cost[e] == graph->cost[erev] )
               {
                  graph->cost[e] = BLOCKED;
                  graph->cost[erev] = BLOCKED;
               }
            }
         }
      }
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
   propdata->lastnodenumber = -1;
   propdata->nfixededges = 0;
   propdata->postrednfixededges = 0;
   propdata->fixingbounds = NULL;
   propdata->propgraph = NULL;
   propdata->propgraphnodenumber = -1;

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
   if( propdata->propgraph != NULL )
      graph_free(scip, &(propdata->propgraph), TRUE);

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */


/** fix a variable (corresponding to an edge) to zero */
SCIP_RETCODE fixedgevar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   int*                  nfixed              /**< counter that is incriminated if variable could be fixed */
   )
{
   assert(scip != NULL);

   if( SCIPvarGetLbLocal(edgevar) < 0.5 && SCIPvarGetUbLocal(edgevar) > 0.5 )
   {
      SCIP_CALL( SCIPchgVarUb(scip, edgevar, 0.0) );
      (*nfixed)++;
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

   return (propdata->nfixededges);
}

/** gets propagator graph  */
void SCIPStpPropGetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data */
   SCIP_Longint*         graphnodenumber     /**< point to b&b node for which graph is valid */
   )
{
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);

   /* get propagator data */
   assert(SCIPfindProp(scip, "stp") != NULL);
   propdata = SCIPpropGetData(SCIPfindProp(scip, "stp"));
   assert(propdata != NULL);

   *graph = (propdata->propgraph);
   *graphnodenumber = propdata->propgraphnodenumber;
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

   return SCIP_OKAY;
}

/**@} */
