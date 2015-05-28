/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_stp.c
 * @brief  propagator for Steiner problems, using the LP reduced cost
 * @author Daniel Rehfeldt
 *
 * This propagator makes use of the reduced cost of an optimally solved LP relaxation to propagate the variables
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "prop_stp.h"


/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "stp"
#define PROP_DESC              "stp propagator"
#define PROP_TIMING             SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY          +1000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

/**@} */


/**@name Default parameter values
 *
 * @{
 */


/**@} */


/*
 * Data structures
 */


/** propagator data */
struct SCIP_PropData
{
   SCIP_Real             maxredcost;         /**< maximum reduced cost of a single binary variable */
};


/**@name Local methods
 *
 * @{
 */

static
SCIP_RETCODE fixedgevar(
   SCIP* scip,
   SCIP_VAR*  edgevar,
   int*  nfixed
   )
{
   if( !(SCIPvarGetLbLocal(edgevar) > 0.5 || SCIPvarGetUbLocal(edgevar) < 0.5) )
   {
      (*nfixed)++;

      SCIP_CALL( SCIPchgVarUb(scip, edgevar, 0.0) );
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

   SCIPfreeMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolStp)
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->maxredcost = 0.0;

   return SCIP_OKAY;
}

/** reduced cost propagation method for an LP solution */
static
SCIP_DECL_PROPEXEC(propExecStp)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_VAR*  edgevar;
   GRAPH* graph;
   PATH* vnoi;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real*  pathdist;
   SCIP_Real offset;
   SCIP_Real redcost;
   SCIP_Real lpobjval;
   SCIP_Real cutoffbound;
   SCIP_Real minpathcost;
   int*    vbase;
   int*    pathedge;
   int k;
   int e;
   int nedges;
   int nnodes;
   int nfixed;

   *result = SCIP_DIDNOTRUN;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* only call propagator if current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   /* only call propagator if optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call propagator if current LP is valid relaxation */
   if( !SCIPisLPRelax(scip) )
      return SCIP_OKAY;

   /* cannot apply reduced cost strengthening if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   cutoffbound = SCIPgetPrimalbound(scip);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   offset = SCIPprobdataGetOffset(scip);
   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   nedges = graph->edges;
   nnodes = graph->knots;

   /* get all variables (corresponding to the edges) */
   vars = SCIPprobdataGetVars(scip);

   if( vars == NULL )
      return SCIP_OKAY;

   /* get LP objective value */
   lpobjval = SCIPgetLPObjval(scip);

   if( SCIPisEQ(scip, lpobjval, 0.0) )
      return SCIP_OKAY;
   *result = SCIP_DIDNOTFIND;
   nfixed = 0;

   /* the required reduced path cost to be surpassed */
   minpathcost = cutoffbound - (lpobjval + offset);
   printf("cutoffbound %f, lpobjval %f offset: %f \n", cutoffbound, lpobjval, offset);
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   for( e = 0; e < nedges; e += 2)
   {
      assert(SCIPvarIsBinary(vars[e]));
      assert(SCIPvarIsBinary(vars[e + 1]));
      cost[e] = SCIPgetVarRedcost(scip, vars[e]);
      cost[e + 1] = SCIPgetVarRedcost(scip, vars[e + 1]);
      costrev[e] = cost[e + 1];
      costrev[e + 1] = cost[e];
   }


   for( k = 0; k < nnodes; k++ )
      graph->mark[k] = (graph->grad[k] > 0);
   /* distance from root to all nodes */
   graph_path_execX(scip, graph, graph->source[0], cost, pathdist, pathedge);

   /* no paths should go back to the root */
   for( e = graph->outbeg[graph->source[0]]; e != EAT_LAST; e = graph->oeat[e] )
      costrev[e] = FARAWAY;

   /* build voronoi diagram */
   voronoi_terms(scip, graph, costrev, vnoi, vbase, graph->path_heap, graph->path_state);

   for( k = 0; k < nnodes; k++ )
   {
      if( graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
         if( Is_term(graph->term[k]) )
            continue;
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
      {
         for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            edgevar = vars[e];
            redcost = SCIPgetVarRedcost(scip, edgevar);

	    /* try to fix edge */
            SCIP_CALL( fixedgevar(scip, vars[e], &nfixed) );

	    /* try to fix reversed edge */
	    SCIP_CALL( fixedgevar(scip, vars[flipedge(e)], &nfixed) );
         }
      }
      else
      {
         for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            edgevar = vars[e];
            redcost = SCIPgetVarRedcost(scip, edgevar);
            if( SCIPisGT(scip, pathdist[k] + redcost + vnoi[graph->head[e]].dist, minpathcost) )
               /* try to fix edge */
               SCIP_CALL( fixedgevar(scip, vars[e], &nfixed) );
         }
      }
   }
   printf("fixed: %d \n", nfixed );

   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &pathdist);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);
   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the stp propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create stp propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecStp, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyStp) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolStp) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeStp) );

   return SCIP_OKAY;
}

/**@} */
