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

/**@file   redcosts.c
 * @brief  Reduced cost based routines for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various routines for reduced costs based computations
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


//#define SCIP_DEBUG
#include "redcosts.h"
#include "portab.h"


/** reduced cost result data */
struct reduce_costs_data
{
   SCIP_Real*            redEdgeCost;        /**< reduced costs */
   SCIP_Real*            rootToNodeDist;     /**< shortest path distances from root  */
   PATH*                 nodeTo3TermsPaths;  /**< paths to three nearest terminals */
   int*                  nodeTo3TermsBases;  /**< three nearest terminals */
   SCIP_Real             cutoff;             /**< reduced cost cutoff value or -1.0 if not used */
   SCIP_Real             dualBound;          /**< dual bound or -1.0 if not used */
   int                   redCostRoot;        /**< graph root for reduced cost calculation */
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
};

/**@name Local methods
 *
 * @{
 */



/**@} */

/**@name Interface methods
 *
 * @{
 */



/** returns number of nodes for which reduced costs are stored */
int redcosts_getNnodes(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   return redcostdata->nnodes;
}


/** returns number of edges for which reduced costs are stored */
int redcosts_getNedges(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   return redcostdata->nedges;
}



/** returns top level reduced costs */
SCIP_Real* redcosts_getEdgeCostsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(redcostdata->redEdgeCost);

   return redcostdata->redEdgeCost;
}


/** returns root to node distances */
SCIP_Real* redcosts_getRootToNodeDistTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(redcostdata->rootToNodeDist);

   return redcostdata->rootToNodeDist;
}


/** returns paths from nodes to closes terms */
PATH* redcosts_getNodeToTermsPathsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(redcostdata->nodeTo3TermsPaths);

   return redcostdata->nodeTo3TermsPaths;
}


/** returns closest terms to nodes */
int* redcosts_getNodeToTermsBasesTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(redcostdata->nodeTo3TermsBases);

   return redcostdata->nodeTo3TermsBases;
}


/** returns cutoff */
SCIP_Real redcosts_getCutoffTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   return redcostdata->cutoff;
}


/** returns dual-bound */
SCIP_Real redcosts_getDualBoundTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   return redcostdata->dualBound;
}


/** returns root used for reduced cost computation */
int redcosts_getRootTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   return redcostdata->redCostRoot;
}



/** sets cutoff */
void redcosts_setCutoffTop(
   SCIP_Real           cutoff,             /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(GE(cutoff, 0.0));

   redcostdata->cutoff = cutoff;
}


/** sets dual-bound */
void redcosts_setDualBoundTop(
   SCIP_Real           dualbound,          /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);

   redcostdata->dualBound = dualbound;
}



/** sets root used for reduced cost computation */
void redcosts_setRootTop(
   int                   root,               /**< the root */
   REDCOST*              redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(root >= 0);


   redcostdata->redCostRoot = root;
}



/* initialize distances from reduced costs */
SCIP_RETCODE redcosts_initializeDistances(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   )
{
   int* pathedge;
   const int daroot = redcostdata->redCostRoot;
   const SCIP_Real* const redcosts = redcosts_getEdgeCostsTop(redcostdata);
   PATH* const vnoi = redcosts_getNodeToTermsPathsTop(redcostdata);
   SCIP_Real* const pathdist = redcosts_getRootToNodeDistTop(redcostdata);
   int* const vbase = redcostdata->nodeTo3TermsBases;
   SCIP_Real* costrev = NULL;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const SCIP_Bool isRpcmw = graph_pc_isRootedPcMw(g);
   const SCIP_Bool directed = (g->stp_type == STP_SAP || g->stp_type == STP_NWSPG);
   int* state;

   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes + 1) );

   /* distance from root to all nodes */
   graph_path_execX(scip, g, daroot, redcosts, pathdist, pathedge);

   for( int e = 0; e < nedges; e++ )
      costrev[e] = redcosts[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = g->outbeg[daroot]; e != EAT_LAST; e = g->oeat[e] )
      costrev[e] = FARAWAY;

   if( isRpcmw )
   {
      if( !g->extended )
         graph_pc_2trans(scip, g);
      else
         graph_mark(g);
   }

   assert(graph_isMarked(g));

   /* build Voronoi diagram */
   if( directed )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );

      assert(!isRpcmw);
      graph_add1stTermPaths(g, costrev, vnoi, vbase, state);
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );

      graph_get3nextTermPaths(g, costrev, costrev, vnoi, vbase, state);

#ifndef NDEBUG
      {
         for( int i = 0; i < nnodes; i++ )
         {
            if( !g->mark[i] )
               continue;

            if( !Is_term(g->term[i]) )
            {
               assert(vbase[i] != daroot || vnoi[i].dist >= FARAWAY);
               assert(vbase[i + nnodes] != daroot || vnoi[i + nnodes].dist >= FARAWAY);
            }
            else
               assert(vbase[i] == i);
         }
      }
#endif
   }

   if( isRpcmw )
      graph_pc_2org(scip, g);

   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &costrev);

   return SCIP_OKAY;
}



/** initializes reduced costs data structure */
SCIP_RETCODE redcosts_init(
   SCIP*                 scip,               /**< SCIP */
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value or -1.0 if not used */
   int                   redCostRoot,        /**< graph root for reduced cost calculation */
   REDCOST**             redcostdata         /**< data to initialize */
)
{
   REDCOST* reddata;
   SCIP_Real* redEdgeCost;
   SCIP_Real* rootToNodeDist;
   PATH* nodeTo3TermsPaths;
   int* nodeTo3TermsBases;

   assert(scip);
   assert(nnodes >= 0);
   assert(nedges >= 0);
   assert(nedges % 2 == 0);
   assert(redCostRoot >= 0 || redCostRoot == UNKNOWN);
   assert(GE(cutoff, 0.0) || EQ(cutoff, -1.0));

   SCIP_CALL( SCIPallocMemory(scip, redcostdata) );
   reddata = *redcostdata;

   SCIP_CALL( SCIPallocMemoryArray(scip, &redEdgeCost, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &rootToNodeDist, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsPaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsBases, 3 * nnodes) );

   reddata->redEdgeCost = redEdgeCost;
   reddata->rootToNodeDist = rootToNodeDist;
   reddata->nodeTo3TermsPaths = nodeTo3TermsPaths;
   reddata->nodeTo3TermsBases = nodeTo3TermsBases;
   reddata->cutoff = cutoff;
   reddata->redCostRoot = redCostRoot;

   reddata->nnodes = nnodes;
   reddata->nedges = nedges;

   return SCIP_OKAY;
}



/** sets cutoff */
void redcosts_setAndReturnCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata,        /**< reduced cost data */
  SCIP_Real*            cutoffbound         /**< cutoff */
)
{
   assert(redcostdata && cutoffbound);

   *cutoffbound = upperbound - redcosts_getDualBoundTop(redcostdata);

   assert(GE_FEAS_EPS(*cutoffbound, 0.0, EPSILON));

   if( *cutoffbound < 0.0 )
      *cutoffbound = 0.0;

   redcostdata->cutoff = *cutoffbound;
}


/** sets cutoff */
void redcosts_setCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata         /**< reduced cost data */
)
{
   assert(redcostdata);

   redcostdata->cutoff = upperbound - redcosts_getDualBoundTop(redcostdata);

   assert(GE_FEAS_EPS(redcostdata->cutoff, 0.0, EPSILON));

   if( redcostdata->cutoff < 0.0 )
      redcostdata->cutoff = 0.0;
}


/** increases reduced cost for deleted arcs */
void redcosts_increaseOnDeletedArcs(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   REDCOST*              redcostdata         /**< reduced cost data */
)
{
   SCIP_Real* redEdgeCost = redcosts_getEdgeCostsTop(redcostdata);
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real offset = 2.0 * redcosts_getCutoffTop(redcostdata) + 1.0;

   assert(GE(offset, 1.0));

   for( int i = 0; i < nedges; i++ )
   {
      if( arcsdeleted[i] )
         redEdgeCost[i] += offset;
   }
}


/** frees */
void redcosts_free(
   SCIP*                 scip,               /**< SCIP */
   REDCOST**             redcostdata         /**< data */
)
{
   REDCOST* reddata;

   assert(scip && redcostdata);

   reddata = *redcostdata;

   SCIPfreeMemoryArray(scip, &(reddata->nodeTo3TermsBases));
   SCIPfreeMemoryArray(scip, &(reddata->nodeTo3TermsPaths));
   SCIPfreeMemoryArray(scip, &(reddata->rootToNodeDist));
   SCIPfreeMemoryArray(scip, &(reddata->redEdgeCost));

   SCIPfreeMemory(scip, redcostdata);
}


/** reduced costs available? */
SCIP_Bool redcosts_forLPareAvailable(
   SCIP*                 scip                /**< SCIP structure */
)
{
   /* only execute if current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
   {
      SCIPdebugMessage("!SCIPhasCurrentNodeLP \n");
      return FALSE;
   }

   /* only execute dualcostVarfixing if optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL (%d) \n", SCIPgetLPSolstat(scip));
      return FALSE;
   }

   /* only execute if current LP is valid relaxation */
   if( !SCIPisLPRelax(scip) )
   {
      SCIPdebugMessage("!SCIPisLPRelax \n");
      return FALSE;
   }

   /* we cannot apply reduced cost strengthening if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
   {
      SCIPdebugMessage("!SCIPisLPSolBasic \n");
      return FALSE;
   }

   /* reduced cost strengthening can only be applied if cutoff is finite */
   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
   {
      SCIPdebugMessage("!SCIPgetCutoffbound \n");
      return FALSE;
   }

   SCIPdebugMessage("reduced costs are available! \n");


   return TRUE;
}


/** initialize reduced costs*/
void redcosts_forLPget(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables (in) */
   const GRAPH*          graph,              /**< graph data */
   SCIP_Real*            redcosts            /**< reduced costs (out) */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(nedges >= 0);
   assert(vars && redcosts && scip);

   for( int e = 0; e < nedges; e++ )
   {
      assert(SCIPvarIsBinary(vars[e]));

      /* variable is already fixed, we must not trust the reduced cost */
      if( SCIPvarGetLbLocal(vars[e]) + 0.5 > SCIPvarGetUbLocal(vars[e]) )
      {
         if( SCIPvarGetLbLocal(vars[e]) > 0.5 )
            redcosts[e] = 0.0;
         else
         {
            assert(SCIPvarGetUbLocal(vars[e]) < 0.5);
            redcosts[e] = FARAWAY;
         }
      }
      else
      {
         if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[e])) )
         {
            assert(!SCIPisDualfeasNegative(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = SCIPgetVarRedcost(scip, vars[e]);
         }
         else
         {
            assert(!SCIPisDualfeasPositive(scip, SCIPgetVarRedcost(scip, vars[e])));
            assert(SCIPisFeasEQ(scip, SCIPgetSolVal(scip, NULL, vars[e]), 1.0) || SCIPisDualfeasZero(scip, SCIPgetVarRedcost(scip, vars[e])));
            redcosts[e] = 0.0;
         }
      }

      if( redcosts[e] < 0.0 )
         redcosts[e] = 0.0;
   }
#ifdef SCIP_DISABLED_CODE
   if( graph_pc_isPcMw(graph) )
   {
      /* we do some clean-up */
      const int nnodes = graph_get_nNodes(graph);
      const int root = graph->source;

      assert(graph->term2edge);

      for( int i = 0; i < nnodes; i++ )
      {
         if( graph_pc_knotIsDummyTerm(graph, i) && i != root )
         {
            const int edge2dummy = flipedge(graph->term2edge[i]);

            assert(edge2dummy >= 0 && graph->head[edge2dummy] == i);
            assert(EQ(graph->cost[edge2dummy], 0.0));
            assert(EQ(redcosts[edge2dummy], 0.0) || EQ(redcosts[edge2dummy], FARAWAY));

            if( EQ(redcosts[edge2dummy], 0.0) )
               redcosts[edge2dummy] = 0.0;
         }
      }
   }
#endif
}
