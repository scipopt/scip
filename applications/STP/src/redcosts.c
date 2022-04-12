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
#include "graph.h"
#include "portab.h"


/** reduced cost result data */
struct reduce_costs_data
{
   SCIP_Real*            lvl_redEdgeCost;        /**< for all levels: reduced costs */
   SCIP_Real*            lvl_rootToNodeDist;     /**< for all levels: shortest path distances from root  */
   PATH*                 lvl_nodeToTermsPaths;   /**< for all levels: paths to nCloseTerms nearest terminals */
   int*                  lvl_nodeToTermsBases;   /**< for all levels: nCloseTerms nearest terminals */
   SCIP_Real*            lvl_cutoff;             /**< for all levels: reduced cost cutoff value or -1.0 if not used */
   SCIP_Real*            lvl_dualBound;          /**< for all levels: dual bound or -1.0 if not used */
   int*                  lvl_redCostRoot;        /**< for all levels: graph root for reduced cost calculation */
   int                   nCloseTerms;        /**< number of close terminals: 1,2, or 3 */
   int                   toplevel;           /**< current top level; 0 <= toplevel <  nLevelsMax     */
   int                   nLevelsMax;         /**< maximum number of levels; >= 1*/
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
};

/**@name Local methods
 *
 * @{
 */


/** initializes reduced costs data structure from given parameter struct */
static
SCIP_RETCODE initFromParams(
   SCIP*                 scip,               /**< SCIP */
   const RCPARAMS*       parameters,         /**< parameters for initialization */
   REDCOST**             redcostdata         /**< data to initialize */
)
{
   REDCOST* rc;
   const int nnodes = parameters->nnodes;
   const int nedges = parameters->nedges;
   const int nCloseTerms = parameters->nCloseTerms;
   const int nLevels = parameters->nLevels;
   const int redCostRoot = parameters->redCostRoot;
   const SCIP_Real cutoff = parameters->cutoff;

   assert(nnodes >= 1 && nedges >= 1);
   assert(nedges % 2 == 0);
   assert(redCostRoot >= 0 || redCostRoot == UNKNOWN);
   assert(GE(cutoff, 0.0) || EQ(cutoff, -1.0));
   assert(nLevels >= 1);
   assert(nCloseTerms == 1 || nCloseTerms == 2 || nCloseTerms == 3);

   SCIP_CALL( SCIPallocMemory(scip, redcostdata) );
   rc = *redcostdata;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_redEdgeCost ), nLevels * nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_rootToNodeDist), nLevels * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_nodeToTermsPaths), nLevels * nCloseTerms * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_nodeToTermsBases), nLevels * nCloseTerms * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_cutoff), nLevels) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_redCostRoot), nLevels) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rc->lvl_dualBound), nLevels) );


   rc->toplevel = 0;
   rc->nLevelsMax = nLevels;
   rc->nCloseTerms = nCloseTerms;
   rc->nnodes = nnodes;
   rc->nedges = nedges;

   for( int i = 0; i < nLevels; i++ )
   {
      rc->lvl_cutoff[i] = -1.0;
      rc->lvl_dualBound[i] = -1.0;
      rc->lvl_redCostRoot[i] = -1;
   }

   rc->lvl_cutoff[0] = cutoff;
   rc->lvl_redCostRoot[0] = redCostRoot;

   return SCIP_OKAY;
}


/** returns top level*/
static inline
int getTopLevel(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(0 <= redcostdata->toplevel && redcostdata->toplevel < redcostdata->nLevelsMax);

   return redcostdata->toplevel;
}


#ifndef NDEBUG
static
SCIP_Bool levelIsValid(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< the level */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return (0 <= level && level <= toplevel);
}
#endif


/** returns start position of given level for lvl_nodeToTermsPaths and lvl_nodeToTermsBases */
static
int getStartPositionCloseTerms(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< the level */
   )
{
   const int nnodes = redcostdata->nnodes;
   const int nCloseTerms = redcostdata->nCloseTerms;

   assert(nnodes >= 1 && nCloseTerms >= 1);
   assert(levelIsValid(redcostdata, level));

   return (nnodes * nCloseTerms * level);
}


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


/** returns reduced costs */
SCIP_Real* redcosts_getEdgeCosts(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   )
{
   const int nedges = redcostdata->nedges;

   assert(levelIsValid(redcostdata, level));
   assert(redcostdata->lvl_redEdgeCost);

   return &(redcostdata->lvl_redEdgeCost[nedges * level]);
}


/** returns top level reduced costs */
SCIP_Real* redcosts_getEdgeCostsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getEdgeCosts(redcostdata, toplevel);
}


/** returns root to node distances */
SCIP_Real* redcosts_getRootToNodeDist(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get distances for */
   )
{
   const int nnodes = redcostdata->nnodes;

   assert(levelIsValid(redcostdata, level));
   assert(redcostdata->lvl_rootToNodeDist);

   return &(redcostdata->lvl_rootToNodeDist[level * nnodes]);
}


/** returns root to node distances for top level */
SCIP_Real* redcosts_getRootToNodeDistTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getRootToNodeDist(redcostdata, toplevel);
}



/** returns paths from nodes to closes terms */
PATH* redcosts_getNodeToTermsPaths(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   )
{
   const int position = getStartPositionCloseTerms(redcostdata, level);
   assert(redcostdata->lvl_nodeToTermsPaths);

   return &(redcostdata->lvl_nodeToTermsPaths[position]);
}


/** returns paths from nodes to closes terms for top level */
PATH* redcosts_getNodeToTermsPathsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getNodeToTermsPaths(redcostdata, toplevel);
}


/** returns closest terminals to nodes */
int* redcosts_getNodeToTermsBases(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to terminals for */
   )
{
   const int position = getStartPositionCloseTerms(redcostdata, level);
   assert(redcostdata->lvl_nodeToTermsBases);

   return &(redcostdata->lvl_nodeToTermsBases[position]);
}


/** returns closest terms to nodes for top level */
int* redcosts_getNodeToTermsBasesTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getNodeToTermsBases(redcostdata, toplevel);
}


/** returns cutoff */
SCIP_Real redcosts_getCutoff(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get cutoff for */
   )
{
   assert(levelIsValid(redcostdata, level));

   return redcostdata->lvl_cutoff[level];
}


/** returns cutoff for top level */
SCIP_Real redcosts_getCutoffTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getCutoff(redcostdata, toplevel);
}


/** returns dual-bound */
SCIP_Real redcosts_getDualBound(
   int                   level,              /**< level */
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   assert(levelIsValid(redcostdata, level));
   assert(redcostdata->lvl_dualBound);

   return redcostdata->lvl_dualBound[level];
}


/** returns dual-bound for top level */
SCIP_Real redcosts_getDualBoundTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getDualBound(toplevel, redcostdata);
}


/** returns root used for reduced cost computation */
int redcosts_getRoot(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get root for */
   )
{
   assert(levelIsValid(redcostdata, level));

   return redcostdata->lvl_redCostRoot[level];
}


/** returns root used for reduced cost computation for to level */
int redcosts_getRootTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return redcosts_getRoot(redcostdata, toplevel);
}


/** returns current (top) level; 0-indexed */
int redcosts_getLevel(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return toplevel;
}


/** returns current number of levels*/
int redcosts_getNlevels(
   const REDCOST*        redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   return toplevel + 1;
}


/** sets cutoff */
void redcosts_setCutoff(
   SCIP_Real           cutoff,             /**< the value */
   int                 level,              /**< level to set cutoff for */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   assert(levelIsValid(redcostdata, level));
   assert(GE(cutoff, 0.0));

   redcostdata->lvl_cutoff[level] = cutoff;
}


/** sets cutoff for top level */
void redcosts_setCutoffTop(
   SCIP_Real           cutoff,             /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   redcosts_setCutoff(cutoff, toplevel, redcostdata);
}


/** sets dual-bound */
void redcosts_setDualBound(
   SCIP_Real           dualbound,          /**< the value */
   int                 level,              /**< level to set dual bound for */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   assert(levelIsValid(redcostdata, level));
   assert(GE(dualbound, 0.0));

   redcostdata->lvl_dualBound[level] = dualbound;
}


/** sets dual-bound */
void redcosts_setDualBoundTop(
   SCIP_Real           dualbound,          /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);
   assert(GE(dualbound, 0.0));

   redcosts_setDualBound(dualbound, toplevel, redcostdata);
}


/** sets root used for reduced cost computation */
void redcosts_setRoot(
   int                   root,               /**< the root */
   int                   level,              /**< level to set dual bound for */
   REDCOST*              redcostdata         /**< reduced costs data */
   )
{
   assert(levelIsValid(redcostdata, level));
   assert(root >= 0);

   redcostdata->lvl_redCostRoot[level] = root;
}


/** sets root used for reduced cost computation */
void redcosts_setRootTop(
   int                   root,               /**< the root */
   REDCOST*              redcostdata         /**< reduced costs data */
   )
{
   const int toplevel = getTopLevel(redcostdata);
   assert(root >= 0);

   redcosts_setRoot(root, toplevel, redcostdata);
}


/** adds a new level */
void redcosts_addLevel(
   REDCOST*              redcostdata         /**< reduced costs data */
   )
{
   assert(redcostdata);
   assert(0 <= redcostdata->toplevel && redcostdata->toplevel < redcostdata->nLevelsMax - 1);

   redcostdata->toplevel++;
}



/* initialize distances from reduced costs */
SCIP_RETCODE redcosts_initializeDistances(
   SCIP*                 scip,               /**< SCIP */
   int                   level,              /**< level to inizialize for*/
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   )
{
   int* pathedge;
   const int daroot = redcosts_getRoot(redcostdata, level);
   const SCIP_Real* const redcosts = redcosts_getEdgeCosts(redcostdata, level);
   PATH* const vnoi = redcosts_getNodeToTermsPaths(redcostdata, level);
   SCIP_Real* const pathdist = redcosts_getRootToNodeDist(redcostdata, level);
   int* const vbase = redcosts_getNodeToTermsBases(redcostdata, level);
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

   /* build Voronoi diagram
    * todo very hacky, should be done properly by the calling method */
   if( directed )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );

      assert(!isRpcmw);
      graph_add1stTermPaths(g, costrev, vnoi, vbase, state);
   }
   else
   {
      const int nCloseTerms = redcostdata->nCloseTerms;
      SCIP_CALL( SCIPallocBufferArray(scip, &state, nCloseTerms * nnodes) );

      if( nCloseTerms == 2 )
      {
         graph_get2nextTermPaths(g, costrev, costrev, vnoi, vbase, state);
      }
      else
      {
         // todo cover case == 1?
         assert(nCloseTerms == 3);
         graph_get3nextTermPaths(g, costrev, costrev, vnoi, vbase, state);
      }

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



/* initialize distances from reduced costs */
SCIP_RETCODE redcosts_initializeDistancesTop(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   )
{
   const int toplevel = getTopLevel(redcostdata);

   assert(g && scip);

   SCIP_CALL( redcosts_initializeDistances(scip, toplevel, g, redcostdata) );

   return SCIP_OKAY;
}


/** initializes reduced costs data structure from given parameter struct */
SCIP_RETCODE redcosts_initFromParams(
   SCIP*                 scip,               /**< SCIP */
   const RCPARAMS*       parameters,         /**< parameters for initialization */
   REDCOST**             redcostdata         /**< data to initialize */
)
{
   assert(scip && parameters && redcostdata);

   SCIP_CALL( initFromParams(scip, parameters, redcostdata) );

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
   RCPARAMS params = { .cutoff = cutoff, .nLevels = 1, .nCloseTerms = 3, .nnodes = nnodes,
                       .nedges = nedges, .redCostRoot = redCostRoot };

   SCIP_CALL( redcosts_initFromParams(scip, &params, redcostdata) );

   return SCIP_OKAY;
}


/** sets cutoff */
void redcosts_setAndReturnCutoffFromBoundTop(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata,        /**< reduced cost data */
  SCIP_Real*            cutoffbound         /**< cutoff */
)
{
   assert(redcostdata && cutoffbound);
   assert(GE(redcosts_getDualBoundTop(redcostdata), 0.0));

   *cutoffbound = upperbound - redcosts_getDualBoundTop(redcostdata);

   assert(GE_FEAS_EPS(*cutoffbound, 0.0, EPSILON));

   if( *cutoffbound < 0.0 )
      *cutoffbound = 0.0;

   redcosts_setCutoffTop(*cutoffbound, redcostdata);
}


/** sets cutoff */
void redcosts_setCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  int                   level,              /**< level */
  REDCOST*              redcostdata         /**< reduced cost data */
)
{
   SCIP_Real cutoff;

   assert(redcostdata);
   assert(levelIsValid(redcostdata, level));
   assert(GE(redcosts_getDualBound(level, redcostdata), 0.0));
   assert(GE(upperbound, 0.0));

   cutoff = upperbound - redcosts_getDualBound(level, redcostdata);

   assert(GE_FEAS_EPS(cutoff, 0.0, EPSILON));

   if( cutoff < 0.0 )
      cutoff = 0.0;

   redcosts_setCutoff(cutoff, level, redcostdata);
}


/** sets cutoff for top level */
void redcosts_setCutoffFromBoundTop(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata         /**< reduced cost data */
)
{
   const int toplevel = getTopLevel(redcostdata);

   redcosts_setCutoffFromBound(upperbound, toplevel, redcostdata);
}


/** increases reduced cost for deleted arcs */
void redcosts_increaseOnDeletedArcs(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   int                   level,              /**< the level */
   REDCOST*              redcostdata         /**< reduced cost data */
)
{
   SCIP_Real* const redEdgeCost = redcosts_getEdgeCosts(redcostdata, level);
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real offset = 2.0 * redcosts_getCutoff(redcostdata, level) + 1.0;

   assert(GE(offset, 1.0));
   assert(nedges == redcostdata->nedges);

   for( int i = 0; i < nedges; i++ )
   {
      if( arcsdeleted[i] )
         redEdgeCost[i] += offset;
   }
}


/** unifies costs  */
void redcosts_unifyBlockedEdgeCosts(
   const GRAPH*          graph,              /**< graph */
   int                   level,              /**< the level */
   REDCOST*              redcostdata         /**< reduced cost data */
)
{
   SCIP_Real* const redEdgeCost = redcosts_getEdgeCosts(redcostdata, level);
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real bound = 2.0 * redcosts_getCutoff(redcostdata, level) + 1.0;

   assert(GE(bound, 1.0));
   assert(LT(bound, FARAWAY));
   assert(nedges == redcostdata->nedges);

   for( int i = 0; i < nedges; i++ )
   {
      if( redEdgeCost[i] > bound )
         redEdgeCost[i] = bound;
   }
}


/** increases reduced cost for deleted arcs for top level */
void redcosts_increaseOnDeletedArcsTop(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   REDCOST*              redcostdata         /**< reduced cost data */
)
{
   const int toplevel = getTopLevel(redcostdata);

   assert(graph && arcsdeleted);

   redcosts_increaseOnDeletedArcs(graph, arcsdeleted, toplevel, redcostdata);
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

   SCIPfreeMemoryArray(scip, &(reddata->lvl_dualBound));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_redCostRoot));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_cutoff));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_nodeToTermsBases));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_nodeToTermsPaths));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_rootToNodeDist));
   SCIPfreeMemoryArray(scip, &(reddata->lvl_redEdgeCost));

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


/** are reduced costs reliable? */
SCIP_Bool redcosts_forLPareReliable(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables (in) */
   const GRAPH*          graph               /**< graph data */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(nedges >= 0);
   assert(vars && scip);

   for( int e = 0; e < nedges; e++ )
   {
      assert(SCIPvarIsBinary(vars[e]));

      /* variable is already fixed? */
      if( SCIPvarGetLbLocal(vars[e]) + 0.5 > SCIPvarGetUbLocal(vars[e]) )
      {
         continue;
      }

      if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[e])) )
      {
         if( SCIPisDualfeasNegative(scip, SCIPgetVarRedcost(scip, vars[e])) )
         {
            return FALSE;
         }
      }
      else
      {
         if( SCIPisDualfeasPositive(scip, SCIPgetVarRedcost(scip, vars[e])) )
         {
            return FALSE;
         }

         if( !(SCIPisFeasEQ(scip, SCIPgetSolVal(scip, NULL, vars[e]), 1.0) || SCIPisDualfeasZero(scip, SCIPgetVarRedcost(scip, vars[e]))) )
         {
            return FALSE;
         }
      }
   }

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

/**@} */
