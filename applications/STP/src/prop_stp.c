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
#include "cons_stp.h"
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
   SCIP_Real*            deg2bounds;         /**< saves largest upper bound to each variable that would allow to set degree 2 constraint */
   SCIP_Bool*            deg2bounded;        /**< maximum degree of vertex is 2? */
   SCIP_Longint          nfails;             /**< number of failures since last successful call */
   SCIP_Longint          ncalls;             /**< number of calls */
   SCIP_Longint          nlastcall;          /**< number of last call */
   SCIP_Longint          nlastnlpiter;       /**< number of last LP iterations */
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

      printf("FIXED TO ONE %d \n", 0);

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


/** gets bound changes specific for PC/MW */
static
void getBoundchangesPcMW(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          propgraph,          /**< graph data structure */
   int*                  nodestate           /**< state */
)
{
   const int nnodes = propgraph->knots;

   assert(graph_pc_isPcMw(propgraph));
   assert(propgraph->extended);

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
            assert(BRANCH_STP_VERTEX_KILLED != nodestate[k]);

            nodestate[k] = BRANCH_STP_VERTEX_TERM;
         }
         /* root->dummyterm edge fixed to 1? Then delete terminal */
         else if( SCIPvarGetLbLocal(vars[root2term]) > 0.5 )
         {
            assert(BRANCH_STP_VERTEX_TERM != nodestate[k]);

            nodestate[k] = BRANCH_STP_VERTEX_KILLED;
         }

         /* term->dummyterm edge fixed to 0? Then delete terminal */
         if( SCIPvarGetUbLocal(vars[pterm2term]) < 0.5 )
         {
            assert(BRANCH_STP_VERTEX_TERM != nodestate[k]);

            nodestate[k] = BRANCH_STP_VERTEX_KILLED;
         }
         /* term->dummyterm edge fixed to 1? Then take terminal */
         else if( SCIPvarGetLbLocal(vars[pterm2term]) > 0.5 )
         {
            assert(BRANCH_STP_VERTEX_KILLED != nodestate[k]);

            nodestate[k] = BRANCH_STP_VERTEX_TERM;
         }
      }
   }
#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      assert(!(BRANCH_STP_VERTEX_TERM == nodestate[k] && graph_pc_knotIsDummyTerm(propgraph, k)));
#endif
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

   for( int k = 0; k < nnodes; k++ )
      nodestate[k] = BRANCH_STP_VERTEX_UNSET;

   for( int e = 0; e < nedges; e++ )
      edgestate[e] = PROP_STP_EDGE_UNSET;

   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      if( pcmw )
      {
         const int tail = propgraph->tail[e];
         const int head = propgraph->head[e];

         if( graph_pc_knotIsDummyTerm(propgraph, tail) || graph_pc_knotIsDummyTerm(propgraph, head) )
            continue;
      }

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
      getBoundchangesPcMW(scip, vars, propgraph, nodestate);
   }

   /* not at root? */
   if( SCIPgetDepth(scip) > 0 )
   {
      SCIP_Bool conflict = FALSE;

      SCIP_CALL( SCIPStpBranchruleGetVertexChgs(scip, nodestate, &conflict) );

      if( conflict )
      {
         *probisinfeas = TRUE;
      }
   }

   return SCIP_OKAY;
}


/** initialize reduced cost distances */
static
SCIP_RETCODE getRedCostDistances(
   SCIP*                 scip,               /**< SCIP structure */
   const SCIP_Real*      redcost,            /**< reduced costs */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**> Voronoi paths  */
   SCIP_Real*            pathdist,           /**< path distance */
   int*                  pathedge,           /**< path edge */
   int*                  vbase,              /**< Voronoi base */
   int*                  state               /**< state  */
)
{
   SCIP_Real* redcostrev;
   const int nedges = g->edges;

   SCIP_CALL( SCIPallocBufferArray(scip, &redcostrev, nedges) );

   graph_mark(g);

   /* distance from root to all nodes */
   graph_path_execX(scip, g, g->source, redcost, pathdist, pathedge);

   for( unsigned int e = 0; e < (unsigned) nedges; e++ )
      redcostrev[e] = redcost[flipedge(e)];

   /* no paths should go back to the root */
   for( int e = g->outbeg[g->source]; e != EAT_LAST; e = g->oeat[e] )
      redcostrev[e] = FARAWAY;

   graph_get3nextTerms(scip, g, redcostrev, redcostrev, vnoi, vbase, g->path_heap, state);

   SCIPfreeBufferArray(scip, &redcostrev);

   return SCIP_OKAY;
}


/** helper method for reduction based variable fixings */
static inline
void setEdgestate(
   const GRAPH*          graph,              /**< graph structure */
   IDX*                  curr,               /**< current ancestor */
   int*                  edgestate           /**< edge state array */
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
   int*                  edgestate           /**< edge state array */
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
   int*                  nfixed              /**< number of fixed edges */
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

         SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nfixed) );
         SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[erev], nfixed) );
      }
   }

   return SCIP_OKAY;
}


/** can edge be excluded? */
static
SCIP_Bool edgestateEdgeIsExcludedPcMw(
   const GRAPH*          graph,              /**< graph structure */
   int                   e,                  /**< edge */
   const int*            nodestate           /**< node state array */
   )
{
   const int tail = graph->tail[e];
   const int head = graph->head[e];

   assert(graph_pc_isPcMw(graph));

   if( graph_pc_knotIsDummyTerm(graph, tail) || graph_pc_knotIsDummyTerm(graph, head) )
      return TRUE;

#if 0
   int dummyterm = -1;

   if( head != graph->source && graph_pc_knotIsDummyTerm(graph, head) )
      dummyterm = head;
   else if( tail != graph->source && graph_pc_knotIsDummyTerm(graph, tail) )
      dummyterm = tail;

   if( dummyterm != -1 )
   {
      const int twin = graph_pc_getTwinTerm(graph, dummyterm);
      assert(graph_pc_knotIsDummyTerm(graph, dummyterm));

      if( BRANCH_STP_VERTEX_TERM == nodestate[twin] || BRANCH_STP_VERTEX_KILLED == nodestate[twin] )
      {
         return TRUE;
      }
   }
#endif

   return FALSE;
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

   *error = FALSE;

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

   fixEdgestate(propgraph, graph_get_fixedges(propgraph), edgestate);

   /* 1-fixed edge to be deleted? */
   for( int e = 0; e < nedges; e++ )
   {
      if( (edgestate[e] == PROP_STP_EDGE_UNSET || edgestate[e] == PROP_STP_EDGE_KILLED) && (SCIPvarGetLbLocal(vars[e]) > 0.5) )
      {
         if( pcmw && edgestateEdgeIsExcludedPcMw(graph, e, nodestate) )
            continue;

         graph_edge_printInfo(propgraph, e);
         printf("1-fixed arc deleted by reduction methods ... can't propagate  \n \n \n");
         *error = TRUE;
         return;
      }
   }

   /* 0-fixed edge been contracted? */
   for( int e = 0; e < nedges; e++ )
   {
      if( edgestate[e] == PROP_STP_EDGE_FIXED && (SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[flipedge(e)]) < 0.5) )
      {
         graph_edge_printInfo(propgraph, e);
         printf("0-fixed arc contracted by reduction methods ... can't propagate  \n \n \n");
         *error = TRUE;
         return;
      }
   }
}


/** updates fixing bounds for reduced cost fixings */
static
void updateEdgeLurkingBounds(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**> shortest path distances  */
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
   const SCIP_Real*      cost,               /**< reduced costs */
   const SCIP_Real*      pathdist,           /**> shortest path distances  */
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
   GRAPH*                propgraph           /**< propagator graph */
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
            graph_pc_knotTofixedTerm(scip, propgraph, node);

            SCIPdebugMessage("...fixed \n");
         }
      }
      else
      {
         graph_pc_enforceNode(scip, propgraph, node);

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

   assert(!pcmw || !graph_pc_knotIsDummyTerm(propgraph, tail));
   assert(!pcmw || !graph_pc_knotIsDummyTerm(propgraph, head));

   assert(SCIPisEQ(scip, propgraph->cost[e], propgraph->cost[erev]));

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
            graph_pc_chgPrize(scip, propgraph, BLOCKED_MINOR, k);
            *hasFixedTerms = TRUE;
         }
      }
   }
}


/** Apply current implications resulting from bound changes to propgraph. */
static
SCIP_RETCODE propgraphApplyImplicationsPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the original graph */
   SCIP_VAR**            vars,               /**< variables */
   GRAPH*                propgraph,          /**< propagator graph */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   SCIP_Real*            offset              /**< pointer to the offset */
   )
{
   const int* verts = SCIPStpGetPcImplVerts(scip);
   const int* start = SCIPStpGetPcImplStarts(scip);

   assert(g && vars && propgraph && nfixed && probisinfeas && offset);
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
                  SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[rootedge], nfixed ) );

                  graph_pc_knotTofixedTerm(scip, propgraph, i);
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

/** Apply current bound changes to propgraph.
 *  Note that the graph type of propgraph might be changed! */
static
SCIP_RETCODE propgraphApplyBoundchanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          g,                  /**< graph data structure */
   int*                  nodestate,          /**< node state (uninitialized) */
   int*                  edgestate,          /**< edge state (uninitialized) */
   GRAPH*                propgraph,          /**< propagator graph */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   SCIP_Bool*            probisinfeas,       /**< is problem infeasible? */
   SCIP_Real*            offset              /**< pointer to the offset */
   )
{
   const int nnodes = graph_get_nNodes(propgraph);
   const int nedges = graph_get_nEdges(propgraph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(propgraph);

   assert(scip && vars && g && edgestate && nodestate && nfixed && probisinfeas && offset);
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

         graph_pc_pcmw2rooted(scip, propgraph, BLOCKED_MINOR, FALSE);

         assert(graph_pc_isRootedPcMw(propgraph));
      }
   }

   SCIP_CALL( graph_path_init(scip, propgraph) );

   for( int e = 0; e < nedges; e += 2 )
   {
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
         propgraphFixNode(scip, k, propgraph);
      }
      else if( BRANCH_STP_VERTEX_KILLED == nodestate[k] )
      {
         propgraphDeleteNode(scip, k, propgraph, offset);
      }
   }

   *probisinfeas = FALSE;


   if( pcmw )
   {
      SCIP_CALL( propgraphApplyImplicationsPcMw(scip, g, vars, propgraph, nfixed, probisinfeas, offset) );
      graph_pc_2trans(scip, propgraph);
   }

   return SCIP_OKAY;
}

/** update the graph from propdata from given graph */
static
SCIP_RETCODE updatePropgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure to use for the update */
   SCIP_PROPDATA*        propdata            /**< propagator data */
)
{
   GRAPH* propgraph = propdata->propgraph;

   assert(propgraph != NULL);

   graph_free_history(scip, propgraph);
   graph_free_historyDeep(scip, propgraph);

   SCIP_CALL( graph_copy_data(scip, graph, propgraph) );

   propgraph->norgmodeledges = propgraph->edges;
   propgraph->norgmodelknots = propgraph->knots;
   propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   assert(!graph_pc_isRootedPcMw(graph) || graph_pc_nFixedTerms(graph) == graph_pc_nFixedTerms(propgraph));

   SCIP_CALL( graph_init_history(scip, propdata->propgraph) );

   return SCIP_OKAY;
}


/** try to make global fixings based on lurking bounds */
static
SCIP_RETCODE fixVarsDualcostLurking(
   SCIP*                 scip,               /**< SCIP structure */
   const SCIP_PROPDATA*  propdata,           /**< propagator data */
   const GRAPH*          graph,              /**< graph structure */
   SCIP_Real             cutoffbound,        /**> cutoff bound  */
   SCIP_VAR**            vars,               /**< variables */
   int*                  nfixededges         /**< points to number of fixed edges */
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
            (*nfixededges)++;
         }
      }
   }

   for( int i = 0; i < nnodes; i++ )
      if( SCIPisLT(scip, cutoffbound, deg2bounds[i]) && !Is_term(graph->term[i]) )
         deg2bounded[i] = TRUE;

   return SCIP_OKAY;
}


/** dual cost based fixing of variables */
static
SCIP_RETCODE fixVarsDualcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   GRAPH*                graph               /**< graph data structure */
)
{
   PATH* vnoi;
   SCIP_Real* redcost;
   SCIP_Real* pathdist;
   const SCIP_Real cutoffbound = SCIPgetCutoffbound(scip);
   const SCIP_Real minpathcost = cutoffbound - lpobjval;
   int* vbase;
   int* state;
   int* pathedge;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(SCIPisGE(scip, minpathcost, 0.0));

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

      /* first call, so we can also fix incoming arcs of root to zero */
      for( int e = graph->inpbeg[graph->source]; e != EAT_LAST; e = graph->ieat[e] )
         SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nfixed) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );

   SCIPStpGetRedcosts(scip, vars, nedges, redcost);
   getRedCostDistances(scip, redcost, graph, vnoi, pathdist, pathedge, vbase, state);

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(graph->term[k]) && SCIPisGT(scip, pathdist[k] + vnoi[k].dist, minpathcost) )
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            /* try to fix edge and reversed one */
            SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nfixed) );
            SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[flipedge(e)], nfixed) );
         }
      }
      else
      {
         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            if( SCIPisGT(scip, pathdist[k] + redcost[e] + vnoi[graph->head[e]].dist, minpathcost) )
               SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], nfixed) );
      }
   }

   /* at root? */
   if( SCIPgetDepth(scip) == 0 )
   {
      updateEdgeLurkingBounds(graph, redcost, pathdist, vnoi, lpobjval, propdata->fixingbounds);
      updateDeg2LurkingBounds(graph, redcost, pathdist, vnoi, lpobjval, propdata->deg2bounds);
   }

   SCIP_CALL( fixVarsDualcostLurking(scip, propdata, graph, cutoffbound, vars, nfixed) );

   SCIPfreeBufferArray(scip, &pathedge);
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
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   GRAPH*                propgraph           /**< graph data structure */
   )
{
   PATH* vnoi;
   SCIP_Real* redcost;
   SCIP_Real* pathdist;
   STP_Bool* marked;
   int* nodearr;
   int* pathedge;
   int* vbase;
   int* state;
   const SCIP_Real cutoffbound = SCIPgetCutoffbound(scip);
   SCIP_Real minpathcost;
   int extnfixed;
   const int nnodes = propgraph->knots;
   const int nedges = propgraph->edges;

   minpathcost = cutoffbound - lpobjval;

   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, nedges) );

   for( int e = 0; e < nedges; e++ )
       if( SCIPvarGetUbLocal(vars[e]) < 0.5 )
          marked[e] = TRUE;
       else
          marked[e] = FALSE;

   SCIPStpGetRedcosts(scip, vars, nedges, redcost);
   getRedCostDistances(scip, redcost, propgraph, vnoi, pathdist, pathedge, vbase, state);

   if( graph_pc_isPcMw(propgraph) )
      graph_pc_2org(scip, propgraph);

   /* reduce graph */
   extnfixed = reduce_extendedEdge(scip, propgraph, vnoi, redcost, pathdist, NULL, minpathcost, propgraph->source, nodearr, marked, TRUE);

   if( graph_pc_isPcMw(propgraph) )
      graph_pc_2trans(scip, propgraph);

   for( int e = 0; e < nedges; e++ )
       if( SCIPvarGetUbLocal(vars[e]) > 0.5 && marked[e] )
       {
          if( propgraph->cost[e] != propgraph->cost[flipedge_Uint(e)] )
          {
             assert(graph_pc_isPcMw(propgraph));
             continue;
          }

          SCIPdebugMessage("fix edge %d to 0 \n", e);

          SCIP_CALL( SCIPStpFixEdgeVar(scip, vars[e], &extnfixed) );
       }

   SCIPdebugMessage("extended fixes: %d \n", extnfixed);

   SCIPfreeBufferArray(scip, &marked);
   SCIPfreeBufferArray(scip, &nodearr);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathdist);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &state);

   assert(extnfixed >= 0);

   return SCIP_OKAY;
}

/** This methods tries to fix edges by performing reductions on the current graph.
 *  To this end, the already 0-fixed edges are (temporarily) removed from the underlying graph to strengthen the reduction methods. */
static
SCIP_RETCODE fixVarsRedbased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             lpobjval,           /**< LP value */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int*                  nfixed,             /**< pointer to number of fixed edges */
   SCIP_Bool*            probisinfeas        /**< is problem infeasible? */
   )
{
   GRAPH* propgraph = NULL;
   int* nodestate = NULL;
   int* edgestate = NULL;
   SCIP_Real offset;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   SCIP_Bool error;

   assert(graph->stp_type != STP_MWCSP);
   assert(graph->stp_type != STP_RMWCSP); // todo implement! especially reduction techniques for rooted

   assert(propdata && scip && vars && probisinfeas && nfixed);
   assert(!graph_pc_isPcMw(graph) || graph->extended);

   SCIPdebugMessage("start redbasedVarfixing, fixings: %d \n", *nfixed);

   *probisinfeas = FALSE;

   /* don't use it for pcmw! */
   offset = -1.0;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodestate, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgestate, nedges) );

   SCIP_CALL( updatePropgraph(scip, graph, propdata) );
   propgraph = propdata->propgraph;

   /* note: graph type might change (from PC/MW to RPC/RMW) */
   SCIP_CALL( propgraphApplyBoundchanges(scip, vars, graph, nodestate, edgestate, propgraph, nfixed, probisinfeas, &offset) );

   if( *probisinfeas )
   {
      SCIPdebugMessage("problem has become infeasible after applying bound changes: terminating! \n");
      goto TERMINATE;
   }

   assert(propgraph->extended);

   if( graph_pc_isRootedPcMw(propgraph) )
      SCIP_CALL( level0RpcRmwInfeas(scip, propgraph, &offset, probisinfeas) );
   else
      SCIP_CALL( level0infeas(scip, propgraph, probisinfeas) );

   if( *probisinfeas )
   {
      SCIPdebugMessage("problem has become infeasible because terminals are not connected: terminating! \n");
      goto TERMINATE;
   }

   if( !graph_valid(scip, propgraph) )
   {
      printf("FAIL: problem in propagation has become invalid! \n");
      return SCIP_ERROR;
   }

   /* now reduce the graph */
   if( graph_pc_isPc(propgraph) )
   {
      SCIP_CALL( fixVarsExtendedRed(scip, lpobjval, vars, propgraph) );
      SCIP_CALL( reducePc(scip, NULL, propgraph, &offset, 2, FALSE, FALSE, FALSE) );
   }
   else
   {
      SCIP_CALL( fixVarsExtendedRed(scip, lpobjval, vars, propgraph) );
      SCIP_CALL( level0(scip, propgraph) );
      SCIP_CALL( reduceStp(scip, propgraph, &offset, 2, FALSE, FALSE, FALSE) );
   }

   assert(graph_valid(scip, propgraph));

   /* mark surviving original edges of propgraph reductions in array edgestate */
   updateEdgestateFromRed(graph, propgraph, vars, nodestate, edgestate, &error);

   if( error )
   {
      return SCIP_ERROR;
   //   goto TERMINATE; // todo check!
   }

#ifdef SCIP_DEBUG
   /* has potentially 0-fixed edge been contracted? */
   for( int e = 0; e < nedges; e++ )
   {
      if( edgestate[e] == PROP_STP_EDGE_FIXED && (SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetLbLocal(vars[flipedge(e)]) < 0.5) )
      {
         graph_edge_printInfo(graph, e);
         printf(" SCIPvarGetLbLocal(vars[flipedge(e)])=%f \n", SCIPvarGetLbLocal(vars[flipedge(e)]));
         printf(" SCIPvarGetUbLocal(vars[flipedge(e)])=%f \n", SCIPvarGetUbLocal(vars[flipedge(e)]));
         printf("potentially 0-fixed arc contracted by reduction methods ... can't propagate  \n \n \n");
         //goto TERMINATE;
         return SCIP_ERROR; // todo is that really an issue?
      }
   }
#endif

   SCIP_CALL( applyEdgestateToProb(scip, graph, vars, edgestate, nfixed) );

   SCIPdebugMessage("reduction based fixings: %d \n", *nfixed);

TERMINATE:
   graph_path_exit(scip, propgraph);
   SCIPfreeBufferArray(scip, &edgestate);
   SCIPfreeBufferArray(scip, &nodestate);

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

   if( !SCIPStpRedcostAvailable(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->ncalls++;

   if( propdata->nfails > 0 && (propdata->nlastcall + propdata->maxnwaitrounds >= propdata->ncalls)
     && (propdata->nlastcall + propdata->nfails > propdata->ncalls) )
      return SCIP_OKAY;

   /* new LP needs to have been solved to avoid conflicts between reduction based and red. cost propagation */
   if( propdata->nlastnlpiter == SCIPgetNLPIterations(scip) )
      return SCIP_OKAY;

   propdata->nlastnlpiter = SCIPgetNLPIterations(scip);
   propdata->nlastcall = propdata->ncalls;

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   nfixed = 0;
   *result = SCIP_DIDNOTFIND;

   lpobjval = SCIPgetLPObjval(scip);

   /* call dual cost based variable fixing */
   SCIP_CALL( fixVarsDualcost(scip, lpobjval, vars, propdata, &nfixed, graph) );

   callreduce = FALSE;

   if( graph_typeIsSpgLike(graph) || graph_pc_isPcMw(graph) )
   {
      const SCIP_Real redratio = ((SCIP_Real) propdata->postrednfixededges ) / (graph->edges);

      if( propdata->propgraph == NULL )
      {
         propdata->propgraphnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
         SCIP_CALL( graph_copy(scip, graph, &(propdata->propgraph)) );
         propdata->propgraph->norgmodeledges = propdata->propgraph->edges;
         propdata->propgraph->norgmodelknots = propdata->propgraph->knots;
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
      /* at root and is ratio of newly fixed edges higher than bound? (== 0 is necessary, could be -1) */
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
      SCIP_Bool probisinfeas = FALSE;

      SCIPdebugMessage("use reduction techniques \n");

      /* call reduced cost based based variable fixing */
      SCIP_CALL( fixVarsRedbased(scip, graph, lpobjval, vars, propdata, &nfixed, &probisinfeas) );

      propdata->postrednfixededges = 0;

      if( probisinfeas )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   if( nfixed > 0 )
   {
      SCIPdebugMessage("newly fixed by STP propagator: %d of (%d) \n", nfixed, propdata->nfixededges);

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
   propdata->nlastnlpiter = 0;
   propdata->lastnodenumber = -1;
   propdata->nfixededges = 0;
   propdata->postrednfixededges = 0;
   propdata->fixingbounds = NULL;
   propdata->deg2bounded = NULL;
   propdata->deg2bounds = NULL;
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
   SCIPfreeMemoryArrayNull(scip, &(propdata->deg2bounded));
   SCIPfreeMemoryArrayNull(scip, &(propdata->deg2bounds));

   if( propdata->propgraph != NULL )
      graph_free(scip, &(propdata->propgraph), TRUE);

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */



/** reduced costs available? */
SCIP_Bool SCIPStpRedcostAvailable(
   SCIP*                 scip                /**< SCIP structure */
)
{
   /* only execute if current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
      return FALSE;

   /* only execute dualcostVarfixing if optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return FALSE;

   /* only execute if current LP is valid relaxation */
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


/** initialize reduced costs*/
void SCIPStpGetRedcosts(
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


/** fix a variable (corresponding to an edge) to zero */
SCIP_RETCODE SCIPStpFixEdgeVar(
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

/** gives propagator graph  */
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

   return SCIP_OKAY;
}

/**@} */
