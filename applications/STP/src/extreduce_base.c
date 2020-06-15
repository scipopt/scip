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

/**@file   extreduce_base.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements interface methods for extended reduction techniques for several Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "solstp.h"
#include "extreduce.h"

#define EXT_PSEUDO_DEGREE_MIN 3
#define EXT_PSEUDO_DEGREE_MAX 5

#ifndef NDEBUG
/** all good with the graph->mark array? */
static
SCIP_Bool graphmarkIsClean(
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const GRAPH*          graph              /**< graph data structure */
)
{
   const int nnodes = graph_get_nNodes(graph);

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph->grad[k] == 0 && k != redcostdata->redCostRoot && !Is_term(graph->term[k]) )
      {
         if( graph->mark[k] )
         {
            graph_knot_printInfo(graph, k);
            return FALSE;
         }
      }
   }

   return TRUE;
}
#endif

/** initialize */
static inline
SCIP_RETCODE extInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   DISTDATA*             distdata,           /**< distance data (out) */
   EXTPERMA*             extpermanent        /**< permanent extension data (out) */
)
{
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, FALSE, distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeletable, extpermanent) );

   return SCIP_OKAY;
}


/** free */
static inline
void extFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   DISTDATA*             distdata,           /**< distance data (in/out) */
   EXTPERMA*             extpermanent        /**< permanent extension data (in/out) */
)
{
   extreduce_extPermaFreeMembers(scip, extpermanent);
   extreduce_distDataFreeMembers(scip, graph, distdata);
   graph_free_dcsr(scip, graph);
}


/** initializes new star */
static inline
SCIP_RETCODE pseudodeleteInitStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure  */
   int                   node,               /**< node */
   STAR*                 stardata,           /**< star */
   int**                 compedges,          /**< to be allocated */
   int**                 extleaves           /**< to be allocated */
)
{
   const int degree = g->grad[node];

   reduce_starReset(g, node, stardata);

   SCIP_CALL( SCIPallocBufferArray(scip, compedges, degree) );
   SCIP_CALL( SCIPallocBufferArray(scip, extleaves, degree - 1) );

   return SCIP_OKAY;
}


/** frees new star */
static inline
void pseudodeleteFreeStar(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 compedges,          /**< to be freed */
   int**                 extleaves           /**< to be freed */
)
{
   SCIPfreeBufferArray(scip, extleaves);
   SCIPfreeBufferArray(scip, compedges);
}


/** gets next star */
static inline
void pseudodeleteGetNextStar(
   const GRAPH*          g,                  /**< graph data structure  */
   STAR*                 stardata,           /**< star */
   EXTCOMP*              extcomp             /**< to be filled */
)
{
   int nstaredges;
   int* const compedges = extcomp->compedges;
   int* const extleaves = extcomp->extleaves;
   const int* staredges = reduce_starGetNext(stardata, &(nstaredges));

   extcomp->ncompedges = nstaredges;
   extcomp->nextleaves = nstaredges - 1;
   compedges[0] = flipedge(staredges[0]);
   extcomp->comproot = g->head[staredges[0]];

   for( int i = 1; i < nstaredges; i++ )
   {
      const int outedge = staredges[i];
      compedges[i] = outedge;
      extleaves[i - 1] = g->head[outedge];
   }

   assert(extcomp->ncompedges >= 3);
   assert(extcomp->comproot >= 0);
}


/** frees new star */
static inline
SCIP_Bool pseudodeleteAllStarsChecked(
   const STAR*           stardata            /**< star */

)
{
   return reduce_starAllAreChecked(stardata);
}


/** is node a good candidate for pseudo deletion? */
static inline
SCIP_Bool pseudodeleteNodeIsPromising(
   const GRAPH*          g,                  /**< graph data structure  */
   int                   node                /**< node */
)
{
   const SCIP_Bool pc = graph_pc_isPc(g);
   int degree = -1;

   assert(node >= 0);
   assert(!g->extended);

   if( pc )
   {
      if( !g->mark[node] || graph_pc_knotIsFixedTerm(g, node) )
         return FALSE;

      if( Is_term(g->term[node]) && !graph_pc_termIsNonLeafTerm(g, node) )
         return FALSE;

      degree = graph_pc_realDegree(g, node, FALSE);

      if( Is_term(g->term[node]) && degree != 3 )
         return FALSE;
   }
   else
   {
      if( Is_term(g->term[node]) )
         return FALSE;

      degree = g->grad[node];
   }

   assert(degree >= 0);

   if( degree < EXT_PSEUDO_DEGREE_MIN || degree > EXT_PSEUDO_DEGREE_MAX  )
      return FALSE;

   return TRUE;
}


/** Extended reduction test for pseudo-eliminating nodes.
 *  Applies pseudo-elimination if pseudoDelNodes != NULL, or just marks them otherwise */
static
SCIP_RETCODE pseudodeleteExecute(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   SCIP_Bool*            pseudoDelNodes,     /**< pseudo deletable nodes or NULL (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Bool* pseudoDeletable = NULL;
   STAR* stardata;
   DISTDATA distdata;
   EXTPERMA extpermanent;
   const SCIP_Bool markPseudoDeletion = (pseudoDelNodes != NULL);
   const SCIP_Bool isExtendedOrg = graph->extended;

   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);
   assert(graph_isMarked(graph));

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   if( graph_pc_isPc(graph) )
      graph_pc_2orgcheck(scip, graph);

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );
   if( !markPseudoDeletion )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pseudoDeletable, nnodes) );
   }

   SCIP_CALL( reduce_starInit(scip, EXT_PSEUDO_DEGREE_MAX, &stardata) );

   /* pseudo-elimination loop */
   for( int i = 0; i < nnodes; ++i )
   {
      SCIP_Bool nodeisDeletable = FALSE;

      /* node already pseudo-deletable? */
      if( markPseudoDeletion && pseudoDelNodes[i] )
         continue;

      if( pseudodeleteNodeIsPromising(graph, i) )
      {
         extpermanent.redcostEqualAllow = (result && !solstp_containsNode(graph, result, i));

         SCIP_CALL( extreduce_checkNode(scip, graph, redcostdata, i, stardata, &distdata, &extpermanent, &nodeisDeletable) );
      }

    //  if( nodeisDeletable )
    //  printf("mark %d deg=%d \n", i, graph->grad[i]);

      if( markPseudoDeletion && nodeisDeletable )
         pseudoDelNodes[i] = TRUE;

      if( !markPseudoDeletion )
         pseudoDeletable[i] = nodeisDeletable;

   }

   // todo: if single edges are ledge, try to eliminate via extended reduction!
   // todo: fill cuttoff value from SD...probably would be good to give both SD and DA dists to method! need extra struct.
   // todo better give a function for checking! can generate SDs and DAs on the fly!  Algorithm pattern!

   /* apply pseudo-eliminations? */
   if( !markPseudoDeletion )
   {
      assert(!pseudoDelNodes);
      assert(nelims);

      SCIP_CALL( reduce_applyPseudoDeletions(scip, redcostdata, pseudoDeletable, graph, offsetp, nelims) );
   }

   reduce_starFree(scip, &stardata);
   SCIPfreeBufferArrayNull(scip, &pseudoDeletable);
   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   if( graph_pc_isPc(graph) && isExtendedOrg != graph->extended )
      graph_pc_2trans(scip, graph);

   return SCIP_OKAY;
}


/** deletes an edge and makes corresponding adaptations */
static inline
void removeEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata            /**< distance data (in/out) */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("removing edge ");
   graph_edge_printInfo(graph, edge);
#endif

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));

   graph_edge_delFull(scip, graph, edge, TRUE);
   extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

   if( graph->grad[tail] == 0 )
   {
      if( Is_term(graph->term[tail])  )
      {
         assert(graph_pc_isPcMw(graph) || tail == graph->source);
      }
      else
      {
         graph->mark[tail] = FALSE;
      }
   }

   if( graph->grad[head] == 0 )
   {
      if( Is_term(graph->term[head]) || head == graph->source )
      {
         assert(graph_pc_isPcMw(graph));
      }
      else
      {
         graph->mark[head] = FALSE;
      }
   }

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));
}


/** Extended reduction test for arcs.
 * This method will also set edgedeletable[a] to TRUE if arc 'a' can be deleted, but its anti-parallel arc not. */
SCIP_RETCODE extreduce_deleteArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         extpermanent.redcostEqualAllow = allowequality;

         if( !edgedeletable[e] )
         {
            SCIP_Bool deletable;
            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, e, &distdata, &extpermanent,
                  &deletable) );

            if( deletable )
               edgedeletable[e] = TRUE;
         }

         if( !edgedeletable[erev] )
         {
            SCIP_Bool erevdeletable = FALSE;

            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, erev, &distdata, &extpermanent,
                  &erevdeletable) );

            if( erevdeletable )
               edgedeletable[erev] = TRUE;
         }

         if( edgedeletable[e] && edgedeletable[erev] )
         {
            assert(edgedeletable[e] && edgedeletable[erev]);

            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** extended reduction test for edges */
SCIP_RETCODE extreduce_deleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         extpermanent.redcostEqualAllow = allowequality;

         // todo acually utilize the allow equality!

         SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, e, &distdata, &extpermanent, &deletable) );

         if( deletable )
         {
            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** extended reduction test for marking pseudo-eliminable nodes */
SCIP_RETCODE extreduce_updatePseudoDeletableNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   SCIP_Bool*            pseudoDelNodes,     /**< array marking the pseudo-deletable nodes (in/out)  */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable       /**< edge array to mark which (directed) edge can be removed (in/out) */
)
{

   assert(scip && redcostdata && graph);

   SCIP_CALL( pseudodeleteExecute(scip, redcostdata, result, graph, pseudoDelNodes, edgedeletable, NULL, NULL) );

   return SCIP_OKAY;
}


/** extended reduction test for pseudo-eliminating nodes */
SCIP_RETCODE extreduce_pseudoDeleteNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   assert(scip && redcostdata && graph && nelims);

   SCIP_CALL( pseudodeleteExecute(scip, redcostdata, result, graph, NULL, edgedeletable, offsetp, nelims) );

   return SCIP_OKAY;
}


/** check (directed) arc */
SCIP_RETCODE extreduce_checkArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* isterm = extpermanent->isterm;
   const int root = redcostdata->redCostRoot;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real* rootdist = redcostdata->rootToNodeDist;
   const PATH* nodeToTermpaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real cutoff = redcostdata->cutoff;
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   const SCIP_Real edgebound = redcost[edge] + rootdist[tail] + nodeToTermpaths[head].dist;
   SCIP_Bool restoreAntiArcDeleted = FALSE;
   STP_Bool* const edgedeleted = extpermanent->edgedeleted;

   assert(scip && graph && redcost && rootdist && nodeToTermpaths && distdata);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);
   assert(graph_isMarked(graph));
   assert(extreduce_extPermaIsClean(graph, extpermanent));

   /* trivial rule-out? */
   if( SCIPisGT(scip, edgebound, cutoff) || (extpermanent->redcostEqualAllow && SCIPisEQ(scip, edgebound, cutoff)) || head == root )
   {
      *edgeIsDeletable = TRUE;
      return SCIP_OKAY;
   }

   if( edgedeleted && !edgedeleted[flipedge(edge)] )
   {
      edgedeleted[flipedge(edge)] = TRUE;
      restoreAntiArcDeleted = TRUE;
   }

   *edgeIsDeletable = FALSE;

   /* can we extend from head of 'edge'? */
   if( extLeafIsExtendable(graph, isterm, head) )
   {
      int comphead = graph->head[edge];
      int compedge = edge;
      EXTCOMP extcomp = { .compedges = &compedge, .extleaves = &(comphead),
         .nextleaves = 1, .ncompedges = 1, .comproot = graph->tail[edge],
         .allowReversion = FALSE };

      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, distdata, extpermanent, edgeIsDeletable) );
   }

   if( restoreAntiArcDeleted )
   {
      assert(edgedeleted);
      edgedeleted[flipedge(edge)] = FALSE;
   }

   return SCIP_OKAY;
}


/** check edge */
SCIP_RETCODE extreduce_checkEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* const isterm = extpermanent->isterm;

   assert(scip && graph && edgeIsDeletable && distdata && extpermanent);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph_isMarked(graph));
   assert(extreduce_extPermaIsClean(graph, extpermanent));

   *edgeIsDeletable = FALSE;

   /* is any extension possible? */
   if( extLeafIsExtendable(graph, isterm, graph->tail[edge]) || extLeafIsExtendable(graph, isterm, graph->head[edge]) )
   {
      int comphead = graph->head[edge];
      int compedge = edge;
      EXTCOMP extcomp = { .compedges = &compedge, .extleaves = &(comphead),
                          .nextleaves = 1, .ncompedges = 1,
                          .comproot = graph->tail[edge], .allowReversion = TRUE };

      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, distdata, extpermanent, edgeIsDeletable) );
   }

   return SCIP_OKAY;
}


/** check node for possible  */
SCIP_RETCODE extreduce_checkNode(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   node,               /**< the node */
   STAR*                 stardata,           /**< star */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            isPseudoDeletable   /**< is node pseudo-deletable? */
)
{
   int* compedges;
   int* extleaves;

   assert(scip && graph && redcostdata && distdata && extpermanent && isPseudoDeletable);
   assert(node >= 0 && node < graph->knots);

   *isPseudoDeletable = TRUE;

   SCIP_CALL( pseudodeleteInitStar(scip, graph, node, stardata, &compedges, &extleaves) );

   /* check all stars of degree >= 3, as long as they can be ruled-out */
   while ( *isPseudoDeletable )
   {
      EXTCOMP extcomp = { .compedges = compedges, .extleaves = extleaves,
                          .nextleaves = -1, .ncompedges = -1,
                          .comproot = -1, .allowReversion = TRUE };

      pseudodeleteGetNextStar(graph, stardata, &extcomp);

      *isPseudoDeletable = FALSE;
      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, distdata, extpermanent, isPseudoDeletable) );

      if( pseudodeleteAllStarsChecked(stardata) )
         break;
   }


   // todo check pseudo-deletable edges!


   pseudodeleteFreeStar(scip, &compedges, &extleaves);


   return SCIP_OKAY;
}


/** deletes an edge and makes corresponding adaptations */
void extreduce_edgeRemove(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata            /**< distance data (in/out) */
)
{
   removeEdge(scip, edge, graph, distdata);
}


/** is the edge valid? */
SCIP_Bool extreduce_edgeIsValid(
   const GRAPH*          graph,              /**< graph data structure */
   int                   e                   /**< edge to be checked */
)
{
   if( EAT_FREE == graph->oeat[e] )
   {
      return FALSE;
   }
   else if( graph_pc_isPcMw(graph) )
   {
      const int tail = graph->tail[e];
      const int head = graph->head[e];

      if( (!graph->mark[tail] || !graph->mark[head]) )
      {
         assert(graph_pc_knotIsDummyTerm(graph, tail) || graph_pc_knotIsDummyTerm(graph, head));

         return FALSE;
      }

      assert(!graph_pc_knotIsDummyTerm(graph, tail));
      assert(!graph_pc_knotIsDummyTerm(graph, head));
   }

   return TRUE;
}


/** recompute costs and reduced costs for current tree */
void extreduce_treeRecompCosts(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
#ifndef NDEBUG
   const int tree_nDelUpArcs = extdata->tree_nDelUpArcs;
#endif
   REDDATA* const reddata = extdata->reddata;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   SCIP_Real tree_cost = 0.0;
   SCIP_Real tree_redcost = 0.0;
   const SCIP_Real* const cost = graph->cost;
   const SCIP_Real* const redcost = reddata->redCosts;
   const int* const tree_edges = extdata->tree_edges;
   const int tree_nedges = extdata->tree_nedges;
   const SCIP_Bool isPc = (graph->prize != NULL);

   extdata->tree_nDelUpArcs = 0;

   assert(!extreduce_treeIsFlawed(scip, graph, extdata));
   assert(isPc == graph_pc_isPc(graph));

   for( int i = 0; i < tree_nedges; i++ )
   {
      const int edge = tree_edges[i];
      const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);

      assert(edge >= 0 && edge < graph->edges);

      tree_cost += cost[edge];

      if( !edgeIsDeleted )
      {
         tree_redcost += redcost[edge];
         assert(LT(tree_redcost, FARAWAY));
      }
      else
      {
         extdata->tree_nDelUpArcs++;
      }
   }

   assert(SCIPisEQ(scip, tree_cost, extdata->tree_cost));
   assert(SCIPisEQ(scip, tree_redcost, extdata->tree_redcost));
   assert(tree_nDelUpArcs == extdata->tree_nDelUpArcs);

   extdata->tree_cost = tree_cost;
   extdata->tree_redcost = tree_redcost;

   if( isPc )
   {
      const int* const innerNodes = extdata->tree_innerNodes;
      const SCIP_Real* const prizes = graph->prize;
      SCIP_Real tree_innerPrize = 0.0;
      const int ninnnerNodes = extdata->tree_ninnerNodes;

      for( int i = 0; i < ninnnerNodes; ++i )
      {
         const int node = innerNodes[i];
         tree_innerPrize += prizes[node];
      }

      assert(EQ(tree_innerPrize, extdata->pcdata->tree_innerPrize));

      extdata->pcdata->tree_innerPrize = tree_innerPrize;
   }
}

/** get maximum allowed stack size */
int extreduce_getMaxStackSize(void)
{
   return STP_EXT_MAXSTACKSIZE;
}


/** get maximum allowed number of components */
int extreduce_getMaxStackNcomponents(
   const GRAPH*          graph               /**< graph data structure */
)
{
	assert(graph);

	return STP_EXT_MAXNCOMPS;
}


/** get maximum allowed depth for extended tree in given graph */
int extreduce_getMaxTreeDepth(
   const GRAPH*          graph               /**< graph data structure */
)
{
   const int maxdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;

   assert(maxdepth > 0);

   return maxdepth;
}
