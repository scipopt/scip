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

#define STP_GENSTAR_MAXDEG 6
#define STP_GENSTAR_MAXENDDEG 4

#define GENSTAR_NODE_OTHER 0
#define GENSTAR_NODE_TAIL  1
#define GENSTAR_NODE_HEAD  2
#define GENSTAR_NODE_COMBI 3

#define EXT_PROFIT_MINRATIO 0.05


enum EXTPSEUDO_MODE { delete_all = 0, delete_profits = 1, delete_nonprofits = 2 };



/** generalized star */
typedef struct general_star
{
   STAR* star;
   STP_Vectype(int) edges_tail;
   STP_Vectype(int) edges_head;
   STP_Vectype(int) edges_all;
   int* nodes_mark;
   int* tmp_compedges;
   int* tmp_extleaves;
   int edge;
   int replacedge_tail;
   int replacedge_head;
} GENSTAR;


/** helper */
typedef struct extension_pseudo_deletion
{
   SCIP_Real*            cutoffs;            /**< cutoffs array of size STP_DELPSEUDO_MAXNEDGES */
   int*                  nodestouches;       /**< touches array on nodes */
   SCIP_Real*            offsetp;            /**< pointer to store offset */
   STP_Bool*             nodeInSol;          /**< solution marker per node or NULL */
   int                   node;               /**< current node */
   int                   nelims_extended;    /**< number of eliminations */
   int                   nelims_simple;      /**< number of eliminations for simple deletion */
   SCIP_Bool             cutoffIsPromising;  /**< is sufficient cutoff possible? */
   SCIP_Bool             cutoffIsComputed;   /**< already computed?? */
   SCIP_Bool             useSolForEq;        /**< use solution for equality? */
   enum EXTPSEUDO_MODE   deletionMode;       /**< what to delete */
} EXTPSEUDO;


#ifndef NDEBUG
/** all good with the graph->mark array? */
static
SCIP_Bool graphmarkIsClean(
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const GRAPH*          graph              /**< graph data structure */
)
{
   const int nnodes = graph_get_nNodes(graph);
   const int redcostroot = redcosts_getRootTop(redcostdata);

   assert(graph_knot_isInRange(graph, redcostroot));

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph->grad[k] == 0 && k != redcostroot && !Is_term(graph->term[k]) )
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


/** replaces edge by a path */
static inline
SCIP_RETCODE replaceEdgeByPath(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to replace */
   const GENSTAR*        genstar,            /**< star */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   REDCOST*              redcostdata,        /**< reduced cost data structures */
   DISTDATA*             distdata,           /**< distance data (in/out) */
   EXTPERMA*             extpermanent        /**< (in/out) */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];
   const int path_edgein = flipedge(genstar->replacedge_tail);
   const int path_edgeout = genstar->replacedge_head;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("replacing edge ");
   graph_edge_printInfo(graph, edge);

   SCIPdebugMessage("by path \n");
   graph_edge_printInfo(graph, path_edgein);
   graph_edge_printInfo(graph, edge);
   graph_edge_printInfo(graph, path_edgeout);
#endif

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));

   distdata->hasPathReplacement = TRUE;

   if( distdata->sdistdata )
   {
      SCIP_CALL_ABORT( reduce_sdRepair(scip, edge, TRUE, graph, distdata->sdistdata) );
   }

   {
       const int csredge = graph->dcsr_storage->id2csredge[edge];
       assert(graph_valid_dcsr(graph, FALSE));
       graph_dcsr_deleteEdgeBi(scip, graph->dcsr_storage, csredge);
       assert(graph_valid_dcsr(graph, FALSE));
    }

   SCIP_CALL( graph_edge_delPseudoPath(scip, graph, edge, path_edgein, path_edgeout, redcosts_getEdgeCostsTop(redcostdata)) );
   extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

   if( extpermanent )
   {
      const int path_tail = graph->tail[path_edgein];
      const int path_head = graph->head[path_edgeout];

      reduce_impliedNodesRepair(scip, graph, path_tail, path_head, extpermanent->nodes_implications);
      reduce_impliedNodesRepair(scip, graph, tail, head, extpermanent->nodes_implications);
      assert(reduce_impliedNodesIsValid(graph, (const STP_Vectype(int)*) extpermanent->nodes_implications));
   }

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


   return SCIP_OKAY;
}


/** deletes an edge and makes corresponding adaptations */
static inline
void removeEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata,           /**< distance data (in/out) */
   EXTPERMA*             extpermanent        /**< (in/out) */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("removing edge ");
   graph_edge_printInfo(graph, edge);
#endif

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));

   if( distdata->sdistdata )
   {
      SCIP_CALL_ABORT( reduce_sdRepair(scip, edge, distdata->hasPathReplacement, graph, distdata->sdistdata) );
   }

   /*
   static int cc = 0;
   char name[1000];
   sprintf(name, "outprev%d.stp", cc);
   graph_writeStpByName(scip, graph, name, 0.0);
      graph_edge_printInfo(graph, edge);

   */

   graph_edge_delFull(scip, graph, edge, TRUE);
   extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

   if( extpermanent )
   {
      reduce_impliedNodesRepair(scip, graph, tail, head, extpermanent->nodes_implications);
      assert(reduce_impliedNodesIsValid(graph, (const STP_Vectype(int)*) extpermanent->nodes_implications));
   }

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

/** initialize */
static inline
SCIP_RETCODE extInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useSd,              /**< use special distance? */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   DISTDATA**            distdata,           /**< distance data (out) */
   EXTPERMA**            extpermanent        /**< permanent extension data (out) */
)
{
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, useSd, FALSE, distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_full, graph, edgedeletable, extpermanent) );

   return SCIP_OKAY;
}


/** free */
static inline
void extFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   DISTDATA**            distdata,           /**< distance data (in/out) */
   EXTPERMA**            extpermanent        /**< permanent extension data (in/out) */
)
{
   extreduce_extPermaFree(scip, extpermanent);
   extreduce_distDataFree(scip, graph, distdata);
   graph_free_dcsr(scip, graph);
}


/** initializes for check */
static inline
void generalStarCheckInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure  */
   GENSTAR*              genstar             /**< general star */
)
{
   int* const nodes_mark = genstar->nodes_mark;
   const int ntails = StpVecGetSize(genstar->edges_tail);
   const int nheads = StpVecGetSize(genstar->edges_head);

   assert(ntails >= 1 && nheads >= 1);
   assert(ntails + nheads >= 3);

   StpVecClear(genstar->edges_all);

   for( int i = 0; i < ntails; i++ )
   {
      const int edge = genstar->edges_tail[i];
      const int head = g->head[edge];
      assert(graph_edge_isInRange(g, edge));
      assert(nodes_mark[head] == GENSTAR_NODE_OTHER);

      StpVecPushBack(scip, genstar->edges_all, edge);
      nodes_mark[head] = GENSTAR_NODE_TAIL;
   }

   for( int i = 0; i < nheads; i++ )
   {
      const int edge = genstar->edges_head[i];
      const int head = g->head[edge];
      assert(graph_edge_isInRange(g, edge));
      assert(nodes_mark[head] == GENSTAR_NODE_OTHER || nodes_mark[head] == GENSTAR_NODE_TAIL);

      StpVecPushBack(scip, genstar->edges_all, edge);

      if( nodes_mark[head] == GENSTAR_NODE_OTHER )
      {
         nodes_mark[head] = GENSTAR_NODE_HEAD;
      }
      else
      {
         assert(nodes_mark[head] == GENSTAR_NODE_TAIL);
         nodes_mark[head] = GENSTAR_NODE_COMBI;
      }
   }

   assert(StpVecGetSize(genstar->edges_all) == ntails + nheads);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("\n all general star edges: \n");
   for( int i = 0; i < StpVecGetSize(genstar->edges_all); i++ )
   {
      graph_edge_printInfo(g, genstar->edges_all[i]);
   }
#endif

   reduce_starResetWithEdges(g, genstar->edges_all, genstar->star);
}


/** exits from check */
static inline
void generalStarCheckExit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure  */
   GENSTAR*              genstar             /**< general star */
)
{
   int* const nodes_mark = genstar->nodes_mark;
   const int degree = StpVecGetSize(genstar->edges_all);

   assert(degree == (StpVecGetSize(genstar->edges_head) + StpVecGetSize(genstar->edges_tail)));

   for( int i = 0; i < degree; i++ )
   {
      const int edge = genstar->edges_all[i];
      const int head = g->head[edge];
      assert(graph_edge_isInRange(g, edge));

      nodes_mark[head] = GENSTAR_NODE_OTHER;
   }
}


/** gets next star */
static inline
void generalStarCheckGetNextStar(
   const GRAPH*          g,                  /**< graph data structure  */
   GENSTAR*              genstar,            /**< general star */
   EXTCOMP*              extcomp,            /**< to be filled */
   SCIP_Bool*            allVisited          /**< all stars visited? */
)
{
   const int center_edge = genstar->edge;
   const int center_tail = g->tail[center_edge];
   int nstaredges;
   int* const compedges = extcomp->compedges;
   int* const extleaves = extcomp->extleaves;
   const int* staredges = reduce_starGetNext(genstar->star, &(nstaredges));
   int* const nodes_mark = genstar->nodes_mark;
   SCIP_Bool hasConflict = FALSE;
   int ntailnodes = 0;
   int nheadnodes = 0;

   *allVisited = FALSE;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("star edges: \n");
   for( int i = 0; i < nstaredges; i++ )
   {
      graph_edge_printInfo(g, staredges[i]);
   }
#endif

   assert(graph_edge_isInRange(g, center_edge));
   assert(nstaredges >= 3);

   /* add component edges and check for conflicts */
   for( int i = 0; i < nstaredges; i++ )
   {
      const int outedge = staredges[i];
      const int head = g->head[outedge];
      assert(graph_edge_isInRange(g, outedge));
      assert(nodes_mark[head] != GENSTAR_NODE_OTHER);

      if( nodes_mark[head] < 0 )
      {
         SCIPdebugMessage("double end node: %d ... \n", head);
         assert(nodes_mark[head] == -GENSTAR_NODE_COMBI);
         hasConflict = TRUE;
         continue;
      }

      if( nodes_mark[head] == GENSTAR_NODE_TAIL )
      {
         ntailnodes++;
         assert(g->tail[outedge] == center_tail);
      }
      else if( nodes_mark[head] == GENSTAR_NODE_HEAD )
      {
         nheadnodes++;
         assert(g->tail[outedge] == g->head[center_edge]);
      }
      else if( nodes_mark[head] == GENSTAR_NODE_COMBI )
      {
         if( g->tail[outedge] == center_tail )
         {
            ntailnodes++;
         }
         else
         {
            assert(g->tail[outedge] == g->head[center_edge]);
            nheadnodes++;
         }
      }

      /* mark as visited */
      nodes_mark[head] *= -1;
   }

   /* clean up */
   for( int i = 0; i < nstaredges; i++ )
   {
      const int outedge = staredges[i];
      const int head = g->head[outedge];
      assert(graph_edge_isInRange(g, outedge));
      assert(nodes_mark[head] != GENSTAR_NODE_OTHER);

      if( nodes_mark[head] < 0 )
         nodes_mark[head] *= -1;
   }

   if( !hasConflict && (ntailnodes == 0 || nheadnodes == 0) )
   {
      SCIPdebugMessage("one sided star... \n");
      hasConflict = TRUE;
   }

   if( hasConflict )
   {
      SCIPdebugMessage("...conflict found! \n");
      if( reduce_starAllAreChecked(genstar->star) )
      {
         *allVisited = TRUE;
      }
      else
      {
         generalStarCheckGetNextStar(g, genstar, extcomp, allVisited);
      }

      return;
   }

   extcomp->ncompedges = nstaredges;
   extcomp->nextleaves = nstaredges - 1;

   assert(g->tail[staredges[0]] == center_tail);
   compedges[0] = flipedge(staredges[0]);

   for( int i = 1; i < nstaredges; i++ )
   {
      const int outedge = staredges[i];

      compedges[i] = outedge;
   }


   for( int i = 1; i < nstaredges; i++ )
   {
      const int outedge = staredges[i];
      extleaves[i - 1] = g->head[outedge];
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("select component edges: \n");
   for( int i = 0; i < nstaredges; i++ )
   {
      graph_edge_printInfo(g, compedges[i]);
   }
#endif

   assert(compedges[0] >= 0);
   extcomp->comproot = g->tail[compedges[0]];

   assert(extcomp->ncompedges >= 3);
   assert(extcomp->comproot >= 0);
}


/** check for elimination */
static inline
SCIP_RETCODE generalStarCheck(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   GENSTAR*              genstar,            /**< general star */
   EXTPERMA*             extpermanent,       /**< data */
   SCIP_Bool*            isDeletable         /**< deletable? */
)
{
   const int* result = extpermanent->result;
   const int edge = genstar->edge;
   const int degree_tail = graph->grad[graph->tail[edge]];
   const int degree_head = graph->grad[graph->head[edge]];
#ifndef NDEBUG
   const int tree_maxdepth_org = extpermanent->tree_maxdepth;
#endif

#ifdef SCIP_DEBUG
   SCIPdebugMessage("checking GENERAL STAR edge ");
   graph_edge_printInfo(graph, edge);
#endif

   *isDeletable = TRUE;

   assert(graph_edge_isInRange(graph, edge));
   assert(degree_tail + degree_head <= STP_GENSTAR_MAXDEG + 2);


   if( degree_tail < 3 && degree_head < 3 )
   {
      SCIPdebugMessage("general-star early rule-out! \n");
      return SCIP_OKAY;
   }

   if( degree_tail == 1 || degree_head == 1 )
   {
      SCIPdebugMessage("general-star early rule-out! \n");
      return SCIP_OKAY;
   }

   if( degree_tail + degree_head == STP_GENSTAR_MAXDEG + 2 )
      extpermanent->tree_maxdepth--;

   generalStarCheckInit(scip, graph, genstar);
   extpermanent->redcostEqualAllow = (extpermanent->solIsValid && result[edge] != CONNECT && result[flipedge(edge)] != CONNECT);

   /* check all general stars of degree >= 3, as long as they can be ruled-out */
   while ( *isDeletable )
   {
      SCIP_Bool allVisited;
      EXTCOMP extcomp = { .compedges = genstar->tmp_compedges, .extleaves = genstar->tmp_extleaves,
                          .nextleaves = -1, .ncompedges = -1,  .genstar_centeredge = genstar->edge,
                          .comproot = -1, .allowReversion = FALSE };

      generalStarCheckGetNextStar(graph, genstar, &extcomp, &allVisited);

      if( allVisited )
         break;

      *isDeletable = FALSE;
      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, extpermanent, isDeletable) );

      if( reduce_starAllAreChecked(genstar->star) )
         break;
   }

   if( degree_tail + degree_head == STP_GENSTAR_MAXDEG + 2 )
      extpermanent->tree_maxdepth++;

   generalStarCheckExit(scip, graph, genstar);

   assert(tree_maxdepth_org == extpermanent->tree_maxdepth);

   return SCIP_OKAY;
}


/** sets up the general star (especially edges_tail and edges_head)
 *  and checks whether deletion attempt makes sense */
static inline
void generalStarSetUp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< inducing edge */
   GENSTAR*              genstar,            /**< general star */
   SCIP_Bool*            isPromising,        /**< promising? */
   DISTDATA*             distdata            /**< distance data */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];
   StpVecClear(genstar->edges_tail);
   StpVecClear(genstar->edges_head);
   genstar->edge = edge;
   genstar->replacedge_tail = -1;
   genstar->replacedge_head = -1;

   *isPromising = FALSE;

   /* NOTE: considered too expensive, since SD tree would need to be rebuilt */
   if( reduce_sdgraphEdgeIsInMst( distdata->sdistdata->sdgraph, edge) )
      return;

   if( (graph->grad[tail] + graph->grad[head]) <= (STP_GENSTAR_MAXDEG + 2) && !Is_term(graph->term[tail]) && !Is_term(graph->term[head]) )
   {
      PSEUDOANS* const pseudoancestors = graph->pseudoancestors;
      const SCIP_Real edgecost = graph->cost[edge];
      const SCIP_Real maxsdcost = reduce_sdgraphGetMaxCost(distdata->sdistdata->sdgraph);
      const STP_Bool* halfedges_isInSdMst = reduce_sdgraphGetMstHalfMark(distdata->sdistdata->sdgraph);
      int* hasharr;
      int ntails;
      int nfails;
      int nheads;
      const int hashsize = graph_pseudoAncestorsGetHashArraySize(pseudoancestors);
      const STP_Vectype(int) edges_tail;
      const STP_Vectype(int) edges_head;

      assert(hashsize > 0);

      *isPromising = TRUE;

      SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &hasharr, hashsize) );

      for( int e = graph->outbeg[tail]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int myhead = graph->head[e];

         if( myhead != head )
         {
            StpVecPushBack(scip, genstar->edges_tail, e);
         }
      }

      for( int e = graph->outbeg[head]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int myhead = graph->head[e];

         if( myhead != tail )
         {
            StpVecPushBack(scip, genstar->edges_head, e);
         }
      }

      ntails = StpVecGetSize(genstar->edges_tail);
      nheads = StpVecGetSize(genstar->edges_head);
      edges_tail = genstar->edges_tail;
      edges_head = genstar->edges_head;
      assert(ntails + nheads <= 6);

      nfails = 0;

      graph_pseudoAncestors_hashEdge(pseudoancestors, edge, hasharr);

      for( int j = 0; j < ntails && *isPromising; j++ )
      {
         SCIP_Bool conflict;
         const int edge_j = edges_tail[j];
         const int node_j = graph->head[edge_j];
         graph_pseudoAncestors_hashEdgeDirty(pseudoancestors, edge_j, TRUE, &conflict, hasharr);

         if( conflict )
            continue;

         for( int k = 0; k < nheads; k++ )
         {
            const int edge_k = edges_head[k];
            const SCIP_Real pathcost = graph->cost[edge_j] + graph->cost[edge_k] + edgecost;
            const int node_k = graph->head[edges_head[k]];

            assert(*isPromising);

            if( node_j == node_k )
               continue;

            graph_pseudoAncestors_hashEdgeDirty(pseudoancestors, edge_k, TRUE, &conflict, hasharr);

            if( conflict )
               continue;

            graph_pseudoAncestors_unhashEdge(pseudoancestors, edge_k, hasharr);

            if( GE(maxsdcost, pathcost) || halfedges_isInSdMst[edge / 2] )
            {
               const SCIP_Real sd = extreduce_distDataGetSdDoubleForbiddenSingle(scip, graph, edge, node_j, node_k, distdata);

               SCIPdebugMessage("%d->%d sd=%f pathcost=%f \n", node_j, node_k, sd, pathcost);

               if( sd < -0.5 || GT(sd, pathcost) )
               {
                  SCIPdebugMessage("...not promising, skip edge \n");

                  if( ++nfails > 1 )
                  {
                     *isPromising = FALSE;
                     break;
                  }

                  genstar->replacedge_tail = edge_j;
                  genstar->replacedge_head = edge_k;
               }
            }
         }

         graph_pseudoAncestors_unhashEdge(pseudoancestors, edge_j, hasharr);
      }

      assert(*isPromising == (nfails <= 1));

      graph_pseudoAncestors_unhashEdge(pseudoancestors, edge, hasharr);

#ifndef NDEBUG
      for( int i = 0; i < hashsize; i++ )
      {
         assert(hasharr[i] == 0);
      }
#endif

      SCIPfreeCleanBufferArray(scip, &hasharr);
   }
}


/** initializes */
static
SCIP_RETCODE generalStarInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   GENSTAR*              genstar             /**< general star */
)
{
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(genstar->nodes_mark), graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(genstar->tmp_compedges), STP_GENSTAR_MAXDEG) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(genstar->tmp_extleaves), STP_GENSTAR_MAXDEG - 1) );
   StpVecReserve(scip, genstar->edges_tail, STP_GENSTAR_MAXENDDEG);
   StpVecReserve(scip, genstar->edges_head, STP_GENSTAR_MAXENDDEG);
   StpVecReserve(scip, genstar->edges_all, STP_GENSTAR_MAXDEG);
   SCIP_CALL( reduce_starInit(scip, STP_GENSTAR_MAXDEG, &(genstar->star)) );

   return SCIP_OKAY;
}


/** exits */
static
void generalStarExit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   GENSTAR*              genstar             /**< general star */
)
{
   StpVecFree(scip, genstar->edges_tail);
   StpVecFree(scip, genstar->edges_head);
   StpVecFree(scip, genstar->edges_all);
   reduce_starFree(scip, &(genstar->star));
   SCIPfreeBufferArray(scip, &(genstar->tmp_compedges));
   SCIPfreeBufferArray(scip, &(genstar->tmp_extleaves));
#ifndef NEBDUG
   for( int i = 0; i < graph->knots; i++ )
      assert(genstar->nodes_mark[i] == 0);
#endif
   SCIPfreeCleanBufferArray(scip,  &(genstar->nodes_mark));
}


/** deletes center edges of general stars */
static
SCIP_RETCODE generalStarDeleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   REDCOST*              redcostdata,        /**< reduced cost data structures */
   EXTPERMA*             extpermanent,       /**< extension data */
   GRAPH*                graph,              /**< graph data structure */
   DISTDATA*             distdata,           /**< distance data */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   GENSTAR genstar = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1, -1, -1 };
   const int* result = extpermanent->result;
   int ncands = 0;
   int npseudoelims = 0;

   assert(!extpermanent->solIsValid || result);

   *nelims = 0;

   SCIPdebugMessage("General-star deletion starts \n");

   SCIP_CALL( generalStarInit(scip, graph, &genstar) );

   for( int i = 0; i < graph->edges; i += 2 )
   {
      SCIP_Bool isPromising;

      if( graph_edge_isDeleted(graph, i) )
         continue;

      generalStarSetUp(scip, graph, i, &genstar, &isPromising, distdata);

      if( isPromising )
      {
         SCIP_Bool isDeletable;
         SCIP_CALL( generalStarCheck(scip, graph, redcostdata, &genstar, extpermanent, &isDeletable ) );
         ncands++;
         if( isDeletable )
         {
            if( extpermanent->solIsValid && (result[i] == CONNECT || result[flipedge(i)] == CONNECT ) )
               extpermanent->solIsValid = FALSE;

            if( genstar.replacedge_head != -1 )
            {
               npseudoelims++;
               SCIP_CALL( replaceEdgeByPath(scip, i, &genstar, graph, redcostdata, distdata, extpermanent) );

               continue;
            }

            removeEdge(scip, i, graph, distdata, extpermanent);
            (*nelims)++;
         }
      }
   }

   generalStarExit(scip, graph, &genstar);

   SCIPdebugMessage("number of general star eliminations=%d \n", *nelims);
   SCIPdebugMessage("number of general star replacements=%d \n", npseudoelims);

   return SCIP_OKAY;
}

/** initializes */
static inline
SCIP_RETCODE pseudodeleteInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            result,             /**< solution array or NULL */
   const GRAPH*          g,                  /**< graph data structure  */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   EXTPSEUDO*            extpseudo           /**< to initialize */
   )
{
   SCIP_Real* cutoffs;
   int* nodestouches;
   const int nnodes = graph_get_nNodes(g);

   assert(extpseudo);

   if( result )
   {
      STP_Bool* nodeInSol;
      assert(solstp_isValid(scip, g, result));

      SCIP_CALL( SCIPallocBufferArray(scip, &nodeInSol, nnodes) );
      solstp_setVertexFromEdge(g, result, nodeInSol);

      extpseudo->nodeInSol = nodeInSol;
      extpseudo->useSolForEq = TRUE;
   }
   else
   {
      extpseudo->nodeInSol = NULL;
      extpseudo->useSolForEq = FALSE;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &cutoffs, STP_DELPSEUDO_MAXNEDGES) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodestouches, nnodes) );

   for( int i = 0; i < nnodes; ++i )
      nodestouches[i] = 0;

   extpseudo->cutoffIsPromising = FALSE;
   extpseudo->cutoffIsComputed = FALSE;
   extpseudo->node = -1;
   extpseudo->offsetp = offsetp;
   extpseudo->nelims_extended = 0;
   extpseudo->nelims_simple = 0;
   extpseudo->nodestouches = nodestouches;
   extpseudo->cutoffs = cutoffs;
   extpseudo->deletionMode = delete_all;

   return SCIP_OKAY;
}


/** frees */
static inline
void pseudodeleteExit(
   SCIP*                 scip,               /**< SCIP data structure */
   EXTPSEUDO*            extpseudo,          /**< to initialize */
   int*                  nelimsp             /**< pointer: number of eliminations (OUT) */
   )
{
   assert(nelimsp);
   assert(*nelimsp == 0);
   assert(extpseudo->nelims_extended >= 0 && extpseudo->nelims_simple >= 0);


   *nelimsp = extpseudo->nelims_extended + extpseudo->nelims_simple;

   SCIPfreeBufferArray(scip, &(extpseudo->nodestouches));
   SCIPfreeBufferArray(scip, &(extpseudo->cutoffs));
   SCIPfreeBufferArrayNull(scip, &(extpseudo->nodeInSol));
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


/** is biased SD extension promising? */
static inline
SCIP_Bool pseudodeleteBiasedIsPromising(
   SCIP*                 scip,               /**< SCIP data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   const GRAPH*          g                   /**< graph data structure  */
)
{
   if( graph_pc_isPc(g) )
   {
      return FALSE;
   }
   else
   {
      SDPROFIT* sdprofit;
      double profitsratio;
      int nprofits = 0;
      const int nnodes = graph_get_nNodes(g);
      const int nterms = g->terms;
      // todo might also compute proper sd first..and give to routine
      SCIP_CALL_ABORT( reduce_sdprofitInit1stOnly(scip, g, g->cost, &sdprofit));

      for( int i = 0; i < nnodes; i++ )
      {
         if( Is_term(g->term[i]) )
            continue;

         if( GT(sdprofit->nodes_bias[i], 0.0) )
         {
            nprofits++;
         }
      }

      profitsratio = (double) nprofits / (double) nterms;
      SCIPdebugMessage("profitsratio=%f \n", profitsratio);

      reduce_sdprofitFree(scip, &sdprofit);

      if( profitsratio <  EXT_PROFIT_MINRATIO )
      {
         SCIPdebugMessage("biased SD for extension not promising: %f < %f \n", profitsratio, EXT_PROFIT_MINRATIO);

         return FALSE;
      }
   }

   return TRUE;
}

/** is node a good candidate for pseudo deletion? */
static inline
SCIP_Bool pseudodeleteNodeIsPromising(
   const GRAPH*          g,                  /**< graph data structure  */
   const EXTPERMA*       extperma,           /**< extension data */
   const EXTPSEUDO*      extpseudo,          /**< pseudo */
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

      assert(degree == g->grad[node]); // otherwise the outside loop does not work..
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

   if( extpseudo->deletionMode != delete_all )
   {
      SCIP_Bool hasProfit;
      assert(extperma->distdata_biased);
      assert(extperma->distdata_biased->sdistdata);
      assert(extperma->distdata_biased->sdistdata->sdprofit);
      assert(!graph_pc_isPcMw(g));

      hasProfit = GT(reduce_sdprofitGetProfit(extperma->distdata_biased->sdistdata->sdprofit, node, -1, -1), 0.0);

      if( extpseudo->deletionMode == delete_profits && !hasProfit )
         return FALSE;

      if( extpseudo->deletionMode == delete_nonprofits && hasProfit )
         return FALSE;
   }

   return TRUE;
}



/** computes cutoffs for pseudo-deletion of given node */
static inline
SCIP_RETCODE pseudodeleteDeleteComputeCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             checkpromising,     /**< check whether promising? */
   SCIP_Bool             abortDeg3,          /**< abort for degree 3? */
   DISTDATA*             distdata,           /**< distance data */
   int                   node,               /**< to be deleted */
   GRAPH*                graph,              /**< graph data structure */
   EXTPSEUDO*            extpseudo           /**< data */
)
{
   int adjedges[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real* cutoffs = extpseudo->cutoffs;
   const int* const nodestouches = extpseudo->nodestouches;
   int edgecount = 0;
   const int degree = graph->grad[node];
   const SCIP_Bool isPc = graph_pc_isPc(graph);
   const SCIP_Real eps = 2.0 * SCIPepsilon(scip);

   extpseudo->node = node;

   if( abortDeg3 && degree <= 3 )
   {
      extpseudo->cutoffIsPromising = TRUE;
      extpseudo->cutoffIsComputed = FALSE;
      return SCIP_OKAY;
   }
   extpseudo->cutoffIsComputed = TRUE;

   if( isPc && graph_pc_knotIsNonLeafTerm(graph, node) )
   {
      for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int head = graph->head[e];
         if( !graph_pc_knotIsDummyTerm(graph, head) )
            adjedges[edgecount++] = e;
      }
   }
   else
   {
      for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
         adjedges[edgecount++] = e;
   }

   edgecount = 0;
   for( int i = 0; i < degree - 1; i++ )
   {
      const int vert = graph->head[adjedges[i]];
      const SCIP_Real edgecost = graph->cost[adjedges[i]];
      for( int i2 = i + 1; i2 < degree; i2++ )
      {
         const int vert2 = graph->head[adjedges[i2]];
         const SCIP_Real edgecost2 = graph->cost[adjedges[i2]];
         const SCIP_Real newedgecost = edgecost + edgecost2;

         assert(edgecount < STP_DELPSEUDO_MAXNEDGES);

         cutoffs[edgecount] = extreduce_distDataGetSdDouble(scip, graph, vert, vert2, distdata);

         if( SCIPisEQ(scip, newedgecost, cutoffs[edgecount]) && nodestouches[node] == 0 && nodestouches[vert] == 0 && nodestouches[vert2] == 0 )
         {
            cutoffs[edgecount] = extreduce_distDataGetSdDoubleForbiddenLast(scip, graph, vert, vert2,
                  adjedges[i2], adjedges[i], distdata);
            assert(SCIPisLE(scip, newedgecost, cutoffs[edgecount]) || EQ(cutoffs[edgecount], -1.0));

            if( SCIPisEQ(scip, newedgecost, cutoffs[edgecount]) )
            {
               cutoffs[edgecount] -= eps;
               assert(SCIPisGT(scip, newedgecost, cutoffs[edgecount]));
            }
         }

         if( LT(cutoffs[edgecount], 0.0) )
         {
            cutoffs[edgecount] = FARAWAY;
         }

         edgecount++;
      }
   }

   extpseudo->cutoffIsPromising = FALSE;

   if( checkpromising )
   {
      SCIP_CALL( graph_knot_delPseudoCheckIfPossible(scip, graph, graph->cost, cutoffs, NULL, node, &(extpseudo->cutoffIsPromising)) );
   }

   return SCIP_OKAY;
}


/** pseudo-deletes single node */
static inline
SCIP_RETCODE pseudodeleteDeleteNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node,               /**< to be deleted */
   REDCOST*              redcostdata,        /**< reduced cost data */
   DISTDATA*             distdata,           /**< distance data */
   GRAPH*                graph,              /**< graph data structure */
   EXTPSEUDO*            extpseudo,          /**< data */
   SCIP_Bool*            success             /**< success? */
)
{
   const SCIP_Real* cutoffs = extpseudo->cutoffs;
   int* const nodestouches = extpseudo->nodestouches;
   SCIP_Real prize = -1.0;
   SCIP_Bool rpc3term = FALSE;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   assert(redcostdata);
   assert(node == extpseudo->node);

   if( !extpseudo->cutoffIsComputed )
   {
      SCIP_CALL( pseudodeleteDeleteComputeCutoffs(scip, FALSE, FALSE, distdata, node, graph, extpseudo) );
   }

   if( isPc && graph_pc_knotIsNonLeafTerm(graph, node) && graph->grad[node] == 3 )
   {
      rpc3term = TRUE;
      prize = graph->prize[node];
   }

   for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
      nodestouches[graph->head[e]]++;

   /* now try to eliminate... */
   SCIP_CALL( graph_knot_delPseudo(scip, graph, graph->cost, cutoffs, NULL, node, redcostdata, success) );

   if( *success )
   {
      graph->mark[node] = FALSE;

      if( extpseudo->useSolForEq )
      {
         assert(extpseudo->nodeInSol);

         if( extpseudo->nodeInSol[node] )
            extpseudo->useSolForEq = FALSE;
      }

      SCIPdebugMessage("deletion successful! \n");

      if( rpc3term )
      {
         assert(isPc);
         assert(extpseudo->offsetp);
         assert(GE(prize, 0.0));

         *(extpseudo->offsetp) += prize;
      }
   }
   else
   {
      for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
         nodestouches[graph->head[e]]--;
   }

   return SCIP_OKAY;
}


/** apply pseudo eliminations for marked nodes
 *  NOTE: bad design, but useful to reuse the SDs... */
static
SCIP_RETCODE pseudodeleteDeleteMarkedNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Bool*      pseudoDelNodes,     /**< node with pseudo deletable nodes */
   REDCOST*              redcostdata,        /**< reduced cost data */
   DISTDATA*             distdata,           /**< distance data */
   GRAPH*                graph,              /**< graph data structure */
   EXTPSEUDO*            extpseudo           /**< data */
)
{
   const int nnodes = graph_get_nNodes(graph);

   assert(pseudoDelNodes);
   assert(extpseudo->nelims_simple == 0);

   for( int degree = 2; degree <= STP_DELPSEUDO_MAXGRAD; degree++ )
   {
      for( int k = 0; k < nnodes; ++k )
      {
         if( pseudoDelNodes[k] && graph->grad[k] == degree  )
         {
            SCIP_Bool success;
            assert(!Is_term(graph->term[k]));

            SCIP_CALL( pseudodeleteDeleteComputeCutoffs(scip, FALSE, FALSE, distdata, k, graph, extpseudo) );
            SCIP_CALL( pseudodeleteDeleteNode(scip, k, redcostdata, distdata, graph, extpseudo, &success) );

            if( success )
               extpseudo->nelims_simple++;
         }
      }
   }

   return SCIP_OKAY;
}


/** Extended reduction test for pseudo-eliminating nodes.
 *  Applies pseudo-elimination if pseudoDelNodes != NULL, or just marks them otherwise */
static
SCIP_RETCODE pseudodeleteExecute(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            result,             /**< solution array or NULL */
   const SCIP_Bool*      pseudoDelNodes,     /**< nodes to pseudo-eliminate already */
   EXTPERMA*             extperma,           /**< extension data */
   EXTPSEUDO*            extpseudo,          /**< pseudo-deletion data */
   GRAPH*                graph               /**< graph data structure (in/out) */
)
{
   const int nnodes = graph_get_nNodes(graph);
   STAR* stardata;
   DISTDATA* const distdata_default = extperma->distdata_default;
   REDCOST* const redcostdata = extperma->redcostdata;

   assert(graph_isMarked(graph));

   extreduce_distDataRecomputeDirtyPaths(scip, graph, distdata_default);

   if( pseudoDelNodes )
   {
      SCIP_CALL( pseudodeleteDeleteMarkedNodes(scip, pseudoDelNodes, redcostdata, distdata_default, graph, extpseudo) );
      SCIPdebugMessage("number of eliminations by simple pseudo-elimination %d \n", extpseudo->nelims_simple);
   }

   SCIP_CALL( reduce_starInit(scip, EXT_PSEUDO_DEGREE_MAX, &stardata) );

   reduce_nodesDeg1(scip, graph);

   for( int degree = EXT_PSEUDO_DEGREE_MIN; degree <= EXT_PSEUDO_DEGREE_MAX; degree++  )
   {
      /* pseudo-elimination loop */
      for( int i = 0; i < nnodes; ++i )
      {
         SCIP_Bool isReplaced;
         SCIP_Bool nodeisDeletable = FALSE;

         if( graph->grad[i] == 1 && !Is_term(graph->term[i]) )
         {
            graph_knot_del(scip, graph, i, TRUE);
            graph->mark[i] = FALSE;
            continue;
         }

         if( graph->grad[i] != degree )
            continue;

         if( pseudodeleteNodeIsPromising(graph, extperma, extpseudo, i) )
         {
            SCIP_CALL( pseudodeleteDeleteComputeCutoffs(scip, TRUE, TRUE, distdata_default, i, graph, extpseudo) );

            if( extpseudo->cutoffIsPromising )
            {
               extperma->redcostEqualAllow = (extpseudo->useSolForEq && !extpseudo->nodeInSol[i]);

               if( extperma->useSdBias && degree == 3 )
                  SCIP_CALL( extreduce_mstbiasedCheck3NodeSimple(scip, graph, i, distdata_default, extperma->distdata_biased, &nodeisDeletable) );

               if( !nodeisDeletable )
                  SCIP_CALL( extreduce_checkNode(scip, graph, redcostdata, i, stardata, extperma, &nodeisDeletable) );
            }
         }

         if( !nodeisDeletable )
            continue;

         SCIP_CALL( pseudodeleteDeleteNode(scip, i, redcostdata,
         //   todo   extperma->useSdBias ? extperma->distdata_biased :
                     distdata_default, graph, extpseudo, &isReplaced) );

         if( isReplaced )
            extpseudo->nelims_extended++;
      }
   }

   reduce_starFree(scip, &stardata);
   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** initializes */
SCIP_RETCODE extreduce_init(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useSd,              /**< use special distance? */
   enum EXTRED_MODE      mode,               /**< mode */
   GRAPH*                graph,              /**< graph data structure */
   REDCOST*              redcostdata,        /**< reduced costs data */
   EXTPERMA**            extpermanent        /**< permanent extension data (out) */
)
{
   EXTPERMA* extperma;

   assert(scip && graph && redcostdata && extpermanent);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   graph_mark(graph);
   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_extPermaInit(scip, mode, graph, NULL, extpermanent) );
   extperma = *extpermanent;

   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, useSd, FALSE, &(extperma->distdata_default)) );
   extperma->redcostdata = redcostdata;

   return SCIP_OKAY;
}


/** frees */
void extreduce_exit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   EXTPERMA**            extpermanent        /**< permanent extension data */
)
{
   EXTPERMA* extperma;

   assert(scip && graph && extpermanent);

   extperma = *extpermanent;

   if( extperma->distdata_biased )
      extreduce_distDataFree(scip, graph, &(extperma->distdata_biased));

   graph_free_dcsr(scip, graph);
   extreduce_distDataFree(scip, graph, &(extperma->distdata_default));
   extreduce_extPermaFree(scip, extpermanent);
}



/** Extended reduction test for arcs.
 * This method will also set edgedeletable[a] to TRUE if arc 'a' can be deleted, but its anti-parallel arc not. */
SCIP_RETCODE extreduce_deleteArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   REDCOST*              redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   const SCIP_Bool useSd = !graph_pc_isPc(graph);
   const int nedges = graph_get_nEdges(graph);
   DISTDATA* distdata;
   EXTPERMA* extpermanent;
   SCIP_Bool withSol = (result != NULL);

   assert(scip && redcostdata && edgedeletable);

   *nelims = 0;

   if( SCIPisZero(scip, redcosts_getCutoffTop(redcostdata)) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, useSd, graph, NULL, &distdata, &extpermanent) );

   if( useSd )
   {
      assert(distdata->sdistdata);
      SCIP_CALL( reduce_sdRepairSetUp(scip, graph, distdata->sdistdata) );
   }

   extpermanent->distdata_default = distdata;
   extpermanent->result = result;

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, redcostdata, e) )
      {
         const int erev = e + 1;
         extpermanent->redcostEqualAllow = (withSol && result[e] != CONNECT && result[erev] != CONNECT);
         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( !edgedeletable[e] )
         {
            SCIP_Bool deletable;
            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, e, extpermanent,
                  &deletable) );

            if( deletable )
            {
               if( withSol && result[e] == CONNECT )
                  withSol = FALSE;

               edgedeletable[e] = TRUE;
            }
         }

         if( !edgedeletable[erev] )
         {
            SCIP_Bool erevdeletable = FALSE;

            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, erev, extpermanent,
                  &erevdeletable) );

            if( erevdeletable )
            {
               if( withSol && result[erev] == CONNECT )
                  withSol = FALSE;

               edgedeletable[erev] = TRUE;
            }
         }

         if( edgedeletable[e] && edgedeletable[erev] )
         {
            assert(edgedeletable[e] && edgedeletable[erev]);

            removeEdge(scip, e, graph, distdata, extpermanent);

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
   EXTPERMA*             extperma,           /**< extension data */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const SCIP_Bool useSd = !graph_pc_isPc(graph);
   const int nedges = graph_get_nEdges(graph);
   REDCOST* const redcostdata = extperma->redcostdata;
   DISTDATA* const distdata = extperma->distdata_default;
   SD* const sddata = distdata->sdistdata;
   const int* result = extperma->result;
   SCIP_Bool withSol;

   assert(scip && redcostdata);
   assert(!graph_pc_isMw(graph) && "not supported yet");
   assert(redcosts_getRootTop(redcostdata) >= 0 && redcosts_getRootTop(redcostdata) < graph->knots);
   assert(graph_isMarked(graph));
   assert((result != NULL) == extperma->solIsValid);
   assert(!extperma->solIsValid || solstp_isValid(scip, graph, result));

   *nelims = 0;

   if( SCIPisZero(scip, redcosts_getCutoffTop(redcostdata)) )
      return SCIP_OKAY;

   SCIPdebugMessage("run extreduce_deleteEdges with depth %d \n", extperma->tree_maxdepth);

   if( useSd )
   {
      assert(sddata);
      SCIP_CALL( reduce_sdRepairSetUp(scip, graph, sddata) );
   }

   if( useSd )
   {
      int npathelims = 0;
      SCIP_CALL( reduce_pathreplaceExt(scip, graph, extperma, &npathelims) );

      (*nelims) += npathelims;
   }

   withSol = extperma->solIsValid;

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, redcostdata, e) )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( useSd && extperma->mode == extred_fast && reduce_sdgraphEdgeIsInMst(sddata->sdgraph, e) )
            continue;

         extperma->redcostEqualAllow = (withSol && result[e] != CONNECT && result[erev] != CONNECT);
         SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, e, extperma, &deletable) );

         if( deletable )
         {
            removeEdge(scip, e, graph, distdata, extperma);

            if( withSol && (result[e] == CONNECT || result[erev] == CONNECT ) )
               withSol = FALSE;

            (*nelims)++;
         }
      }
   }

   extperma->solIsValid = withSol;

   if( graph_typeIsSpgLike(graph) )
   {
      int sepanelims = 0;
      SCIP_CALL( reduce_termsepaDaWithExperma(scip, graph, extperma, NULL, &sepanelims) );
      *nelims += sepanelims;
      graph_mark(graph);

      if( sepanelims > 0 && extperma->solIsValid )
         extperma->solIsValid = solstp_isUnreduced(scip, graph, result);
   }

   if( extperma->mode == extred_full && !graph_pc_isPc(graph) )
   {
      int ngenstarelims = 0;
      SCIP_CALL( generalStarDeleteEdges(scip, redcostdata, extperma, graph, distdata, &ngenstarelims) );
      *nelims += ngenstarelims;
   }

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}

/** extended reduction test for pseudo-eliminating nodes */
SCIP_RETCODE extreduce_pseudoDeleteNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Bool*      pseudoDelNodes,     /**< nodes to pseudo-eliminate already */
   EXTPERMA*             extperma,           /**< extension data */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   EXTPSEUDO extpseudo;
   /* NOTE: hacky, needs to be kept because sub-methods do not use solIsValid flag */
   const int* result = extperma->solIsValid ? extperma->result : NULL;
   SCIP_Bool isExtendedOrg;

   assert(scip && extperma && extperma->redcostdata && graph && nelims);
   assert(!extperma->useSdBias);
   assert(!extperma->solIsValid || solstp_isValid(scip, graph, result));

   *nelims = 0;

   if( SCIPisZero(scip, redcosts_getCutoffTop(extperma->redcostdata)) )
      return SCIP_OKAY;

   SCIPdebugMessage("run extreduce_pseudoDeleteNodes with depth %d \n", extperma->tree_maxdepth);

   isExtendedOrg = graph->extended;
   if( graph_pc_isPc(graph) )
   {
      graph_pc_2orgcheck(scip, graph);
      reduce_removeDeg0NonLeafTerms(scip, graph, offsetp);
   }
   else
   {
      graph_mark(graph);
   }

   SCIP_CALL( pseudodeleteInit(scip, result, graph, offsetp, &extpseudo) );

   if( pseudodeleteBiasedIsPromising(scip, extperma, graph) )
   {
      assert(!extperma->distdata_biased);
      extperma->useSdBias = TRUE;
      extpseudo.deletionMode = delete_nonprofits;

      SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, TRUE, TRUE, &(extperma->distdata_biased)) );
      extperma->distdata_biased->hasPathReplacement = extperma->distdata_default->hasPathReplacement;
      SCIP_CALL( extreduce_extPermaAddMLdistsbiased(scip, extperma) );
   }
   assert(extpseudo.deletionMode == delete_nonprofits || extpseudo.deletionMode == delete_all);

   /* call actual method */
   SCIP_CALL( pseudodeleteExecute(scip, result, NULL, extperma, &extpseudo, graph) );
   SCIPdebugMessage("number of eliminations after extended pseudo-elimination %d \n", extpseudo.nelims_extended);

   if( extperma->useSdBias )
   {
      assert(extperma->distdata_biased);

      extperma->useSdBias = FALSE;
      extpseudo.deletionMode = delete_profits;
      SCIP_CALL( pseudodeleteExecute(scip, result, pseudoDelNodes, extperma, &extpseudo, graph) );

      SCIPdebugMessage("number of eliminations after extended pseudo-elimination on profits %d \n", extpseudo.nelims_extended);
   }

   extperma->solIsValid = extpseudo.useSolForEq;

   /* NOTE: also sets "nelims" pointer*/
   pseudodeleteExit(scip, &extpseudo, nelims);

   /* todo otherwise (I015) we get isolated vertices...not sure whether this is a bug or normal behavior */
   if( *nelims > 0 )
      SCIP_CALL( reduce_unconnected(scip, graph));

   if( graph_pc_isPc(graph) && isExtendedOrg != graph->extended )
      graph_pc_2trans(scip, graph);

   return SCIP_OKAY;
}


/** deletes center edges of general stars */
SCIP_RETCODE extreduce_deleteGeneralStars(
   SCIP*                 scip,               /**< SCIP data structure */
   REDCOST*              redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const SCIP_Bool useSd = !graph_pc_isPc(graph);
   DISTDATA* distdata;
   EXTPERMA* extpermanent;

   assert(scip && redcostdata);
   assert(redcosts_getRootTop(redcostdata) >= 0 && redcosts_getRootTop(redcostdata) < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcosts_getCutoffTop(redcostdata)) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, useSd, graph, edgedeletable, &distdata, &extpermanent) );

   if( useSd )
   {
      assert(distdata->sdistdata);
      SCIP_CALL( reduce_sdRepairSetUp(scip, graph, distdata->sdistdata) );
   }

   extpermanent->distdata_default = distdata;
   extpermanent->result = result;

   SCIP_CALL( generalStarDeleteEdges(scip, redcostdata, extpermanent, graph, distdata, nelims) );

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;

}


/** check (directed) arc */
SCIP_RETCODE extreduce_checkArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   REDCOST*              redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* isterm = extpermanent->isterm;
   const int root = redcosts_getRootTop(redcostdata);
   SCIP_Real* const redcost = redcosts_getEdgeCostsTop(redcostdata);
   const SCIP_Real* rootdist = redcosts_getRootToNodeDistTop(redcostdata);
   const PATH* nodeToTermpaths = redcosts_getNodeToTermsPathsTop(redcostdata);
   const SCIP_Real cutoff = redcosts_getCutoffTop(redcostdata);
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   const int edge_anti = flipedge(edge);
   const SCIP_Real edgebound = redcost[edge] + rootdist[tail] + nodeToTermpaths[head].dist;

   assert(scip && graph && redcost && rootdist && nodeToTermpaths);
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

   *edgeIsDeletable = FALSE;

   /* can we extend from head of 'edge'? */
   if( extLeafIsExtendable(graph, isterm, head) )
   {
      const SCIP_Real restoreCost = redcost[edge_anti];
      int comphead = graph->head[edge];
      int compedge = edge;
      EXTCOMP extcomp = { .compedges = &compedge, .extleaves = &(comphead),
         .nextleaves = 1, .ncompedges = 1, .comproot = graph->tail[edge],
         .genstar_centeredge = -1,
         .allowReversion = FALSE };

      redcost[edge_anti] += 2.0 * cutoff;

      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, extpermanent, edgeIsDeletable) );

      redcost[edge_anti] = restoreCost;
   }


   return SCIP_OKAY;
}


/** check edge */
SCIP_RETCODE extreduce_checkEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* const isterm = extpermanent->isterm;

   assert(scip && graph && edgeIsDeletable && extpermanent);
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
                          .nextleaves = 1, .ncompedges = 1, .genstar_centeredge = -1,
                          .comproot = graph->tail[edge], .allowReversion = TRUE };

      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, extpermanent, edgeIsDeletable) );
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
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            isPseudoDeletable   /**< is node pseudo-deletable? */
)
{
   int* compedges;
   int* extleaves;

   assert(scip && graph && redcostdata && extpermanent && isPseudoDeletable);
   assert(node >= 0 && node < graph->knots);

   *isPseudoDeletable = TRUE;

   SCIP_CALL( pseudodeleteInitStar(scip, graph, node, stardata, &compedges, &extleaves) );

   if( graph->grad[node] == 3 )
      extpermanent->tree_maxdepth++;

   /* check all stars of degree >= 3, as long as they can be ruled-out */
   while ( *isPseudoDeletable )
   {
      EXTCOMP extcomp = { .compedges = compedges, .extleaves = extleaves,
                          .nextleaves = -1, .ncompedges = -1, .genstar_centeredge = -1,
                          .comproot = -1, .allowReversion = TRUE };

      pseudodeleteGetNextStar(graph, stardata, &extcomp);

      *isPseudoDeletable = FALSE;
      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, extpermanent, isPseudoDeletable) );

      if( !(*isPseudoDeletable) && extcomp.ncompedges == 3 && graph->grad[node] <= 4 )
      {
         const SCIP_Bool allowEquality = (graph->grad[node] == 3);
         SCIP_CALL( extreduce_spgCheck3ComponentSimple(scip, graph, node, &extcomp, allowEquality,
               extpermanent->distdata_default, isPseudoDeletable) );
      }

      if( pseudodeleteAllStarsChecked(stardata) )
         break;
   }

   if( graph->grad[node] == 3 )
        extpermanent->tree_maxdepth--;

   pseudodeleteFreeStar(scip, &compedges, &extleaves);

   return SCIP_OKAY;
}


/** deletes an edge and makes corresponding adaptations */
void extreduce_edgeRemove(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata,           /**< distance data (in/out) */
   EXTPERMA*             extpermanent        /**< (in/out) can be NULL */
)
{
   removeEdge(scip, edge, graph, distdata, extpermanent);
}


/** is the edge valid? */
SCIP_Bool extreduce_edgeIsValid(
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   int                   e                   /**< edge to be checked */
)
{
   assert(graph && redcostdata);
   assert(graph_edge_isInRange(graph, e));

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

   if( EQ(0.0, redcosts_getEdgeCostsTop(redcostdata)[e]) && EQ(0.0, redcosts_getEdgeCostsTop(redcostdata)[flipedge(e)]) )
   {
      return FALSE;
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
   const SCIP_Real* const cost = graph->cost;
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

      if( edgeIsDeleted )
      {
         extdata->tree_nDelUpArcs++;
      }
   }

   assert(SCIPisEQ(scip, tree_cost, extdata->tree_cost));
   assert(tree_nDelUpArcs == extdata->tree_nDelUpArcs);

   extreduce_redcostTreeRecompute(scip, graph, extdata);

   extdata->tree_cost = tree_cost;

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
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extpermanent         /**< extension data */
)
{
   int nedges;
   int maxdepth;

   graph_get_nVET(graph, NULL, &nedges, NULL);

   if( nedges > STP_EXT_DEPTH_MAXNEDGES )
   {
      maxdepth = STP_EXT_MINDFSDEPTH;
   }
   else if( nedges > STP_EXT_DEPTH_MIDNEDGES )
   {
      maxdepth = STP_EXT_MIDDFSDEPTH;
   }
   else
   {
      assert(nedges <= STP_EXT_DEPTH_MIDNEDGES);
      maxdepth = STP_EXT_MAXDFSDEPTH;
   }

   if( extpermanent->mode == extred_fast )
   {
      maxdepth--;
   }

   assert(maxdepth > 0);

   return maxdepth;
}
