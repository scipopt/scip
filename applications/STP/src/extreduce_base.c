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

#define STP_GENSTAR_MAXDEG 6
#define STP_GENSTAR_MAXENDDEG 4

#define GENSTAR_NODE_OTHER 0
#define GENSTAR_NODE_TAIL  1
#define GENSTAR_NODE_HEAD  2
#define GENSTAR_NODE_COMBI 3


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
   int*                  nelimsp;            /**< pointer number of eliminations */
   int                   node;               /**< current node */
   SCIP_Bool             cutoffIsPromising;  /**< is sufficient cutoff possible? */
   SCIP_Bool             cutoffIsComputed;   /**< already computed?? */
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
   SCIP_Bool             useSd,              /**< use special distance? */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   DISTDATA*             distdata,           /**< distance data (out) */
   EXTPERMA*             extpermanent        /**< permanent extension data (out) */
)
{
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, useSd, distdata) );
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

   assert(EQ(degree, StpVecGetSize(genstar->edges_head) + StpVecGetSize(genstar->edges_tail)));

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
   SCIP_Bool*            allVisited
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
   DISTDATA*             distdata,           /**< distance data */
   SCIP_Bool*            isDeletable         /**< deletable? */
)
{
   const int edge = genstar->edge;
   *isDeletable = TRUE;

   assert(graph_edge_isInRange(graph, edge));

#ifdef SCIP_DEBUG
   SCIPdebugMessage("checking GENERAL STAR edge ");
   graph_edge_printInfo(graph, edge);
#endif

   if( graph->grad[graph->tail[edge]] < 3 && graph->grad[graph->head[edge]] < 3 )
   {
      SCIPdebugMessage("general-star early rule-out! \n");
      return SCIP_OKAY;
   }

   if( graph->grad[graph->tail[edge]] == 1 || graph->grad[graph->head[edge]] == 1 )
   {
      SCIPdebugMessage("general-star early rule-out! \n");
      return SCIP_OKAY;
   }

   generalStarCheckInit(scip, graph, genstar);

   /* check all general stars of degree >= 3, as long as they can be ruled-out */
   while ( *isDeletable )
   {
      SCIP_Bool allVisited;
      // todo allow reversion
      EXTCOMP extcomp = { .compedges = genstar->tmp_compedges, .extleaves = genstar->tmp_extleaves,
                          .nextleaves = -1, .ncompedges = -1,  .genstar_centeredge = genstar->edge,
                          .comproot = -1, .allowReversion = FALSE };

      generalStarCheckGetNextStar(graph, genstar, &extcomp, &allVisited);

      if( allVisited )
         break;

      *isDeletable = FALSE;
      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, distdata, extpermanent, isDeletable) );

      if( reduce_starAllAreChecked(genstar->star) )
         break;
   }

   generalStarCheckExit(scip, graph, genstar);

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

#if 0
      printf("next: \n");
      graph_edge_printInfo(graph, edge);
      graph_knot_printInfo(graph, tail);
      graph_knot_printInfo(graph, head);
#endif

      for( int e = graph->outbeg[tail]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int myhead = graph->head[e];

         if( myhead != head )
         {
          //  graph_edge_printInfo(graph, e);
            StpVecPushBack(scip, genstar->edges_tail, e);
         }
      }

      for( int e = graph->outbeg[head]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int myhead = graph->head[e];

         if( myhead != tail )
         {
          //  graph_edge_printInfo(graph, e);
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

                  if( nfails++ > 1 )
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

      graph_pseudoAncestors_unhashEdge(pseudoancestors, edge, hasharr);


      for( int i = 0; i < hashsize; i++ )
      {
         assert(hasharr[i] == 0);
      }
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
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   EXTPERMA*             extpermanent,       /**< extension data */
   GRAPH*                graph,              /**< graph data structure */
   DISTDATA*             distdata,           /**< distance data */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   GENSTAR genstar = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, -1, -1, -1 };
   *nelims = 0;
int ncands = 0;
int npseudoelims = 0;
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
         SCIP_CALL( generalStarCheck(scip, graph, redcostdata, &genstar, extpermanent, distdata, &isDeletable ) );
         ncands++;
         if( isDeletable )
         {
            if( genstar.replacedge_head != -1 )
            {
               npseudoelims++;
               continue;
            }

            extreduce_edgeRemove(scip, i, graph, distdata, extpermanent);
            (*nelims)++;
         }
      }
   }

   generalStarExit(scip, graph, &genstar);

   printf("ncands=%d \n", ncands);
   printf("npseudoelims=%d \n", npseudoelims);

   printf("number of general star eliminations=%d \n", *nelims);
   return SCIP_OKAY;
}

/** initializes */
static inline
SCIP_RETCODE pseudodeleteInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure  */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims,             /**< number of eliminations (out) */
   EXTPSEUDO*            extpseudo           /**< to initialize */
   )
{
   SCIP_Real* cutoffs;
   int* nodestouches;
   const int nnodes = graph_get_nNodes(g);

   assert(nelims && extpseudo);

   *nelims = 0;


   SCIP_CALL( SCIPallocBufferArray(scip, &cutoffs, STP_DELPSEUDO_MAXNEDGES) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodestouches, nnodes) );

   for( int i = 0; i < nnodes; ++i )
      nodestouches[i] = 0;

   extpseudo->cutoffIsPromising = FALSE;
   extpseudo->cutoffIsComputed = FALSE;
   extpseudo->node = -1;
   extpseudo->offsetp = offsetp;
   extpseudo->nelimsp = nelims;
   extpseudo->nodestouches = nodestouches;
   extpseudo->cutoffs = cutoffs;

   return SCIP_OKAY;
}


/** frees */
static inline
void pseudodeleteExit(
   SCIP*                 scip,               /**< SCIP data structure */
   EXTPSEUDO*            extpseudo           /**< to initialize */
   )
{
   SCIPfreeBufferArray(scip, &(extpseudo->nodestouches));
   SCIPfreeBufferArray(scip, &(extpseudo->cutoffs));
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



/** computes cutoffs for pseudo-deletion of given node */
static inline
SCIP_RETCODE pseudodeleteDeleteComputeCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             checkpromising,     /**< check whether promising? */
   SCIP_Bool             abortDeg3,          /**< abort for degree 3? */
   DISTDATA*             distdata,           /**< distance data */
   int                   node,               /**> to be deleted */
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
               //const SCIP_Real dist_dbg = extreduce_distComputeRestrictedDist(scip, graph, node, distdata, vert, vert2);
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
   int                   node,               /**> to be deleted */
   REDCOST*              redcostdata,        /**< reduced cost data */
   DISTDATA*             distdata,           /**< distance data */
   GRAPH*                graph,              /**< graph data structure */
   EXTPSEUDO*            extpseudo           /**< data */
)
{
   const SCIP_Real* cutoffs = extpseudo->cutoffs;
   int* const nodestouches = extpseudo->nodestouches;
   SCIP_Real* redcost = redcostdata->redEdgeCost;
   SCIP_Real prize = -1.0;
   SCIP_Bool rpc3term = FALSE;
   SCIP_Bool success;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   assert(redcost);
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
   SCIP_CALL( graph_knot_delPseudo(scip, graph, graph->cost, cutoffs, NULL, node, redcost, &success) );

   if( success )
   {
      (*(extpseudo->nelimsp))++;
      graph->mark[node] = FALSE;

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

   for( int k = 0; k < nnodes; ++k )
   {
      if( pseudoDelNodes[k] && 2 <= graph->grad[k] && graph->grad[k] <= STP_DELPSEUDO_MAXGRAD )
      {
         assert(!Is_term(graph->term[k]));

         SCIP_CALL( pseudodeleteDeleteComputeCutoffs(scip, FALSE, FALSE, distdata, k, graph, extpseudo) );
         SCIP_CALL( pseudodeleteDeleteNode(scip, k, redcostdata, distdata, graph, extpseudo) );
      }
   }

   return SCIP_OKAY;
}


/** Extended reduction test for pseudo-eliminating nodes.
 *  Applies pseudo-elimination if pseudoDelNodes != NULL, or just marks them otherwise */
static
SCIP_RETCODE pseudodeleteExecute(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Bool*      pseudoDelNodes,     /**< nodes to pseudo-eliminate already */
   REDCOST*              redcostdata,        /**< reduced cost data */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nnodes = graph_get_nNodes(graph);
   STAR* stardata;
   DISTDATA distdata;
   EXTPERMA extpermanent;
   EXTPSEUDO extpseudo;

   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);
   assert(graph_isMarked(graph));

   extpermanent.redcostEqualAllow = FALSE;

   SCIP_CALL( pseudodeleteInit(scip, graph, offsetp, nelims, &extpseudo) );
   SCIP_CALL( extInit(scip, !graph_pc_isPc(graph), graph, NULL, &distdata, &extpermanent) );

   if( pseudoDelNodes )
   {
      SCIP_CALL( pseudodeleteDeleteMarkedNodes(scip, pseudoDelNodes, redcostdata, &distdata, graph, &extpseudo) );
      SCIPdebugMessage("number of eliminations after initial pseudo-elimination %d \n", *nelims);
   }

   SCIP_CALL( reduce_starInit(scip, EXT_PSEUDO_DEGREE_MAX, &stardata) );

   /* pseudo-elimination loop */
   for( int i = 0; i < nnodes; ++i )
   {
      SCIP_Bool nodeisDeletable = FALSE;

      if( pseudodeleteNodeIsPromising(graph, i) )
      {
         SCIP_CALL( pseudodeleteDeleteComputeCutoffs(scip, TRUE, TRUE, &distdata, i, graph, &extpseudo) );

         if( extpseudo.cutoffIsPromising )
         {
            SCIP_CALL( extreduce_checkNode(scip, graph, redcostdata, i, stardata, &distdata, &extpermanent, &nodeisDeletable) );
         }
      }

      if( !nodeisDeletable )
         continue;

      SCIP_CALL( pseudodeleteDeleteNode(scip, i, redcostdata, &distdata, graph, &extpseudo) );
   }

   reduce_starFree(scip, &stardata);
   extFree(scip, graph, &distdata, &extpermanent);
   pseudodeleteExit(scip, &extpseudo);

   SCIPdebugMessage("number of eliminations after extended pseudo-elimination %d \n", *nelims);

   /* todo in I015 we get isolated vertices...not sure whether this is a bug or normal behaviour */
   SCIP_CALL( reduceLevel0(scip, graph));

   assert(graphmarkIsClean(redcostdata, graph));

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
      SCIP_CALL_ABORT( reduce_sdRepair(scip, edge, graph, distdata->sdistdata) );
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
   DISTDATA distdata;
   EXTPERMA extpermanent;
   SCIP_Bool withSol = (result != NULL);

   assert(scip && redcostdata && edgedeletable);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, FALSE, graph, NULL, &distdata, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, redcostdata, e) )
      {
         const int erev = e + 1;
         extpermanent.redcostEqualAllow = (withSol && result[e] != CONNECT && result[erev] != CONNECT);
         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( !edgedeletable[e] )
         {
            SCIP_Bool deletable;
            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, e, &distdata, &extpermanent,
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

            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, erev, &distdata, &extpermanent,
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

            removeEdge(scip, e, graph, &distdata, &extpermanent);

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
   const SCIP_Bool useSd = !graph_pc_isPc(graph);
   const int nedges = graph_get_nEdges(graph);
   DISTDATA distdata;
   EXTPERMA extpermanent;
   SCIP_Bool withSol = (result != NULL);

   assert(scip && redcostdata);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, useSd, graph, edgedeletable, &distdata, &extpermanent) );

   if( useSd )
   {
      assert(distdata.sdistdata);
      SCIP_CALL( reduce_sdRepairSetUp(scip, graph, distdata.sdistdata) );
   }

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, redcostdata, e) )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         extpermanent.redcostEqualAllow = (withSol && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, e, &distdata, &extpermanent, &deletable) );

         if( deletable )
         {
            removeEdge(scip, e, graph, &distdata, &extpermanent);

            if( withSol && (result[e] == CONNECT || result[erev] == CONNECT ) )
               withSol = FALSE;

            (*nelims)++;
         }
      }
   }

  // printf("number of extended edge eliminations=%d \n", *nelims);

/*
   if( 1 )
   {
      int ngenstarelims = 0;
      SCIP_CALL( generalStarDeleteEdges(scip, redcostdata, &extpermanent, graph, &distdata, &ngenstarelims) );
      *nelims += ngenstarelims;
   }
*/


   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}

/** extended reduction test for pseudo-eliminating nodes */
SCIP_RETCODE extreduce_pseudoDeleteNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Bool*      pseudoDelNodes,     /**< nodes to pseudo-eliminate already */
   REDCOST*              redcostdata,        /**< reduced cost data */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   SCIP_Real*            offsetp,            /**< pointer to store offset */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   SCIP_Bool isExtendedOrg;

   assert(scip && redcostdata && graph && nelims);

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   isExtendedOrg = graph->extended;
   if( graph_pc_isPc(graph) )
      graph_pc_2orgcheck(scip, graph);

   /* call actual method */
   SCIP_CALL( pseudodeleteExecute(scip, pseudoDelNodes, redcostdata, graph, offsetp, nelims) );

   if( graph_pc_isPc(graph) && isExtendedOrg != graph->extended )
      graph_pc_2trans(scip, graph);

   return SCIP_OKAY;
}


/** deletes center edges of general stars */
SCIP_RETCODE extreduce_deleteGeneralStars(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const SCIP_Bool useSd = !graph_pc_isPc(graph);
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, useSd, graph, edgedeletable, &distdata, &extpermanent) );

   if( useSd )
   {
      assert(distdata.sdistdata);
      SCIP_CALL( reduce_sdRepairSetUp(scip, graph, distdata.sdistdata) );
   }

   SCIP_CALL( generalStarDeleteEdges(scip, redcostdata, &extpermanent, graph, &distdata, nelims) );

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

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
         .genstar_centeredge = -1,
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
                          .nextleaves = 1, .ncompedges = 1, .genstar_centeredge = -1,
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
                          .nextleaves = -1, .ncompedges = -1, .genstar_centeredge = -1,
                          .comproot = -1, .allowReversion = TRUE };

      pseudodeleteGetNextStar(graph, stardata, &extcomp);

      *isPseudoDeletable = FALSE;
      SCIP_CALL( extreduce_checkComponent(scip, graph, redcostdata, &extcomp, distdata, extpermanent, isPseudoDeletable) );

      if( pseudodeleteAllStarsChecked(stardata) )
         break;
   }

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

   if( EQ(0.0, redcostdata->redEdgeCost[e]) && EQ(0.0, redcostdata->redEdgeCost[flipedge(e)]) )
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

   assert(maxdepth > 0);

   return maxdepth;
}
