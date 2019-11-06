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

/**@file   reduce_ext.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements extended reduction techniques for several Steiner problems.
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
#include "reduce.h"

#define EXT_ANCESTORS_MAX  16
#define EXT_STATE_NONE     0
#define EXT_STATE_EXPANDED 1
#define EXT_STATE_MARKED   2
#define EXT_REDCOST_NRECOMP 10

#define STP_EXT_MAXDFSDEPTH 6
#define STP_EXT_MINDFSDEPTH 4
#define STP_EXT_MAXGRAD 8
#define STP_EXT_MINGRAD 6
#define STP_EXT_EDGELIMIT 50000


/** deletes an edge and makes corresponding adaptations */
static
void removeEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph               /**< graph data structure */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];
   graph_edge_del(scip, graph, edge, TRUE);

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
}

/** mark ancestors of given edge */
static
SCIP_Bool markAncestorsConflict(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      const unsigned idx = ((unsigned) curr->index) / 2;

      assert(curr->index >= 0 && idx < (unsigned) (MAX(graph->edges, graph->orgedges) / 2));

      if( ancestormark[idx] )
         return TRUE;

      ancestormark[idx] = 1;
   }

   return FALSE;
}

/** mark ancestors of given edge */
static
void markAncestors(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      assert((ancestormark[((unsigned) curr->index) / 2]) == 0);

      ancestormark[((unsigned) curr->index) / 2] = 1;
   }
}

/** unmark ancestors of given edge */
static
void unmarkAncestorsConflict(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      ancestormark[((unsigned) curr->index) / 2] = 0;
   }
}

/** unmark ancestors of given edge */
static
void unmarkAncestors(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      const unsigned idx = ((unsigned) curr->index) / 2;

      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      assert(ancestormark[idx] == 1);

      ancestormark[idx] = 0;
   }
}

/** finalize subtree computations (clean up, update global bound)  */
static
void finalizeSubtree(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            edgeends,           /**< heads or tail of edge */
   const int*            treeedges,          /**< tree edges */
   int                   dfsdepth,           /**< dfs depth */
   SCIP_Bool             stopped,            /**< has extension been stopped? */
   SCIP_Real             localbound,         /**< local bound */
   SCIP_Real*            globalbound,        /**< points to global bound */
   SCIP_Bool*            ruleout,            /**< rule out? */
   int*                  nodepos,            /**< node position in tree */
   int*                  ancestormark        /**< ancestor mark array */
)
{

   if( localbound > *globalbound )
        *globalbound = localbound;

   if( !stopped )
      *ruleout = TRUE;

   if( !(*ruleout) )
   {
      for( int etree = 0; etree < dfsdepth; etree++ )
      {
         const int node = edgeends[treeedges[etree]];
         assert(nodepos[node]);
         nodepos[node] = 0;
         unmarkAncestors(graph, treeedges[etree], ancestormark);
      }
   }
   else
   {
      assert(dfsdepth == 0);
   }
}

/** should we truncate subtree? */
static
SCIP_Bool truncateSubtree(
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             extendedcost,       /**< cost of subtree */
   int                   root,               /**< DA root */
   int                   currhead,           /**< latest node */
   int                   maxgrad,            /**< maximum allowed degree */
   int                   dfsdepth,           /**< dfs depth */
   int                   maxdfsdepth,        /**< max dfs depth */
   SCIP_Real*            minbound,           /**< bound */
   SCIP_Bool*            stopped             /**< real truncation? */
)
{
   assert(graph->mark[currhead]);

   if( Is_term(graph->term[currhead]) || graph->grad[currhead] > maxgrad || dfsdepth >= maxdfsdepth )
   {
      /* real truncation? */
      if( root != currhead )
      {
         *stopped = TRUE;
         if( extendedcost < *minbound )
            *minbound = extendedcost;
      }
      return TRUE;
   }
   return FALSE;
}

/** remove latest edge from subtree and returns new tree cost */
static
SCIP_Real shortenSubtree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      redcost,            /**< reduced costs */
   const int*            treeedges,          /**< edges of tree */
   SCIP_Real             treecostold,        /**< old tree cost */
   SCIP_Real             treecostoffset,     /**< tree cost offset */
   int                   dfsdepth,           /**< DFS depth */
   int                   lastnode,           /**< last node */
   int*                  nodepos,            /**< array to mark node position*/
   int*                  ancestormark,       /**< ancestor mark array */
   unsigned int*         nstallrounds        /**< rounds since latest update */
)
{
   const int lastedge = treeedges[dfsdepth - 1];

   assert(dfsdepth >= 1 && lastnode >= 0);
   assert(SCIPisGE(scip, treecostold, 0.0) && treecostoffset >= 0.0);

   nodepos[lastnode] = 0;
   unmarkAncestors(graph, lastedge, ancestormark);

   /* recompute tree cost? */
   if( (*nstallrounds)++ >= 9 )
   {
      SCIP_Real treecost = treecostoffset;

      *nstallrounds = 0;

      for( int i = 0; i < dfsdepth - 1; i++ )
         treecost += redcost[treeedges[i]];

#if 0 // todo deleteeme
      if( !SCIPisEQ(scip, treecost, treecostold - redcost[lastedge]) )
      {
         printf("%.15f %.15f \n", treecost, treecostold - redcost[lastedge]);
         assert(0);
         exit(1);
      }
#endif

      assert(SCIPisEQ(scip, treecost, treecostold - redcost[lastedge]));
   }

   return (treecostold - redcost[lastedge]);
}

/** extend subtree and return new nadded_edges */
static
int extendSubtreeHead(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      redcost,            /**< reduced costs */
   int                   curredge,           /**< added edges */
   int                   currhead,           /**< latest node */
   int                   dfsdepth,           /**< DFS depth*/
   int                   nadded_edges,       /**< added edges */
   SCIP_Real*            treecost,           /**< pointer to treecost */
   SCIP_Real*            treecosts,          /**< edge costs of tree */
   int*                  nodepos,            /**< node position in tree */
   int*                  treeedges,          /**< edges of tree */
   int*                  edgestack,          /**< stack */
   int*                  ancestormark        /**< ancestor mark array */
)
{
   int n = nadded_edges + 1;
   nodepos[currhead] = dfsdepth + 1;
   *treecost += redcost[curredge];
   treeedges[dfsdepth] = curredge;
   treecosts[dfsdepth] = graph->cost[curredge];

   markAncestors(graph, curredge, ancestormark);
   assert(graph->mark[currhead]);

   for( int e = graph->outbeg[currhead]; e != EAT_LAST; e = graph->oeat[e] )
   {
      const int head = graph->head[e];
      if( !nodepos[head] && graph->mark[head] )
         edgestack[n++] = e;
   }

   return n;
}

/** extend subtree and return new nadded_edges */
static
int extendSubtreeTail(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      redcost,            /**< reduced costs */
   int                   curredge,           /**< added edges */
   int                   currtail,           /**< latest node */
   int                   dfsdepth,           /**< DFS depth*/
   int                   nadded_edges,       /**< added edges */
   SCIP_Real*            treecost,           /**< pointer to treecost */
   SCIP_Real*            treecosts,          /**< edge costs of tree */
   int*                  nodepos,            /**< node position in tree */
   int*                  treeedges,          /**< edges of tree */
   int*                  edgestack,          /**< stack */
   int*                  ancestormark        /**< ancestor mark array */
)
{
   int n = nadded_edges + 1;
   nodepos[currtail] = dfsdepth + 1;
   *treecost += redcost[curredge];
   treeedges[dfsdepth] = curredge;
   treecosts[dfsdepth] = graph->cost[curredge];

   markAncestors(graph, curredge, ancestormark);
   assert(graph->mark[currtail]);

   for( int e = graph->inpbeg[currtail]; e != EAT_LAST; e = graph->ieat[e] )
   {
      const int tail = graph->tail[e];
      if( !nodepos[tail] && graph->mark[tail] )
         edgestack[n++] = e;
   }

   return n;
}


/** small helper */
static
void updateEqArrays(
   int                   edge,               /**< the edge */
   unsigned int*         eqstack,            /**< stores edges that were used for equality comparison */
   int*                  eqstack_size,       /**< pointer to size of eqstack */
   SCIP_Bool*            eqmark              /**< marks edges that were used for equality comparison */
)
{
   const unsigned int halfedge = ((unsigned int) edge) / 2;

   assert(edge >= 0);

   if( !eqmark[halfedge]  )
   {
      eqmark[halfedge] = TRUE;
      eqstack[*eqstack_size] = halfedge;

      *eqstack_size = *eqstack_size + 1;
   }
}


/** compare with tree bottleneck */
static
SCIP_Bool bottleneckRuleOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      treecosts,          /**< costs of edges of the tree */
   const SCIP_Real*      basebottlenecks,    /**< bottleneck costs for innode and basenode */
   const int*            nodepos,            /**< node position in tree */
   const int*            treeedges,          /**< edges in tree (corresponding to treecosts) */
   SCIP_Real             orgedgecost,        /**< cost of original edge */
   SCIP_Real             extedgecost,        /**< cost of short-cut edge */
   int                   orgedge,            /**< original edge */
   int                   extedge,            /**< short-cut edge */
   int                   dfsdepth,           /**< dfs depth */
   SCIP_Bool             allow_eq,           /**< test for equality? */
   unsigned int*         eqstack,            /**< stores edges that were used for equality comparison */
   int*                  eqstack_size,       /**< pointer to size of eqstack */
   SCIP_Bool*            eqmark              /**< marks edges that were used for equality comparison */
)
{
   int start;

   assert(!allow_eq ||(eqstack != NULL && eqstack_size != NULL && eqmark != NULL));

   if( allow_eq )
   {
      if( SCIPisLE(scip, extedgecost, orgedgecost) )
      {
         if( SCIPisEQ(scip, extedgecost, orgedgecost) )
            updateEqArrays(orgedge, eqstack, eqstack_size, eqmark);

         return TRUE;
      }
   }
   else if( SCIPisLT(scip, extedgecost, orgedgecost) )
      return TRUE;

   if( nodepos[graph->tail[extedge]] > graph->knots )
   {
      start = nodepos[graph->tail[extedge]] - 1 - graph->knots;
      assert(start >= 0 && start <= 3);

      if( start <= 2 )
      {
         assert(basebottlenecks[start] > 0);

         if( SCIPisLT(scip, extedgecost, basebottlenecks[start]) )
            return TRUE;
      }
      start = 0;
   }
   else
   {
      start = nodepos[graph->tail[extedge]];  /* not -1! We save the incoming bottlenecks */
      assert(start >= 1 && start <= dfsdepth);
      assert(start < dfsdepth || graph->tail[orgedge] == graph->tail[extedge]);

   }

   for( int i = start; i < dfsdepth; i++ )
   {
      assert(treecosts[i] >= 0.0);

      if( allow_eq )
      {
         if( SCIPisLE(scip, extedgecost, treecosts[i]) )
         {
            if( SCIPisEQ(scip, extedgecost, treecosts[i]) )
               updateEqArrays(treeedges[i], eqstack, eqstack_size, eqmark);

            return TRUE;
         }
      }
      else if( SCIPisLT(scip, extedgecost, treecosts[i]) )
            return TRUE;
   }

   return FALSE;
}

/** can subtree be ruled out? */
static
SCIP_Bool ruleOutSubtree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      treecosts,          /**< costs of edges of the tree */
   const SCIP_Real*      basebottlenecks,    /**< bottleneck costs for innode and basenode */
   const int*            nodepos,            /**< node position in tree */
   const int*            treeedges,          /**< edges in tree (corresponding to treecosts) */
   SCIP_Real             cutoff,             /**< cut-off bound */
   SCIP_Real             extendedcost,       /**< cost of subtree */
   int                   dfsdepth,           /**< dfs depth */
   int                   curredge,           /**< latest edge */
   SCIP_Bool             allowequality,      /**< rule out also in case of equality */
   const int*            ancestormark,       /**< ancestor mark array */
   unsigned int*         eqstack,            /**< stores edges that were used for equality comparison */
   int*                  eqstack_size,       /**< pointer to size of eqstack */
   SCIP_Bool*            eqmark              /**< marks edges that were used for equality comparison */
)
{

   if( allowequality ? (extendedcost >= cutoff) : SCIPisGT(scip, extendedcost, cutoff) )
   {
      return TRUE;
   }
   else
   {
      const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
      SCIP_Real currcost;
      int currhead;

      if( ancestormark != NULL )
      {
         int count = 0;
         for( IDX* curr = graph->ancestors[curredge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
         {
            assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
            if( ancestormark[((unsigned) curr->index) / 2] )
               return TRUE;
         }
      }

      currcost = graph->cost[curredge];
      currhead = graph->head[curredge];

      for( int e = graph->inpbeg[currhead]; e != EAT_LAST; e = graph->ieat[e] )
      {
         int tail;
         SCIP_Real ecost;

         if( e == curredge )
            continue;

         assert(flipedge(e) != curredge);

         tail = graph->tail[e];
         ecost = graph->cost[e];

         if( nodepos[tail] )
         {
            const SCIP_Bool allow_eq = (eqmark != NULL && eqmark[((unsigned) e) / 2] == FALSE);
            assert(graph->mark[tail]);

            if( bottleneckRuleOut(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, currcost,
                  ecost, curredge, e, dfsdepth, allow_eq, eqstack, eqstack_size, eqmark) )
               return TRUE;
         }
         else
         {
            const SCIP_Bool is_term = Is_term(graph->term[tail]);

            for( int e2 = graph->outbeg[tail], count = 0; e2 != EAT_LAST && count < STP_EXT_MAXGRAD;
                  e2 = graph->oeat[e2], count++ )
            {
               int head2;

               if( e2 == e )
                  continue;

               assert(flipedge(e2) != e && e2 != curredge && e2 != flipedge(curredge) );

               head2 = graph->head[e2];
               assert(head2 != tail);

               if( nodepos[head2] )
               {
                  SCIP_Real extcost = is_term ? MAX(ecost, graph->cost[e2]) : (ecost + graph->cost[e2]);
                  const SCIP_Bool allow_eq = (eqmark != NULL && eqmark[((unsigned) e) / 2] == FALSE && eqmark[((unsigned) e2) / 2] == FALSE);
                  assert(graph->mark[head2]);

                  if( pcmw && is_term )
                     extcost = MAX(ecost + graph->cost[e2] - graph->prize[tail], extcost);

                  assert(head2 != currhead);

                  if( bottleneckRuleOut(scip, graph, treecosts, basebottlenecks,
                        nodepos, treeedges, currcost, extcost, curredge, flipedge(e2), dfsdepth, allow_eq, eqstack, eqstack_size, eqmark) )
                  {
                      return TRUE;
                  }
               }
            }
         }
      }
   }

   return FALSE;
}

/** check (directed) edge */
static
SCIP_RETCODE reduceCheckEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   root,               /**< graph root from dual ascent */
   const SCIP_Real*      redcost,            /**< reduced costs */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   SCIP_Real             cutoff,             /**< cutoff value */
   int                   edge,               /**< (directed) edge to be checked */
   SCIP_Bool             equality,           /**< allow equality? */
   int*                  edgestack,          /**< array of size nodes for internal computations */
   SCIP_Bool*            deletable,          /**< is edge deletable? */
   unsigned int*         eqstack,            /**< stores edges that were used for equality comparison */
   int*                  eqstack_size,       /**< pointer to size of eqstack */
   SCIP_Bool*            eqmark              /**< marks edges that were used for equality comparison */
)
{
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   const int maxgrad = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINGRAD : STP_EXT_MAXGRAD;
   SCIP_Real edgebound = redcost[edge] + pathdist[tail] + vnoi[head].dist;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(redcost != NULL);
   assert(pathdist != NULL);
   assert(vnoi != NULL);
   assert(edge >= 0 && edge < graph->edges);

#ifndef NDEBUG
   if( !graph_pc_isPcMw(graph) )
      for( int k = 0; k < graph->knots; k++ )
         assert(graph->mark[k] == (graph->grad[k] > 0));

   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);
#endif

   *deletable = FALSE;

   /* trivial rule-out? */
   if( SCIPisGT(scip, edgebound, cutoff) || (equality && SCIPisEQ(scip, edgebound, cutoff)) || head == root )
   {
      *deletable = TRUE;
   }
   else if( (graph->grad[tail] <= maxgrad && !Is_term(graph->term[tail]))
         || (graph->grad[head] <= maxgrad && !Is_term(graph->term[head])) )
   {
      SCIP_Real treecosts[STP_EXT_MAXDFSDEPTH + 1];
      int treeedges[STP_EXT_MAXDFSDEPTH + 1];
      SCIP_Real basebottlenecks[3];
      int* nodepos;
      int* ancestormark;
      const int nnodes = graph->knots;
      const int maxdfsdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;

      /* allocate clean arrays */
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &nodepos, nnodes) );
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &ancestormark, (MAX(graph->edges, graph->orgedges) / 2)) );

      basebottlenecks[0] = graph->cost[edge];
#ifndef NEBUG
      basebottlenecks[1] = -1.0;
      basebottlenecks[2] = -1.0;

      for( int i = 0; i < STP_EXT_MAXDFSDEPTH + 1; i++ )
      {
         treecosts[i] = -1.0;
         treeedges[i] = -1;
      }
#endif

      markAncestors(graph, edge, ancestormark);

      /* check whether to extend from head */
      if( !Is_term(graph->term[head]) && graph->grad[head] <= maxgrad )
      {
         const SCIP_Real treecostoffset = redcost[edge] + pathdist[tail];
         SCIP_Real minbound = FARAWAY;
         SCIP_Real treecost = treecostoffset;
         unsigned int rounds = 0;
         int dfsdepth = 0;
         int nadded_edges = 0;
         SCIP_Bool stopped = FALSE;

         nodepos[tail] = nnodes + 1;
         nodepos[head] = nnodes + 4;

         for( int e = graph->outbeg[head]; e != EAT_LAST; e = graph->oeat[e] )
         {
            const int head2 = graph->head[e];
            if( head2 != tail && graph->mark[head2] )
               edgestack[nadded_edges++] = e;
         }

         /* limited DFS */
         while( nadded_edges > 0 )
         {
            const int curredge = edgestack[--nadded_edges];
            const int currhead = graph->head[curredge];

            assert(graph->mark[currhead]);

            /*  subtree already processed? */
            if( nodepos[currhead] )
            {
               assert(dfsdepth >= 1 && treeedges[dfsdepth - 1] == curredge);

               treecost = shortenSubtree(scip, graph, redcost, treeedges, treecost, treecostoffset, dfsdepth,
                     currhead, nodepos, ancestormark, &rounds);

               dfsdepth--;
            }
            else
            {
               const SCIP_Real extendedcost = treecost + redcost[curredge] + vnoi[currhead].dist;

               if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, cutoff, extendedcost, dfsdepth, curredge, equality,
                     ancestormark, eqstack, eqstack_size, eqmark) )
                  continue;

               if( truncateSubtree(graph, extendedcost, root, currhead, maxgrad, dfsdepth, maxdfsdepth, &minbound, &stopped) )
               {
                  /* stopped and no further improvement of bound possible? */
                 if( stopped && minbound <= edgebound )
                    break;

                  continue;
               }

               nadded_edges = extendSubtreeHead(scip, graph, redcost, curredge, currhead, dfsdepth, nadded_edges, &treecost, treecosts, nodepos,
                     treeedges, edgestack, ancestormark);
               dfsdepth++;
            }
         } /* DFS loop */

         assert(stopped || minbound == FARAWAY);
         assert(SCIPisGE(scip, minbound, redcost[edge] + pathdist[tail] + vnoi[head].dist));

         finalizeSubtree(graph, graph->head, treeedges, dfsdepth, stopped, minbound, &edgebound, deletable, nodepos, ancestormark);
      } /* extend from head */

      /* check whether to extend from tail */
      if( !(*deletable) && !Is_term(graph->term[tail]) && graph->grad[tail] <= maxgrad )
      {
         const SCIP_Real treecostoffset = edgebound - pathdist[tail];
         SCIP_Real minbound = FARAWAY;
         SCIP_Real treecost = treecostoffset;
         int dfsdepth = 0;
         int nadded_edges = 0;
         unsigned int rounds = 0;
         SCIP_Bool stopped = FALSE;

         nodepos[head] = nnodes + 1;
         nodepos[tail] = nnodes + 4;

         for( int e = graph->inpbeg[tail]; e != EAT_LAST; e = graph->ieat[e] )
         {
            const int tail2 = graph->tail[e];
            if( tail2 != head && graph->mark[tail2] )
               edgestack[nadded_edges++] = e;
         }

         /* limited DFS */
         while( nadded_edges > 0 )
         {
            const int curredge = edgestack[--nadded_edges];
            const int currtail = graph->tail[curredge];
            assert(graph->mark[currtail]);

            /*  subtree already processed? */
            if( nodepos[currtail] )
            {
               assert(dfsdepth >= 1 && treeedges[dfsdepth - 1] == curredge);

               treecost = shortenSubtree(scip, graph, redcost, treeedges, treecost, treecostoffset, dfsdepth,
                     currtail, nodepos, ancestormark, &rounds);

               dfsdepth--;
            }
            else
            {
               const SCIP_Real extendedcost = treecost + redcost[curredge] + pathdist[currtail];

               if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, cutoff, extendedcost, dfsdepth, flipedge((unsigned)curredge), equality,
                     ancestormark, eqstack, eqstack_size, eqmark) )
                  continue;

               if( truncateSubtree(graph, extendedcost, -1, currtail, maxgrad, dfsdepth, maxdfsdepth, &minbound, &stopped) )
               {
                  if( stopped )
                     break;
                  continue;
               }

               nadded_edges = extendSubtreeTail(graph, redcost, curredge, currtail, dfsdepth, nadded_edges, &treecost, treecosts, nodepos, treeedges,
                     edgestack, ancestormark);
               dfsdepth++;
            }
         } /* DFS loop */

         assert(stopped || minbound == FARAWAY);
         assert(SCIPisGE(scip, minbound, redcost[edge] + pathdist[tail] + vnoi[head].dist));

         finalizeSubtree(graph, graph->tail, treeedges, dfsdepth, stopped, minbound, &edgebound, deletable, nodepos, ancestormark);
      } /* extend from tail */

      /* clean arrays */
      nodepos[head] = 0;
      nodepos[tail] = 0;
      unmarkAncestors(graph, edge, ancestormark);

      /* free memory */
      for( int i = 0; i < nnodes; i++ )
         assert(nodepos[i] == 0);
      for( int i = 0; i < MAX(graph->edges, graph->orgedges) / 2; i++ )
         assert(ancestormark[i] == 0);

      SCIPfreeCleanBufferArray(scip, &ancestormark);
      SCIPfreeCleanBufferArray(scip, &nodepos);
   }

   return SCIP_OKAY;
}


/** todo don't use this function, but check for conflicts in misc_stp.c and handle them */
SCIP_RETCODE reduce_deleteConflictEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   int* edgemark;
   const int nedges = g->edges;
   const int nancestors = MAX(g->edges, g->orgedges) / 2;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);

   assert(scip != NULL && g != NULL);
   assert(g->ancestors != NULL);
   assert(!graph_pc_isPcMw(g) || !g->extended);

   SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, nancestors) );

   for( int e = 0; e < nancestors; e++ )
      edgemark[e] = FALSE;

   for( int e = 0; e < nedges; e += 2 )
   {
      SCIP_Bool conflict;

      if( g->oeat[e] == EAT_FREE )
         continue;

      if( pcmw && g->cost[e] != g->cost[e + 1] )
         continue;

      conflict = markAncestorsConflict(g, e, edgemark);

      unmarkAncestorsConflict(g, e, edgemark);
      if( conflict )
      {
         conflict = markAncestorsConflict(g, e + 1, edgemark);
         unmarkAncestorsConflict(g, e + 1, edgemark);

         if( conflict )
         {
            assert(!graph_pc_isPcMw(g) || (!Is_pseudoTerm(g->term[g->head[e]]) && !Is_pseudoTerm(g->term[g->tail[e]])));
            graph_edge_del(scip, g, e, TRUE);
         }
      }
   }

   SCIPfreeBufferArray(scip, &edgemark);

   return SCIP_OKAY;
}

/** check 3-tree */
SCIP_RETCODE reduce_extendedCheck3Tree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   root,               /**< graph root from dual ascent */
   const SCIP_Real*      redcost,            /**< reduced costs */
   const SCIP_Real*      pathdist,           /**< shortest path distances  */
   const PATH*           vnoi,               /**< Voronoi paths  */
   const int*            vbase,              /**< bases to Voronoi paths */
   SCIP_Real             cutoff,             /**< cutoff value */
   const int*            outedges,           /**< two outgoing edge */
   int                   inedge,             /**< incoming edge */
   int*                  edgestack,          /**< array of size nodes for internal computations */
   SCIP_Real*            treebound,          /**< to store a lower bound for the tree */
   SCIP_Bool*            ruleout,             /**< could tree be ruled out? */
   unsigned int*         eqstack,            /**< stores edges that were used for equality comparison */
   int*                  eqstack_size,       /**< pointer to size of eqstack */
   SCIP_Bool*            eqmark              /**< marks edges that were used for equality comparison */
)
{
#ifndef NDEBUG
   const SCIP_Real orgtreebound = *treebound;
#endif
   assert(scip != NULL);
   assert(graph != NULL);
   assert(redcost != NULL);
   assert(ruleout != NULL);
   assert(pathdist != NULL);
   assert(vnoi != NULL);
   assert(vbase != NULL);
   assert(outedges != NULL);
   assert(inedge >= 0 && outedges[0] >= 0 && outedges[1] >= 0);
   assert(inedge < graph->edges && outedges[0] < graph->edges && outedges[1] < graph->edges);

   /* trivial rule-out? */
   if( SCIPisGT(scip, *treebound, cutoff) )
   {
      *ruleout = TRUE;
   }
   else
   {
      SCIP_Real treecosts[STP_EXT_MAXDFSDEPTH + 1];
      int treeedges[STP_EXT_MAXDFSDEPTH + 1];
      SCIP_Real basebottlenecks[3];
      const SCIP_Real basecost = redcost[inedge] + redcost[outedges[0]] + redcost[outedges[1]] + pathdist[graph->tail[inedge]];
      int* nodepos;
      int* ancestormark;
      const int nnodes = graph->knots;
      const int innode = graph->tail[inedge];
      const int basenode = graph->head[inedge];
      const int maxdfsdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;
      const int maxgrad = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINGRAD : STP_EXT_MAXGRAD;

      assert(basenode == graph->tail[outedges[0]] && basenode == graph->tail[outedges[1]]);

      /* allocate clean array */
      SCIP_CALL(SCIPallocCleanBufferArray(scip, &nodepos, nnodes));
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &ancestormark, (MAX(graph->edges, graph->orgedges) / 2)) );

      nodepos[basenode] = nnodes + 1;            /* basenode */
      nodepos[innode] = nnodes + 2;              /* innode */

#ifndef NDEBUG
      for( int i = 0; i < STP_EXT_MAXDFSDEPTH + 1; i++ )
      {
         treecosts[i] = -1.0;
         treeedges[i] = -1;
      }
#endif

      *ruleout = FALSE;

      /* can we rule out subgraph beforehand? */
      if( markAncestorsConflict(graph, inedge, ancestormark)
        || markAncestorsConflict(graph, outedges[0], ancestormark)
        || markAncestorsConflict(graph, outedges[1], ancestormark) )
      {
         *ruleout = TRUE;
      }

      /* main loop: extend tree and try to rule it out */
      for( int i = 0; i < 2 && !(*ruleout); i++ )
      {
         SCIP_Real minbound = FARAWAY;
         SCIP_Real treecost = basecost;
         const int startedge = outedges[i];
         const int costartedge = outedges[(i + 1) % 2];
         const int startnode = graph->head[startedge];
         const int costartnode = graph->head[costartedge];
         int dfsdepth = 0;
         int nadded_edges = 0;
         SCIP_Bool stopped = FALSE;
         unsigned int rounds = 0;

         /* try to extend the tree from startnode */

         assert(startnode != root && costartnode != root);

         /* for basenode */
         basebottlenecks[0] = graph->cost[startedge];

         /* for innode */
         basebottlenecks[1] = MAX(graph->cost[startedge], graph->cost[inedge]);

         /* for costartnode */
         basebottlenecks[2] = MAX(graph->cost[startedge], graph->cost[costartedge]);

         nodepos[costartnode] = nnodes + 3;
         nodepos[startnode] = nnodes + 4;

         /* can we rule out entire subtree already? */
         if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, FARAWAY, 0.0, 0, startedge, FALSE,
               NULL, NULL, NULL, NULL) )
         {
            *ruleout = TRUE;
            break;
         }

         /* cannot extend over terminals or vertices of high degree */
         if( Is_term(graph->term[startnode]) || graph->grad[startnode] > maxgrad )
            continue;

         assert(nodepos[startnode]);

         for( int e = graph->outbeg[startnode]; e != EAT_LAST; e = graph->oeat[e] )
            if( !nodepos[graph->head[e]] )
               edgestack[nadded_edges++] = e;

         /* limited DFS starting from startnode */
         while ( nadded_edges > 0 )
         {
            const int curredge = edgestack[--nadded_edges];
            const int currhead = graph->head[curredge];

            /*  subtree already processed? */
            if( nodepos[currhead] )
            {
               assert(dfsdepth >= 1 && treeedges[dfsdepth - 1] == curredge);

               treecost = shortenSubtree(scip, graph, redcost, treeedges, treecost, basecost, dfsdepth,
                     currhead, nodepos, ancestormark, &rounds);

               dfsdepth--;
            }
            else
            {
               SCIP_Real extendedcost = treecost + redcost[curredge];

               if( vbase[costartnode] == vbase[currhead] )
               {
                  /* also covers the case that leaf is a terminal */
                  const SCIP_Real costartfar = vnoi[costartnode + nnodes].dist;
                  const SCIP_Real currentfar = vnoi[currhead + nnodes].dist;

                  assert(vbase[costartnode + nnodes] >= 0 || costartfar == FARAWAY);
                  assert(vbase[currhead + nnodes] >= 0 || currentfar == FARAWAY);
                  extendedcost += MIN(costartfar + vnoi[currhead].dist, vnoi[costartnode].dist + currentfar);
               }
               else
                  extendedcost += vnoi[costartnode].dist + vnoi[currhead].dist;

               /* can we rule out subtree? */
               if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, cutoff, extendedcost, dfsdepth, curredge, FALSE,
                     ancestormark, eqstack, eqstack_size, eqmark) )
                  continue;

               /* do we need to stop extension? */
               if( truncateSubtree(graph, extendedcost, root, currhead, maxgrad, dfsdepth, maxdfsdepth, &minbound, &stopped) )
               {
                  if( stopped && minbound <= (*treebound) )
                     break;

                  continue;
               }

               nadded_edges = extendSubtreeHead(scip, graph, redcost, curredge, currhead, dfsdepth, nadded_edges, &treecost, treecosts, nodepos,
                     treeedges, edgestack, ancestormark);
               dfsdepth++;
            }
         } /* DFS loop */

         assert(stopped || minbound == FARAWAY);
#ifndef NDEBUG
         assert(SCIPisGE(scip, minbound, orgtreebound));
#endif

         finalizeSubtree(graph, graph->head, treeedges, dfsdepth, stopped, minbound, treebound, ruleout, nodepos, ancestormark);
      } /* main loop for outgoing edges */

#ifndef NDEBUG
      assert(*treebound >= orgtreebound);
#endif

      if( !(*ruleout) && !Is_term(graph->term[innode]) && graph->grad[innode] <= maxgrad )
      {
         /* move down the incoming edge */
         const SCIP_Real treecostoffset = *treebound - pathdist[innode];
         SCIP_Real minbound = FARAWAY;
         SCIP_Real treecost = treecostoffset;
         int dfsdepth = 0;
         int nadded_edges = 0;
         SCIP_Bool stopped = FALSE;
         unsigned int rounds = 0;

         assert(treecost >= 0.0);
         assert(nodepos[innode] == nnodes + 2);
         assert(nodepos[basenode] == nnodes + 1);

         nodepos[graph->head[outedges[0]]] = nnodes + 2;
         nodepos[graph->head[outedges[1]]] = nnodes + 3;
         nodepos[innode] = nnodes + 4;

         /* for basenode */
         basebottlenecks[0] = graph->cost[inedge];

         /* for outedges[0] */
         basebottlenecks[1] = MAX(graph->cost[outedges[0]], graph->cost[inedge]);

         /* for outedges[1] */
         basebottlenecks[2] = MAX(graph->cost[outedges[1]], graph->cost[inedge]);

         /* can we rule out entire subtree already? */
         if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, FARAWAY, 0.0, 0, flipedge(inedge), FALSE,
               NULL, NULL, NULL, NULL) )
         {
            *ruleout = TRUE;
         }
         else
         {
            for( int e = graph->inpbeg[innode]; e != EAT_LAST; e = graph->ieat[e] )
               if( !nodepos[graph->tail[e]] )
                  edgestack[nadded_edges++] = e;
         }

         /* limited DFS starting from inedge */
         while( nadded_edges > 0 )
         {
            const int curredge = edgestack[--nadded_edges];
            const int currtail = graph->tail[curredge];

            /*  subtree already processed? */
            if( nodepos[currtail] )
            {
               assert(dfsdepth >= 1 && treeedges[dfsdepth - 1] == curredge);

               treecost = shortenSubtree(scip, graph, redcost, treeedges, treecost, treecostoffset, dfsdepth,
                     currtail, nodepos, ancestormark, &rounds);

               dfsdepth--;
            }
            else
            {
               const SCIP_Real extendedcost = treecost + redcost[curredge] + pathdist[currtail];

               if( ruleOutSubtree(scip, graph, treecosts, basebottlenecks, nodepos, treeedges, cutoff, extendedcost, dfsdepth, flipedge((unsigned)curredge), FALSE,
                     ancestormark, eqstack, eqstack_size, eqmark) )
                  continue;

               if( truncateSubtree(graph, extendedcost, -1, currtail, maxgrad, dfsdepth, maxdfsdepth, &minbound, &stopped) )
               {
                  if( stopped && minbound <= (*treebound) )
                     break;

                  continue;
               }

               nadded_edges = extendSubtreeTail(graph, redcost, curredge, currtail, dfsdepth, nadded_edges, &treecost, treecosts, nodepos,
                     treeedges, edgestack, ancestormark);
               dfsdepth++;
            }
         } /* DFS loop */

#ifndef NDEBUG
         assert(stopped || minbound == FARAWAY);
         assert(SCIPisGE(scip, minbound, orgtreebound));
#endif

         finalizeSubtree(graph, graph->tail, treeedges, dfsdepth, stopped, minbound, treebound, ruleout, nodepos, ancestormark);
      } /* inedge */

      /* clean nodepos array */
      nodepos[basenode] = 0;
      nodepos[innode] = 0;

      unmarkAncestorsConflict(graph, inedge, ancestormark);

      for( int i = 0; i < 2; i++ )
      {
         nodepos[graph->head[outedges[i]]] = 0;
         unmarkAncestorsConflict(graph, outedges[i], ancestormark);
      }

#ifndef NDEBUG
      assert(*treebound >= orgtreebound);

      for( int i = 0; i < nnodes; i++ )
         assert(nodepos[i] == 0);

      for( int i = 0; i < MAX(graph->edges, graph->orgedges) / 2; i++ )
         assert(ancestormark[i] == 0);
#endif

      SCIPfreeCleanBufferArray(scip, &ancestormark);
      SCIPfreeCleanBufferArray(scip, &nodepos);
   }

   return SCIP_OKAY;
}

/** reduce SPG graph based on reduced cost information and given upper bound todo to be deleted and replaced by reduce_extendedEdge2  */
int reduce_extendedEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const PATH*           vnoi,               /**< Voronoi data structure */
   const SCIP_Real*      cost,               /**< dual ascent costs */
   const SCIP_Real*      pathdist,           /**< distance array from shortest path calculations */
   const int*            result,             /**< sol int array */
   SCIP_Real             minpathcost,        /**< the required reduced path cost to be surpassed */
   int                   root,               /**< the root */
   int*                  nodearr,            /**< for internal stuff */
   STP_Bool*             marked,             /**< edge array to mark which (directed) edge can be removed */
   SCIP_Bool             markdirected        /**< try to also mark edge if anti-parallel is not marked */
)
{
   unsigned int* eqstack;
   SCIP_Bool* eqmark;
   int nfixed = 0;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   const int halfnedges = graph->edges / 2;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(marked != NULL);
   assert(!pcmw || !graph->extended);

   if( SCIPisZero(scip, minpathcost) )
      return 0;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &eqmark, halfnedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eqstack, halfnedges) );

   if( !pcmw )
      for( int k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int erev = e + 1;

      if( pcmw && !SCIPisEQ(scip, graph->cost[e], graph->cost[erev]) )
         continue;

      if( graph->oeat[e] != EAT_FREE )
      {
         int eqstack_size = 0;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(graph->oeat[erev] != EAT_FREE);

         if( SCIPisZero(scip, cost[e]) || SCIPisZero(scip, cost[erev]) )
            continue;

         if( !marked[e] )
         {
            SCIP_CALL_ABORT(reduceCheckEdge(scip, graph, root, cost, pathdist, vnoi, minpathcost, e, allowequality, nodearr,
                  &deletable, eqstack, &eqstack_size, eqmark));

            if( deletable )
               marked[e] = TRUE;
         }

         if( !marked[erev] && (deletable || markdirected) )
         {
            SCIP_Bool erevdeletable = TRUE;
            SCIP_CALL_ABORT(reduceCheckEdge(scip, graph, root, cost, pathdist, vnoi, minpathcost, erev, allowequality,
                  nodearr, &erevdeletable, eqstack, &eqstack_size, eqmark));

            if( erevdeletable )
               marked[erev] = TRUE;

            deletable = (deletable && erevdeletable);
         }

         for( int i = 0; i < eqstack_size; i++ )
            eqmark[eqstack[i]] = FALSE;

         for( int i = 0; i < halfnedges; i++ )
            assert(eqmark[i] == FALSE);

         if( deletable )
         {
            removeEdge(scip, e, graph);

            nfixed++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &eqstack);
   SCIPfreeCleanBufferArray(scip, &eqmark);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      if( graph->grad[k] == 0 && k != root )
         assert(!graph->mark[k]|| (Is_term(graph->term[k]) && graph_pc_isPcMw(graph)));
#endif

   return nfixed;
}
