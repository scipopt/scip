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
 * A list of all interface methods can be found in graph.h.
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
#define STP_EXT_MAXEDGES 500
#define STP_EXT_MAXTREESIZE 20
#define STP_EXT_MAXNLEAVES 20
#define STP_EXT_EDGELIMIT 50000

#define EXEDGE_FREE 0
#define EXEDGE_FIXED 1
#define EXEDGE_KILLED 2

/** reduction data */
typedef struct reduction_data
{
   const SCIP_Real* const reducedcosts;
   const SCIP_Real* const rootdist;
   const PATH* const termpaths;
   const STP_Bool* const edgedeleted;
   int*  const ancestormark;
   const SCIP_Real cutoff;
   const SCIP_Real treeredcostoffset;
   const SCIP_Bool equality;
} REDDATA;

/** distance data */
typedef struct distance_data
{
   SCIP_Bool* const nodepaths_dirty;
   RANGE* const closenodes_range;
   int* const closenodes_index;
   SCIP_Real* const closenodes_dist;
   RANGE* const pathroots_range;   /* of size nedges / 2*/
   int* const pathroots;
} DISTDATA;

/** extension data */
typedef struct extension_data
{
   int* const extstack_data;
   int* const extstack_start;
   int* const extstack_state;
   int* const tree_leaves;
   int* const tree_edges;
   int* const tree_deg;        /**< array of size nnodes */
   SCIP_Real tree_redcost;
   int tree_nleaves;
   int tree_size;
   int tree_depth;
   int extstack_size;
   const int extstack_maxsize;
   const int extstack_maxedges;
   const int tree_maxnleaves;
   const int tree_maxdepth;
   const int tree_maxsize;
   REDDATA* const reddata;
   DISTDATA* const distdata;
} EXTDATA;

/** initializes distance data */
static
SCIP_RETCODE distanceData_initMembers(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = graph->knots;
   int nnodes_real; /* number of not yet deleted nodes */
   int nedges_real;

   assert(distdata && graph && scip);
   assert(distdata->nodepaths_dirty == NULL && distdata->closenodes_dist == NULL && distdata->closenodes_range == NULL);
   assert(distdata->closenodes_index == NULL && distdata->pathroots == NULL && distdata->pathroots_range == NULL);



   return SCIP_OKAY;
}

/** frees members of distance data */
static
void distanceData_freeMembers(
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA*             distdata            /**< to be freed */
)
{
   return;
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



#ifndef NDEBUG
static
SCIP_Bool extTreeIsFlawed(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   int* edgecount;
   int* degreecount;
   const int* const tree_edges = extdata->tree_edges;
   const int* const tree_deg = extdata->tree_deg;
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const int treesize = extdata->tree_size;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   int leavescount;

   SCIP_Bool flawed = FALSE;

   assert(nleaves >= 1);

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &edgecount, nedges) );
   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &degreecount, nnodes) );

   for( int i = 0; i < treesize; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];
      const int tail = graph->tail[e];

      assert(e >= 0 && e < nedges);

      if( edgecount[e] > 0 || tree_deg[tail] <= 0 || tree_deg[head] <= 0 )
         flawed = TRUE;

      degreecount[tail]++;
      degreecount[head]++;

      edgecount[e]++;
   }


   leavescount = 1; /* for tail of initial edge */

   /* degree check */
   for( int i = 0; i < treesize && !flawed; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];

      if( degreecount[head] == 1 )
         leavescount++;

      if( degreecount[head] != extdata->tree_deg[head] )
         flawed = TRUE;
   }

   /* leaves check */
   if( !flawed && leavescount != nleaves )
      flawed = TRUE;

   for( int i = 0; i < nleaves && !flawed; i++ )
   {
      const int leaf = tree_leaves[i];
      if( degreecount[leaf] != 1 )
         flawed = TRUE;
   }

   /* clean-up */
   for( int i = 0; i < treesize; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];
      const int tail = graph->tail[e];

      edgecount[e] = 0;
      degreecount[tail] = 0;
      degreecount[head] = 0;
   }

   SCIPfreeCleanBufferArray(scip, &degreecount);
   SCIPfreeCleanBufferArray(scip, &edgecount);

   return flawed;
}
#endif


/** should we truncate from current component? */
static
void printStack(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
#ifdef SCIP_DEBUG
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extdata->extstack_size - 1;

   for( int j = 0; j <= stackpos; j++ )
   {
      if( extdata->extstack_state[j] == EXT_STATE_NONE )
         printf("pos=%d state=NONE \n", j);
      else if( extdata->extstack_state[j] == EXT_STATE_EXPANDED )
         printf("pos=%d state=EXPANDED \n", j);
      else
      {
         assert(extdata->extstack_state[j] == EXT_STATE_MARKED);

         printf("pos=%d state=MARKED \n", j);
      }

      /* check all leaves of current component */
      for( int i = extstack_start[j]; i < extstack_start[j + 1]; i++ )
      {
         const int edge = extstack_data[i];
         assert(edge >= 0 && edge < graph->edges);

         printf("  ");
         graph_edge_printInfo(graph, edge);
      }
   }
#endif
}


inline static
SCIP_Bool extTreeEdgeAncestorConflict(
   const REDDATA*        reddata,            /**< reduction data */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge                /**< edge to check for conflict */
)
{
   const int* const ancestormark = reddata->ancestormark;
   int count = 0;
   for( const IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      if( ancestormark[((unsigned) curr->index) / 2] )
      {
#ifdef SCIP_DEBUG
         printf("conflict found for ");
         graph_edge_printInfo(graph, edge);
#endif
         return TRUE;
      }
   }
   return FALSE;
}

inline static
SCIP_Bool extLeafIsExtendable(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   int                   leaf                /**< the leaf */
)
{
   assert(graph && isterm);
   assert(leaf >= 0 && leaf < graph->knots);

   // todo if not a terminal, check whether number of neigbhors not contained in current tree < STP_EXT_MAXGRAD

   return (!isterm[leaf] && graph->grad[leaf] <= STP_EXT_MAXGRAD);
}

/** finds position of given leaf in leaves data */
static
int extFindLeafPos(
   const EXTDATA*        extdata,            /**< extension data */
   int                   leaf,               /**< leaf to find */
   int                   startpos            /**< position to start from (going backwards) */
)
{
   int i;
   const int* const tree_leaves = extdata->tree_leaves;

   assert(extdata && tree_leaves);
   assert(startpos > 0);
   assert(leaf >= 0 && extdata->tree_deg[leaf] >= 1);

   for( i = startpos; i >= 0; i-- )
   {
      const int currleaf = tree_leaves[i];

      if( currleaf == leaf )
         break;
   }

   return i;
}

/** adds top component of stack to tree */
static
void extTreeAddStackTop(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   int* const tree_edges = extdata->tree_edges;
   int* const tree_leaves = extdata->tree_leaves;
   int* const tree_deg = extdata->tree_deg;
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real* const redcost = reddata->reducedcosts;
   const int stackpos = extdata->extstack_size - 1;
   const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];

   assert(!(*conflict));
   assert(stackpos >= 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   /* not the initial edge? */
   if( tree_deg[comproot] != graph->knots )
   {
      /* update tree leaves array */

      int comprootpos;

      assert(tree_deg[comproot] == 1);

      /* switch last leaf and root component */
      extdata->tree_nleaves--;
      assert(extdata->tree_nleaves > 0);

      comprootpos = extFindLeafPos(extdata, comproot, extdata->tree_nleaves);
      assert(comprootpos > 0);

      tree_leaves[comprootpos] = tree_leaves[extdata->tree_nleaves];
   }

   /* add top expanded component to tree data */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      assert(extdata->tree_size < extdata->extstack_maxsize);
      assert(edge >= 0 && edge < graph->edges);
      assert(tree_deg[graph->head[edge]] == 0 && tree_deg[graph->tail[edge]] > 0);

      extdata->tree_redcost += redcost[edge];

      tree_edges[(extdata->tree_size)++] = edge;
      tree_leaves[(extdata->tree_nleaves)++] = head;
      tree_deg[head] = 1;

      /* not the initial edge? */
      if( stackpos > 0 )
         tree_deg[graph->tail[edge]]++;

      // todo break and use unmark conflict in extBacktrack
      if( markAncestorsConflict(graph, edge, reddata->ancestormark) )
         *conflict = TRUE;
   }

   extdata->tree_depth++;

   assert(!extTreeIsFlawed(scip, graph, extdata));
}

/** some updates */
static
void extTreeSyncWithStack(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int*                  nupdatestalls,      /**< update stalls counter */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int stackposition = extdata->extstack_size - 1;
   REDDATA* const reddata = extdata->reddata;

   assert(scip && graph && extdata && reddata && nupdatestalls && conflict);
   assert(!(*conflict));

   printStack(graph, extdata);

   /* is current component expanded? */
   if( extdata->extstack_state[stackposition] == EXT_STATE_EXPANDED )
      extTreeAddStackTop(scip, graph, extdata, conflict); /* add component to tree */

   /* recompute reduced costs? */
   if( ++(*nupdatestalls) > EXT_REDCOST_NRECOMP )
   {
      SCIP_Real treecost = reddata->treeredcostoffset;
      const SCIP_Real* const redcost = reddata->reducedcosts;
      const int* const tree_edges = extdata->tree_edges;
      const int tree_size = extdata->tree_size;

      *nupdatestalls = 0;

      assert(!extTreeIsFlawed(scip, graph, extdata));

      for( int i = 0; i < tree_size; i++ )
      {
         const int edge = tree_edges[i];
         assert(edge >= 0 && edge < graph->edges);

         treecost += redcost[edge];
      }

      assert(SCIPisEQ(scip, treecost, extdata->tree_redcost));

      extdata->tree_redcost = treecost;
   }
}

/** should we truncate from current component? */
static
SCIP_Bool extTruncate(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extdata->extstack_size - 1;

   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   if( extdata->tree_depth >= extdata->tree_maxdepth )
   {
      SCIPdebugMessage("truncate (depth) \n");
      return TRUE;
   }

   if( extdata->tree_size >= extdata->tree_maxsize )
   {
      SCIPdebugMessage("truncate (tree size) \n");
      return TRUE;
   }

   if( extdata->tree_nleaves >= extdata->tree_maxnleaves )
   {
      SCIPdebugMessage("truncate (number of leaves) \n");
      return TRUE;
   }

   if( extstack_start[stackpos] >= extdata->extstack_maxedges )
   {
      SCIPdebugMessage("truncate (edges on stack) \n");
      return TRUE;
   }

   /* check whether at least one leaf is extendable */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int leaf = graph->head[edge];

      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[leaf] > 0);

      if( extLeafIsExtendable(graph, isterm, leaf) )
         return FALSE;
   }

   SCIPdebugMessage("truncate (non-promising) \n");
   return TRUE;
}

/** get peripheral reduced cost of current tree including  */
static
SCIP_Real extTreeGetRedcosts(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   extedge             /**< edge for extension or -1 */
)
{
   REDDATA* const reddata = extdata->reddata;
   const int* const tree_leaves = extdata->tree_leaves;
   const PATH* const termpaths = reddata->termpaths;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Real tree_redcost = extdata->tree_redcost; /* includes reduced costs to initial tail */

   assert(graph && reddata && extdata);
   assert(nleaves > 1);

   for( int i = 1; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      assert(extdata->tree_deg[leaf] == 1);

      tree_redcost += termpaths[leaf].dist;
   }

   // todo for the general case we also must consider the case that tail[edge] is the root!
   if( extedge != -1 )
   {
      const int base = graph->tail[extedge];
      const int extvert = graph->head[extedge];
      const SCIP_Real* const redcost = reddata->reducedcosts;

      assert(extedge >= 0 && extedge < graph->edges);
      assert(extdata->tree_deg[base] == 1);

      tree_redcost += redcost[extedge] + termpaths[extvert].dist - termpaths[base].dist;
   }

   return tree_redcost;
}

/** can any extension via edge be ruled out? */
static
SCIP_Bool extRuleOutSimple(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   edge                /**< edge to be tested */
)
{
   REDDATA* const reddata = extdata->reddata;
   const int* const tree_deg = extdata->tree_deg;
   const int extvert = graph->head[edge];

   assert(scip && graph && reddata && extdata);
   assert(edge >= 0 && edge < graph->edges);

   if( tree_deg[extvert] != 0 )
      return TRUE;

   if( reddata->edgedeleted && reddata->edgedeleted[edge] )
      return TRUE;
   else
   {
      const SCIP_Real tree_redcost = extTreeGetRedcosts(graph, extdata, edge);
      const SCIP_Real cutoff = reddata->cutoff;

      if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
         return TRUE;

      if( extTreeEdgeAncestorConflict(reddata, graph, edge) )
         return TRUE;
   }

   return FALSE;
}

/** can subtree be peripherally ruled out? */
static
SCIP_Bool extRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcosts(graph, extdata, -1);
   const SCIP_Real cutoff = reddata->cutoff;

   if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
   {
      SCIPdebugMessage("Rule-out periph (red. cost) \n");
      return TRUE;
   }
   else
   {
      // todo do tree bottleneck test

      // todo do MST test

#ifndef NDEBUG
      const int stackpos = extdata->extstack_size - 1;
      const int* const extstack_data = extdata->extstack_data;
      const int* const extstack_start = extdata->extstack_start;
      const int* const ancestormark = reddata->ancestormark;

      for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
      {
         const int curredge = extstack_data[i];
         int count = 0;
         for( IDX* curr = graph->ancestors[curredge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
         {
            assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
            assert(ancestormark[((unsigned) curr->index) / 2] == 1);
         }
      }
#endif
   }

   return FALSE;
}


/** top component is rebuilt, and
 *  if success == TRUE: goes back to first marked component
 *  if success == FALSE: goes back to first marked or non-expanded component
 *   */
static
void extBacktrack(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             success,            /**< backtrack from success? */
   SCIP_Bool             ancestor_conflict,  /**< backtrack triggered by ancestor conflict? */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int* const tree_deg = extdata->tree_deg;
   int stackpos = extdata->extstack_size - 1;
   const int stackstart = extstack_start[stackpos];

   assert(graph && reddata && extdata);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);

   /* top component already expanded? */
   if( extstack_state[stackpos] != EXT_STATE_NONE )
   {
      const SCIP_Real* const redcost = reddata->reducedcosts;
      int* const tree_leaves = extdata->tree_leaves;
      const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];
      const int compsize = extstack_start[stackpos + 1] - extstack_start[stackpos];

      assert(compsize > 0);
      assert(tree_deg[comproot] > 1);

      /* remove top component */
      for( int i = stackstart; i < extstack_start[stackpos + 1]; i++ )
      {
         const int edge = extstack_data[i];
         const int head = graph->head[edge];
         const int tail = graph->tail[edge];

         assert(edge >= 0 && edge < graph->edges);
         assert(tree_deg[head] == 1 && tree_deg[tail] > 1);

         (extdata->tree_redcost) -= redcost[edge];
         tree_deg[head] = 0;
         tree_deg[tail]--;

         if( ancestor_conflict )
            unmarkAncestorsConflict(graph, edge, reddata->ancestormark);
         else
            unmarkAncestors(graph, edge, reddata->ancestormark);
      }

      (extdata->tree_size) -= compsize;
#if 0
      for( int k = extdata->tree_nleaves - 1; k >= extdata->tree_nleaves - compsize; k-- )
         printf("bt remove leaf %d \n", tree_leaves[k]);
#endif
      (extdata->tree_nleaves) -= compsize;
      (extdata->tree_depth)--;

      /* add component root to leaves array and remove current one */
      assert(tree_deg[comproot] == 1);
      tree_leaves[extdata->tree_nleaves++] = comproot;

      assert(extdata->tree_size >= 0 && extdata->tree_depth >= 0);
   }

   stackpos--;

   /* backtrack */
   if( success )
   {
      while( extstack_state[stackpos] == EXT_STATE_NONE )
      {
         stackpos--;
         assert(stackpos >= 0);
      }

      SCIPdebugMessage("backtrack SUCCESS \n");
      assert(extstack_state[stackpos] == EXT_STATE_EXPANDED || extstack_state[stackpos] == EXT_STATE_MARKED);
   }
   else
   {
      while( extstack_state[stackpos] == EXT_STATE_EXPANDED )
      {
         stackpos--;
         assert(stackpos >= 0);
      }

      SCIPdebugMessage("backtrack FAILURE \n");
      assert(extstack_state[stackpos] == EXT_STATE_NONE || extstack_state[stackpos] == EXT_STATE_MARKED);
   }

   extdata->extstack_size = stackpos + 1;

   assert(!extTreeIsFlawed(scip, graph, extdata));
}

/** expands top component of stack (backtracks if stack is full) */
static
void extStackExpand(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extdata->extstack_size - 1;
   int datasize = extstack_start[stackpos];
   const int setsize = extstack_start[stackpos + 1] - extstack_start[stackpos];
   const uint32_t powsize = pow(2, setsize);

   assert(extdata && scip && graph && isterm && success);
   assert(setsize <= STP_EXT_MAXGRAD);
   assert(setsize > 0 && setsize <= 32);
   assert(stackpos >= 1);
   assert(extstack_state[stackpos] == EXT_STATE_NONE);

   /* stack too full? */
   if( (datasize + (int) pow(2, setsize)) > extdata->extstack_maxsize )
   {
      *success = FALSE;
      extBacktrack(scip, graph, *success, FALSE, extdata);

      return;
   }

   /* collect edges for new component and find conflicts */
   for( int i = extstack_start[stackpos], j = 0; i < extstack_start[stackpos + 1]; i++, j++ )
   {
      const int edge = extstack_data[i];
      assert(j < STP_EXT_MAXGRAD);
      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[graph->head[edge]] == 0);

      // todo find excluding pairs, use ancestormark and bottleneck
      extedges[j] = edge;
   }

   /* compute and add components (overwrite previous, non-expanded component) */
   for( uint32_t counter = 1; counter < powsize; counter++ )
   {
      for( int j = 0; j < setsize; j++ )
      {
         /* Check if jth bit in counter is set */
         if( counter & (1 << j) )
         {
            extstack_data[datasize++] = extedges[j];
            SCIPdebugMessage("  %d \n", graph->head[extedges[j]]);
         }
      }

      SCIPdebugMessage("... added \n");
      assert(stackpos < extdata->extstack_maxsize - 1);

      extstack_state[stackpos] = EXT_STATE_EXPANDED;
      extstack_start[++stackpos] = datasize;

      assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
   }

   assert(stackpos >= extdata->extstack_size);

   extdata->extstack_size = stackpos;
}

/** extend top component of stack (backtracks if stack is full) */
static
void extExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD * STP_EXT_MAXGRAD];
   int extedgesstart[STP_EXT_MAXGRAD + 1];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extdata->extstack_size - 1;
   int nfullextensions;
   int nsingleextensions;

#ifndef NDEBUG
   assert(stackpos >= 0);
   assert(extstack_state[stackpos] == EXT_STATE_EXPANDED );
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] <= STP_EXT_MAXGRAD);

   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedgesstart[i] = -1;

   for( int i = 0; i < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;
#endif

   extstack_state[stackpos] = EXT_STATE_MARKED;

   nfullextensions = 0;
   nsingleextensions = 0;
   extedgesstart[0] = 0;

   /* loop over all leaves of extension */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int leaf = graph->head[extstack_data[i]];

      assert(extstack_data[i] >= 0 && extstack_data[i] < graph->edges);

      /* extensions from leaf not possible? */
      if( !extLeafIsExtendable(graph, isterm, leaf) )
         continue;

      /* assemble feasible single edge extensions from leaf */
      for( int e = graph->outbeg[leaf]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( !extRuleOutSimple(scip, graph, extdata, e) )
         {
            assert(nsingleextensions < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD);
            extedges[nsingleextensions++] = e;
         }
#ifdef SCIP_DEBUG
         else
         {
            printf("simple rule out: ");
            graph_edge_printInfo(graph, e);
         }
#endif
      }

      extedgesstart[++nfullextensions] = nsingleextensions;
   }

   assert(nfullextensions <= STP_EXT_MAXGRAD);

   /* found no valid extensions? */
   if( nfullextensions == 0 )
   {
      *success = FALSE;
   }
   /* found valid extensions, but all ruled out already? */
   else if( nsingleextensions == 0 )
   {
      *success = TRUE;
   }
   /* found non-empty valid extensions */
   else
   {
      int datasize = extstack_start[stackpos + 1];
      int extsize[STP_EXT_MAXGRAD];
      int extindex[STP_EXT_MAXGRAD];

      assert(extedgesstart[nfullextensions] - extedgesstart[0] > 0);

      /* stack too small? */
      if( datasize + (extedgesstart[nfullextensions] - extedgesstart[0]) > extdata->extstack_maxsize )
      {
         *success = FALSE;
         extBacktrack(scip, graph, *success, FALSE, extdata);

         return;
      }

      for( int i = 0; i < nfullextensions; i++ )
      {
         assert(extedgesstart[i + 1] >= 0);

         extsize[i] = extedgesstart[i + 1] - extedgesstart[i];
         extindex[i] = i;
         assert(extsize[i] >= 0);
      }

      SCIPsortDownIntInt(extsize, extindex, nfullextensions);

      /* put the non-empty extensions on the stack, with smallest last */
      for( int i = 0; i < nfullextensions; i++ )
      {
         const int index = extindex[i];

         if( extsize[i] == 0 )
         {
            assert(i > 0);
            assert(extedgesstart[index + 1] - extedgesstart[index] == 0);

            for( int j = i; j < nfullextensions; j++ )
               assert(extsize[j] == 0);

            break;
         }

         for( int j = extedgesstart[index]; j < extedgesstart[index + 1]; j++ )
         {
            assert(extedges[j] >= 0);
            extstack_data[datasize++] = extedges[j];
         }

         assert(stackpos < extdata->extstack_maxsize - 2);

         extstack_state[++stackpos] = EXT_STATE_NONE;
         extstack_start[stackpos + 1] = datasize;
      }

#ifdef SCIP_DEBUG
      printf("added extending edges:  \n");

      for( int i = extstack_start[extdata->extstack_size]; i < extstack_start[stackpos + 1]; i++ )
         graph_edge_printInfo(graph, extstack_data[i]);
#endif

      extdata->extstack_size = stackpos + 1;

      *success = TRUE;

      /* try to expand last (smallest) component */
      extStackExpand(scip, graph, isterm, extdata, success);
   }
}


/** check (directed) arc */
SCIP_RETCODE reduce_extendedCheckArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   root,               /**< graph root from dual ascent */
   const SCIP_Real*      redcost,            /**< reduced costs */
   const SCIP_Real*      rootdist,           /**< shortest path distances  */
   const PATH*           termpaths,          /**< paths to nearest terminals  */
   const STP_Bool*       edgedeleted,        /**< edge array to mark which directed edge can be removed or NULL */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value */
   int                   edge,               /**< directed edge to be checked */
   SCIP_Bool             equality,           /**< allow equality? */
   int*                  tree_deg,           /**< -1 for forbidden nodes (e.g. PC terminals), 0 otherwise; in method: position ( > 0) for nodes in tree */
   SCIP_Bool*            deletable           /**< is edge deletable? */
)
{
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   SCIP_Real edgebound = redcost[edge] + rootdist[tail] + termpaths[head].dist;

#ifndef NDEBUG
   assert(scip && graph && redcost && rootdist && termpaths && deletable && tree_deg);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);

   if( !graph_pc_isPcMw(graph) )
      for( int k = 0; k < graph->knots; k++ )
         assert(graph->mark[k] == (graph->grad[k] > 0));
#endif

   /* trivial rule-out? */
   if( SCIPisGT(scip, edgebound, cutoff) || (equality && SCIPisEQ(scip, edgebound, cutoff)) || head == root )
   {
      *deletable = TRUE;
      return SCIP_OKAY;
   }

   *deletable = FALSE;

   /* can we extend from 'edge'? */
   if( extLeafIsExtendable(graph, isterm, tail) || extLeafIsExtendable(graph, isterm, head) )
   {
      int* extstack_data;
      int* extstack_start;
      int* extstack_state;
      int* tree_edges;
      int* tree_leaves;
      int* ancestormark;
      const int nnodes = graph->knots;
      const int maxdfsdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;
      const int maxstackedges = MIN(nnodes / 2, STP_EXT_MAXEDGES);

      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_data, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_state, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_edges, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_leaves, nnodes) );
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &ancestormark, (MAX(graph->edges, graph->orgedges) / 2)) );

      tree_deg[root] = -1;

      /* can we extend from head? */
      if( extLeafIsExtendable(graph, isterm, head) )
      {
         const SCIP_Real treeredcostoffset = rootdist[tail];
         int nupdatestalls = 0;
         SCIP_Bool success = TRUE;
         SCIP_Bool conflict = FALSE;
         DISTDATA distdata = {NULL, NULL, NULL, NULL, NULL, NULL};
         REDDATA reddata = {redcost, rootdist, termpaths, edgedeleted, ancestormark, cutoff, treeredcostoffset, equality};
         EXTDATA extdata = {extstack_data, extstack_start, extstack_state, tree_leaves, tree_edges,
            tree_deg, 0.0, 0, 0, 0, 0, nnodes - 1, maxstackedges, STP_EXT_MAXNLEAVES, maxdfsdepth,
            STP_EXT_MAXTREESIZE, &reddata, &distdata};

         extdata.tree_redcost = treeredcostoffset;
         extdata.tree_depth = 0;
         tree_deg[tail] = nnodes;

         /* put 'edge' on the stack */
         extdata.extstack_size = 1;
         extstack_start[0] = 0;
         extstack_start[1] = 1;
         extstack_data[0] = edge;
         extstack_state[0] = EXT_STATE_EXPANDED;
         tree_leaves[0] = tail;
         extdata.tree_nleaves = 1;

         extTreeSyncWithStack(scip, graph, &extdata, &nupdatestalls, &conflict);

         assert(!conflict);

         extExtend(scip, graph, isterm, &extdata, &success);

         assert(extstack_state[0] == EXT_STATE_MARKED);
         assert(success || 1 == extdata.extstack_size);

         /* limited DFS backtracking; stops once back at 'edge' */
         while( extdata.extstack_size > 1 )
         {
            const int stackposition = extdata.extstack_size - 1;
            conflict = FALSE;

            extTreeSyncWithStack(scip, graph, &extdata, &nupdatestalls, &conflict);

            /* has current component already been extended? */
            if( extstack_state[stackposition] == EXT_STATE_MARKED )
            {
               extBacktrack(scip, graph, success, FALSE, &extdata);
               continue;
            }

            /* component not expanded yet? */
            if( extstack_state[stackposition] != EXT_STATE_EXPANDED )
            {
               assert(extstack_state[stackposition] == EXT_STATE_NONE);

               extStackExpand(scip, graph, isterm, &extdata, &success);
               continue;
            }

            assert(extstack_state[stackposition] == EXT_STATE_EXPANDED);

            if( conflict || extRuleOutPeriph(scip, graph, &extdata) )
            {
               success = TRUE;
               extBacktrack(scip, graph, success, conflict, &extdata);
               continue;
            }

            if( extTruncate(graph, isterm, &extdata) )
            {
               success = FALSE;
               extBacktrack(scip, graph, success, FALSE, &extdata);
               continue;
            }

            /* neither ruled out nor truncated, so extend */
            extExtend(scip, graph, isterm, &extdata, &success);

         } /* DFS loop */

         *deletable = success;
         assert(tree_deg[head] == 1 && tree_deg[tail] == nnodes);

         tree_deg[head] = 0;
         tree_deg[tail] = 0;

         unmarkAncestors(graph, edge, ancestormark);
      } /* extend from head */


      /* finalize arrays */

      tree_deg[root] = 0;

#ifndef NDEBUG
      for( int i = 0; i < MAX(graph->edges, graph->orgedges) / 2; i++ )
         assert(ancestormark[i] == 0);

      for( int i = 0; i < nnodes; i++ )
         assert(tree_deg[i] == 0 || tree_deg[i] == -1);
#endif
      SCIPfreeCleanBufferArray(scip, &ancestormark);
      SCIPfreeBufferArray(scip, &tree_leaves);
      SCIPfreeBufferArray(scip, &tree_edges);
      SCIPfreeBufferArray(scip, &extstack_state);
      SCIPfreeBufferArray(scip, &extstack_start);
      SCIPfreeBufferArray(scip, &extstack_data);
   }

   return SCIP_OKAY;
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
            assert(!graph_pc_isPcMw(g) || (!Is_pterm(g->term[g->head[e]]) && !Is_pterm(g->term[g->tail[e]])));
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


/** extended reduction test for edges */
SCIP_RETCODE reduce_extendedEdge2(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const PATH*           termpaths,          /**< paths to nearest terminals */
   const SCIP_Real*      redcost,            /**< dual ascent costs */
   const SCIP_Real*      rootdist,           /**< shortest path distance (w.r.t reduced costs) from root to any node */
   const int*            result,             /**< solution array */
   SCIP_Real             minpathcost,        /**< the required reduced path cost to be surpassed */
   int                   root,               /**< the root */
   SCIP_Bool             markdirected,       /**< try to also mark edge if anti-parallel is not marked */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   SCIP_Bool* isterm;
   int* tree_deg;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(scip && graph && redcost);
   assert(!pcmw || !graph->extended);
   assert(root >= 0 && root < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, minpathcost) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_deg, nnodes) );

   graph_get_isTerm(graph, isterm);

   if( !pcmw )
      for( int k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph->mark[k] )
         tree_deg[k] = 0;
      else
         tree_deg[k] = -1;
   }

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int tail = graph->tail[e];
      const int head = graph->head[e];

      if( pcmw && (!graph->mark[tail] || !graph->mark[head]) )
         continue;

      if( graph->oeat[e] != EAT_FREE )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev);
         assert(SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) || SCIPisZero(scip, redcost[erev]) )
            continue;

         if( !edgedeletable[e] )
         {
            SCIP_CALL( reduce_extendedCheckArc(scip, graph, root, redcost, rootdist, termpaths, edgedeletable,
                  isterm, minpathcost, e, allowequality, tree_deg, &deletable) );

            if( deletable )
               edgedeletable[e] = TRUE;
         }

         if( !edgedeletable[erev] && (deletable || markdirected) )
         {
            SCIP_Bool erevdeletable = TRUE;

            SCIP_CALL( reduce_extendedCheckArc(scip, graph, root, redcost, rootdist, termpaths, edgedeletable,
                  isterm, minpathcost, erev, allowequality, tree_deg, &erevdeletable) );

            if( erevdeletable )
               edgedeletable[erev] = TRUE;

            deletable = (deletable && erevdeletable);
         }

         if( deletable )
         {
            graph_edge_del(scip, graph, e, TRUE);

            // todo also delete CSR

            if( graph->grad[tail] == 0 )
               graph->mark[tail] = FALSE;

            if( graph->grad[head] == 0 )
               graph->mark[head] = FALSE;

            (*nelims)++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &tree_deg);
   SCIPfreeBufferArray(scip, &isterm);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      if( graph->grad[k] == 0 && k != root )
         assert(!graph->mark[k]);
#endif

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
            const int tail = graph->tail[e];
            const int head = graph->head[e];
            graph_edge_del(scip, graph, e, TRUE);

            if( graph->grad[tail] == 0 )
               graph->mark[tail] = FALSE;

            if( graph->grad[head] == 0 )
               graph->mark[head] = FALSE;

            nfixed++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &eqstack);
   SCIPfreeCleanBufferArray(scip, &eqmark);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      if( graph->grad[k] == 0 && k != root )
         assert(!graph->mark[k]);
#endif

   return nfixed;
}
