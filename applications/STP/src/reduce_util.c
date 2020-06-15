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

/**@file   reduce_util.c
 * @brief  utility methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem reduction techniques.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "reduce.h"
#include "portab.h"


/** storage for edge on complete graph */
typedef struct complete_edge
{
   int                   tail;              /**< tail vertex */
   int                   head;              /**< head vertex */
   SCIP_Real             cost;              /**< edge cost */
} CEDGE;


/** lightweight minimum spanning tree structure that allows to add vertices to given MST on complete graph (in CSR format) */
struct dynamic_complete_minimum_spanning_tree
{
   CEDGE*                edgestore;         /**< storage for edges (of size maxnnodes) */
   SCIP_Real*            adjcost_buffer;    /**< distances buffer (of size maxnnodes) */
   SCIP_Bool*            nodemark;          /**< array for marking nodes (of size maxnnodes) */
   int                   maxnnodes;         /**< maximum number of nodes that can be handled */
};


/** see reduce.h */
struct node_one_hop_star
{
   int*                  edgeId;            /**< IDs for each adjacent edge of current node (of size maxNodeDegree) */
   int*                  edgesSelected;     /**< list of currently selected edges (of size maxNodeDegree) */
   int*                  edgesSelectedPos;  /**< list of position of currently selected edges w.r.t. edgeId (of size maxNodeDegree) */
   int                   nodeDegree;        /**< degree of current node */
   int                   starDegree;        /**< degree of current star */
   int                   maxNodeDegree;     /**< maximum allowed node degree */
   int                   starcenter;        /**< node for which the star is created */
   SCIP_Bool             allStarsChecked;   /**< have all stars been checked? */
};


/** bottleneck link-cut tree node */
typedef struct bottleneck_link_cut_node
{
   SCIP_Real             dist_edgetail;      /**< distance for edge tail */
   SCIP_Real             dist_edgehead;      /**< distance for edge tail */
   SCIP_Real             edgecost;           /**< edge cost */
   SCIP_Real             edgebottleneck;     /**< restricted bottleneck of edge */
   int                   head;               /**< head of node in tree */
   int                   edge;               /**< edge of node (w.r.t original graph used for initialization) */
} BLCNODE;


/** see reduce.h */
struct bottleneck_link_cut_tree
{
   BLCNODE*              mst;               /**< the actual tree, represented by its nodes */
   int*                  nodes_curr2org;    /**< map */
   int*                  nodes_org2curr;    /**< map */
   SCIP_Bool*            mstedges_isLink;   /**< is the edge on a fundamental path between terminals? */
   int                   root;              /**< root of the tree */
   int                   nnodes_org;        /**< original number of graph nodes (before reductions) */
   int                   nnodes_curr;       /**< current number of graph nodes */
};


/*
 * local methods
 */


/** builds mapping */
static inline
SCIP_RETCODE blctreeBuildNodeMap(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< the graph */
   BLCTREE* RESTRICT     blctree             /**< the tree */
)
{
   const int* nodes_isMarked = graph->mark;
   int* RESTRICT nodes_curr2org;
   int* RESTRICT nodes_org2curr;
   int nodecount_curr;
   const int nnodes_curr = blctree->nnodes_curr;
   const int nnodes_org = blctree->nnodes_org;

   assert(0 < nnodes_curr && nnodes_curr <= nnodes_org);
   assert(nnodes_org == graph->knots);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(blctree->nodes_curr2org), nnodes_curr) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(blctree->nodes_org2curr), nnodes_org) );

   nodes_curr2org = blctree->nodes_curr2org;
   nodes_org2curr = blctree->nodes_org2curr;
   nodecount_curr = 0;

   for( int i = 0; i < nnodes_org; ++i )
   {
      if( nodes_isMarked[i] )
      {
         assert(nodecount_curr < nnodes_curr);
         nodes_org2curr[i] = nodecount_curr;
         nodes_curr2org[nodecount_curr++] = i;
      }
      else
      {
         assert(graph_pc_isPc(graph) || graph->grad[i] == 0);
         nodes_org2curr[i] = UNKNOWN;
      }
   }

   assert(nodecount_curr <= nnodes_curr);
   assert(graph_pc_isPc(graph) || nodecount_curr == nnodes_curr);

   if( graph_pc_isPc(graph) )
   {
      blctree->nnodes_curr = nodecount_curr;
   }

   return SCIP_OKAY;
}


/** allocates BLC tree memory and builds mapping */
static
SCIP_RETCODE blctreeInitPrimitives(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< the graph */
   BLCTREE**             blctree             /**< to be initialized */
)
{
   BLCTREE *tree;
   int nnodes_curr;
   const int nnodes = graph_get_nNodes(graph);

   graph_get_nVET(graph, &(nnodes_curr), NULL, NULL);

   SCIP_CALL( SCIPallocMemory(scip, blctree) );
   tree = *blctree;

   tree->root = UNKNOWN;
   tree->nnodes_org = nnodes;
   tree->nnodes_curr = nnodes_curr;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tree->mst), nnodes_curr) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tree->mstedges_isLink), nnodes_curr) );
   SCIP_CALL( blctreeBuildNodeMap(scip, graph, tree) );

   return SCIP_OKAY;
}


/** resets given node */
static inline
void blctreeResetNode(
   int                   i,               /**< the node */
   BLCTREE* RESTRICT     blctree          /**< the tree */
)
{
   BLCNODE* RESTRICT tree = blctree->mst;
   assert(0 <= i && i < blctree->nnodes_curr);

   tree[i].dist_edgetail = 0.0;
   tree[i].dist_edgehead = 0.0;
   tree[i].edgebottleneck = FARAWAY;
   tree[i].edgecost = -FARAWAY;
   tree[i].edge = UNKNOWN;
   tree[i].head = UNKNOWN;
}


/** sets new root*/
static inline
void blctreeEvert(
   const GRAPH*          graph,              /**< the graph */
   int                   newroot,            /**< new root */
   BLCTREE*              blctree             /**< tree */
)
{
   BLCNODE blcnode_org;
   BLCNODE* RESTRICT tree = blctree->mst;
   int node;
   const int oldroot = blctree->root;

   assert(tree);
   assert(newroot != oldroot);
   assert(0 <= newroot && newroot < blctree->nnodes_org);
   assert(0 <= tree[newroot].head && newroot < blctree->nnodes_org);

   node = newroot;
   blcnode_org = tree[newroot];
   while( node != oldroot )
   {
      const SCIP_Real dist_tail = blcnode_org.dist_edgetail;
      const SCIP_Real dist_head = blcnode_org.dist_edgehead;
      const SCIP_Real edgecost = blcnode_org.edgecost;
      const SCIP_Real edgebottleneck = blcnode_org.edgebottleneck;
      const int head = blcnode_org.head;
      const int edge = blcnode_org.edge;

      assert(edge != UNKNOWN);
      assert(head != UNKNOWN);
      assert(graph->tail[edge] == blctree->nodes_curr2org[node]);

      blcnode_org = tree[head];

      /* swap! */
      tree[head].dist_edgetail = dist_head;
      tree[head].dist_edgehead = dist_tail;
      tree[head].edgecost = edgecost;
      tree[head].edgebottleneck = edgebottleneck;
      tree[head].head = node;
      tree[head].edge = flipedge_Uint(edge);

      node = head;
   }

   blctreeResetNode(newroot, blctree);
   blctree->root = newroot;
}


/** gets cost of path */
static inline
SCIP_Real blctreeGetRootPathCost(
   int                   startnode,          /**< node to start from */
   const BLCTREE*        blctree             /**< tree */
)
{
   const BLCNODE* const tree = blctree->mst;
   SCIP_Real pathcost = 0.0;
   int node = startnode;
   const int root = blctree->root;

   assert(0 <= startnode && startnode < blctree->nnodes_org);

   while( node != root )
   {
      const int head = tree[node].head;
      assert(head != UNKNOWN);

      SCIPdebugMessage("%d  ", node);

      pathcost += tree[node].edgecost;
      node = head;
   }

   return pathcost;
}


/** updates path to root */
static inline
void blctreeUpdateRootPath(
   int                   startnode,          /**< node to start from */
   SCIP_Real             bottleneck,         /**< bottleneck */
   BLCTREE*              blctree             /**< tree */
)
{
   BLCNODE* RESTRICT tree = blctree->mst;
   SCIP_Real nodedist_tail = 0.0;
   const SCIP_Real pathcost = blctreeGetRootPathCost(startnode, blctree);
   int node = startnode;
   const int root = blctree->root;

   SCIPdebugMessage(" ....  blctreeUpdateRootPath %d->%d  pathcost=%f \n", startnode, root, pathcost);

   assert(0 <= startnode && startnode < blctree->nnodes_org);

   while( node != root )
   {
      SCIP_Real nodedist_head;
      const int head = tree[node].head;
      assert(head != UNKNOWN);

      /* take the MIN */
      if( tree[node].edgebottleneck > bottleneck )
      {
         tree[node].edgebottleneck = bottleneck;
      }

      /* take the MAX */
      if( tree[node].dist_edgetail < nodedist_tail )
      {
         tree[node].dist_edgetail = nodedist_tail;
      }

      nodedist_tail += tree[node].edgecost;

      nodedist_head = pathcost - nodedist_tail;
      assert(GE(nodedist_head, 0.0));

      /* take the MAX */
      if( tree[node].dist_edgehead < nodedist_head )
      {
         tree[node].dist_edgehead = nodedist_head;
      }

      SCIPdebugMessage("node=%d dist_tail=%f  dist_head=%f edge %d->%d: \n", node, tree[node].dist_edgetail,
            tree[node].dist_edgehead, node, head);


      node = head;
   }
}


/** computes bottlenecks */
static
void blctreeComputeBottlenecks(
   const GRAPH*          graph,              /**< the graph */
   BLCTREE*              blctree             /**< to be built */
)
{
   BLCNODE* RESTRICT mst = blctree->mst;
   const int* nodes_curr2org = blctree->nodes_curr2org;
   const int* nodes_org2curr = blctree->nodes_org2curr;
   const int* nodes_isMarked = graph->mark;
   const int nnodes_curr = blctree->nnodes_curr;

   assert(mst);
   assert(0 <= blctree->root && blctree->root < nnodes_curr);
   assert(graph_isMarked(graph));

   for( int node = 0; node < nnodes_curr; node++ )
   {
      const int node_org = nodes_curr2org[node];

      assert(graph_knot_isInRange(graph, node_org));
      assert(graph->grad[node_org] > 0);
      assert(!graph_pc_isPc(graph) || !graph_pc_knotIsDummyTerm(graph, node_org));

      if( blctree->root != node )
      {
         /* make i the new root */
         blctreeEvert(graph, node, blctree);
      }

      for( int e = graph->outbeg[node_org]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int head_org = graph->head[e];

         if( nodes_isMarked[head_org] )
         {
            const int head = nodes_org2curr[head_org];
            assert(0 <= head && head < nnodes_curr);

            /* is the edge in the tree or already checked? */
            if( mst[head].head == node || head < node )
               continue;

            assert(mst[head].edge != e || mst[head].edge != flipedge(e));

            SCIPdebugMessage("---CHECK EDGE: %d->%d  root=%d \n", node_org, head_org, node);

            blctreeUpdateRootPath(head, graph->cost[e], blctree);
         }
      }
   }

   /* NOTE necessary for correct computation of edge state! */
   if( blctree->root != nodes_org2curr[graph->source] )
   {
      blctreeEvert(graph, nodes_org2curr[graph->source], blctree);
   }
}



/** sets stage of MST edges */
static
void blctreeComputeEdgesState(
   const GRAPH*          graph,              /**< the graph */
   BLCTREE*              blctree             /**< to be built */
)
{
   BLCNODE* RESTRICT mst = blctree->mst;
   SCIP_Bool* RESTRICT isTermLink = blctree->mstedges_isLink;
   const int* nodes_curr2org = blctree->nodes_curr2org;
   const int nnodes_curr = blctree->nnodes_curr;
   const int blcroot = blctree->root;

   assert(mst);
   assert(0 <= blctree->root && blctree->root < nnodes_curr);
   assert(graph_isMarked(graph));
   assert(blctree->root == blctree->nodes_org2curr[graph->source]);

   for( int i = 0; i < nnodes_curr; i++ )
      isTermLink[i] = FALSE;

   /* mark all edges from terminals to root */
   for( int i = 0; i < nnodes_curr; i++ )
   {
      const int node_org = nodes_curr2org[i];
      assert(graph_knot_isInRange(graph, node_org));
      assert(!graph_pc_isPc(graph) || !graph_pc_knotIsDummyTerm(graph, node_org));

      if( Is_term(graph->term[node_org]) )
      {
         assert(graph->grad[node_org] > 0);

         for( int node = i; node != blcroot; node = mst[node].head )
         {
            assert(0 <= node && node < nnodes_curr);

            isTermLink[node] = TRUE;
         }
      }
   }

   assert(!isTermLink[blcroot]);

}

/** builds MST and sets BLC tree accordingly */
static
SCIP_RETCODE blctreeBuildMst(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                graph,              /**< the graph */
   BLCTREE*              blctree             /**< to be built */
)
{
   PATH* pmst;
   BLCNODE* RESTRICT tree = blctree->mst;
   const int* nodes_curr2org = blctree->nodes_curr2org;
   const int* nodes_org2curr = blctree->nodes_org2curr;
   const int nnodes_org = graph_get_nNodes(graph);
   const int nnodes_curr = blctree->nnodes_curr;
   blctree->root = nodes_org2curr[graph->source];

   assert(tree);
   assert(nnodes_org == blctree->nnodes_org);
   assert(nnodes_curr > 0);
   assert(graph->grad[graph->source] > 0);
   assert(!graph_pc_isMw(graph));
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(!graph_pc_isPcMw(graph) || graph_pc_isRootedPcMw(graph));

   SCIP_CALL( SCIPallocBufferArray(scip, &pmst, nnodes_org) );
   graph_path_exec(scip, graph, MST_MODE, graph->source, graph->cost, pmst);

   /* fill the tree */
   for( int k_curr = 0; k_curr < nnodes_curr; k_curr++ )
   {
      const int k_org = nodes_curr2org[k_curr];
      const int mstedge = pmst[k_org].edge;

      assert(graph_knot_isInRange(graph, k_org));

      if( mstedge >= 0 )
      {
         const int revedge = flipedge(mstedge);
         const int head_curr = nodes_org2curr[graph->head[revedge]];

         SCIPdebugMessage("MST edge %d->%d \n", k_org, graph->head[revedge]);

         tree[k_curr].edge = revedge;
         tree[k_curr].edgecost = graph->cost[revedge];
         tree[k_curr].head = head_curr;
         tree[k_curr].dist_edgetail = 0.0;
         tree[k_curr].dist_edgehead = 0.0;
         tree[k_curr].edgebottleneck = FARAWAY;

         assert(k_org != graph->source);
         assert(k_curr != blctree->root);
         assert(0 <= head_curr && head_curr < nnodes_curr);
         assert(graph->tail[tree[k_curr].edge] == k_org);
      }
#ifndef NDEBUG
      else
      {
         assert(k_curr == blctree->root);
         assert(k_org == graph->source);
      }
#endif
   }

   blctreeResetNode(blctree->root, blctree);

   SCIPfreeBufferArray(scip, &pmst);

   return SCIP_OKAY;
}


/** builds BLC tree */
static
SCIP_RETCODE blctreeInitBottlenecks(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                graph,              /**< the graph */
   BLCTREE*              blctree             /**< to be built */
)
{
   SCIP_CALL( blctreeBuildMst(scip, graph, blctree) );
   blctreeComputeBottlenecks(graph, blctree);
   blctreeComputeEdgesState(graph, blctree);

   return SCIP_OKAY;
}


/** recursive method for adding node to MST */
static
void dcmstInsert(
   const CSR*            org_mst,            /**< the base MST */
   const SCIP_Real       adjcosts[],         /**< (undirected) adjacency costs for new node */
   int                   root,               /**< the current root */
   CEDGE                 new_mst[],          /**< new MST */
   SCIP_Bool             new_nodemarked[],   /**< array */
   CEDGE*                max_path_edge,      /**< pointer to maximum edge on path to new node */
   int*                  new_nedges          /**< pointer to current number of edges */
)
{
   CEDGE root2new = { .tail = root, .head = org_mst->nnodes, .cost = adjcosts[root] };
   const int* const org_start = org_mst->start;
   const int* const org_head = org_mst->head;
   const SCIP_Real* const org_cost = org_mst->cost;

   assert(new_nodemarked[root]);

   /* visit all neighbors or root in the original MST */
   for( int i = org_start[root]; i != org_start[root + 1]; ++i )
   {
      const int w = org_head[i];

      /* node not visited yet? */
      if( !new_nodemarked[w] )
      {
         const SCIP_Real costroot2w = org_cost[i];

         new_nodemarked[w] = TRUE;
         dcmstInsert(org_mst, adjcosts, w, new_mst, new_nodemarked, max_path_edge, new_nedges);

         assert(max_path_edge->tail >= 0);
         assert(*new_nedges >= 0 && *new_nedges < org_mst->nnodes);

         if( max_path_edge->cost < costroot2w )
         {
            new_mst[(*new_nedges)++] = *max_path_edge;

            if( costroot2w < root2new.cost )
            {
               root2new.tail = root;
               root2new.head = w;
               root2new.cost = costroot2w;
            }
         }
         else
         {
            const int nedges = (*new_nedges);

            new_mst[nedges].tail = root;
            new_mst[nedges].head = w;
            new_mst[nedges].cost = costroot2w;

            (*new_nedges)++;

            if( max_path_edge->cost < root2new.cost )
            {
               root2new = *max_path_edge;
            }
         }
      }
   }

   *max_path_edge = root2new;
}


/** add node to MST */
static inline
void dcmstAddNode(
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst                /**< underlying structure */
)
{
   CEDGE max_path_edge = { .tail = -1, .head = -1, .cost = -1.0 };
   CEDGE* const edgestore = dmst->edgestore;
   SCIP_Bool* const nodemark = dmst->nodemark;
   int nedges_new = 0;
   const int nnodes_in = mst_in->nnodes;

   assert(nnodes_in >= 1);

   nodemark[0] = TRUE;

   for( int i = 1; i < nnodes_in; ++i )
      nodemark[i] = FALSE;

   dcmstInsert(mst_in, adjcosts, 0, edgestore, nodemark, &max_path_edge, &nedges_new);

   assert(nedges_new == nnodes_in - 1);

   edgestore[nedges_new] = max_path_edge;
}


/** transforms edge-store to CSR  */
static inline
void dcmstGetCSRfromStore(
   const DCMST*          dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   const CEDGE* const edgestore = dmst->edgestore;
   int* const mst_start = mst_out->start;
   int* const mst_head = mst_out->head;
   SCIP_Real* const mst_cost = mst_out->cost;
   const int mst_nnodes = mst_out->nnodes;

   /* undirected edges */
   const int mst_nedges = mst_nnodes - 1;

   assert(mst_nnodes <= dmst->maxnnodes);
   assert(2 * mst_nedges == mst_out->nedges_max);

   BMSclearMemoryArray(mst_start, mst_nnodes + 1);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;

      assert(v1 >= 0 && v1 < mst_nnodes);
      assert(v2 >= 0 && v2 < mst_nnodes);

      mst_start[v1]++;
      mst_start[v2]++;
   }

   assert(mst_start[mst_nnodes] == 0);

   for( int i = 1; i <= mst_nnodes; ++i )
   {
      mst_start[i] += mst_start[i - 1];
   }

   assert(mst_start[mst_nnodes] == mst_out->nedges_max);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;
      const SCIP_Real cost = edgestore[i].cost;

      assert(mst_start[v1] >= 1);
      assert(mst_start[v2] >= 1);

      mst_head[--mst_start[v1]] = v2;
      mst_cost[mst_start[v1]] = cost;

      mst_head[--mst_start[v2]] = v1;
      mst_cost[mst_start[v2]] = cost;
   }
}


/** gets weight of MST from DCMST store */
static inline
SCIP_Real dcmstGetWeightFromStore(
   SCIP*                 scip,               /**< SCIP */
   int                   mst_nedges,         /**< number of edges */
   const DCMST*          dmst                /**< underlying structure */
)
{
   const CEDGE* const edgestore = dmst->edgestore;
   SCIP_Real weight = 0.0;

   assert(mst_nedges < dmst->maxnnodes);

   for( int i = 0; i < mst_nedges; ++i )
   {
#ifndef NDEBUG
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;

      assert(v1 >= 0 && v1 < dmst->maxnnodes);
      assert(v2 >= 0 && v2 < dmst->maxnnodes);
      assert(GE(dmst->edgestore[i].cost, 0.0));
      assert(LE(dmst->edgestore[i].cost, FARAWAY));
#endif

      weight += edgestore[i].cost;
   }

   return weight;
}


/** sets star position array to initial setting for current star degree */
static inline
void starSelectedPositionsReset(
   STAR*                 star                /**< the star */
)
{
   int* const edgesSelectedPos =  star->edgesSelectedPos;
   const int starDegree = star->starDegree;

   for( int i = 0; i < starDegree; i++ )
   {
      edgesSelectedPos[i] = i;
   }
}


/** fills array star->edgesSelected by using the current positions */
static inline
void starSelectedEdgesUpdate(
   STAR*                 star                /**< the star (in/out) */
)
{
   int* const edgesSelected = star->edgesSelected;
   const int* const edgesSelectedPos = star->edgesSelectedPos;
   const int* const edgeId = star->edgeId;
   const int starDegree = star->starDegree;

   assert(starDegree >= 3);

   for( int i = 0; i < starDegree; i++ )
   {
      const int pos = edgesSelectedPos[i];
      edgesSelected[i] = edgeId[pos];
   }
}

/** copies selected positions into given array */
static inline
void starSelectedPositionsCopy(
   const STAR*           star,               /**< the star (in/out) */
   int*                  posStorage          /**< to copy into */
)
{
   const int* const edgesSelectedPos = star->edgesSelectedPos;
   const int starDegree = star->starDegree;

   assert(starDegree >= 3);

   for( int i = 0; i < starDegree; i++ )
   {
      const int pos = edgesSelectedPos[i];
      assert(0 <= pos && pos < star->nodeDegree);

      posStorage[i] = pos;
   }
}


/** moves to next star */
static inline
void starSelectedPositionsSetNext(
   STAR*                 star                /**< the star (in/out) */
)
{
   int pos;
   const int nodeDegree = star->nodeDegree;
   const int starDegree = star->starDegree;
   int* const edgesSelectedPos = star->edgesSelectedPos;

   assert(3 <= starDegree && starDegree <= nodeDegree);

   /* all current positions are stored in edgesSelectedPos[0,...,starDegree-1] */

   /* check for each position, bottom-up, whether it can be increased without it hitting the border */
   for( pos = starDegree - 1; pos >= 0; pos-- )
   {
      const SCIP_Bool isLastPos = (pos == (starDegree - 1));
      const int border = isLastPos ? nodeDegree : edgesSelectedPos[pos + 1];

      /* still space? */
      if( edgesSelectedPos[pos] < border - 1 )
      {
         break;
      }
   }

   if( pos >= 0 )
   {
      edgesSelectedPos[pos]++;

      /* adapt all following positions */
      for( int i = pos + 1; i < starDegree; i++ )
      {
         edgesSelectedPos[i] = edgesSelectedPos[i - 1] + 1;
      }
   }
   else
   {
      assert(pos == -1);
      star->starDegree--;
      starSelectedPositionsReset(star);
   }
}


/** have we just move past the last star? */
static inline
SCIP_Bool starIsDeg2(
   const STAR*           star                /**< the star (in/out) */
)
{
   if( star->starDegree <= 2 )
   {
      assert(star->starDegree == 2);
      return TRUE;
   }

   return FALSE;
}


/*
 * Interface methods
 */

/** apply pseudo eliminations provided */
SCIP_RETCODE reduce_applyPseudoDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const SCIP_Bool*      pseudoDelNodes,     /**< node with pseudo deletable nodes */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real*            offsetp,            /**< offset pointer (for PC) */
   int*                  nelims              /**< number of eliminations */
)
{
   int adjvert[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real cutoffs[STP_DELPSEUDO_MAXNEDGES];
   SCIP_Real cutoffsrev[STP_DELPSEUDO_MAXNEDGES];
   const PATH* nodeTo3TermsPaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real* rootToNodeDist = redcostdata->rootToNodeDist;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real cutoffbound = redcostdata->cutoff;
   int* nodetouchcount;
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Bool success;
   const SCIP_Bool isPc = graph_pc_isPc(graph);
   const SCIP_Bool isExtendedOrg = graph->extended;

   assert(GE(cutoffbound, 0.0));
   assert(nodeTo3TermsPaths && rootToNodeDist && redcost);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodetouchcount, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      nodetouchcount[k] = 0;

   if( isPc )
      graph_pc_2orgcheck(scip, graph);

   *nelims = 0;

   for( int degree = 3; degree <= STP_DELPSEUDO_MAXGRAD; degree++ )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         SCIP_Real prize = -1.0;
         int edgecount = 0;
         SCIP_Bool rpc3term = FALSE;

         if( !pseudoDelNodes[k] || nodetouchcount[k] > 0 )
            continue;

         if( isPc && degree == 3 && graph_pc_knotIsNonLeafTerm(graph, k) && graph->grad[k] == 3 )
         {
            rpc3term = TRUE;
            prize = graph->prize[k];
         }
         else if( (degree != graph->grad[k] || Is_anyTerm(graph->term[k])) )
         {
            continue;
         }

         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            nodetouchcount[graph->head[e]]++;

         if( rpc3term )
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            {
               const int head = graph->head[e];
               if( !graph_pc_knotIsDummyTerm(graph, head) )
                  adjvert[edgecount++] = head;
            }
         }
         else
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               adjvert[edgecount++] = graph->head[e];
         }

         assert(edgecount == degree);

         edgecount = 0;
         for( int i = 0; i < degree - 1; i++ )
         {
            const int vert = adjvert[i];
            for( int i2 = i + 1; i2 < degree; i2++ )
            {
               const int vert2 = adjvert[i2];

               assert(edgecount < STP_DELPSEUDO_MAXNEDGES);

               cutoffs[edgecount] = cutoffbound - (rootToNodeDist[vert] + nodeTo3TermsPaths[vert2].dist);
               cutoffsrev[edgecount] = cutoffbound - (rootToNodeDist[vert2] + nodeTo3TermsPaths[vert].dist);

               edgecount++;
            }
         }

         assert(edgecount > 0);

#ifdef SCIP_DEBUG
         SCIPdebugMessage("try pseudo-deletion of ");
         graph_knot_printInfo(graph, k);
#endif

         /* now try to eliminate */
         SCIP_CALL( graph_knot_delPseudo(scip, graph, redcost, cutoffs, cutoffsrev, k, &success) );

         if( success )
         {
            (*nelims)++;
            graph->mark[k] = FALSE;

            SCIPdebugMessage("deletion successful! \n");

            if( rpc3term )
            {
               assert(isPc);
               assert(offsetp);
               assert(GE(prize, 0.0));

               *offsetp += prize;
            }
         }
         else
         {
            for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
            {
               nodetouchcount[graph->head[e]]--;
               assert(nodetouchcount[graph->head[e]] >= 0);
            }
         }
      }
   }

   assert(graph_valid(scip, graph));

   if( isPc && isExtendedOrg != graph->extended )
      graph_pc_2trans(scip, graph);


   SCIPfreeBufferArray(scip, &nodetouchcount);

   return SCIP_OKAY;
}


/** initializes BLC tree */
SCIP_RETCODE reduce_blctreeInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                graph,              /**< the graph */
   BLCTREE**             blctree             /**< to be initialized */
)
{
   assert(scip && graph && blctree);

   graph_mark(graph);

   SCIP_CALL( blctreeInitPrimitives(scip, graph, blctree) );
   SCIP_CALL( blctreeInitBottlenecks(scip, graph, *blctree) );

   return SCIP_OKAY;
}


/** rebuilds BLC tree (needs to be initializes) */
SCIP_RETCODE reduce_blctreeRebuild(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                graph,              /**< the graph */
   BLCTREE*              blctree             /**< to be initialized */
)
{
   assert(scip && graph && blctree);

   SCIP_CALL( blctreeInitBottlenecks(scip, graph, blctree) );

   return SCIP_OKAY;
}


/** frees BLC tree */
void reduce_blctreeFree(
   SCIP*                 scip,               /**< SCIP */
   BLCTREE**             blctree             /**< to be freed */
)
{
   assert(scip && blctree);

   SCIPfreeMemoryArray(scip, &((*blctree)->nodes_org2curr));
   SCIPfreeMemoryArray(scip, &((*blctree)->nodes_curr2org));
   SCIPfreeMemoryArray(scip, &((*blctree)->mstedges_isLink));
   SCIPfreeMemoryArray(scip, &((*blctree)->mst));

   SCIPfreeMemory(scip, blctree);
}


/** gets number of BLC MST edges */
int reduce_blctreeGetMstNedges(
   const BLCTREE*        blctree             /**< BLC tree */
)
{
   assert(blctree);
   assert(blctree->nnodes_curr > 0);

   return blctree->nnodes_curr - 1;
}


/** gets BLC MST edges */
void reduce_blctreeGetMstEdges(
   const GRAPH*          graph,              /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   int*                  edgelist            /**< of size nodes - 1 */
)
{
   const BLCNODE* mst;
   int nodecount;
   const int nnodes_curr = blctree->nnodes_curr;

   assert(edgelist);
   assert(graph_get_nNodes(graph) == blctree->nnodes_org);

   mst = blctree->mst;
   nodecount = 0;

   for( int i = 0; i < nnodes_curr; ++i )
   {
      const int edge = mst[i].edge;

      if( edge >= 0 )
      {
         edgelist[nodecount++] = edge;
      }
      else
      {
         assert(i == blctree->root);
         assert(edge == UNKNOWN);
      }
   }

   assert(nodecount == nnodes_curr - 1);
}


/** gets BLC MST edges (maximum) distances from tail to cut */
void reduce_blctreeGetMstEdgesToCutDist(
   const GRAPH*          graph,              /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   SCIP_Real* RESTRICT   tails2CutDist,      /**< of size nodes - 1 */
   SCIP_Real* RESTRICT   heads2CutDist       /**< of size nodes - 1 */
)
{
   const BLCNODE* mst;
   int nodecount;
   const int nnodes_curr = blctree->nnodes_curr;

   assert(tails2CutDist && heads2CutDist);
   assert(graph_get_nNodes(graph) == blctree->nnodes_org);

   mst = blctree->mst;
   nodecount = 0;

   for( int i = 0; i < nnodes_curr; ++i )
   {
      const int edge = mst[i].edge;

      if( edge >= 0 )
      {
         const SCIP_Real taildist = mst[i].dist_edgetail;
         const SCIP_Real headdist = mst[i].dist_edgehead;
         assert(GE(taildist, 0.0));
         assert(GE(headdist, 0.0));

#ifdef SCIP_DEBUG
         graph_edge_printInfo(graph, edge);
         printf("taildist=%f headdist=%f \n", taildist, headdist);
#endif

         tails2CutDist[nodecount] = taildist;
         heads2CutDist[nodecount++] = headdist;
      }
      else
      {
         assert(i == blctree->root);
         assert(edge == UNKNOWN);
      }
   }

   assert(nodecount == nnodes_curr - 1);
}


/** returns BLC MST edges state.
 *  I.e. whether the ede is between two terminal or not */
void reduce_blctreeGetMstEdgesState(
   const GRAPH*          graph,              /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   SCIP_Bool*            statelist           /**< of size nodes - 1 */
)
{
   const BLCNODE* mst = blctree->mst;
   int nodecount = 0;
   const int nnodes_curr = blctree->nnodes_curr;

   assert(graph && blctree && statelist);
   assert(graph_get_nNodes(graph) == blctree->nnodes_org);
   assert(blctree->mstedges_isLink);

   for( int i = 0; i < nnodes_curr; ++i )
   {
      const int edge = mst[i].edge;

      if( edge >= 0 )
      {
         statelist[nodecount++] = blctree->mstedges_isLink[i];
      }
      else
      {
         assert(i == blctree->root);
         assert(edge == UNKNOWN);
      }
   }

   assert(nodecount == nnodes_curr - 1);
}


/** gets BLC MST bottleneck costs */
void reduce_blctreeGetMstBottlenecks(
   const GRAPH*          graph,              /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   SCIP_Real*            costlist            /**< of size nodes - 1 */
)
{
   const BLCNODE* mst;
   int nodecount;
   const int nnodes_curr = blctree->nnodes_curr;

   assert(blctree && costlist);
   assert(graph_get_nNodes(graph) == blctree->nnodes_org);

   mst = blctree->mst;
   nodecount = 0;

   for( int i = 0; i < nnodes_curr; ++i )
   {
      const int edge = mst[i].edge;

      if( edge >= 0 )
      {
         costlist[nodecount++] = mst[i].edgebottleneck;
      }
      else
      {
         assert(i == blctree->root);
         assert(edge == UNKNOWN);
      }
   }

   assert(nodecount == nnodes_curr - 1);
}


/** initializes dynamic MST structure */
SCIP_RETCODE reduce_dcmstInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnnodes,          /**< maximum number of nodes that can be handled */
   DCMST**               dcmst               /**< to be initialized */
)
{
   DCMST* mst;

   assert(scip && dcmst);
   assert(maxnnodes >= 1);

   SCIP_CALL( SCIPallocMemory(scip, dcmst) );

   mst = *dcmst;

   mst->maxnnodes = maxnnodes;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->edgestore), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->adjcost_buffer), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->nodemark), maxnnodes) );


   return SCIP_OKAY;
}


/** adds node to CSR "mst_in" and saves result in "mst_out" */
void reduce_dcmstAddNode(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   assert(mst_in && adjcosts && dmst && mst_out);

   assert(reduce_dcmstMstIsValid(scip, mst_in));

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(mst_out->nedges_max == mst_in->nedges_max + 2);
   assert(mst_in->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst_in, adjcosts, dmst);

   dcmstGetCSRfromStore(dmst, mst_out);

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(reduce_dcmstMstIsValid(scip, mst_out));
}


/** Adds node to CSR "mst".
 *  NOTE: There needs to be enough space in CSR arrays for one more node! */
void reduce_dcmstAddNodeInplace(
   SCIP*                 scip,               /**< SCIP */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst                 /**< source/target */
)
{
   assert(mst && adjcosts && dmst);

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(mst->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst, adjcosts, dmst);

   mst->nnodes += 1;
   mst->nedges_max += 2;

   dcmstGetCSRfromStore(dmst, mst);

   assert(reduce_dcmstMstIsValid(scip, mst));
}

/** computes MST on 0 node */
void reduce_dcmstGet0NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 0);
   assert(mst->nedges_max == 0);

   start[0] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 1 node */
void reduce_dcmstGet1NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 1);
   assert(mst->nedges_max == 0);

   start[0] = 0;
   start[1] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 2 nodes */
void reduce_dcmstGet2NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost,           /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   SCIP_Real* const cost = mst->cost;
   int* const start = mst->start;
   int* const head = mst->head;

   assert(edgecost > 0.0);
   assert(mst->nnodes == 2);
   assert(mst->nedges_max == 2);

   start[0] = 0;
   start[1] = 1;
   start[2] = 2;

   head[0] = 1;
   head[1] = 0;

   cost[0] = edgecost;
   cost[1] = edgecost;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(edgecost, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 3 nodes */
void reduce_dcmstGet3NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost01,         /**< edge cost */
   SCIP_Real             edgecost02,         /**< edge cost */
   SCIP_Real             edgecost12,         /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   assert(0 && "implement me");
}


/** gets weight of MST extended along given vertex */
SCIP_Real reduce_dcmstGetExtWeight(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst,                /**< MST for which to compute extended weight */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst                /**< underlying structure */
)
{
   SCIP_Real weight;
#ifndef NDEBUG
   /* since we have a tree, |E_{ext}| = |E| + 1 = |V| */
   const int nedges_ext = mst->nnodes;
#endif

   assert(scip && adjcosts && dmst);
   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(mst->nedges_max % 2 == 0);
   assert((mst->nedges_max / 2) + 1 == nedges_ext);

   dcmstAddNode(mst, adjcosts, dmst);

   weight = dcmstGetWeightFromStore(scip, mst->nnodes, dmst);

   if( GT(weight, FARAWAY) )
      weight = FARAWAY;

   assert(GE(weight, 0.0));

   return weight;
}


/** gets weight of MST */
SCIP_Real reduce_dcmstGetWeight(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in              /**< source */
)
{
   SCIP_Real weight = 0.0;
   const int nedges = mst_in->nedges_max;
   const SCIP_Real* cost = mst_in->cost;

   assert(scip);
   assert(reduce_dcmstMstIsValid(scip, mst_in));

   for( int i = 0; i < nedges; i++ )
   {
      assert(cost[i] >= 0.0);

      weight += cost[i];
   }

   weight /= 2.0;

   if( GT(weight, FARAWAY) )
      weight = FARAWAY;

   assert(GE(weight, 0.0));

   return weight;
}


/** returns maximum number of nodes */
int reduce_dcmstGetMaxnnodes(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);

   return dmst->maxnnodes;
}


/** Returns buffer of size 'reduce_dcmstGetMaxnnodes'.
  * NOTE: buffer is never used within any other function, apart from allocation and freeing.
  * NOTE: in debug mode the array is initialized to -1.0 */
SCIP_Real* reduce_dcmstGetAdjcostBuffer(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);
   assert(dmst->adjcost_buffer);

#ifndef NDEBUG
   for( int i = 0; i < dmst->maxnnodes; i++ )
      dmst->adjcost_buffer[i] = -1.0;
#endif

   return dmst->adjcost_buffer;
}


/** frees dynamic MST structure */
void reduce_dcmstFree(
   SCIP*                 scip,               /**< SCIP */
   DCMST**               dcmst               /**< to be initialized */
)
{
   assert(scip && dcmst);

   SCIPfreeMemoryArray(scip, &((*dcmst)->nodemark));
   SCIPfreeMemoryArray(scip, &((*dcmst)->adjcost_buffer));
   SCIPfreeMemoryArray(scip, &((*dcmst)->edgestore));

   SCIPfreeMemory(scip, dcmst);
}


/** is the CSR a valid MST on any underlying graph (with number of nodes and edges of the CSR)? */
SCIP_Bool reduce_dcmstMstIsValid(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            cmst                /**< the MST candidate */
)
{
   SCIP_Bool* visited;
   const int* const head_csr = cmst->head;
   const int nnodes = cmst->nnodes;
   SCIP_Bool isValid = TRUE;

#ifndef NDEBUG
   const int* const start_csr = cmst->start;
#endif

   if( nnodes == 0 )
   {
      assert(cmst->nedges_max == 0);
      assert(start_csr[0] == 0);

      return TRUE;
   }

   assert(nnodes >= 1);
   assert(cmst->nedges_max % 2 == 0);
   assert(start_csr[0] == 0);

   if( !graph_csr_isValid(cmst, FALSE) )
   {
      SCIPdebugMessage("CSR is broken! \n");
      return FALSE;
   }

   if( cmst->nnodes != (cmst->nedges_max / 2) + 1 )
   {
      SCIPdebugMessage("wrong nodes/edges ratio \n");
      return FALSE;
   }

   if( nnodes == 1 )
   {
      return TRUE;
   }

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &visited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      visited[i] = FALSE;

   for( int i = 0; i < cmst->nedges_max; i++ )
   {
      const int head = head_csr[i];

      assert(head >= 0 && head < nnodes);

      visited[head] = TRUE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !visited[i] )
      {
         SCIPdebugMessage("mst does not contain node %d \n", i);

         isValid = FALSE;
         break;
      }
   }

   SCIPfreeMemoryArray(scip, &visited);

   return isValid;
}


/** initializes STAR structure */
SCIP_RETCODE reduce_starInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxdegree,          /**< maximum node degree that can be handled */
   STAR**                star                /**< the star */
)
{
   STAR* s;

   assert(scip && star);
   assert(maxdegree >= 3);

   SCIP_CALL( SCIPallocMemory(scip, star) );

   s = *star;

   s->starcenter = -1;
   s->nodeDegree = -1;
   s->starDegree = -1;
   s->maxNodeDegree = maxdegree;
   s->allStarsChecked = FALSE;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgeId), maxdegree) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgesSelected), maxdegree) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(s->edgesSelectedPos), maxdegree) );

   return SCIP_OKAY;
}


/** frees STAR structure */
void reduce_starFree(
   SCIP*                 scip,               /**< SCIP */
   STAR**                star                /**< the star */
)
{
   STAR* s;
   assert(scip && star);

   s = *star;
   assert(s);

   SCIPfreeMemoryArray(scip, &(s->edgesSelectedPos));
   SCIPfreeMemoryArray(scip, &(s->edgesSelected));
   SCIPfreeMemoryArray(scip, &(s->edgeId));

   SCIPfreeMemory(scip, star);
}


/** resets star data structure with new node data */
void reduce_starReset(
   const GRAPH*          g,                  /**< graph */
   int                   node,               /**< the node (degree <= STP_DELPSEUDO_MAXGRAD) */
   STAR*                 star                /**< the star */
)
{
   assert(g && star);
   assert(0 <= node && node < g->knots);
   assert(g->grad[node] <= star->maxNodeDegree);

   star->nodeDegree = g->grad[node];
   star->starDegree = star->nodeDegree;
   star->allStarsChecked = FALSE;
   star->starcenter = node;

   for( int e = g->outbeg[node], i = 0; e != EAT_LAST; e = g->oeat[e], i++ )
   {
      assert(i < star->nodeDegree);
      star->edgeId[i] = e;
   }

   /* initially, select the entire star */
   starSelectedPositionsReset(star);
}


/** gets center */
int reduce_starGetCenter(
   const STAR*           star                /**< the star (in/out) */
)
{
   assert(star);

   return star->starcenter;
}

/** gets next star */
const int* reduce_starGetNext(
   STAR*                 star,               /**< the star (in/out) */
   int*                  nedges              /**< number of edges of next star (out) */
)
{
   assert(star && nedges);
   assert(!reduce_starAllAreChecked(star));

   *nedges = star->starDegree;
   starSelectedEdgesUpdate(star);
   starSelectedPositionsSetNext(star);

   /* just finished? */
   if( starIsDeg2(star) )
      star->allStarsChecked = TRUE;

   assert(3 <= *nedges && *nedges <= star->maxNodeDegree);

   return star->edgesSelected;
}


/** gets next star */
const int* reduce_starGetNextAndPosition(
   STAR*                 star,               /**< the star (in/out) */
   int*                  position,           /**< array to store the positions */
   int*                  nedges              /**< number of edges of next star (out) */
)
{
   assert(star);
   assert(!reduce_starAllAreChecked(star));

   if( nedges )
   {
      *nedges = star->starDegree;
   }

   starSelectedEdgesUpdate(star);
   starSelectedPositionsCopy(star, position);
   starSelectedPositionsSetNext(star);

   /* just finished? */
   if( starIsDeg2(star) )
      star->allStarsChecked = TRUE;

   if( nedges )
   {
      assert(3 <= *nedges && *nedges <= star->maxNodeDegree);
   }

   return star->edgesSelected;
}



/** gets ruled out edges after termination */
const int* reduce_starGetRuledOutEdges(
   STAR*                 star,               /**< the star */
   int*                  nedges              /**< number of edges of next star (out) */
)
{
   assert(star);
   assert(reduce_starAllAreChecked(star));

   // todo fill later
   // todo also add a method to abort early if no ruled-out edges can be found anymore

   return NULL;
}


/** sets current star to ruled-out */
void reduce_starSetRuledOut(
   STAR*                 star                /**< the star */
)
{
   assert(star);
   assert(!reduce_starAllAreChecked(star));


   // todo fill later for edge rule out
}


/** sets current star to failed */
void reduce_starSetFailed(
   STAR*                 star                /**< the star */
)
{
   assert(star);
   assert(!reduce_starAllAreChecked(star));

   // todo fill later for edge rule out
}


/** have all stars been checked? */
SCIP_Bool reduce_starAllAreChecked(
   const STAR*           star                /**< the star */
)
{
   assert(star);

   return star->allStarsChecked;
}


/** builds reduced costs data structure and returns it.
 * NOTE: memory needs to be freed again! */
SCIP_RETCODE reduce_redcostdataInit(
   SCIP*                 scip,               /**< SCIP */
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value or -1.0 if not used */
   int                   redCostRoot,        /**< graph root for reduced cost calculation */
   REDCOST*              redcostdata         /**< data to initialize */
)
{
   SCIP_Real* redEdgeCost;
   SCIP_Real* rootToNodeDist;
   PATH* nodeTo3TermsPaths;
   int* nodeTo3TermsBases;

   assert(scip);
   assert(nnodes >= 0);
   assert(nedges >= 0);
   assert(nedges % 2 == 0);
   assert(redCostRoot >= 0);
   assert(GE(cutoff, 0.0) || EQ(cutoff, 0.0));
   assert(redcostdata);

   SCIP_CALL( SCIPallocMemoryArray(scip, &redEdgeCost, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &rootToNodeDist, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsPaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeTo3TermsBases, 3 * nnodes) );

   redcostdata->redEdgeCost = redEdgeCost;
   redcostdata->rootToNodeDist = rootToNodeDist;
   redcostdata->nodeTo3TermsPaths = nodeTo3TermsPaths;
   redcostdata->nodeTo3TermsBases = nodeTo3TermsBases;
   redcostdata->cutoff = cutoff;
   redcostdata->redCostRoot = redCostRoot;

#ifndef NDEBUG
   redcostdata->nnodes = nnodes;
   redcostdata->nedges = nedges;
#endif

   return SCIP_OKAY;
}


/** frees member arrays */
void reduce_redcostdataFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   REDCOST*              redcostdata         /**< data */
)
{
   assert(scip && redcostdata);

   SCIPfreeMemoryArray(scip, &(redcostdata->nodeTo3TermsBases));
   SCIPfreeMemoryArray(scip, &(redcostdata->nodeTo3TermsPaths));
   SCIPfreeMemoryArray(scip, &(redcostdata->rootToNodeDist));
   SCIPfreeMemoryArray(scip, &(redcostdata->redEdgeCost));
}
