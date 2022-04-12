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

/**@file   reduce_sdgraph.c
 * @brief  bottleneck distance graph methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements methods for Steiner tree problem special distance (bottleneck Steiner distance) graph.
 * This graph is the (complete) distance graph on the terminal vertex set.
 * Note that the complete graph is in general not built.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "reduce.h"
#include "misc_stp.h"
#include "portab.h"

#define STP_SDQUERYFULL_MAXTERMS 200

//#define SDQUERYFULL_HALFMATRIX


/*
 * STRUCTS
 */


/** see reduce.h */
struct special_distance_graph
{
   GRAPH*                distgraph;          /**< (complete) distance graph */
   PATH*                 sdmst;              /**< MST on sdgraph */
   SCIP_Real*            mstcosts;           /**< maximum MST edge costs in descending order */
   SCIP_Real*            mstsdist;           /**< helper array; each entry needs to -1.0; of size nnodesorg */
   int*                  nodemapOrgToDist;   /**< node mapping from original graph to distance graph */
   STP_Bool*             halfedge_isInMst;   /**< signifies whether edge of original graph is part of MST
                                                  NOTE: operates on edges / 2! */
   SCIP_Real             mstmaxcost;         /**< maximum edge cost */
   int                   nnodesorg;          /**< number of nodes of original graph */
   int                   nedgesorg;          /**< number of edges of original graph */
   SCIP_Bool             mstcostsReady;      /**< are the mstcosts already available? */
   SCIP_Bool             edgemarkReady;      /**< edge mark available? */
   SCIP_Bool             usingRMQ;           /**< is RMQ being used? */
   /* Full query data: */
   SCIP_Real*            fullq_dists;        /**< of size |T| * |T - 1| or NULL */
   int*                  fullq_nodeToIdx;    /**< maps nodes of original graphs to full index; of size |V| */
#ifdef SDQUERYFULL_HALFMATRIX
   int*                  fullq_IdxToStart;   /**< start for dists */
#endif
   int                   fullq_size;         /**< |T - 1| * |T - 1| or  |T - 1| * |T - 1|  / 2  */
   int                   fullq_dimension;    /**< |T - 1| */
   /* RMQ query data: */
   int*                  rmq_nodeToRmqEntry; /**< of size |V| */
   SCIP_Real*            rmq_sparseTable;    /**< of size rmq_log * 2 |T| */
   SCIP_Real*            rmq_edgecosts;      /**< of size 2|T| - 3 */
   int                   rmq_loglength;      /**< log2(2|T| - 3) */
};



/** Simple binary tree node.
 *  Value of UNKNOWN signifies that child does not exist. */
typedef struct binary_tree_node
{
   int                   child1;             /**< first child */
   int                   child2;             /**< second child */
} BINARYNODE;


/** data needed for building the RMQ based SD query structures */
typedef struct lowest_common_ancestor_tree_builder
{
   BINARYNODE*           lcatree;            /**< binary tree used for LCA computation; of size |T| - 1 */
   SCIP_Real*            lcatree_costs;      /**< cost per node of binary tree used for LCA computation; of size 2|T| - 3 */
   int*                  termToLcatreeNode;  /**< node mapping from terminals (SD graph nodes) to binary LCA tree nodes */
   int*                  lcatreeNodeToInds;  /**< mapping from binary LCA tree nodes to RMQ/Full indices; of size |T| - 1 */
} LCABUILDER;


/*
 * local methods
 */

static const int deBruijnBits[32] =
{
   0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
   8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
};

/** returns the floor of log2(number) */
static inline
int log2floor(
   uint32_t              number              /**< number */
   )
{
   uint32_t v = number;
   v |= v >> 1;
   v |= v >> 2;
   v |= v >> 4;
   v |= v >> 8;
   v |= v >> 16;

   assert(deBruijnBits[(uint32_t)(v * 0x07C4ACDDU) >> 27] == (int) log2(number));

   return deBruijnBits[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}

#ifdef SCIP_DEBUG
/* prints the tree */
static
void lcabuilderPrintLcaTree(
   const SDGRAPH*        sdgraph,            /**< the SD graph */
   const LCABUILDER*     lcabuilder          /**< the builder */
)
{
   const BINARYNODE* const lcatree = lcabuilder->lcatree;
   const SCIP_Real* const lcatree_costs = lcabuilder->lcatree_costs;
   const int* const termToLcatreeNode = lcabuilder->termToLcatreeNode;
   const int nterms = graph_get_nNodes(sdgraph->distgraph);

   assert(termToLcatreeNode && lcatree && lcatree_costs && lcabuilder->lcatree && lcabuilder->termToLcatreeNode);

   printf("\n SD QUERY LCA TREE: \n");

   for( int i = 0; i < nterms - 1; i++ )
   {
      const SCIP_Real cost = lcatree_costs[i];
      printf("node %d: \n", i);

      printf("...cost=%f \n", cost);
      printf("...child1=%d, child2=%d \n", lcatree[i].child1, lcatree[i].child2);
   }

   printf("\n");

   for( int i = 0; i < nterms; i++ )
   {
      printf("term %d attached to tree node %d \n", i, termToLcatreeNode[i]);
   }
}
#endif



/** Gets special distance (i.e. bottleneck distance) from graph.
 *  Corresponds to bottleneck length of path between term1 and term2 on distance graph */
static inline
SCIP_Real sdgraphGetSd(
   int                    term1,             /**< terminal 1 */
   int                    term2,             /**< terminal 2 */
   SDGRAPH*               sdgraph            /**< the SD graph */
)
{
   SCIP_Real* RESTRICT mstsdist = sdgraph->mstsdist;
   const GRAPH* distgraph = sdgraph->distgraph;
   const PATH* sdmst = sdgraph->sdmst;
   const int* nodesid = sdgraph->nodemapOrgToDist;
   const int sdnode1 = nodesid[term1];
   const int sdnode2 = nodesid[term2];
   SCIP_Real sdist = 0.0;
   int tempnode = sdnode1;

   assert(sdnode1 != sdnode2);
   assert(distgraph->source == 0);

   mstsdist[tempnode] = 0.0;

   /* not at root? */
   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;

      assert(distgraph->head[ne] == tempnode);
      tempnode = distgraph->tail[ne];

      if( distgraph->cost[ne] > sdist )
         sdist = distgraph->cost[ne];

      mstsdist[tempnode] = sdist;
      if( tempnode == sdnode2 )
         break;
   }

   /* already finished? */
   if( tempnode == sdnode2 )
   {
      tempnode = 0;
   }
   else
   {
      tempnode = sdnode2;
      sdist = 0.0;
   }

   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;
      tempnode = distgraph->tail[ne];

      if( distgraph->cost[ne] > sdist )
         sdist = distgraph->cost[ne];

      assert(GE(mstsdist[tempnode], 0.0) == (mstsdist[tempnode] > -0.5));

      /* already visited? */
      if( mstsdist[tempnode] > -0.5 )
      {
         if( mstsdist[tempnode] > sdist )
            sdist = mstsdist[tempnode];
         break;
      }

#ifndef NDEBUG
      assert(EQ(mstsdist[tempnode], -1.0));

      if( tempnode == 0 )
      {
         assert(sdnode1 == 0);
      }
#endif
   }

   /* restore mstsdist */
   tempnode = sdnode1;
   mstsdist[tempnode] = -1.0;
   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;
      tempnode = distgraph->tail[ne];
      mstsdist[tempnode] = -1.0;
      if( tempnode == sdnode2 )
         break;
   }

   assert(GE(sdist, 0.0));

   return sdist;
}

/** gets length of RMQ array */
static inline
int sdqueryGetRmqLength(
   const SDGRAPH*        sdgraph              /**< the SD graph */
   )
{
   const int nterms = graph_get_nNodes(sdgraph->distgraph);
   assert(nterms >= 2);
   assert((2 * nterms - 3) >= 1);

   return (2 * nterms - 3);
}

/** initializes LCA builder */
static
SCIP_RETCODE sdqueryLcaBuilderInit(
   SCIP*                 scip,               /**< SCIP */
   const SDGRAPH*        sdgraph,            /**< the SD graph */
   LCABUILDER**          lcabuilder          /**< the builder */
)
{
   LCABUILDER* rmqb;
   const int nterms = graph_get_nNodes(sdgraph->distgraph);

   assert(nterms >= 2);

   SCIP_CALL( SCIPallocMemory(scip, lcabuilder) );
   rmqb = *lcabuilder;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(rmqb->lcatree), nterms - 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rmqb->lcatree_costs), nterms - 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rmqb->lcatreeNodeToInds), nterms - 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(rmqb->termToLcatreeNode), nterms) );

   for( int i = 0; i < nterms - 1; i++ )
   {
      rmqb->lcatree[i].child1 = UNKNOWN;
      rmqb->lcatree[i].child2 = UNKNOWN;
   }

#ifndef NDEBUG
   for( int i = 0; i < nterms - 1; i++ )
   {
      rmqb->lcatreeNodeToInds[i] = -1;
   }

   for( int i = 0; i < nterms; i++ )
   {
      rmqb->termToLcatreeNode[i] = -1;
   }
#endif

   return SCIP_OKAY;
}


/** frees LCA builder */
static
void sdqueryLcaBuilderFree(
   SCIP*                 scip,               /**< SCIP */
   LCABUILDER**          lcabuilder          /**< the builder */
)
{
   LCABUILDER* rmqb = *lcabuilder;

   SCIPfreeMemoryArray(scip, &(rmqb->termToLcatreeNode));
   SCIPfreeMemoryArray(scip, &(rmqb->lcatreeNodeToInds));
   SCIPfreeMemoryArray(scip, &(rmqb->lcatree_costs));
   SCIPfreeMemoryArray(scip, &(rmqb->lcatree));

   SCIPfreeMemory(scip, lcabuilder);
}


/* attaches node */
static inline
void sdqueryAttachBinaryTreeNode(
   const GRAPH*          distgraph,          /**< SD distance graph */
   int                   parentposition,     /**< position of parent (edge) in binary tree array */
   int                   graphnode,          /**< graph node to attach */
   LCABUILDER*           lcabuilder,         /**< the builder */
   UF*                   uf                  /**< union-find data structure */
)
{
   BINARYNODE* RESTRICT lcatree = lcabuilder->lcatree;
   int* RESTRICT termToLcatreeNode = lcabuilder->termToLcatreeNode;
   const int childidentifier = SCIPStpunionfindFind(uf, graphnode);
   const int nterms = graph_get_nNodes(distgraph);

   graph_knot_isInRange(distgraph, graphnode);
   assert(0 <= parentposition && parentposition < nterms - 1);

   if( childidentifier == graphnode )
   {
      assert(graphnode < nterms);

      /* just need to set a pointer in this case */
      termToLcatreeNode[graphnode] = parentposition;
   }
   else
   {
      /* child identifier marks an edge position */
      const int childosition = childidentifier - nterms;
      assert(0 <= childosition && childosition < nterms - 1);

      if( lcatree[parentposition].child1 == UNKNOWN )
      {
         lcatree[parentposition].child1 = childosition;
      }
      else
      {
         assert(lcatree[parentposition].child2 == UNKNOWN);
         lcatree[parentposition].child2 = childosition;
      }
   }

   SCIPStpunionfindUnion(uf, nterms + parentposition, childidentifier, FALSE);
}


/** builds binary tree for LCA computation */
static
SCIP_RETCODE sdqueryBuildBinaryTree(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH*              sdgraph,            /**< the SD graph */
   LCABUILDER*           lcabuilder          /**< the builder */
)
{
   const GRAPH* const distgraph = sdgraph->distgraph;
   const PATH* const sdmst = sdgraph->sdmst;
   int* RESTRICT mstedges_id;
   SCIP_Real* RESTRICT lcatree_costs = lcabuilder->lcatree_costs;
   const int nterms = graph_get_nNodes(sdgraph->distgraph);
   int edgecount = 0;

   /* Union-find: entries 0,...,|T|-1 are used for terminals,
    *             entries |T|,...,2 |T| - 2 are used for the ordered MST edges */
   UF uf;

   SCIP_CALL( SCIPStpunionfindInit(scip, &uf, 2 * nterms - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mstedges_id, nterms - 1) );

   /* fill the MST edge arrays */
   for( int i = 0; i < nterms; i++ )
   {
      const int edge = sdmst[i].edge;

      if( edge >= 0 )
      {
         assert(graph_edge_isInRange(sdgraph->distgraph, edge));
         assert(EQ(sdmst[i].dist, sdgraph->distgraph->cost[edge]));

         mstedges_id[edgecount] = edge;
         lcatree_costs[edgecount++] = sdmst[i].dist;
      }
   }
   assert(edgecount == nterms - 1);
   assert(edgecount >= 1);

   /* build the actual tree */

   SCIPsortDownRealInt(lcatree_costs, mstedges_id, edgecount);

   for( int i = edgecount - 1; i >= 0; i-- )
   {
      const int edge = mstedges_id[i];
      const int tail = distgraph->tail[edge];
      const int head = distgraph->head[edge];

      sdqueryAttachBinaryTreeNode(distgraph, i, tail, lcabuilder, &uf);
      sdqueryAttachBinaryTreeNode(distgraph, i, head, lcabuilder, &uf);
   }

   SCIPfreeBufferArray(scip, &mstedges_id);
   SCIPStpunionfindFreeMembers(scip, &uf);

#ifdef SCIP_DEBUG
   lcabuilderPrintLcaTree(sdgraph, lcabuilder);
#endif

   return SCIP_OKAY;
}


/** Builds RMQ by DFS.
 *  NOTE: recursive method. */
static
void sdqueryRmqDfs(
   int                   root,               /**< root node (LCA tree position) */
   const BINARYNODE*     lcatree,            /**< tree */
   const SCIP_Real*      lcatree_costs,      /**< LCA tree costs */
   SCIP_Real*            rmq_edgecosts,      /**< edge costs for RMQ */
   int*                  lcatreeNodeToRmq,   /**< mapping from binary LCA tree nodes to RMQ indices */
   int*                  rmq_count           /**< current position */
)
{
   const BINARYNODE bnode = lcatree[root];

   assert(lcatreeNodeToRmq[root] == -1);
   lcatreeNodeToRmq[root] = *rmq_count;
   rmq_edgecosts[(*rmq_count)++] = lcatree_costs[root];

   SCIPdebugMessage("checking node %d \n", root);

   if( bnode.child1 == UNKNOWN )
   {
      assert(bnode.child2 == UNKNOWN );
      return;
   }

   sdqueryRmqDfs(bnode.child1, lcatree, lcatree_costs, rmq_edgecosts, lcatreeNodeToRmq, rmq_count);
   rmq_edgecosts[(*rmq_count)++] = lcatree_costs[root];

   if( bnode.child2 == UNKNOWN )
   {
      return;
   }

   sdqueryRmqDfs(bnode.child2, lcatree, lcatree_costs, rmq_edgecosts, lcatreeNodeToRmq, rmq_count);
   rmq_edgecosts[(*rmq_count)++] = lcatree_costs[root];
}


/** builds RMQ sparse table */
static
void sdqueryBuildRmqSparseTable(
   SDGRAPH*              sdgraph              /**< the SD graph */
   )
{
   const SCIP_Real* const rmq_edgecosts = sdgraph->rmq_edgecosts;
   SCIP_Real* const rmq_sparsetable = sdgraph->rmq_sparseTable;
   const int rmq_length = sdqueryGetRmqLength(sdgraph);
   const int logt = sdgraph->rmq_loglength;
   const int rowlength = logt + 1;

   SCIPdebugMessage("\n building sparse table with floor(log2(|T|))=%d \n", logt);

   /* base computation */
   for( int i = 0; i < rmq_length; i++ )
   {
      const int start_i = i * rowlength;
      rmq_sparsetable[start_i] = rmq_edgecosts[i];
      SCIPdebugMessage("max[%d, %d)=%f \n", i, i + 1, rmq_sparsetable[start_i]);
   }

   /* main (dynamic programming) computation */
   for( int j = 1; j <= logt; j++ )
   {
      const int shift = 1 << (unsigned int) (j - 1);

      SCIPdebugMessage("round %d (2^%d), shift=%d \n", j, j, shift);

      for( int i = 0; i < rmq_length; i++ )
      {
         const int start_i = i * rowlength;

         if( i + shift < rmq_length )
         {
            /* position for interval [i, i + 2^(j-1)) */
            const int pos_iprev = start_i + j - 1;
            const int start_shift = (i + shift) * rowlength;
            /* position for interval [(i+2^(j-1)), (i 2^(j-1)) + 2^(j-1)) */
            const int pos_shift = start_shift + j - 1 ;

            assert(pos_shift < rmq_length * rowlength);

            rmq_sparsetable[start_i + j] = MAX(rmq_sparsetable[pos_iprev], rmq_sparsetable[pos_shift]);
         }
         else
         {
            /* NOTE: basically debug marker */
            rmq_sparsetable[start_i + j] = FARAWAY;
         }

         SCIPdebugMessage("max[%d, 2^%d)=%f \n", i, j, rmq_sparsetable[start_i + j]);
      }
   }
}


/** builds RMQ, including sparse table representation */
static
void sdqueryBuildRmq(
   SDGRAPH*              sdgraph,            /**< the SD graph */
   LCABUILDER*           lcabuilder          /**< the builder */
)
{
   const BINARYNODE* const lcatree = lcabuilder->lcatree;
   const SCIP_Real* const lcatree_costs = lcabuilder->lcatree_costs;
   SCIP_Real* const rmq_edgecosts = sdgraph->rmq_edgecosts;
   int* lcatreeNodeToRmq = lcabuilder->lcatreeNodeToInds;
   int rmq_length = 0;

   /* NOTE: because we have sorted in ascending edge cost order, the root is 0 */

   sdqueryRmqDfs(0, lcatree, lcatree_costs, rmq_edgecosts, lcatreeNodeToRmq, &rmq_length);
   assert(rmq_length == sdqueryGetRmqLength(sdgraph));

#ifdef SCIP_DEBUG
   for( int i = 0; i < rmq_length; i++ )
   {
      printf("i=%d ocst=%f \n", i, rmq_edgecosts[i]);
   }
#endif

   sdqueryBuildRmqSparseTable(sdgraph);
}


/** builds RMQ map */
static
void sdqueryBuildNodesToRmqMap(
   const GRAPH*          g,                  /**< graph */
   const LCABUILDER*     lcabuilder,         /**< the builder */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   const int* const termToLcatreeNode = lcabuilder->termToLcatreeNode;
   const int* const nodemapOrgToDist = sdgraph->nodemapOrgToDist;
   const int* const lcatreeNodeToRmq = lcabuilder->lcatreeNodeToInds;
   int* const rmq_nodeToRmqEntry = sdgraph->rmq_nodeToRmqEntry;
   const int nnodes = graph_get_nNodes(g);

   assert(termToLcatreeNode && nodemapOrgToDist && lcatreeNodeToRmq && rmq_nodeToRmqEntry);

   for( int i = 0; i < nnodes; i++ )
   {
      assert(rmq_nodeToRmqEntry[i] == UNKNOWN);

      if( Is_term(g->term[i]) )
      {
         const int node_distgraph = nodemapOrgToDist[i];
         const int node_lcatree = termToLcatreeNode[node_distgraph];
         const int node_rmq = lcatreeNodeToRmq[node_lcatree];

         assert(graph_knot_isInRange(sdgraph->distgraph, node_distgraph));
         rmq_nodeToRmqEntry[i] = node_rmq;

         SCIPdebugMessage("node %d to %d \n", i, node_rmq);
      }
   }
}



/** initializes SD constant time query data */
static
SCIP_RETCODE sdqueryRmqInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int rmqlength = sdqueryGetRmqLength(sdgraph);
   const int logrmqlength = (int) log2(rmqlength);
   LCABUILDER* lcabuilder;
#ifndef NDEBUG
   const int nterms = graph_get_nNodes(sdgraph->distgraph);

   assert(nterms >= 2);
   assert(nterms == g->terms);
   assert(logrmqlength >= 0);
   assert(!sdgraph->rmq_edgecosts);
   assert(!sdgraph->rmq_sparseTable);
   assert(sdgraph->usingRMQ);
#endif

   sdgraph->rmq_loglength = logrmqlength;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->rmq_nodeToRmqEntry), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->rmq_sparseTable), (logrmqlength + 1) * rmqlength) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->rmq_edgecosts), rmqlength) );

   for( int i = 0; i < nnodes; i++ )
      sdgraph->rmq_nodeToRmqEntry[i] = UNKNOWN;

   SCIP_CALL( sdqueryLcaBuilderInit(scip, sdgraph, &lcabuilder) );
   SCIP_CALL( sdqueryBuildBinaryTree(scip, sdgraph, lcabuilder) );
   sdqueryBuildRmq(sdgraph, lcabuilder);
   sdqueryBuildNodesToRmqMap(g, lcabuilder, sdgraph);

   sdqueryLcaBuilderFree(scip, &lcabuilder);

   return SCIP_OKAY;
}





/** builds map */
static
void sdqueryBuildNodesToFullMap(
   const GRAPH*          g,                  /**< graph */
   const LCABUILDER*     lcabuilder,         /**< the builder */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   const int* const termToLcatreeNode = lcabuilder->termToLcatreeNode;
   const int* const nodemapOrgToDist = sdgraph->nodemapOrgToDist;
   int* const fullq_nodeToIdx = sdgraph->fullq_nodeToIdx;
   const int nnodes = graph_get_nNodes(g);

   assert(termToLcatreeNode && nodemapOrgToDist && fullq_nodeToIdx);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         const int node_distgraph = nodemapOrgToDist[i];
         const int node_lcatree = termToLcatreeNode[node_distgraph];

         assert(graph_knot_isInRange(sdgraph->distgraph, node_distgraph));
         fullq_nodeToIdx[i] = node_lcatree;

         SCIPdebugMessage("FullMap: node %d to %d \n", i, node_lcatree);
      }
      else
      {
         fullq_nodeToIdx[i] = UNKNOWN;
      }
   }
}


/** Gets special distance (i.e. bottleneck distance) from graph.
 *  Corresponds to bottleneck length of path between term1 and term2 on distance graph */
static inline
SCIP_Real sdqueryGetSd(
   int                    term1,             /**< terminal 1 */
   int                    term2,             /**< terminal 2 */
   const SDGRAPH*         sdgraph            /**< the SD graph */
)
{
   const int* const rmq_nodeToRmqEntry = sdgraph->rmq_nodeToRmqEntry;
   SCIP_Real sd;
   int rmqindex_min;
   int rmqindex_max;

   assert(term1 != term2);

   SCIPdebugMessage("query %d--%d \n", term1, term2);

   if( rmq_nodeToRmqEntry[term1] <= rmq_nodeToRmqEntry[term2] )
   {
      rmqindex_min = rmq_nodeToRmqEntry[term1];
      rmqindex_max = rmq_nodeToRmqEntry[term2];
   }
   else
   {
      rmqindex_min = rmq_nodeToRmqEntry[term2];
      rmqindex_max = rmq_nodeToRmqEntry[term1];
   }

   assert(rmqindex_min <= rmqindex_max);
   assert(0 <= rmqindex_min);
   assert(rmqindex_max < sdqueryGetRmqLength(sdgraph));

   if( rmqindex_min == rmqindex_max )
   {
      sd = sdgraph->rmq_edgecosts[rmqindex_min];
      SCIPdebugMessage("...mapped to same RMQ index \n");
   }
   else
   {
      const SCIP_Real* const rmq_sparsetable = sdgraph->rmq_sparseTable;
      const int rowlength = sdgraph->rmq_loglength + 1;
      const unsigned int msbit = (unsigned int) log2floor(rmqindex_max - rmqindex_min);
      const int pos_lower = rmqindex_min * rowlength + msbit;
      const int pos_upper = (rmqindex_max - (int) (1 << msbit)) * rowlength + msbit;

      SCIPdebugMessage("MSB=%u \n", msbit);
      SCIPdebugMessage("RMQ indices: %d, %d \n", rmqindex_min, rmqindex_max);
      SCIPdebugMessage("start indices: %d, %d  \n", rmqindex_min, (rmqindex_max - (int) (1 << msbit)));
      SCIPdebugMessage("max{%f, %f} \n", rmq_sparsetable[pos_lower], rmq_sparsetable[pos_upper]);

      sd = MAX(rmq_sparsetable[pos_lower], rmq_sparsetable[pos_upper]);
   }

   assert(LT(sd, FARAWAY));
   //printf("%f vs %f \n", sd, sdgraphGetSd(term1, term2, (SDGRAPH*) sdgraph));

   assert(EQ(sd, sdgraphGetSd(term1, term2, (SDGRAPH*) sdgraph)));

   return sd;
}


/** Gets special distance (i.e. bottleneck distance) from graph.
 *  Corresponds to bottleneck length of path between term1 and term2 on distance graph */
static inline
SCIP_Real sdqueryFullGetSd(
   int                    term1,             /**< terminal 1 */
   int                    term2,             /**< terminal 2 */
   const SDGRAPH*         sdgraph            /**< the SD graph */
)
{
   int pos;
   const int index1 = sdgraph->fullq_nodeToIdx[term1];
   const int index2 = sdgraph->fullq_nodeToIdx[term2];

   SCIPdebugMessage("query %d--%d \n", term1, term2);
   assert(term1 != term2);
   assert(sdgraph->fullq_dimension >= 1);
   assert(0 <= index1 && index1 < sdgraph->fullq_dimension);
   assert(0 <= index2 && index2 < sdgraph->fullq_dimension);

#ifdef SDQUERYFULL_HALFMATRIX
   if( index1 <= index2 )
   {
      pos = sdgraph->fullq_IdxToStart[index2] + index1;
      //pos = (index2 * (index2 + 1)) / 2 + index1;
   }
   else
   {
      pos = sdgraph->fullq_IdxToStart[index1] + index2;
     // pos = (index1 * (index1 + 1)) / 2 + index2;
   }
#else
   pos = index1 * sdgraph->fullq_dimension + index2;
#endif

   assert(pos < sdgraph->fullq_size);
   return sdgraph->fullq_dists[pos];
}


/** initializes full SD query data structure by DFS */
static
void sdqueryFullDfs(
   int                   root,               /**< root node (LCA tree position) */
   const BINARYNODE*     lcatree,            /**< tree */
   int                   nlcanodes,          /**< number of nodes */
   const SCIP_Real*      lcatree_costs,      /**< LCA tree costs */
   SCIP_Bool* RESTRICT   nodes_isVisited,    /**< per node: is visited? */
   UF*                   uf,                 /**< union-find data structure */
   SCIP_Real* RESTRICT   fullq_dists         /**< distances to be computed */
)
{
   const BINARYNODE bnode = lcatree[root];
#ifdef SDQUERYFULL_HALFMATRIX
   const int offset_root = (root * (root + 1)) / 2;
#endif

   assert(!nodes_isVisited[root]);
   SCIPdebugMessage("checking node %d (cost=%f)\n", root, lcatree_costs[root]);

   nodes_isVisited[root] = TRUE;

   if( bnode.child1 != UNKNOWN )
   {
      if( !nodes_isVisited[bnode.child1] )
      {
         sdqueryFullDfs(bnode.child1, lcatree, nlcanodes, lcatree_costs, nodes_isVisited, uf, fullq_dists);
         /* NOTE: root is also made the identifier of the set */
         SCIPStpunionfindUnion(uf, root, bnode.child1, FALSE);
      }
   }

   if( bnode.child2 != UNKNOWN )
   {
      if( !nodes_isVisited[bnode.child2] )
      {
         sdqueryFullDfs(bnode.child2, lcatree, nlcanodes, lcatree_costs, nodes_isVisited, uf, fullq_dists);
         SCIPStpunionfindUnion(uf, root, bnode.child2, FALSE);
      }
   }

   /* 2. fill distances, also from root to root! */

#ifdef SDQUERYFULL_HALFMATRIX
   for( int i = 0; i <= root; i++ )
   {
      if( nodes_isVisited[i] )
      {
         const int ancestor = SCIPStpunionfindFind(uf, i);
         const int pos = offset_root + i;

         assert(0 <= ancestor && ancestor < nlcanodes);

         SCIPdebugMessage("add %d->%d=%f \n", root, i, lcatree_costs[ancestor]);
         fullq_dists[pos] = lcatree_costs[ancestor];
      }
   }

   for( int i = root + 1; i < nlcanodes; i++ )
   {
      if( nodes_isVisited[i] )
      {
         const int ancestor = SCIPStpunionfindFind(uf, i);
         int pos = (i * (i + 1)) / 2 + root;

         assert(0 <= ancestor && ancestor < nlcanodes);

         SCIPdebugMessage("add %d->%d=%f \n", i, root, lcatree_costs[ancestor]);
         fullq_dists[pos] = lcatree_costs[ancestor];
      }
   }
#else
   for( int i = 0; i < nlcanodes; i++ )
   {
      if( nodes_isVisited[i] )
      {
         const int ancestor = SCIPStpunionfindFind(uf, i);

         assert(0 <= ancestor && ancestor < nlcanodes);

         SCIPdebugMessage("double add %d->%d=%f \n", i, root, lcatree_costs[ancestor]);
         fullq_dists[root * nlcanodes + i] = lcatree_costs[ancestor];
         fullq_dists[i * nlcanodes + root] = lcatree_costs[ancestor];
      }
   }
#endif
}


/** builds full representation */
static
SCIP_RETCODE sdqueryFullBuild(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH*              sdgraph,            /**< the SD graph */
   LCABUILDER*           lcabuilder          /**< the builder */
)
{
   UF uf;
   const BINARYNODE* const lcatree = lcabuilder->lcatree;
   const SCIP_Real* const lcatree_costs = lcabuilder->lcatree_costs;
   SCIP_Real* const fullq_dists = sdgraph->fullq_dists;
   SCIP_Bool* nodes_isVisited;
   const int nlcanodes = sdgraph->distgraph->knots - 1;

   assert(fullq_dists);
   assert(nlcanodes >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &(nodes_isVisited), nlcanodes) );

   SCIP_CALL( SCIPStpunionfindInit(scip, &uf, nlcanodes) );

   for( int i = 0; i < nlcanodes; i++ )
      nodes_isVisited[i] = FALSE;

#ifdef SDQUERYFULL_HALFMATRIX
   for( int i = 0; i < nlcanodes; i++ )
   {
      const int pos = (i * (i + 1)) / 2;
      sdgraph->fullq_IdxToStart[i] = pos;
   }
#endif

   /* NOTE: because we have sorted in ascending edge cost order, the root is 0 */
   sdqueryFullDfs(0, lcatree, nlcanodes, lcatree_costs, nodes_isVisited, &uf, fullq_dists);

   SCIPStpunionfindFreeMembers(scip, &uf);
   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}


/** initializes full SD constant time query data */
static
SCIP_RETCODE sdqueryFullInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   const int nterms = graph_get_nNodes(sdgraph->distgraph);
   const int fullq_dimension = (nterms - 1);
#ifdef SDQUERYFULL_HALFMATRIX
   const int fullq_size = (fullq_dimension * fullq_dimension) / 2 + fullq_dimension;
#else
   const int fullq_size = (fullq_dimension * fullq_dimension);
#endif
   LCABUILDER* lcabuilder;

   assert(!sdgraph->fullq_dists);
   assert(!sdgraph->usingRMQ);
   assert(nterms == g->terms);
   assert(g->terms >= 2 && g->terms <= 1000);
   SCIPdebugMessage("fullq_dimension=%d, fullq_size=%d \n", fullq_dimension, fullq_size);

   sdgraph->fullq_dimension = fullq_dimension;
   sdgraph->fullq_size = fullq_size;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->fullq_dists), fullq_size) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->fullq_nodeToIdx), g->knots) );
#ifdef SDQUERYFULL_HALFMATRIX
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdgraph->fullq_IdxToStart), g->knots) );
#endif

   SCIP_CALL( sdqueryLcaBuilderInit(scip, sdgraph, &lcabuilder) );
   SCIP_CALL( sdqueryBuildBinaryTree(scip, sdgraph, lcabuilder) );
   SCIP_CALL( sdqueryFullBuild(scip, sdgraph, lcabuilder) );
   sdqueryBuildNodesToFullMap(g, lcabuilder, sdgraph);

   sdqueryLcaBuilderFree(scip, &lcabuilder);

   return SCIP_OKAY;
}


/** initializes SD constant time query data */
static
SCIP_RETCODE sdqueryInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   if( g->terms > STP_SDQUERYFULL_MAXTERMS )
   {
      sdgraph->usingRMQ = TRUE;
      SCIP_CALL( sdqueryRmqInit(scip, g, sdgraph) );
   }
   else
   {
      sdgraph->usingRMQ = FALSE;
      SCIP_CALL( sdqueryFullInit(scip, g, sdgraph) );
   }

   return SCIP_OKAY;
}


/** frees SD constant time query data */
static
void sdqueryFree(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   if( sdgraph->usingRMQ )
   {
      assert(sdgraph->rmq_edgecosts && sdgraph->rmq_sparseTable && sdgraph->rmq_nodeToRmqEntry);
      assert(!sdgraph->fullq_dists && !sdgraph->fullq_nodeToIdx);

      SCIPfreeMemoryArray(scip, &(sdgraph->rmq_edgecosts));
      SCIPfreeMemoryArray(scip, &(sdgraph->rmq_sparseTable));
      SCIPfreeMemoryArray(scip, &(sdgraph->rmq_nodeToRmqEntry));
   }
   else
   {
      assert(sdgraph->fullq_dists);
      assert(!sdgraph->rmq_edgecosts && !sdgraph->rmq_sparseTable && !sdgraph->rmq_nodeToRmqEntry);
#ifdef SDQUERYFULL_HALFMATRIX
      SCIPfreeMemoryArray(scip, &(sdgraph->fullq_IdxToStart));
#endif
      SCIPfreeMemoryArray(scip, &(sdgraph->fullq_nodeToIdx));
      SCIPfreeMemoryArray(scip, &(sdgraph->fullq_dists));
   }
}


/** gets number of edges */
static
int distgraphGetNedges(
   const GRAPH*          g                   /**< graph to initialize from */
   )
{
   int maxnedges;
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   const SCIP_Longint terms2 = (SCIP_Longint) (nterms - 1) * nterms;

   if( nedges >= terms2 )
   {
      assert(terms2 <= INT_MAX);
      maxnedges = terms2;
   }
   else
   {
      maxnedges = nedges;
   }

   return maxnedges;
}


/** adds nodes to distance graph */
static
void distgraphAddNodes(
   const GRAPH*          g,                  /**< graph to initialize from */
   int* RESTRICT         distnodes_id,       /**< IDs of nodes */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nnodes_distgraph = 0;

   /* add the nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         graph_knot_add(distgraph, STP_TERM_NONE);
         distnodes_id[k] = nnodes_distgraph++;
      }
      else
      {
         distnodes_id[k] = UNKNOWN;
      }
   }

   assert(distgraph->knots == nnodes_distgraph);
   assert(distgraph->knots == g->terms);

   graph_knot_chg(distgraph, 0, STP_TERM);
   distgraph->source = 0;
}


/** gets SD for boundary edge */
static inline
SCIP_Real distgraphGetBoundaryEdgeDist(
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   int                   vbase_tail,         /**< base */
   int                   vbase_head,         /**< base */
   SCIP_Real             edgecost,           /**< cost */
   const SCIP_Real*      nodes_vdist,        /**< distance */
   const SDPROFIT*       sdprofit            /**< profit */
   )
{
   const SCIP_Real profit_tail = reduce_sdprofitGetProfit(sdprofit, tail, vbase_head, vbase_tail);
   const SCIP_Real profit_head = reduce_sdprofitGetProfit(sdprofit, head, vbase_head, vbase_tail);
   SCIP_Real distSimple = edgecost + nodes_vdist[tail] + nodes_vdist[head];
   const SCIP_Real distAll = distSimple - profit_tail - profit_head;
   const SCIP_Real distTailHead = edgecost + nodes_vdist[tail] - profit_tail;
   const SCIP_Real distHeadTail = edgecost + nodes_vdist[head] - profit_head;

   const SCIP_Real dist = miscstp_maxReal((SCIP_Real[])
               { distAll, distTailHead, distHeadTail, edgecost,
                 nodes_vdist[tail], nodes_vdist[head] },
                 6);
   return dist;
}


/** gets SD for boundary edge */
static inline
SCIP_Real distgraphGetBoundaryEdgeDist2(
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   int                   vbase_tail,         /**< base */
   int                   vbase_head,         /**< base */
   SCIP_Real             edgecost,           /**< cost */
   SCIP_Real             dist_tail,          /**< distance */
   SCIP_Real             dist_head,          /**< distance */
   const SDPROFIT*       sdprofit            /**< profit */
   )
{
   const SCIP_Real profit_tail = reduce_sdprofitGetProfit(sdprofit, tail, vbase_head, vbase_tail);
   const SCIP_Real profit_head = reduce_sdprofitGetProfit(sdprofit, head, vbase_head, vbase_tail);
   SCIP_Real distSimple = edgecost + dist_tail + dist_head;
   const SCIP_Real distAll = distSimple - profit_tail - profit_head;
   const SCIP_Real distTailHead = edgecost + dist_tail - profit_tail;
   const SCIP_Real distHeadTail = edgecost + dist_head - profit_head;

   const SCIP_Real dist = miscstp_maxReal((SCIP_Real[])
               { distAll, distTailHead, distHeadTail, edgecost,
                   dist_tail, dist_head},
                 6);
   return dist;
}


/** gets SD for boundary edge ... choose among nearest terminals (w.r.t. implied SD) */
static inline
SCIP_Real distgraphGetBoundaryEdgeDistBest(
   const GRAPH*          g,                  /**< graph  */
   const TPATHS*         tpaths,             /**< terminal paths */
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   SCIP_Real             edgecost,           /**< cost */
   const SDPROFIT*       sdprofit,           /**< profit */
   int*                  base_tail,          /**< base of tail */
   int*                  base_head           /**< base of head */
   )
{
   SCIP_Real dists_tail[3];
   SCIP_Real dists_head[3];
   int terms_tail[3];
   int terms_head[3];
   int nterms_tail;
   int nterms_head;
   SCIP_Real dist = FARAWAY;

   assert(g && tpaths && sdprofit);

   *base_tail = -1;
   *base_head = -1;

   graph_tpathsGet3CloseTerms(g, tpaths, tail, FARAWAY, terms_tail, NULL, dists_tail, &nterms_tail);
   graph_tpathsGet3CloseTerms(g, tpaths, head, FARAWAY, terms_head, NULL, dists_head, &nterms_head);

   for( int i = 0; i < nterms_tail; ++i )
   {
      for( int j = 0; j < nterms_head; ++j )
      {
         if( terms_tail[i] != terms_head[j] )
         {
            const SCIP_Real distnew =
            distgraphGetBoundaryEdgeDist2(tail, head, terms_tail[i], terms_head[j], edgecost, dists_tail[i], dists_head[j], sdprofit);

            if( distnew < dist )
            {
               *base_tail = terms_tail[i];
               *base_head = terms_head[j];
               dist = distnew;
            }
         }
      }
   }

   return dist;
}


/** inserts new edge */
static inline
void distgraphInsertEdge(
   SCIP*                  scip,               /**< SCIP data structure */
   int                    sdnode1,           /**< end node 1 */
   int                    sdnode2,           /**< end node 2 */
   SCIP_Real              edgecost,          /**< cost */
   int                    edgeid,            /**< ID or -1 */
   int* RESTRICT          edgeorg,           /**< IDs of edges or NULL */
   GRAPH*                 distgraph,         /**< the SD graph */
   SCIP_Bool*             success            /**< could the edge be added? */
)
{
   int ne;
   *success = TRUE;

#ifndef NDEBUG
   assert(sdnode1 != sdnode2);
   assert(0 <= sdnode1 && sdnode1 < distgraph->knots);
   assert(0 <= sdnode2 && sdnode2 < distgraph->knots);
   assert(GE(edgecost, 0.0));
   assert(edgeorg != NULL || edgeid == -1);
#endif

   /* find the corresponding edge in the distance graph */
   for( ne = distgraph->outbeg[sdnode1]; ne != EAT_LAST; ne = distgraph->oeat[ne] )
   {
      if( distgraph->head[ne] == sdnode2 )
         break;
   }

   /* edge exists already? */
   if( ne != EAT_LAST )
   {
      assert(ne >= 0);
      assert(distgraph->head[ne] == sdnode2);
      assert(distgraph->tail[ne] == sdnode1);

      if( distgraph->cost[ne] > edgecost )
      {
         distgraph->cost[ne]            = edgecost;
         distgraph->cost[Edge_anti(ne)] = edgecost;

         if( edgeorg != NULL )
            edgeorg[ne / 2] = edgeid;
      }
   }
   else
   {
      assert(distgraph->edges <= distgraph->esize);

      if( distgraph->edges == distgraph->esize  )
      {
         *success = FALSE;
      }
      else
      {
         if( edgeorg != NULL )
            edgeorg[distgraph->edges / 2] = edgeid;

         graph_edge_add(scip, distgraph, sdnode1, sdnode2, edgecost, edgecost);
      }
   }
}



/** adds edges to distance graph, given a Voronoi diagram */
static
void distgraphAddEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph to initialize from */
   const int*            distnodes_id,       /**< IDs of nodes */
   const VNOI*           vnoi,               /**< Voronoi */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   int* RESTRICT         edgeorg,            /**< IDs of edges */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const SCIP_Real* nodes_vdist = vnoi->nodes_dist;
   const int* RESTRICT nodes_vbase = vnoi->nodes_base;
   const SCIP_Bool useProfit = (sdprofit != NULL);

   for( int e = 0; e < nedges / 2; e++ )
      edgeorg[e] = UNKNOWN;

   /* add the edges */
   for( int tail = 0; tail < nnodes; tail++ )
   {
      const int vbase_tail = nodes_vbase[tail];

      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         const int vbase_head = nodes_vbase[head];

         assert(tail == g->tail[e]);

         if( vbase_tail != vbase_head )
         {
            SCIP_Bool success;
            const SCIP_Real distance = useProfit ?
               distgraphGetBoundaryEdgeDist(tail, head, vbase_tail, vbase_head, g->cost[e], nodes_vdist, sdprofit)
               : (g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

           // if( LT(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]) )
         //   printf("distance: %f < %f \n", distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

            assert(LE(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]));
            assert(Is_term(g->term[vbase_tail]));
            assert(Is_term(g->term[vbase_head]));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, e,
                  edgeorg, distgraph, &success);

            assert(success);
         }
      }
   }
}



/** helper */
static
void sdgraphSetDefaults(
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   SCIP_Real* mstsdist = g_sd->mstsdist;
   const int nedges = graph_get_nEdges(g);

   for( int i = 0; i < nnodes; i++ )
   {
      mstsdist[i] = -1.0;
   }

   g_sd->nnodesorg = nnodes;
   g_sd->nedgesorg = nedges;
   g_sd->mstcostsReady = FALSE;
   g_sd->edgemarkReady = TRUE;
   g_sd->usingRMQ = TRUE;
   g_sd->fullq_dists = NULL;
   g_sd->fullq_nodeToIdx = NULL;
   g_sd->fullq_size = -1;
   g_sd->fullq_dimension = -1;
   g_sd->rmq_edgecosts = NULL;
   g_sd->rmq_sparseTable = NULL;
   g_sd->rmq_nodeToRmqEntry = NULL;
   g_sd->rmq_loglength = -1;
}


/** allocates memory */
static
SCIP_RETCODE sdgraphAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   SDGRAPH* g_sd;

   SCIP_CALL( SCIPallocMemory(scip, sdgraph) );
   g_sd = *sdgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->sdmst), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstcosts), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstsdist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->nodemapOrgToDist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->halfedge_isInMst), nedges / 2) );

   sdgraphSetDefaults(g, g_sd);

   return SCIP_OKAY;
}



/** allocates memory */
static
SCIP_RETCODE sdgraphAllocRestricted(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nterms = g->terms;
   SDGRAPH* g_sd;

   SCIP_CALL( SCIPallocMemory(scip, sdgraph) );
   g_sd = *sdgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->sdmst), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstcosts), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstsdist), nnodes) );

   g_sd->nodemapOrgToDist = NULL;
   g_sd->halfedge_isInMst = NULL;

   sdgraphSetDefaults(g, g_sd);

   return SCIP_OKAY;
}

/** adds edges to distance graph, given terminal paths */
static
void distgraphAddEdgesFromTpaths(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph to initialize from */
   const int*            distnodes_id,       /**< IDs of nodes */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool useProfit = (sdprofit != NULL);

   /* add the edges */
   for( int tail = 0; tail < nnodes; tail++ )
   {
      SCIP_Real dist_tail;
      int vbase_tail;
      graph_tpathsGetClosestTerm(g, tpaths, tail, &vbase_tail, NULL, &dist_tail);

      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         SCIP_Real dist_head;
         int vbase_head;

         graph_tpathsGetClosestTerm(g, tpaths, head, &vbase_head, NULL, &dist_head);

         if( vbase_tail != vbase_head )
         {
            SCIP_Bool success;
            const SCIP_Real distance = useProfit ?
               distgraphGetBoundaryEdgeDist2(tail, head, vbase_tail, vbase_head, g->cost[e], dist_tail, dist_head, sdprofit)
               : (g->cost[e] + dist_tail + dist_head);

            assert(LE(distance, g->cost[e] + dist_tail + dist_head));
            assert(Is_term(g->term[vbase_tail]) && Is_term(g->term[vbase_head]));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, -1, NULL, distgraph, &success);
            assert(success);
         }
      }
   }
}


/** builds distance graph */
static
SCIP_RETCODE sdgraphBuildDistgraph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   SDGRAPH*              g_sd,               /**< the SD graph */
   VNOI**                vnoi,               /**< Voronoi */
   int**                 distedge2org        /**< array of size nedges / 2 */
)
{
   GRAPH* distgraph;
   int* RESTRICT distnodes_id = g_sd->nodemapOrgToDist;
   const int nedges = graph_get_nEdges(g);
   const int nedges_distgraph = distgraphGetNedges(g);

   assert(g->terms >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, distedge2org, nedges / 2) );

   /* build biased Voronoi diagram */
   SCIP_CALL( graph_vnoiInit(scip, g, TRUE, vnoi) );

   if( sdprofit )
      graph_vnoiComputeImplied(scip, g, sdprofit, *vnoi);
   else
      graph_vnoiCompute(scip, g, *vnoi);

   /* build distance graph from Voronoi diagram */
   SCIP_CALL( graph_init(scip, &(g_sd->distgraph), g->terms, nedges_distgraph, 1) );
   distgraph = g_sd->distgraph;
   distgraphAddNodes(g, distnodes_id, distgraph);
   distgraphAddEdges(scip, g, distnodes_id, *vnoi, sdprofit, *distedge2org, distgraph);

   assert(graph_valid(scip, distgraph));

   return SCIP_OKAY;
}


/** builds distance graph */
static
SCIP_RETCODE sdgraphBuildDistgraphFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   GRAPH* distgraph;
   int* RESTRICT distnodes_id = g_sd->nodemapOrgToDist;
   const int nedges_distgraph = distgraphGetNedges(g);

   assert(g->terms >= 1);

   SCIP_CALL( graph_init(scip, &(g_sd->distgraph), g->terms, nedges_distgraph, 1) );
   distgraph = g_sd->distgraph;
   distgraphAddNodes(g, distnodes_id, distgraph);
   distgraphAddEdgesFromTpaths(scip, g, distnodes_id, sdprofit, tpaths, distgraph);

   assert(graph_valid(scip, distgraph));

   return SCIP_OKAY;
}


/** updates distance graph */
static
SCIP_RETCODE sdgraphUpdateDistgraphFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   int* hasharr;
   STP_Vectype(int) profitnodes_tail = NULL;
   STP_Vectype(int) profitnodes_head = NULL;

   GRAPH* RESTRICT distgraph = g_sd->distgraph;
   const int* distnodes_id = g_sd->nodemapOrgToDist;
   const int nnodes = graph_get_nNodes(g);
   const int edgelimit = MIN(2 * distgraph->edges, distgraph->esize);

   assert(LE(distgraph->edges, distgraph->esize));

   StpVecReserve(scip, profitnodes_tail, nnodes);
   StpVecReserve(scip, profitnodes_head, nnodes);
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, nnodes) );

   for( int tail = 0; tail < nnodes && distgraph->edges < edgelimit; tail++ )
   {
      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         int vbase_tail;
         int vbase_head;
         SCIP_Real distance;

         if( head < tail )
            continue;

         distance = distgraphGetBoundaryEdgeDistBest(g, tpaths, tail, head, g->cost[e], sdprofit, &vbase_tail, &vbase_head);

         /* NOTE: we should not take the fast query method here, because sdgraphGetSd at least partly reflects
          * the change of the distance graph */
         if( LT(distance, FARAWAY) && LT(distance, sdgraphGetSd(vbase_tail, vbase_head, g_sd)) )
         {
            SCIP_Bool success = TRUE;
            assert(vbase_tail >= 0 && vbase_head >= 0);

            graph_tpathsGetProfitNodes(scip, g, tpaths, sdprofit, tail, vbase_tail, profitnodes_tail);
            for( int k = 0; k < StpVecGetSize(profitnodes_tail); k++ )
            {
               assert(!hasharr[profitnodes_tail[k]]);
               hasharr[profitnodes_tail[k]] = 1;
            }
            graph_tpathsGetProfitNodes(scip, g, tpaths, sdprofit, head, vbase_head, profitnodes_head);

            for( int k = 0; k < StpVecGetSize(profitnodes_head); k++ )
            {
               if( hasharr[profitnodes_head[k]] )
               {
                  success = FALSE;
                //  printf("shared node: %d \n", profitnodes_head[k]);
                  break;
               }
            }

            for( int k = 0; k < StpVecGetSize(profitnodes_tail); k++ )
            {
               assert(1 == hasharr[profitnodes_tail[k]]);
               hasharr[profitnodes_tail[k]] = 0;
            }

            if( !success )
            {
              // printf("fail for %d %d \n", vbase_tail, vbase_head);
               continue;
            }
            //SCIPdebugMessage("add biased MST edge %d->%d (%f<%f) \n", vbase_tail, vbase_head, distance, sdgraphGetSd(vbase_tail, vbase_head, g_sd));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, -1, NULL, distgraph, &success);
            assert(success);

            if( distgraph->edges >= edgelimit )
               break;
         }
      }
   }
#ifndef NDEBUG
   for( int i = 0; i < nnodes; i++ )
      assert(hasharr[i] == 0);
#endif

   SCIPfreeCleanBufferArray(scip, &hasharr);
   StpVecFree(scip, profitnodes_head);
   StpVecFree(scip, profitnodes_tail);

   return SCIP_OKAY;
}


/** builds MST costs (ordered) for distance graph */
static
void sdgraphMstSortCosts(
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   const GRAPH* const distgraph = g_sd->distgraph;
   const PATH* const mst = g_sd->sdmst;
   SCIP_Real* RESTRICT mstcosts = g_sd->mstcosts;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);

   assert(mst[0].edge == UNKNOWN);
   assert(nnodes_distgraph >= 1);

   for( int k = 1; k < nnodes_distgraph; k++ )
   {
#ifndef NDEBUG
      const int e = mst[k].edge;
      assert(e >= 0);
      assert(GE(distgraph->cost[e], 0.0));
      assert(EQ(mst[k].dist, distgraph->cost[e]));
#endif

      mstcosts[k - 1] = mst[k].dist;
   }

   SCIPsortDownReal(mstcosts, nnodes_distgraph - 1);

   /* debug sentinel */
   if( nnodes_distgraph > 1 )
   {
      mstcosts[nnodes_distgraph - 1] = -FARAWAY;
   }
   else
   {
      mstcosts[0] = 0.0;
   }
}


/** marks original edges corresponding to MST */
static
void sdgraphMstMarkOrgEdges(
   const GRAPH*          g,                  /**< graph to initialize from */
   const VNOI*           vnoi,               /**< Voronoi */
   const int*            distedge2org,       /**< array of size nedges / 2 */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   GRAPH* distgraph = g_sd->distgraph;
   PATH* RESTRICT mst = g_sd->sdmst;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);
   const int* const nodes_vbase = vnoi->nodes_base;
   const int* const nodes_vpred = vnoi->nodes_predEdge;
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;

#ifndef NDEBUG
   const int nedges = graph_get_nEdges(g);

   for( int e = 0; e < nedges / 2; e++ )
      assert(!orgedges_isInMst[e]);
#endif

   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const int ne = distedge2org[e / 2];

      assert(e >= 0);
      assert(e < nedges);

      orgedges_isInMst[ne / 2] = TRUE;

      for( int v = g->head[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v] / 2] = TRUE;
      }

      for( int v = g->tail[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v]/ 2] = TRUE;
      }

      assert(e != EAT_LAST);
   }
}

/** builds MST on distance graph */
static
SCIP_RETCODE sdgraphMstBuild(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   PATH* RESTRICT mst = g_sd->sdmst;
   GRAPH* distgraph = g_sd->distgraph;
   SCIP_Real maxcost;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;
   const int nedges = graph_get_nEdges(g);
   const SCIP_Bool distgraphIsInit = (distgraph->path_heap != NULL);

   if( orgedges_isInMst)
   {
      for( int e = 0; e < nedges / 2; e++ )
         orgedges_isInMst[e] = FALSE;
   }

   for( int k = 0; k < nnodes_distgraph; k++ )
      distgraph->mark[k] = TRUE;

   if( !distgraphIsInit )
   {
      SCIP_CALL( graph_path_init(scip, distgraph) );
   }
   graph_path_exec(scip, distgraph, MST_MODE, distgraph->source, distgraph->cost, mst);

   if( !distgraphIsInit )
   {
      graph_path_exit(scip, distgraph);
   }

   assert(mst[0].edge == -1);

   maxcost = 0.0;
   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const SCIP_Real cost = distgraph->cost[e];

      assert(graph_edge_isInRange(distgraph, e));

      if( cost > maxcost )
         maxcost = cost;
   }

   g_sd->mstmaxcost = maxcost;

   return SCIP_OKAY;
}


/** finalizes distance graph */
static
void sdgraphFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   VNOI**                vnoi,               /**< Voronoi data structure */
   int**                 edgeorg             /**< array of size nedges / 2 */
)
{
   graph_vnoiFree(scip, vnoi);
   SCIPfreeBufferArray(scip, edgeorg);
}


/*
 * Interface methods
 */


/** initializes SD graph */
SCIP_RETCODE reduce_sdgraphInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   VNOI* vnoi;
   int* edgeorg;
   assert(scip && g && sdgraph);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraph(scip, g, NULL, *sdgraph, &vnoi, &edgeorg) );
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   sdgraphMstMarkOrgEdges(g, vnoi, edgeorg, *sdgraph);

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   SCIP_CALL( sdqueryInit(scip, g, *sdgraph) );

   return SCIP_OKAY;
}


/** initializes SD graph */
SCIP_RETCODE reduce_sdgraphInitFromDistGraph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   GRAPH*                distgraph,          /**< distance graph */
   int*                  node2dist,          /**< map */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   assert(scip && g && sdgraph && distgraph && node2dist);

   SCIP_CALL( sdgraphAllocRestricted(scip, g, sdgraph) );
   (*sdgraph)->distgraph = distgraph;
   (*sdgraph)->nodemapOrgToDist = node2dist;

   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );

   SCIP_CALL( sdqueryInit(scip, g, *sdgraph) );

   return SCIP_OKAY;
}


/** initializes biased SD graph */
SCIP_RETCODE reduce_sdgraphInitBiased(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   VNOI* vnoi;
   int* edgeorg;
   assert(scip && g && sdgraph);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraph(scip, g, sdprofit, *sdgraph, &vnoi, &edgeorg) );
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   sdgraphMstMarkOrgEdges(g, vnoi, edgeorg, *sdgraph);

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   SCIP_CALL( sdqueryInit(scip, g, *sdgraph) );

   return SCIP_OKAY;
}


/** initializes biased SD graph from given terminal paths */
SCIP_RETCODE reduce_sdgraphInitBiasedFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   assert(scip && g && sdgraph && tpaths);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraphFromTpaths(scip, g, sdprofit, tpaths, *sdgraph) );
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );

   if( sdprofit )
   {
      SCIP_CALL( sdgraphUpdateDistgraphFromTpaths(scip, g, sdprofit, tpaths, *sdgraph) );
      SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   }

   /* NOTE probably we never need that...for extending reductions we anyway should only take biased paths */
   (*sdgraph)->edgemarkReady = FALSE;
   SCIPfreeMemoryArray(scip, &(*sdgraph)->halfedge_isInMst);

   SCIP_CALL( sdqueryInit(scip, g, *sdgraph) );

   return SCIP_OKAY;
}


/** return maximum MST edge cost */
SCIP_Real reduce_sdgraphGetMaxCost(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(GE(sdgraph->mstmaxcost, 0.0));

   return sdgraph->mstmaxcost;
}


/** returns edge mark */
const STP_Bool* reduce_sdgraphGetMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasMstHalfMark(sdgraph));

   return sdgraph->halfedge_isInMst;
}


/** has edge mark? */
SCIP_Bool reduce_sdgraphHasMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   if( sdgraph->edgemarkReady )
   {
      assert(sdgraph->halfedge_isInMst);
      return TRUE;
   }

   assert(!sdgraph->halfedge_isInMst);

   return FALSE;
}


/** is edge in current SD MST? */
SCIP_Bool reduce_sdgraphEdgeIsInMst(
   const SDGRAPH*        sdgraph,            /**< the SD graph */
   int                   edge                /**< edge */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasMstHalfMark(sdgraph));
   assert(sdgraph->halfedge_isInMst);
   assert(edge >= 0);

   return sdgraph->halfedge_isInMst[edge / 2];
}


/** MST costs in descending order available? */
SCIP_Bool reduce_sdgraphHasOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   return sdgraph->mstcostsReady;
}




/** initializes all MST costs in descending order */
void reduce_sdgraphInitOrderedMstCosts(
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   if( !sdgraph->mstcostsReady )
   {
      sdgraphMstSortCosts(sdgraph);
      sdgraph->mstcostsReady = TRUE;
   }

   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);
}


/** Gets special distance (e.g. bottleneck distance) from distance graph.
 *  Only works if both nodes are terminals!  */
SCIP_Real reduce_sdgraphGetSd(
   int                    term1,             /**< node 1 */
   int                    term2,             /**< node 2 */
   SDGRAPH*               sdgraph            /**< the SD graph */
)
{
#ifndef NDEBUG
   assert(sdgraph);
   assert(term1 != term2);
   assert(sdgraph->mstcosts);
   assert(sdgraph->mstsdist);
   assert(0 <= term1 && term1 < sdgraph->nnodesorg);
   assert(0 <= term2 && term2 < sdgraph->nnodesorg);
   assert(sdgraph->nodemapOrgToDist[term1] != UNKNOWN);
   assert(sdgraph->nodemapOrgToDist[term2] != UNKNOWN);

   for( int i = 0; i < sdgraph->nnodesorg; i++ )
      assert(EQ(sdgraph->mstsdist[i], -1.0));
#endif

   if( sdgraph->usingRMQ )
   {
      return sdqueryGetSd(term1, term2, sdgraph);
   }
   else
   {
      return sdqueryFullGetSd(term1, term2, sdgraph);
   }
}

/** returns all MST costs in descending order */
const SCIP_Real* reduce_sdgraphGetOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);

   return sdgraph->mstcosts;
}

/** returns mapping from original nodes to distance graph nodes */
const int* reduce_sdgraphGetOrgnodesToSdMap(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);

   return sdgraph->nodemapOrgToDist;
}




/** Inserts new edge.
 *  NOTE: just a wrapper, should only be used by other reduce_sd* methods */
void reduce_sdgraphInsertEdge(
   SCIP*                  scip,              /**< SCIP data structure */
   int                    sdnode1,           /**< end node 1 */
   int                    sdnode2,           /**< end node 2 */
   SCIP_Real              edgecost,          /**< cost */
   int                    edgeid,            /**< ID or -1 */
   int* RESTRICT          edgeorg,           /**< IDs of edges or NULL */
   SDGRAPH*               sdgraph,           /**< the SD graph */
   SCIP_Bool*             success            /**< could the edge be added? */
)
{
   distgraphInsertEdge(scip, sdnode1, sdnode2, edgecost, edgeid, edgeorg, sdgraph->distgraph, success);
}


/** Builds MST on distance graph.
 *  NOTE: just a wrapper, should only be used by other reduce_sd* methods */
SCIP_RETCODE reduce_sdgraphMstBuild(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   SCIP_CALL( sdgraphMstBuild(scip, g, g_sd) );

   return SCIP_OKAY;
}


/** Builds MST costs (ordered) for distance graph
*  NOTE: just a wrapper, should only be used by other reduce_sd* methods */
void reduce_sdgraphMstSortCosts(
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   sdgraphMstSortCosts(g_sd);
}


/** frees SD graph */
void reduce_sdgraphFree(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   SDGRAPH* g_sd;
   assert(scip && sdgraph);

   g_sd = *sdgraph;
   assert(g_sd);

   SCIPfreeMemoryArrayNull(scip, &(g_sd->halfedge_isInMst));
   SCIPfreeMemoryArray(scip, &(g_sd->nodemapOrgToDist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstsdist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstcosts));
   SCIPfreeMemoryArray(scip, &(g_sd->sdmst));

   graph_free(scip, &(g_sd->distgraph), TRUE);
   sdqueryFree(scip, g_sd);

   SCIPfreeMemory(scip, sdgraph);
}


/** frees SD graph, but does not free actual graph and node-map (assumed to be non-owned) */
void reduce_sdgraphFreeFromDistGraph(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   SDGRAPH* g_sd;
   assert(scip && sdgraph);

   g_sd = *sdgraph;
   assert(g_sd);
   assert(!(g_sd->halfedge_isInMst));

   SCIPfreeMemoryArray(scip, &(g_sd->mstsdist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstcosts));
   SCIPfreeMemoryArray(scip, &(g_sd->sdmst));

   sdqueryFree(scip, g_sd);

   SCIPfreeMemory(scip, sdgraph);
}
