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

/**@file   bidecomposition.c
 * @brief  several decomposition methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements bi-connected components decomposition methods for several Steiner problems.
 *
 * A list of all interface methods can be found in bidecomposition.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bidecomposition.h"
#include "stpvector.h"
#include "portab.h"

/*
 * Data structures
 */

struct biconnected_stack_node {
   int                   id;                 /**< node ID in the underlying graph */
   int                   parent;             /**< parent in the underlying graph */
   int                   nextEdge;           /**< next edge in the underlying graph */
   int                   biconnStart;        /**< start of current component in separate stack */
   int                   lowpoint;           /**< low-point of this node */
   SCIP_Bool             isArtPoint;         /**< is articulation point? */
};


/*
 * Local methods
 */


/** sets root  */
static inline
void cutNodesSetDfsRoot(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   if( !graph_pc_isPcMw(g) )
   {
      cutnodes->dfsroot = g->source;
   }
   else
   {
      int i;
      const int nnodes = graph_get_nNodes(g);
      const int pseudoroot = graph_pc_isRootedPcMw(g) ? -1 : g->source;

      assert(!g->extended);

      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(g->term[i]) && i != pseudoroot)
         {
            assert(!graph_pc_knotIsDummyTerm(g, i));
            cutnodes->dfsroot = i;
            break;
         }
      }

      assert(i < nnodes);
   }
}


/** processes bi-connected component */
static
void cutNodesProcessComponent(
   int                   root,               /**< root of component */
   int                   stack_start,        /**< start of component in the stack */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   const int ncomps = ++(cutnodes->biconn_ncomps);
   int* biconn_nodesmark = cutnodes->biconn_nodesmark;
   int* biconn_comproot = cutnodes->biconn_comproots;
   STP_Vectype(int) biconn_stack = cutnodes->biconn_stack;

   assert(stack_start >= 0);
   assert(StpVecGetSize(biconn_stack) > stack_start);
   assert(ncomps > 0);
   assert(biconn_comproot[ncomps] == -1);
   assert(ncomps == 1 || biconn_comproot[ncomps - 1] != -1);

   biconn_comproot[ncomps] = root;

   for( int size = StpVecGetSize(biconn_stack); size != stack_start; size-- )
   {
      const int compnode = biconn_stack[size - 1];
      SCIPdebugMessage("%d \n", compnode);

      assert(size >= 1);
      assert(size == StpVecGetSize(biconn_stack));
      assert(biconn_nodesmark[compnode] == 0);
      assert(compnode != root);

      biconn_nodesmark[compnode] = ncomps;
      StpVecPopBack(biconn_stack);
   }
}


//#define USE_RECURSIVE_DFS
#ifdef USE_RECURSIVE_DFS
/** recursive DFS */
static
void cutNodesComputeRecursive(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< vertex */
   int                   parent,             /**< parent of node */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
#ifdef CUTTREE_PRINT_STATISTICS
   int* childcount_terms = cutnodes->childcount_terms;
   int* childcount_nodes = cutnodes->childcount_nodes;
#endif

   int* const nodes_hittime = cutnodes->nodes_hittime;
   const int myhittime = cutnodes->curr_hittime;
   int mylowpoint = myhittime;
   int nchildren = 0;
   SCIP_Bool isCutNode = FALSE;
   const SCIP_Bool nodeIsRoot = (parent == -1);

   nodes_hittime[node] = myhittime;
   (cutnodes->curr_hittime)++;
   StpVecPushBack(cutnodes->scip, cutnodes->biconn_stack, node);

   for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int head = g->head[e];
      if( !g->mark[head] )
         continue;

      /* not visited? */
      if( nodes_hittime[head] < 0 )
      {
         const int stack_start = StpVecGetSize(cutnodes->biconn_stack);
         assert(head != parent);

         cutNodesComputeRecursive(g, head, node, cutnodes);

#ifdef CUTTREE_PRINT_STATISTICS
         childcount_nodes[node] += childcount_nodes[head] + 1;
         childcount_terms[node] += childcount_terms[head] + (Is_term(g->term[head]) ? 1 : 0);
#endif
         if( nodeIsRoot )
         {
            nchildren++;

            if( nchildren > 1 )
            {
               SCIPdebugMessage("mark new bi-connected component from root \n");
               cutNodesProcessComponent(node, stack_start, cutnodes);
            }
         }
         else if( cutnodes->curr_lowpoint >= myhittime )
         {
            isCutNode = TRUE;

            SCIPdebugMessage("mark new bi-connected component from %d \n", node);
            cutNodesProcessComponent(node, stack_start, cutnodes);
         }

         if( mylowpoint > cutnodes->curr_lowpoint )
            mylowpoint = cutnodes->curr_lowpoint;
      }
      else if( head != parent )
      {
         assert(nodes_hittime[head] >= 0);
         if( mylowpoint > nodes_hittime[head] )
            mylowpoint = nodes_hittime[head];
      }
   }

   if( nodeIsRoot && nchildren > 1 )
   {
      assert(!isCutNode);
      isCutNode = TRUE;
      SCIPdebugMessage("found parent cut node: %d \n", node);
   }

   if( isCutNode )
   {
      StpVecPushBack(cutnodes->scip, cutnodes->artpoints, node);
   }

   cutnodes->curr_lowpoint = mylowpoint;
}
#else
/** non-recursive DFS-based biconnected components helper */
static
void cutNodesProcessNext(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   const int stack_top = cutnodes->stack_size - 1;
   STACK_NODE* stack_top_data = &(cutnodes->stack_nodes[stack_top]);
   const int node = stack_top_data->id;
   int* const nodes_hittime = cutnodes->nodes_hittime;
   int e;
   SCIP_Bool childWasAdded = FALSE;
   const SCIP_Bool isFirstNodeVisit = (stack_top_data->nextEdge == g->outbeg[node]);
   SCIPdebugMessage("processing node %d, global hit-time=%d \n", node, cutnodes->curr_hittime);

   if( isFirstNodeVisit )
   {
      assert(StpVecGetSize(cutnodes->biconn_stack) < StpVecGetcapacity(cutnodes->biconn_stack));
      assert(stack_top_data->lowpoint == -1);

      nodes_hittime[node] = (cutnodes->curr_hittime)++;
      stack_top_data->lowpoint = nodes_hittime[node];
      StpVecPushBack(cutnodes->scip, cutnodes->biconn_stack, node);
      assert(!stack_top_data->isArtPoint);
   }

   for( e = stack_top_data->nextEdge; e != EAT_LAST; e = g->oeat[e] )
   {
      const int head = g->head[e];
      if( !g->mark[head] )
         continue;

      /* not visited? */
      if( nodes_hittime[head] < 0 )
      {
         assert(head != stack_top_data->parent);
         assert(cutnodes->stack_size < g->knots);

         /* add child node to stack */
         cutnodes->stack_nodes[cutnodes->stack_size++] =
               (STACK_NODE) { .id = head, .parent = node, .nextEdge = g->outbeg[head],
                              .biconnStart = -1, .lowpoint = -1, .isArtPoint = FALSE };
         childWasAdded = TRUE;
         break;
      }
      else if( head != stack_top_data->parent )
      {
         assert(nodes_hittime[head] >= 0);
         if( stack_top_data->lowpoint > nodes_hittime[head] ) {
            SCIPdebugMessage("update my lowpoint to %d (from %d) \n", nodes_hittime[head], head);
            stack_top_data->lowpoint = nodes_hittime[head];
         }
      }
   }
   stack_top_data->nextEdge = (e != EAT_LAST) ? g->oeat[e] : EAT_LAST;

   /* are we finished with this node? */
   if( !childWasAdded )
   {
      cutnodes->stack_size--;

      if( cutnodes->stack_size > 0 ) {
         /* update hit-time of parent */
         const int parent_pos = cutnodes->stack_size - 1;
         if( cutnodes->stack_nodes[parent_pos].lowpoint > stack_top_data->lowpoint ) {
            SCIPdebugMessage("update parent(%d) current lowpoint to %d \n", cutnodes->stack_nodes[parent_pos].id, stack_top_data->lowpoint);
            cutnodes->stack_nodes[parent_pos].lowpoint = stack_top_data->lowpoint;
         }
      }
      SCIPdebugMessage("finished node %d \n", node);
   }

   if( !isFirstNodeVisit ) {
      const SCIP_Bool nodeIsRoot = (node == cutnodes->dfsroot);
      assert(nodeIsRoot == (stack_top_data->parent == -1));
      SCIPdebugMessage("return to node %d, cutnodes->curr_lowpoint=%d \n", node, cutnodes->curr_lowpoint);

      if( nodeIsRoot )
      {
         assert(cutnodes->nrootcomps >= 0);
         cutnodes->nrootcomps++;

         /* is root an articulation point? */
         if( cutnodes->nrootcomps > 1 )
         {
            SCIPdebugMessage("mark new bi-connected component from root \n");
            cutNodesProcessComponent(node, stack_top_data->biconnStart, cutnodes);

            /* have we detected for the first time that root is an articulation point? */
            if( cutnodes->nrootcomps == 2 ) {
               assert(!stack_top_data->isArtPoint);
               stack_top_data->isArtPoint = TRUE;
            }
         }
      }
      else if( cutnodes->curr_lowpoint >= nodes_hittime[node]  )
      {
         SCIPdebugMessage("mark new bi-connected component from %d \n", node);
         cutNodesProcessComponent(node, stack_top_data->biconnStart, cutnodes);
         stack_top_data->isArtPoint = TRUE;
      }
   }

   if( !childWasAdded && stack_top_data->isArtPoint ) {
      StpVecPushBack(cutnodes->scip, cutnodes->artpoints, node);
   }

   stack_top_data->biconnStart = StpVecGetSize(cutnodes->biconn_stack);
   cutnodes->curr_lowpoint = stack_top_data->lowpoint;
}


/** non-recursive DFS */
static
void cutNodesCompute(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   const int root = cutnodes->dfsroot;
   assert(graph_knot_isInRange(g, root));
   cutnodes->stack_nodes[0] = (STACK_NODE) { .id = root, .parent = -1, .nextEdge = g->outbeg[root], .biconnStart = -1, .lowpoint = -1, .isArtPoint = FALSE };
   cutnodes->stack_size = 1;
   SCIPdebugMessage("start bidecomposition check with node %d \n", root);

   while( cutnodes->stack_size > 0 )
   {
      cutNodesProcessNext(g, cutnodes);
   }
}
#endif

/** post-process */
static
void cutNodesComputePostProcess(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   const int lastroot = cutnodes->artpoints[StpVecGetSize(cutnodes->artpoints) - 1];
   assert(cutnodes->biconn_comproots[0] == -1);

   cutnodes->biconn_comproots[0] = lastroot;
   cutnodes->biconn_ncomps++;
}

#ifndef NDEBUG
/** builds CSR like arrays for biconnected components */
static
SCIP_Bool decomposeCsrIsValid(
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g,                  /**< graph data structure */
   const BIDECOMP*       bidecomp            /**< bi-decomposition structure */
   )
{
   const int* const nodes = bidecomp->nodes;
   const int* const starts = bidecomp->starts;
   const int* const biconn_nodesmark = cutnodes->biconn_nodesmark;
   const int* const biconn_comproots = cutnodes->biconn_comproots;
   const int ncomps = cutnodes->biconn_ncomps;

   for( int i = 0; i < ncomps; i++ )
   {
      const int comp_start = starts[i];
      const int comp_end = starts[i + 1];
      const int comp_root = biconn_comproots[i];

      assert(comp_start <= comp_end);

      SCIPdebugMessage("component %d is of size %d (root=%d); nodes: \n", i, comp_end - comp_start,
            comp_root);

      for( int j = comp_start; j != comp_end; j++ )
      {
         const int comp_node = nodes[j];
         assert(graph_knot_isInRange(g, comp_node));

#ifdef SCIP_DEBUG
         graph_knot_printInfo(g, comp_node);
#endif

         if( biconn_nodesmark[comp_node] != i )
         {
            int k = 0;
            SCIP_Bool isCompRoot = FALSE;
            for( k = 0; k < ncomps; k++ )
            {
               if( comp_node == biconn_comproots[k] )
               {
                  isCompRoot = TRUE;
                  break;
               }
            }

            if( !isCompRoot )
            {
               printf("ERROR: no component root found \n");
               return FALSE;
            }
         }
      }
   }

   return TRUE;
}

#endif


/** builds CSR like arrays for biconnected components */
static
void decomposeBuildCsr(
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g,                  /**< graph data structure */
   BIDECOMP*             bidecomp            /**< bidecomposition data structure */
   )
{
   int* RESTRICT nodes = bidecomp->nodes;
   int* RESTRICT starts = bidecomp->starts;
   const int* const biconn_nodesmark = cutnodes->biconn_nodesmark;
   const int* const biconn_comproots = cutnodes->biconn_comproots;
   const int ncomps = cutnodes->biconn_ncomps;
   const int nnodes = graph_get_nNodes(g);

   for( int i = 0; i <= ncomps; i++ )
      starts[i] = 0;

   for( int i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && g->grad[i] > 0 )
      {
         const int compid = biconn_nodesmark[i];
         assert(0 <= compid && compid < ncomps);

         starts[compid]++;
      }
   }

   /* we also need to count the component roots */
   for( int i = 0; i < ncomps; i++ )
   {
      const int comproot = biconn_comproots[i];
      if( biconn_nodesmark[comproot] != i && g->grad[comproot] > 0 )
         starts[i]++;
   }

   for( int i = 1; i <= ncomps; i++ )
      starts[i] += starts[i - 1];

   assert(starts[ncomps] == starts[ncomps - 1]);

   /* now fill the values in */
   for( int i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && g->grad[i] > 0 )
      {
         const int compid = biconn_nodesmark[i];
         assert(0 <= compid && compid < ncomps);

         starts[compid]--;
         nodes[starts[compid]] = i;

         assert(compid == 0 || starts[compid - 1] <= starts[compid]);
      }
   }

   for( int i = 0; i < ncomps; i++ )
   {
      const int comproot = biconn_comproots[i];
      if( biconn_nodesmark[comproot] != i && g->grad[comproot] > 0 )
      {
         starts[i]--;
         nodes[starts[i]] = comproot;
      }
   }

   assert(starts[0] == 0);
   assert(starts[ncomps] <= nnodes + ncomps);

   assert(decomposeCsrIsValid(cutnodes, g, bidecomp));
}


/** gets first marked component */
static
SCIP_RETCODE decomposeGetFirstMarked(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< graph data structure */
   int*                  subroot             /**< the new root (out) */
   )
{
   const int* const gMark = orggraph->mark;
   STP_Bool* nodes_isVisited;
   STP_Vectype(int) stack = NULL;
   const int nnodes = graph_get_nNodes(orggraph);
   const int root = orggraph->source;
   SCIP_Bool markedIsFound = FALSE;

   if( gMark[root] )
   {
      *subroot = root;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_isVisited, nnodes) );

    for( int i = 0; i < nnodes; i++ )
       nodes_isVisited[i] = FALSE;

   nodes_isVisited[root] = TRUE;
   StpVecPushBack(scip, stack, root);

   /* DFS loop */
   while( !StpVecIsEmpty(stack) && !markedIsFound )
   {
      const int node = stack[StpVecGetSize(stack) - 1];
      assert(StpVecGetSize(stack) >= 1);
      StpVecPopBack(stack);

      assert(graph_knot_isInRange(orggraph, node));

      for( int a = orggraph->outbeg[node]; a != EAT_LAST; a = orggraph->oeat[a] )
      {
         const int head = orggraph->head[a];

         if( !nodes_isVisited[head] )
         {
            if( gMark[head] )
            {
               assert(!markedIsFound);

               markedIsFound = TRUE;
               *subroot = head;
               break;
            }
            StpVecPushBack(scip, stack, head);
            nodes_isVisited[head] = TRUE;
         }
      }
   }

   assert(markedIsFound);

   StpVecFree(scip, stack);
   SCIPfreeBufferArray(scip, &nodes_isVisited);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE bidecomposition_cutnodesInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES**            cutnodes            /**< cut nodes */
   )
{
   CUTNODES* cn;
   int* nodes_hittime;
   int* biconn_nodesmark;
   int* biconn_comproots;
   const int nnodes = graph_get_nNodes(g);

   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, cutnodes) );
   cn = *cutnodes;

#ifdef CUTTREE_PRINT_STATISTICS
   {
      int* childcount_nodes;
      int* childcount_terms;

      SCIP_CALL( SCIPallocMemoryArray(scip, &childcount_nodes, nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &childcount_terms, nnodes) );

      for( int k = 0; k < nnodes; k++ )
         childcount_nodes[k] = 0;

      for( int k = 0; k < nnodes; k++ )
         childcount_terms[k] = 0;

      cn->childcount_nodes = childcount_nodes;
      cn->childcount_terms = childcount_terms;
   }
#endif

   SCIP_CALL( SCIPallocMemoryArray(scip, &biconn_comproots, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &biconn_nodesmark, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodes_hittime, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cn->stack_nodes, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      nodes_hittime[k] = -1;

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      biconn_comproots[k] = -1;
#endif

   for( int k = 0; k < nnodes; k++ )
      biconn_nodesmark[k] = 0;

   cn->scip = scip;
   cn->artpoints = NULL;
   cn->biconn_stack = NULL;
   cn->biconn_comproots = biconn_comproots;
   cn->biconn_nodesmark = biconn_nodesmark;
   cn->nodes_hittime = nodes_hittime;
   cn->stack_size = 0;
   cn->biconn_ncomps = 0;
   cn->dfsroot = -1;
   cn->curr_lowpoint = -1;
   cn->curr_hittime = -1;
   cn->nrootcomps = -1;

   cutNodesSetDfsRoot(g, cn);

   StpVecReserve(scip, cn->biconn_stack, nnodes);

   return SCIP_OKAY;
}


/** exits */
void bidecomposition_cutnodesFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CUTNODES**            cutnodes            /**< cut nodes */
   )
{
   CUTNODES* cn;

   assert(scip && cutnodes);
   cn = *cutnodes;
   assert(cn);

   StpVecFree(scip, cn->artpoints);
   StpVecFree(scip, cn->biconn_stack);

   SCIPfreeMemoryArray(scip, &(cn->stack_nodes));
   SCIPfreeMemoryArray(scip, &(cn->biconn_nodesmark));
   SCIPfreeMemoryArray(scip, &(cn->biconn_comproots));
   SCIPfreeMemoryArray(scip, &(cn->nodes_hittime));

#ifdef CUTTREE_PRINT_STATISTICS
   SCIPfreeMemoryArray(scip, &(cn->childcount_terms));
   SCIPfreeMemoryArray(scip, &(cn->childcount_nodes));
#endif

   SCIPfreeMemory(scip, cutnodes);
}


/** computes cut-nodes and (implicitly) bi-connected components */
void bidecomposition_cutnodesCompute(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   assert(cutnodes);
   assert(StpVecGetSize(cutnodes->biconn_stack) == 0);
   assert(StpVecGetSize(cutnodes->artpoints) == 0);
   assert(cutnodes->curr_lowpoint == -1);

   cutnodes->curr_hittime = 0;
   cutnodes->nrootcomps = 0;
   /* NOTE: we assume the graph to be connected, so we only do the DFS once */
   /* todo make it non-recursive, otherwise it might crash for big graphs! */
   SCIPdebugMessage("starting DFS from %d \n", cutnodes->dfsroot);
#ifdef USE_RECURSIVE_DFS
   cutNodesComputeRecursive(g, cutnodes->dfsroot, -1, cutnodes);
#else
   cutNodesCompute(g, cutnodes);
#endif

   SCIPdebugMessage("number of cut nodes: %d \n", StpVecGetSize(cutnodes->artpoints));
   assert(cutnodes->biconn_ncomps >= StpVecGetSize(cutnodes->artpoints));

   if( cutnodes->biconn_ncomps > 0 )
   {
      assert(StpVecGetSize(cutnodes->artpoints) > 0);
      cutNodesComputePostProcess(g, cutnodes);

      SCIPdebugMessage("%d bi-connected components found! \n", cutnodes->biconn_ncomps);
   }
}

/** initializes */
SCIP_RETCODE bidecomposition_init(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g,                  /**< graph data structure */
   BIDECOMP**            bidecomposition     /**< bidecomposition data structure */
   )
{
   BIDECOMP* bidecomp;
   int* nodes;
   int* starts;
   const int nnodes = graph_get_nNodes(g);
   const int ncomps = cutnodes->biconn_ncomps;

   assert(scip && cutnodes);
   assert(ncomps >= 2);

   SCIP_CALL( SCIPallocMemory(scip, bidecomposition) );
   bidecomp = *bidecomposition;

   SCIP_CALL( SCIPallocMemoryArray(scip, &nodes, nnodes + ncomps) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &starts, ncomps + 1) );
   bidecomp->subinout = NULL;
   bidecomp->nodes = nodes;
   bidecomp->starts = starts;
   bidecomp->nbicomps = ncomps;

   decomposeBuildCsr(cutnodes, g, bidecomp);

   return SCIP_OKAY;
}


/** initializes */
SCIP_RETCODE bidecomposition_initSubInOut(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   BIDECOMP*             bidecomposition     /**< bidecomposition data structure */
   )
{
   assert(scip && g && bidecomposition);
   assert(!bidecomposition->subinout);

   SCIP_CALL( graph_subinoutInit(scip, g, &(bidecomposition->subinout)) );

   return SCIP_OKAY;
}


/** frees */
void bidecomposition_free(
   SCIP*                 scip,               /**< SCIP data structure */
   BIDECOMP**            bidecomposition     /**< bidecomposition data structure */
   )
{
   BIDECOMP* bidecomp;

   assert(scip && bidecomposition);

   bidecomp = *bidecomposition;
   assert(bidecomp);

   if( bidecomp->subinout )
      graph_subinoutFree(scip, &(bidecomp->subinout));

   SCIPfreeMemoryArray(scip, &(bidecomp->starts));
   SCIPfreeMemoryArray(scip, &(bidecomp->nodes));

   SCIPfreeMemory(scip, bidecomposition);
}


/** marks subgraph of given component index */
void bidecomposition_markSub(
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                g                   /**< graph data structure */
   )
{
   int* const gMark = g->mark;
   const SUBINOUT* const subinout = bidecomp->subinout;
   const int* const contractionRecord = graph_subinoutGetContractionRecord(subinout);
   const int* const compnodes = bidecomp->nodes;
   const int compstart = bidecomp->starts[compindex];
   const int compend = bidecomp->starts[compindex + 1];
   const int nnodes = graph_get_nNodes(g);

   for( int i = 0; i < nnodes; i++ )
      gMark[i] = FALSE;

   SCIPdebugMessage("marking subgraph: \n");

   for( int i = compstart; i != compend; i++ )
   {
      const int compnode = compnodes[i];
      int realnode;

      assert(graph_knot_isInRange(g, compnode));

      if( contractionRecord[compnode] != -1 )
      {
         realnode = graph_knot_getContractionRecordAncestor(compnode, subinout);

         SCIPdebugMessage("(taking contracted node %d instead of %d:) \n", contractionRecord[compnode], compnode);
      }
      else
      {
         realnode = compnode;
      }
#ifdef SCIP_DEBUG
      graph_knot_printInfo(g, realnode);
#endif

      assert(graph_knot_isInRange(g, realnode));
      gMark[realnode] = TRUE;
   }
}


/** gets root of marked sub-component
 *  bidecomposition_markSub needs to be called before! */
SCIP_RETCODE bidecomposition_getMarkedSubRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   const GRAPH*          orggraph,           /**< graph data structure */
   const GRAPH*          subgraph,           /**< graph data structure */
   int*                  subroot             /**< the new root (out) */
   )
{
   const int* nodemap_OrgToSub;

   assert(bidecomp && orggraph && subroot);
   assert(bidecomp->subinout);

   *subroot = -1;

   SCIP_CALL( decomposeGetFirstMarked(scip, orggraph, subroot) );

   assert(graph_knot_isInRange(orggraph, *subroot));
   assert(Is_term(orggraph->term[*subroot]));

   nodemap_OrgToSub = graph_subinoutGetOrgToSubNodeMap(bidecomp->subinout);
   *subroot = nodemap_OrgToSub[*subroot];

   assert(graph_knot_isInRange(subgraph, *subroot));
   assert(Is_term(subgraph->term[*subroot]));

   return SCIP_OKAY;
}


/** component consisting of at most one node? */
SCIP_Bool bidecomposition_componentIsTrivial(
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex           /**< component index */
   )
{
   assert(bidecomp);
   assert(0 <= compindex && compindex < bidecomp->nbicomps);

   {
      const int compstart = bidecomp->starts[compindex];
      const int compend = bidecomp->starts[compindex + 1];

      assert(compstart <= compend);

      if( compend - compstart <= 1 )
      {
         SCIPdebugMessage("component %d is of size %d, SKIP! \n", compindex, compend - compstart);
         return TRUE;
      }
   }

   return FALSE;
}


/** checks whether bidecomposition check is possible */
SCIP_Bool bidecomposition_isPossible(
   const GRAPH*          g                   /**< graph data structure */
   )
{
#ifdef USE_RECURSIVE_DFS
   int nnodes_real = 0;
   const int nnodes = graph_get_nNodes(g);
   const int* const isMarked = g->mark;

   assert(graph_isMarked(g));

   for( int i = 0; i < nnodes; i++ )
   {
      if( isMarked[i] )
         nnodes_real++;
   }

   return (nnodes_real < 75000);
#else
   return TRUE;
#endif
}


/** returns nodes ratio of component and the remaining graph */
SCIP_Real bidecomposition_getCompNodeRatio(
   const BIDECOMP*       bidecomp,           /**< bidecomposition data structure */
   int                   compindex           /**< component index */
   )
{
   const int* const starts = bidecomp->starts;
   SCIP_Real ratio;
   const int ncomps = bidecomp->nbicomps;
   const int nallnodes = starts[ncomps] - starts[0];
   int compnnodes;

   assert(bidecomp);
   assert(0 <= compindex && compindex < ncomps);
   assert(nallnodes > 0);

   compnnodes = starts[compindex + 1] - starts[compindex];
   ratio = (SCIP_Real) compnnodes / (SCIP_Real) nallnodes;

   SCIPdebugMessage("component nodes ratio: %f \n", ratio);
   assert(GE(ratio, 0.0));

   return (ratio);
}


/** returns ratio of nodes of maximum component and the remaining graph */
SCIP_Real bidecomposition_getMaxcompNodeRatio(
   const BIDECOMP*       bidecomp            /**< bidecomposition data structure */
   )
{
   const int* const starts = bidecomp->starts;
   SCIP_Real maxratio;
   const int ncomps = bidecomp->nbicomps;
   const int ncompnodes = starts[ncomps] - starts[0];
   int maxcompnnodes = 0;

   assert(bidecomp);
   assert(0 < ncompnodes);

   SCIPdebugMessage("all component nodes=%d \n", ncompnodes);

   for( int i = 0; i < ncomps; i++ )
   {
      const int compnnodes = starts[i + 1] - starts[i];
      if( maxcompnnodes < compnnodes )
         maxcompnnodes = compnnodes;
   }

   maxratio = (SCIP_Real) maxcompnnodes / (SCIP_Real) ncompnodes;

   SCIPdebugMessage("max. component number of nodes=%d \n", maxcompnnodes);
   SCIPdebugMessage("maxratio=%f \n", maxratio);

   assert(GT(maxratio, 0.0));

   return (maxratio);
}
