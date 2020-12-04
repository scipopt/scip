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

/**@file   reduce_sepa.c
 * @brief  several node-separator based reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements node-separator based reduction techniques for several Steiner problems.
 * It contains rather simple test with articulation points, and more involved ones with terminal-separators.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"

#define BIDECOMP_MINRED_MULTIPLIER 2
#define BIDECOMP_MINCOMPRATIO_FIRST     0.95
#define BIDECOMP_MINCOMPRATIO           0.80
//#define CUTTREE_PRINT_STATISTICS


/** decompose */
typedef struct biconnected_component_decomposition
{
   SUBINOUT*             subinout;           /**< helper */
   int*                  nodes;              /**< nodes */
   int*                  starts;             /**< starts array per component */
   int                   nbicomps;           /**< number of components */
} BIDECOMP;


/** cut nodes/ articulation points */
typedef struct cut_nodes
{
   SCIP*                 scip;               /**< SCIP data structure */
#ifdef CUTTREE_PRINT_STATISTICS
   int*                  childcount_nodes;   /**< number of nodes below each node */
   int*                  childcount_terms;   /**< number of terminals below each node */
#endif
   STP_Vectype(int)      biconn_stack;       /**< stack for marking bi-connected component */
   int*                  biconn_nodesmark;   /**< marks in which component each node is 0, 1,.., biconn_ncomps - 1 */
   int*                  biconn_comproots;   /**< root of each component with index 0,1,...,biconn_ncomps - 1 */
   STP_Vectype(int)      artpoints;          /**< cut nodes */
   int*                  nodes_hittime;      /**< hit time 0,1,... */
   int                   biconn_ncomps;      /**< number of components */
   int                   dfsroot;            /**< root */
   int                   curr_lowpoint;      /**< current low-point */
   int                   curr_hittime;       /**< current hit time */
} CUTNODES;


/** Steiner tree based bi-connected component reduction */
typedef struct cut_tree_data
{
   int*                  nodes_prednode;     /**< predecessor per node */
   SCIP_Bool*            comps_isHit;        /**< of size ncomps */
   SCIP_Bool*            nodes_isTree;       /**< of size |V| */
   SCIP_Bool             cutnode0isNeeded;   /**< special treatment */
} CUTTREE;


/** helper */
static inline
int cutNodesGetLastCutnode(
   const CUTNODES*       cutnodes            /**< cut nodes */
   )
{
   const int lastcutnode = cutnodes->artpoints[StpVecGetSize(cutnodes->artpoints) - 1];
   assert(StpVecGetSize(cutnodes->artpoints) >= 1);

   assert(lastcutnode >= 0);
   assert(cutnodes->biconn_nodesmark[lastcutnode] == 0);

   return lastcutnode;
}

#ifdef CUTTREE_PRINT_STATISTICS
/** helper */
static inline
void cutNodesTraverseSub(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,
   SCIP_Bool*            isVisited,          /**< */
   int*                  nterms,
   int*                  ncompnodes
   )
{
   STP_Vectype(int) stack = NULL;
   const int nnodes = graph_get_nNodes(g);

   *nterms = 0;
   *ncompnodes = 0;
   assert(!isVisited[node]);

   StpVecReserve(scip, stack, nnodes);
   StpVecPushBack(scip, stack, node);
   isVisited[node] = TRUE;

   while( StpVecGetSize(stack)  )
   {
      const int k = stack[StpVecGetSize(stack) - 1];
      StpVecPopBack(stack);

      if( Is_term(g->term[k]) )
      {
         (*nterms)++;
      }
      (*ncompnodes)++;

      assert(isVisited[k]);

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         if( !isVisited[head] )
         {
            isVisited[head] = TRUE;
            StpVecPushBack(scip, stack, head);
         }
      }
   }

   StpVecFree(scip, stack);
}


/** traverses tree from cut-node */
static inline
SCIP_RETCODE cutNodesTraverseFromCutNode(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,
   const CUTNODES*       cutnodes            /**< cut nodes */
   )
{
   SCIP_Bool* isVisited;
   const int nnodes = graph_get_nNodes(g);
   STP_Vectype(int) comps_id = NULL;
   STP_Vectype(int) comps_nterms = NULL;

   assert(StpVecGetSize(cutnodes->artpoints) >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &isVisited, nnodes) );

   printf("check cut-node ");
   graph_knot_printInfo(g, node);

   for( int i = 0; i < nnodes; ++i )
      isVisited[i] = FALSE;

   isVisited[node] = TRUE;

   printf(" belowterms=%d, belownodes=%d  ... \n", cutnodes->childcount_terms[node], cutnodes->childcount_nodes[node]);

   for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
   {
      int nterms;
      int ncompnodes;
      const int base = g->head[e];

      if( isVisited[base] )
         continue;

      cutNodesTraverseSub(scip, g, base, isVisited, &nterms, &ncompnodes);
      StpVecPushBack(scip, comps_id, cutnodes->biconn_nodesmark[base]);
      StpVecPushBack(scip, comps_nterms, nterms);

      printf("component=%d: ", cutnodes->biconn_nodesmark[base]);
      printf("nterms=%d, ", nterms);
      printf("ncompnodes=%d \n", ncompnodes);
   }



   StpVecFree(scip, comps_nterms);
   StpVecFree(scip, comps_id);
   SCIPfreeBufferArray(scip, &isVisited);

   return SCIP_OKAY;
}
#endif

#ifdef XXXXX
/** checks bi-connected leaf components */
static
SCIP_RETCODE cutNodesTreeCheckLeaveComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   CUTTREE*              cuttree             /**< cut tree data */
   )
{
   STP_Vectype(int) stack = NULL;
   int* RESTRICT nodes_pred = cuttree->nodes_prednode;
   SCIP_Bool* RESTRICT nodes_isVisited = cutnodes->nodes_isVisited;
   SCIP_Bool* RESTRICT nodes_isTree = cuttree->nodes_isTree;
   const int* const biconn_nodesmark = biconn_nodesmark;
   const int nnodes = graph_get_nNodes(g);
   const int root = cutnodes->dfsroot;
   const int lastcutnode = cutNodesGetLastCutnode(cutnodes);
   int termscount = g->terms;

   /* save terminals per component, set to -1 if cut node is encountered */

   return SCIP_OKAY;
}
#endif


/** helper */
static inline
void cutNodesTreeAddNode(
   int                   node,
   const CUTNODES*       cutnodes,           /**< cut nodes */
   int                   lastcutnode,
   CUTTREE*              cuttree             /**< cut tree data */
   )
{
   const int nodecomp = cutnodes->biconn_nodesmark[node];
   assert(0 <= nodecomp && nodecomp < cutnodes->biconn_ncomps);

   cuttree->nodes_isTree[node] = TRUE;

   if( node != lastcutnode )
      cuttree->comps_isHit[nodecomp] = TRUE;

   /* NOTE: lastcutnode needs special treatment, because it is considered as part of the
    * bi-connected component that it induces */
   if( nodecomp != 0 && cutnodes->biconn_comproots[nodecomp] == lastcutnode )
   {
      cuttree->cutnode0isNeeded = TRUE;
      assert(node != lastcutnode);
   }
}


/** builds Steiner tree */
static
void cutNodesTreeBuildSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   CUTTREE*              cuttree             /**< cut tree data */
   )
{
   STP_Vectype(int) stack = NULL;
   int* RESTRICT nodes_pred = cuttree->nodes_prednode;
   SCIP_Bool* RESTRICT nodes_isTree = cuttree->nodes_isTree;
   const int* const biconn_nodesmark = biconn_nodesmark;
   const int nnodes = graph_get_nNodes(g);
   const int root = cutnodes->dfsroot;
   const int lastcutnode = cutNodesGetLastCutnode(cutnodes);
   const SCIP_Bool isPcMw = graph_pc_isPcMw(g);
   int termscount = isPcMw ? graph_pc_nNonLeafTerms(g) : g->terms;

   assert(nodes_pred && nodes_isTree && cuttree->comps_isHit);
   assert(!nodes_isTree[root] && !cuttree->comps_isHit[cutnodes->biconn_nodesmark[root]]);
   assert(!cuttree->cutnode0isNeeded);

   StpVecReserve(scip, stack, nnodes);

   nodes_pred[root] = root;
   cutNodesTreeAddNode(root, cutnodes, lastcutnode, cuttree);
   StpVecPushBack(scip, stack, root);
   termscount--;

   /* do a DFS until all terminals have been hit */
   while( StpVecGetSize(stack) > 0 && termscount > 0 )
   {
      const int node = stack[StpVecGetSize(stack) - 1];
      StpVecPopBack(stack);

      for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( nodes_pred[head] < 0 && g->mark[head] )
         {
            nodes_pred[head] = node;

            StpVecPushBack(scip, stack, head);

            if( Is_term(g->term[head]) )
            {
               assert(!nodes_isTree[head]);

               if( isPcMw && graph_pc_termIsNonLeafTerm(g, head) )
                  continue;

               termscount--;

               for( int pred = head; !nodes_isTree[pred]; pred = nodes_pred[pred] )
               {
                  assert(graph_knot_isInRange(g, pred));
                  cutNodesTreeAddNode(pred, cutnodes, lastcutnode, cuttree);
               }
            }
         }
      }
   }

   cuttree->cutnode0isNeeded = (cuttree->cutnode0isNeeded && (nodes_isTree[lastcutnode]));

   assert(termscount == 0 || isPcMw);
   StpVecFree(scip, stack);
}


/** deletes unvisited components */
static
void cutNodesTreeDeleteComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const CUTTREE*        cuttree,            /**< cut tree data */
   GRAPH*                g,                  /**< graph */
   SCIP_Real*            fixedp,             /**< pointer to offset value */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   const SCIP_Bool* const comps_isHit = cuttree->comps_isHit;
   const int* const biconn_nodesmark = cutnodes->biconn_nodesmark;
   const int lastcutnode = cutNodesGetLastCutnode(cutnodes);
   const int nnodes = graph_get_nNodes(g);

   for( int i = 0; i < nnodes; i++ )
   {
      int comp;

      if( !g->mark[i] )
         continue;

      comp = biconn_nodesmark[i];

      if( !comps_isHit[comp] && i != lastcutnode )
      {

#ifdef SCIP_DEBUG
         SCIPdebugMessage("delete: ");
         graph_knot_printInfo(g, i);
#endif
         assert(!cuttree->nodes_isTree[i]);

         if( Is_term(g->term[i]) )
         {
            assert(graph_pc_isPcMw(g));
            assert(graph_pc_termIsNonLeafTerm(g, i));

            graph_pc_deleteTerm(scip, g, i, fixedp);
         }
         else
         {
            graph_knot_del(scip, g, i, TRUE);
         }

         (*nelims)++;
      }
   }
}


#ifndef NDEBUG
/** debug checker: makes sure that all remaining proper cut nodes are terminals */
static
SCIP_Bool cutNodesTreeMakeTermsIsComplete(
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g                   /**< graph */
   )
{
   const int nnodes = graph_get_nNodes(g);

   // todo probably wrong
#ifdef SCIP_DISABLED
   const int* const biconn_nodesmark = cutnodes->biconn_nodesmark;
   const int* biconn_comproot = cutnodes->biconn_comproots;

   /* 1. make sure that the components of non-terminal cut-nodes are empty */

   for( int i = 0; i < cutnodes->biconn_ncomps; i++ )
   {
      const int cutnode = biconn_comproot[i];
      assert(graph_knot_isInRange(g, cutnode));

      if( g->grad[cutnode] == 0 || Is_term(g->term[cutnode]) )
         continue;

      for( int k = 0; k < nnodes; k++ )
      {
         if( biconn_nodesmark[k] != i  )
            continue;

         if( g->grad[k] > 0 && g->mark[k] && k != cutnode )
         {
            printf("issue for \n");
            graph_knot_printInfo(g, k);
            return FALSE;
         }
      }
   }
#endif

   /* 2. make sure that for any nonterminal all incident edges are in the same component! */

   for( int i = 0; i < nnodes; i++ )
   {
      const int compid = cutnodes->biconn_nodesmark[i];

      if( g->grad[i] == 0 || Is_term(g->term[i]) )
         continue;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( cutnodes->biconn_nodesmark[head] != compid && !Is_term(g->term[head]) )
         {
            printf("issue for compid=%d, %d \n", compid, cutnodes->biconn_nodesmark[head] );
            graph_knot_printInfo(g, i);
            graph_knot_printInfo(g, head);

            return FALSE;
         }
      }
   }


   return TRUE;
}
#endif


/** changes required cut-nodes from non-terminals to terminals */
static
void cutNodesTreeMakeTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const CUTTREE*        cuttree,            /**< cut tree data */
   GRAPH*                g                   /**< graph */
   )
{
   const int* const biconn_nodesmark = cutnodes->biconn_nodesmark;
   const SCIP_Bool* const nodes_isTree = cuttree->nodes_isTree;
   const int ncutnodes = StpVecGetSize(cutnodes->artpoints);

   assert(nodes_isTree);
   assert(ncutnodes > 0);

   for( int i = 0; i < ncutnodes; i++ )
   {
      const int cutnode = cutnodes->artpoints[i];
      assert(graph_knot_isInRange(g, cutnode));

      if( g->grad[cutnode] == 0 )
         continue;

      if( i == ncutnodes - 1 && !cuttree->cutnode0isNeeded )
         continue;

      if( Is_term(g->term[cutnode]) )
      {
         assert(nodes_isTree[cutnode]);
         continue;
      }

      if( nodes_isTree[cutnode] )
      {
         const int compid = biconn_nodesmark[g->head[g->outbeg[cutnode]]];
         SCIP_Bool isTerm = FALSE;

         for( int e = g->outbeg[cutnode]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            if( biconn_nodesmark[head] != compid )
            {
               isTerm = TRUE;
               break;
            }
         }

         if( isTerm )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage("cut node to terminal: ");
            graph_knot_printInfo(g, cutnode);
#endif
            graph_knot_chg(g, cutnode, STP_TERM);
         }
      }
   }

   assert(cutNodesTreeMakeTermsIsComplete(cutnodes, g));
}


/** initializes */
static
SCIP_RETCODE cutNodesTreeInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   CUTTREE*              cuttree             /**< cut tree data */
   )
{
   int* nodes_prednode;
   SCIP_Bool* comps_isHit;
   SCIP_Bool* nodes_isTree;
   const int nnodes = graph_get_nNodes(g);
   const int ncomps = cutnodes->biconn_ncomps;

   assert(ncomps > 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_prednode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_isTree, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comps_isHit, ncomps) );

   for( int i = 0; i < nnodes; i++ )
   {
      nodes_isTree[i] = FALSE;
      nodes_prednode[i] = -1;
   }

   for( int i = 0; i < ncomps; i++ )
      comps_isHit[i] = FALSE;

   cuttree->nodes_prednode = nodes_prednode;
   cuttree->nodes_isTree = nodes_isTree;
   cuttree->comps_isHit = comps_isHit;

   return SCIP_OKAY;
}


/** exits */
static
void cutNodesTreeExit(
   SCIP*                 scip,               /**< SCIP data structure */
   CUTTREE*              cuttree             /**< cut tree data */
   )
{
   SCIPfreeBuffer(scip, &(cuttree->comps_isHit));
   SCIPfreeBuffer(scip, &(cuttree->nodes_isTree));
   SCIPfreeBuffer(scip, &(cuttree->nodes_prednode));
}



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


/** initializes */
static
SCIP_RETCODE cutNodesInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   const int nnodes = g->knots;

   int* nodes_hittime;
   int* biconn_nodesmark;
   int* biconn_comproots;

   assert(cutnodes->artpoints == NULL);
   assert(cutnodes->biconn_stack == NULL);
   assert(cutnodes->biconn_ncomps == 0);

#ifdef CUTTREE_PRINT_STATISTICS
   {
      int* childcount_nodes;
      int* childcount_terms;

      SCIP_CALL( SCIPallocBufferArray(scip, &childcount_nodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &childcount_terms, nnodes) );

      for( int k = 0; k < nnodes; k++ )
         childcount_nodes[k] = 0;

      for( int k = 0; k < nnodes; k++ )
         childcount_terms[k] = 0;

      cutnodes->childcount_nodes = childcount_nodes;
      cutnodes->childcount_terms = childcount_terms;
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &biconn_comproots, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &biconn_nodesmark, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_hittime, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      nodes_hittime[k] = -1;

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      biconn_comproots[k] = -1;
#endif

   for( int k = 0; k < nnodes; k++ )
      biconn_nodesmark[k] = 0;

   cutnodes->scip = scip;

   cutnodes->biconn_comproots = biconn_comproots;
   cutnodes->biconn_nodesmark = biconn_nodesmark;
   cutnodes->nodes_hittime = nodes_hittime;
   cutnodes->curr_lowpoint = -1;
   cutnodes->curr_hittime = -1;

   cutNodesSetDfsRoot(g, cutnodes);

   StpVecReserve(scip, cutnodes->biconn_stack, nnodes);

   return SCIP_OKAY;
}


/** exits */
static
void cutNodesExit(
   SCIP*                 scip,               /**< SCIP data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   StpVecFree(scip, cutnodes->artpoints);
   StpVecFree(scip, cutnodes->biconn_stack);

   SCIPfreeBufferArray(scip, &(cutnodes->nodes_hittime));
   SCIPfreeBufferArray(scip, &(cutnodes->biconn_nodesmark));
   SCIPfreeBufferArray(scip, &(cutnodes->biconn_comproots));

#ifdef CUTTREE_PRINT_STATISTICS
   SCIPfreeBufferArray(scip, &(cutnodes->childcount_terms));
   SCIPfreeBufferArray(scip, &(cutnodes->childcount_nodes));
#endif
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


/** recursive DFS */
static
void cutNodesComputeDfs(
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

         cutNodesComputeDfs(g, head, node, cutnodes);

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

#ifdef XXX_XXX
    for( int i = 0; i < g->knots; i++ )
    {
       printf("%d in comp %d \n", i, cutnodes->biconn_nodesmark[i]);
    }

    for( int i = 0; i < cutnodes->biconn_ncomps; i++ )
    {
       printf("comp-root[%d]=%d \n", i, cutnodes->biconn_comproots[i]);
    }
#endif

}

/** computes cut-nodes and (implicitly) bi-connected components */
static inline
void cutNodesCompute(
   const GRAPH*          g,                  /**< graph data structure */
   CUTNODES*             cutnodes            /**< cut nodes */
   )
{
   assert(StpVecGetSize(cutnodes->biconn_stack) == 0);
   assert(StpVecGetSize(cutnodes->artpoints) == 0);

   cutnodes->curr_hittime = 0;
   /* NOTE: we assume the graph to be connected, so we only do the DFS once */
   /* todo make it non-recursive, otherwise it might crash for big graphs! */
   SCIPdebugMessage("starting DFS from %d \n", cutnodes->dfsroot);
   cutNodesComputeDfs(g, cutnodes->dfsroot, -1, cutnodes);

   SCIPdebugMessage("number of cut nodes: %d \n", StpVecGetSize(cutnodes->artpoints));
   assert(cutnodes->biconn_ncomps >= StpVecGetSize(cutnodes->artpoints));

   if( cutnodes->biconn_ncomps > 0 )
   {
      assert(StpVecGetSize(cutnodes->artpoints) > 0);
      cutNodesComputePostProcess(g, cutnodes);

      SCIPdebugMessage("%d bi-connected components found! \n", cutnodes->biconn_ncomps);
   }
}


/** removes non-necessary bi-connected components and creates new terminals from cut-nodes */
static
SCIP_RETCODE cutNodesReduceWithTree(
   SCIP*                 scip,               /**< SCIP data structure */
   CUTNODES*             cutnodes,           /**< cut nodes */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixedp,             /**< pointer to offset value */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   CUTTREE cuttree = { NULL, NULL, NULL, FALSE };

   SCIP_CALL( cutNodesTreeInit(scip, g, cutnodes, &cuttree) );
   cutNodesTreeBuildSteinerTree(scip, g, cutnodes, &cuttree);

#ifdef CUTTREE_PRINT_STATISTICS
   printf("\n STATS:  \n");

   for( int i = 0; i < StpVecGetSize(cutnodes->artpoints); i++ )
   {
      const int cutnode = cutnodes->artpoints[i];
      SCIP_CALL( cutNodesTraverseFromCutNode(scip, g, cutnode, cutnodes) );
   }
#endif

   cutNodesTreeDeleteComponents(scip, cutnodes, &cuttree, g, fixedp, nelims);

   if( !graph_pc_isPcMw(g) )
      cutNodesTreeMakeTerms(scip, cutnodes, &cuttree, g);

   cutNodesTreeExit(scip, &cuttree);

   return SCIP_OKAY;
}


/** builds CSR like arrays for biconnected components */
static
void decomposeBuildCsr(
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g,                  /**< graph data structure */
   BIDECOMP*             bidecomp
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


#ifndef NDEBUG
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

            assert(isCompRoot);
         }
      }
   }
#endif

}


/** initalizes */
static
SCIP_RETCODE decomposeInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   const GRAPH*          g,                  /**< graph data structure */
   BIDECOMP*             bidecomp
   )
{
   int* nodes;
   int* starts;
   const int nnodes = graph_get_nNodes(g);
   const int ncomps = cutnodes->biconn_ncomps;

   assert(ncomps >= 2);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nnodes + ncomps) );
   SCIP_CALL( SCIPallocBufferArray(scip, &starts, ncomps + 1) );
   SCIP_CALL( graph_subinoutInit(scip, g, &(bidecomp->subinout)) );

   bidecomp->nodes = nodes;
   bidecomp->starts = starts;
   bidecomp->nbicomps = ncomps;

   decomposeBuildCsr(cutnodes, g, bidecomp);

   return SCIP_OKAY;
}


/** frees */
static
void decomposeFreeMembers(
   SCIP*                 scip,               /**< SCIP data structure */
   BIDECOMP*             bidecomp
   )
{
   assert(scip && bidecomp);

   graph_subinoutFree(scip, &(bidecomp->subinout));
   SCIPfreeBufferArray(scip, &(bidecomp->starts));
   SCIPfreeBufferArray(scip, &(bidecomp->nodes));
}


/** reduces subproblem */
static
SCIP_RETCODE decomposeReduceSub(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase             /**< reduction stuff */
   )
{
   GRAPH* subgraph;
   const SUBINOUT* const subinout = bidecomp->subinout;
   int* const gMark = g->mark;
   const int* const compnodes = bidecomp->nodes;
   const int* const contractionRecord = graph_subinoutGetContractionRecord(subinout);
   const int nnodes = graph_get_nNodes(g);
   const int compstart = bidecomp->starts[compindex];
   const int compend = bidecomp->starts[compindex + 1];

   assert(compstart <= compend);
   assert(redbase && redbase->bidecompparams);

   if( compend - compstart <= 1 )
   {
      SCIPdebugMessage("component %d is of size %d, SKIP! \n", compindex, compend - compstart);
      return SCIP_OKAY;
   }

   SCIPdebugMessage("(depth %d) reduce component %d: \n", redbase->bidecompparams->depth, compindex);

   for( int i = 0; i < nnodes; i++ )
      gMark[i] = FALSE;

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

   /** extract, reduce, re-insert: */

   SCIP_CALL(graph_subgraphExtract(scip, g, bidecomp->subinout, &subgraph));
#ifdef SCIP_DEBUG
   SCIPdebugMessage("subgraph before reduction: ");
   graph_printInfoReduced(subgraph);
#endif

   {
      const int reductbound_org = redbase->redparameters->reductbound;
      redbase->redparameters->reductbound =
            BIDECOMP_MINRED_MULTIPLIER * reduce_getMinNreductions(subgraph, redbase->redparameters->reductbound_min);

      SCIPdebugMessage("subgraph: reductbound_min=%d reductbound=%d \n",
            redbase->redparameters->reductbound_min,redbase->redparameters->reductbound );

      redbase->bidecompparams->newLevelStarted = TRUE;

      SCIP_CALL(redLoopStp(scip, subgraph, redbase));

      redbase->redparameters->reductbound = reductbound_org;
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("subgraph after reduction: ");
   graph_printInfoReduced(subgraph);
#endif
   SCIP_CALL(graph_subgraphReinsert(scip, bidecomp->subinout, g, &subgraph));

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}




/** is promising? */
static
SCIP_Bool decomposeIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const BIDECPARAMS*    bidecompparams,     /**< bidecomposition */
   const BIDECOMP*       bidecomp
   )
{
   const int* starts = bidecomp->starts;
   SCIP_Real maxratio;
   const int ncomps = bidecomp->nbicomps;
   const int ncompnodes = starts[ncomps] - starts[0];
   int maxcompnnodes = 0;
   const SCIP_Real mincompratio = bidecompparams->depth == 0 ? BIDECOMP_MINCOMPRATIO_FIRST : BIDECOMP_MINCOMPRATIO;

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

   return (maxratio < mincompratio);
}

/** solves biconnected components separately */
static
SCIP_RETCODE decomposeExec(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction stuff */
   SCIP_Bool*            wasDecomposed       /**< performed recursive reduction? */

   )
{
   BIDECOMP bidecomp;

   /* NOTE: does not work because we do node contractions when reinserting the subgraph,
    * and because order of nodes is changed in subbgraph */
   assert(redbase->solnode == NULL && "not supported");
   assert(redbase->bidecompparams);
   assert(redbase->bidecompparams->depth < redbase->bidecompparams->maxdepth);
   assert(graph_valid(scip, g));
   assert(*wasDecomposed == FALSE);

   SCIP_CALL( decomposeInit(scip, cutnodes, g, &bidecomp) );

   if( decomposeIsPromising(g, redbase->bidecompparams, &bidecomp) )
   {
      redbase->bidecompparams->depth++;

      /* reduce each biconnected component individually */
      for( int i = 0; i < bidecomp.nbicomps; i++ )
      {
         SCIP_CALL( decomposeReduceSub(scip, &bidecomp, i, g, redbase) );
      }

      *wasDecomposed = TRUE;

      redbase->bidecompparams->depth--;
   }

   decomposeFreeMembers(scip, &bidecomp);

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/*
 * Interface methods
 */



/** decomposition into biconnected components and recursive reduction */
SCIP_RETCODE reduce_bidecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   SCIP_Bool*            wasDecomposed       /**< performed recursive reduction? */
   )
{
#ifdef CUTTREE_PRINT_STATISTICS
   CUTNODES cutnodes = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, -1, -1, -1 };
#else
   CUTNODES cutnodes = { NULL, NULL, NULL, NULL, NULL, NULL, 0, -1, -1, -1 };
#endif

   assert(scip && g && redbase && wasDecomposed);
   assert(graph_typeIsSpgLike(g) && "only SPG decomposition supported yet");

   *wasDecomposed = FALSE;

   if( g->terms == 1 )
      return SCIP_OKAY;

   graph_mark(g);

   SCIP_CALL( cutNodesInit(scip, g, &cutnodes) );
   cutNodesCompute(g, &cutnodes);

   if( cutnodes.biconn_ncomps > 0 )
   {
      int dummy = 0;
      /* get rid of non-required biconnected components (without terminals) */
      SCIP_CALL( cutNodesReduceWithTree(scip, &cutnodes, g, redbase->fixed, &dummy) );

      /* decompose and reduce recursively? */
      SCIP_CALL( decomposeExec(scip, &cutnodes, g, redbase, wasDecomposed) );
   }

   cutNodesExit(scip, &cutnodes);

   return SCIP_OKAY;
}


/** articulation points based, simple reduction */
SCIP_RETCODE reduce_articulations(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixedp,             /**< pointer to offset value */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
#ifdef CUTTREE_PRINT_STATISTICS
   CUTNODES cutnodes = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, -1, -1, -1 };
#else
   CUTNODES cutnodes = { NULL, NULL, NULL, NULL, NULL, NULL, 0, -1, -1, -1 };
#endif

   assert(scip && g && nelims);
   graph_mark(g);

   *nelims = 0;

   SCIP_CALL( cutNodesInit(scip, g, &cutnodes) );
   cutNodesCompute(g, &cutnodes);

   if( cutnodes.biconn_ncomps > 0 )
   {
      SCIP_CALL( cutNodesReduceWithTree(scip, &cutnodes, g, fixedp, nelims) );
   }

   cutNodesExit(scip, &cutnodes);

   return SCIP_OKAY;
}
