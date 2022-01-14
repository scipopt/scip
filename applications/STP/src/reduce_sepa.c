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
//#define BENCH_SEPA_PARTIAL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "substpsolver.h"
#include "bidecomposition.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"

#ifdef BENCH_SEPA_PARTIAL
#include "time.h"
#endif

#define BIDECOMP_MINRED_MULTIPLIER 2
#define BIDECOMP_MINMAXCOMPRATIO_FIRST     0.95
#define BIDECOMP_MINMAXCOMPRATIO           0.80
#define BIDECOMP_MINMAXCOMPRATIO_PARTIAL   0.99
#define BIDECOMP_MAXCOMPRATIO_PARTIAL      0.05
#define BIDECOMP_MAXCOMPRATIO_AGGRESSIVE   0.5
#define BIDECOMP_MINNODES                  10
#define BIDECOMP_MINNODES_PARTIAL          10


//#define CUTTREE_PRINT_STATISTICS


/** Steiner tree based bi-connected component reduction */
typedef struct cut_tree_data
{
   int*                  nodes_prednode;     /**< predecessor per node */
   SCIP_Bool*            comps_isHit;        /**< of size ncomps */
   SCIP_Bool*            nodes_isTree;       /**< of size |V| */
   SCIP_Bool             cutnode0isNeeded;   /**< special treatment */
} CUTTREE;


/*
 * Local methods
 */

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
   int                   node,               /**< node to add  */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   int                   lastcutnode,        /**< last cut node */
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
         const int compid = biconn_nodesmark[cutnode];
         SCIP_Bool isTerm = FALSE;

         for( int e = g->outbeg[cutnode]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            if( biconn_nodesmark[head] != compid && head != cutnodes->biconn_comproots[compid] )
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

            // todo remove me
                 if(g->grad[cutnode] == 1 )
                 {
                    printf("error in cut nodes terms");
                    exit(1);
                 }

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


/** helper that extracts, reduces, re-inserts */
static
SCIP_RETCODE decomposeReduceSubDoIt(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                graph,              /**< graph data structure */
   REDBASE*              redbase             /**< reduction stuff */
   )
{
   REDSOL* redsol = redbase->redsol;
   SUBINOUT* subinout = bidecomp->subinout;
   GRAPH* subgraph;

   assert(graph_valid(scip, graph));

   SCIP_CALL( graph_subgraphExtract(scip, graph, subinout, &subgraph) );

   SCIP_CALL( reduce_solLevelTopUpdate(scip, subgraph, redsol) );
   SCIP_CALL( reduce_solLevelTopTransferSolTo(graph_subinoutGetOrgToSubNodeMap(subinout), redsol) );

#ifdef SCIP_DEBUG
   SCIPdebugMessage("subgraph before reduction: ");
   graph_printInfoReduced(subgraph);
#endif

   {
      RPARAMS* redparams = redbase->redparameters;
      const int reductbound_org = redparams->reductbound;
      const SCIP_Real offset_org = reduce_solGetOffset(redsol);
      redparams->reductbound = BIDECOMP_MINRED_MULTIPLIER * reduce_getMinNreductions(subgraph, redparams->reductbound_min);
      SCIPdebugMessage("subgraph: reductbound_min=%d reductbound=%d \n",
            redparams->reductbound_min,redparams->reductbound );

      reduce_solSetOffset(0.0, redsol);
      redbase->bidecompparams->newLevelStarted = TRUE;

      SCIP_CALL(reduce_redLoopStp(scip, subgraph, redbase));

      redparams->reductbound = reductbound_org;
      reduce_solSetOffset(reduce_solGetOffset(redsol) + offset_org, redsol);

      assert(GE(reduce_solGetOffset(redsol), offset_org));
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("subgraph after reduction: ");
   graph_printInfoReduced(subgraph);
#endif
   SCIP_CALL(graph_subgraphReinsert(scip, subinout, graph, &subgraph));

   reduce_solLevelTopTransferSolBack(graph_subinoutGetSubToOrgNodeMap(subinout), redsol);
   reduce_solLevelTopClean(scip, redsol);

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
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
   assert(scip && g && bidecomp && redbase);
   assert(redbase->bidecompparams);

   if( bidecomposition_componentIsTrivial(bidecomp, compindex) )
   {
      return SCIP_OKAY;
   }

   SCIPdebugMessage("(depth %d) reduce component %d: \n", redbase->bidecompparams->depth, compindex);

   bidecomposition_markSub(bidecomp, compindex, g);
   SCIP_CALL( decomposeReduceSubDoIt(scip, bidecomp, compindex, g, redbase) );

   return SCIP_OKAY;
}


/** fixes original edges */
static
SCIP_RETCODE decomposeExactFixSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const SUBINOUT*       subinout,           /**< helper for problem mapping */
   SUBSTP*               substp,             /**< sub-problem */
   GRAPH*                orggraph,           /**< original graph */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  nelims              /**< number of eliminations or NULL */
   )
{
   SCIP_Real* offset = reduce_solGetOffsetPointer(redbase->redsol);
   int* subedges_sol;
   const int* const subedgesToOrg = graph_subinoutGetSubToOrgEdgeMap(subinout);
   const int nsubedges = substpsolver_getNsubedges(substp);

   assert(nsubedges > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &subedges_sol, nsubedges) );
   SCIP_CALL( substpsolver_getSolution(substp, subedges_sol) );

   for( int i = 0; i < nsubedges; i += 2 )
   {
      const int orgedge = subedgesToOrg[i];
      assert(graph_edge_isInRange(orggraph, orgedge));

      if( subedges_sol[i] == CONNECT || subedges_sol[i + 1] == CONNECT  )
      {
         SCIPdebugMessage("fix edge %d \n", orgedge);

         *offset += orggraph->cost[orgedge];

         /* NOTE: edge will be automatically contracted later; we avoid trouble with other
          * connected components markers in this way */
         orggraph->cost[orgedge] = 0.0;
         orggraph->cost[flipedge(orgedge)] = 0.0;
         graph_knot_chg(orggraph, orggraph->tail[orgedge], STP_TERM);
         graph_knot_chg(orggraph, orggraph->head[orgedge], STP_TERM);
      }
      else
      {
         SCIPdebugMessage("delete edge %d \n", orgedge);

         if( nelims )
            (*nelims)++;

         assert(subedges_sol[i] == UNKNOWN && subedges_sol[i + 1] == UNKNOWN);
         graph_edge_del(scip, orggraph, orgedge, TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &subedges_sol);

   return SCIP_OKAY;
}



/** solves subproblems */
static
SCIP_RETCODE decomposeExactSubDoIt(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   GRAPH*                orggraph,           /**< graph data structure */
   GRAPH*                subgraph,           /**< sub-graph */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  nelims              /**< number of eliminations or NULL */
   )
{
   SUBINOUT* subinout = bidecomp->subinout;
   SUBSTP* substp;
   SCIP_Bool success;

   /* NOTE: subgraph will be moved into substp! */
   SCIP_CALL( substpsolver_init(scip, subgraph, &substp) );
   SCIP_CALL( substpsolver_transferHistory(graph_subinoutGetSubToOrgEdgeMap(subinout),
      orggraph, substp) );

   SCIP_CALL( substpsolver_setMute(substp) );
   //SCIP_CALL( substpsolver_setProbIsIndependent(substp) );
   SCIP_CALL( substpsolver_setProbFullPresolve(substp) );

#ifdef SCIP_DEBUG
   printf("subgraph: ");
   graph_printInfo(subgraph);
#endif

   SCIP_CALL( substpsolver_solve(scip, substp, &success) );
   assert(success);

   /* fix solution in original graph */
   SCIP_CALL( decomposeExactFixSol(scip, subinout, substp, orggraph, redbase, nelims) );

   graph_subinoutClean(scip, subinout);
   substpsolver_free(scip, &substp);

   return SCIP_OKAY;
}


/** tries to solve subproblem */
static
SCIP_RETCODE decomposeExactSubTry(
   SCIP*                 scip,               /**< SCIP data structure */
   const BIDECOMP*       bidecomp,           /**< all-components storage */
   int                   compindex,          /**< component index */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  nelims              /**< number of eliminations or NULL */
   )
{
   GRAPH* subgraph;
   SUBINOUT* subinout = bidecomp->subinout;
   int subroot;

   assert(redbase->bidecompparams);

   if( bidecomposition_componentIsTrivial(bidecomp, compindex) )
   {
      SCIPdebugMessage("component is trivial, returning \n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("(depth %d) solve component %d to optimality \n", redbase->bidecompparams->depth, compindex);

   bidecomposition_markSub(bidecomp, compindex, g);
   SCIP_CALL( graph_subgraphExtract(scip, g, subinout, &subgraph) );

#ifdef SCIP_DEBUG
   SCIPdebugMessage("original nodes of connected components: \n");

   for( int i = 0; i < orggraph->knots; i++ )
   {
      if( orggraph->mark[i] )
         graph_knot_printInfo(orggraph, i);
   }
   #endif

   SCIP_CALL( bidecomposition_getMarkedSubRoot(scip, bidecomp, g, subgraph, &subroot) );
   subgraph->source = subroot;

   SCIP_CALL( decomposeExactSubDoIt(scip, bidecomp, g, subgraph, redbase, nelims) );

   return SCIP_OKAY;
}


/** is promising? */
static
SCIP_Bool decomposeIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const BIDECPARAMS*    bidecompparams,     /**< bidecomposition parameters */
   const BIDECOMP*       bidecomp            /**< bidecomposition data structure */
   )
{
   const SCIP_Real mincompratio = bidecompparams->depth == 0 ? BIDECOMP_MINMAXCOMPRATIO_FIRST : BIDECOMP_MINMAXCOMPRATIO;
   const SCIP_Real maxratio = bidecomposition_getMaxcompNodeRatio(bidecomp);

   assert(GT(maxratio, 0.0));

   if( bidecompparams->depth > 1 && g->knots < BIDECOMP_MINNODES )
      return FALSE;

   return (maxratio < mincompratio);
}


/** is decomposition with (exact) solution of small components promising? */
static
SCIP_Bool decomposePartialIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const REDBASE*        redbase,            /**< reduction stuff */
   const BIDECOMP*       bidecomp            /**< bidecomposition data structure */
   )
{
   SCIP_Real maxratio;

   /* NOTE: we don't want to solve exactly in recombination heuristic etc. */
   if( !redbase->redparameters->userec )
      return FALSE;

   maxratio = bidecomposition_getMaxcompNodeRatio(bidecomp);
   assert(GT(maxratio, 0.0));

   if( g->knots < BIDECOMP_MINNODES_PARTIAL )
      return FALSE;

   return (maxratio < BIDECOMP_MINMAXCOMPRATIO_PARTIAL);
}


/** reduced biconnected components separately */
static
SCIP_RETCODE decomposeReduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   BIDECOMP*             bidecomp,           /**< bidecomposition data structure */
   REDBASE*              redbase             /**< reduction stuff */
   )
{
   REDSOL* redsol = redbase->redsol;

   SCIP_CALL( bidecomposition_initSubInOut(scip, g, bidecomp) );
   SCIP_CALL( reduce_solLevelAdd(scip, g, redsol) );
   redbase->bidecompparams->depth++;

   /* reduce each biconnected component individually */
   for( int i = 0; i < bidecomp->nbicomps; i++ )
   {
      SCIP_CALL( decomposeReduceSub(scip, bidecomp, i, g, redbase) );
   }

   redbase->bidecompparams->depth--;
   /* NOTE: also removes level */
   reduce_solLevelTopFinalize(scip, g, redsol);

   return SCIP_OKAY;
}


/** solves smaller biconnected components to optimality */
static
SCIP_RETCODE decomposePartialExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             maxcompratio,       /**< max nodes ratio */
   GRAPH*                g,                  /**< graph data structure */
   BIDECOMP*             bidecomp,           /**< bi-decomposition */
   int*                  solnode,            /**< solution nodes or NULL */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  nelims              /**< number of eliminations or NULL */
   )
{
#ifdef BENCH_SEPA_PARTIAL
   static double totalTime = 0.0;
   const double startTime = (double) clock() / (double) CLOCKS_PER_SEC;
   double endTime;
#endif

   SCIP_CALL( bidecomposition_initSubInOut(scip, g, bidecomp) );
   SCIP_CALL( graph_subinoutActivateEdgeMap(g, bidecomp->subinout) );
   graph_subinoutActivateNewHistory(bidecomp->subinout);

   SCIPdebugMessage("solving problem by partial exact decomposition at depth %d (%d components) \n",
         redbase->bidecompparams->depth, bidecomp->nbicomps);

   /* solve small biconnected component individually */
   for( int i = 0; i < bidecomp->nbicomps; i++ )
   {
      const SCIP_Real nodesratio = bidecomposition_getCompNodeRatio(bidecomp, i);
      SCIPdebugMessage("nodes ratio of component %d: %f \n", i, nodesratio);

      if( nodesratio > maxcompratio )
      {
         SCIPdebugMessage("component is too large, skipping component \n");
         continue;
      }

      SCIP_CALL( decomposeExactSubTry(scip, bidecomp, i, g, redbase, nelims) );
   }

   /* NOTE: solution edges have been fixed to 0 before */
   SCIP_CALL( reduce_contract0Edges(scip, g, solnode, TRUE) );

#ifdef BENCH_SEPA_PARTIAL
   endTime = (double) clock() / (double) CLOCKS_PER_SEC;
   totalTime += endTime - startTime;
   printf("TIME FOR PARTIAL EXACT DECOMPISITION: %f (total: %f) \n", endTime - startTime, totalTime);
#endif

   return SCIP_OKAY;
}


/** tries to solve or reduce biconnected components separately */
static
SCIP_RETCODE decomposeExec(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  solnode,            /**< solution nodes or NULL */
   SCIP_Bool*            wasDecomposed       /**< performed recursive reduction? */
   )
{
   BIDECOMP* bidecomp;

   /* NOTE: does not work because we do node contractions when reinserting the subgraph,
    * and because order of nodes is changed in subbgraph */
   assert(redbase->solnode == NULL && "not supported");
   assert(redbase->bidecompparams);
   assert(redbase->bidecompparams->depth < redbase->bidecompparams->maxdepth);
   assert(graph_valid(scip, g));
   assert(*wasDecomposed == FALSE);

   SCIP_CALL( bidecomposition_init(scip, cutnodes, g, &bidecomp) );

   if( decomposeIsPromising(g, redbase->bidecompparams, bidecomp) )
   {
      SCIP_CALL( decomposeReduce(scip, g, bidecomp, redbase) );
      *wasDecomposed = TRUE;
   }
   else if( decomposePartialIsPromising(g, redbase, bidecomp) )
   {
      SCIP_CALL( decomposePartialExact(scip, BIDECOMP_MAXCOMPRATIO_PARTIAL, g, bidecomp, solnode, redbase, NULL) );
   }

   bidecomposition_free(scip, &bidecomp);
   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** solves (not too large) biconnected component separately */
static
SCIP_RETCODE decomposeExecExact(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction stuff */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims              /**< number of eliminations */
   )
{
   BIDECOMP* bidecomp;

   /* NOTE: does not work because we do node contractions when reinserting the subgraph,
    * and because order of nodes is changed in subbgraph */
   assert(redbase->solnode == NULL && "not supported");
   assert(redbase->bidecompparams);
   assert(graph_valid(scip, g));

   SCIP_CALL( bidecomposition_init(scip, cutnodes, g, &bidecomp) );

   if( decomposePartialIsPromising(g, redbase, bidecomp) )
   {
      SCIP_CALL( decomposePartialExact(scip, BIDECOMP_MAXCOMPRATIO_AGGRESSIVE, g, bidecomp, solnode, redbase, nelims) );
   }

   bidecomposition_free(scip, &bidecomp);
   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** removes non-necessary bi-connected components and creates new terminals from cut-nodes */
SCIP_RETCODE reduce_nonTerminalComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   const CUTNODES*       cutnodes,           /**< cut nodes */
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


/** decomposition into biconnected components and recursive reduction */
SCIP_RETCODE reduce_bidecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   int*                  solnode,            /**< solution nodes or NULL */
   SCIP_Bool*            wasDecomposed       /**< performed recursive reduction? */
   )
{
   CUTNODES* cutnodes;

   assert(scip && g && redbase && wasDecomposed);
   assert(graph_typeIsSpgLike(g) && "only SPG decomposition supported yet");

   *wasDecomposed = FALSE;

   if( g->terms == 1 )
      return SCIP_OKAY;

   graph_mark(g);

   if( !bidecomposition_isPossible(g) )
   {
      SCIPdebugMessage("graph is too large...don't decompose \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( bidecomposition_cutnodesInit(scip, g, &cutnodes) );
   bidecomposition_cutnodesCompute(g, cutnodes);

   if( cutnodes->biconn_ncomps > 0 )
   {
      int dummy = 0;
      /* get rid of non-required biconnected components (without terminals) */
      SCIP_CALL( reduce_nonTerminalComponents(scip, cutnodes, g, reduce_solGetOffsetPointer(redbase->redsol), &dummy) );

      /* decompose and reduce recursively? */
      SCIP_CALL( decomposeExec(scip, cutnodes, g, redbase, solnode, wasDecomposed) );
   }

   bidecomposition_cutnodesFree(scip, &cutnodes);

   return SCIP_OKAY;
}


/** solves smaller connected components to optimality  */
SCIP_RETCODE reduce_bidecompositionExact(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims              /**< number of eliminations */
   )
{
   CUTNODES* cutnodes;

   assert(scip && g && redbase && nelims);
   assert(graph_typeIsSpgLike(g) && "only SPG decomposition supported yet");

   *nelims = 0;

   if( g->terms == 1 )
      return SCIP_OKAY;

   graph_mark(g);

   if( !bidecomposition_isPossible(g) )
   {
      SCIPdebugMessage("graph is too large...don't decompose \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( bidecomposition_cutnodesInit(scip, g, &cutnodes) );
   bidecomposition_cutnodesCompute(g, cutnodes);

   if( cutnodes->biconn_ncomps > 0 )
   {
      /* get rid of non-required biconnected components (without terminals) */
      SCIP_CALL( reduce_nonTerminalComponents(scip, cutnodes, g, reduce_solGetOffsetPointer(redbase->redsol), nelims) );

      SCIP_CALL( decomposeExecExact(scip, cutnodes, g, redbase, solnode, nelims) );
   }

   bidecomposition_cutnodesFree(scip, &cutnodes);

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
   CUTNODES* cutnodes;

   assert(scip && g && nelims);
   graph_mark(g);

   *nelims = 0;

   if( !bidecomposition_isPossible(g) )
   {
      SCIPdebugMessage("graph is too large...don't decompose \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( bidecomposition_cutnodesInit(scip, g, &cutnodes) );
   bidecomposition_cutnodesCompute(g, cutnodes);

   if( cutnodes->biconn_ncomps > 0 )
   {
      SCIP_CALL( reduce_nonTerminalComponents(scip, cutnodes, g, fixedp, nelims) );
   }

   bidecomposition_cutnodesFree(scip, &cutnodes);

   return SCIP_OKAY;
}
