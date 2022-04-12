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

/**@file   dpterms_core.c
 * @brief  Core of dynamic programming solver for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Contains core methods.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include "scip/scipdefplugins.h"
#include "scip/rbtree.h"
#include "dpterms.h"
#include "dptermsinterns.h"
#include "stpbitset.h"
#include "stpvector.h"
#include "stpprioqueue.h"


/*
 * Data structures
 */

/** helper data */
typedef struct trace_triplet
{
   SCIP_Real             bdist_local;        /**< temporary tree bottleneck distance */
   SCIP_Real             bdist_global;       /**< tree bottleneck distance */
   int                   index;              /**< index */
} TTRIPLET;


/** saves some data updated in every iteration */
typedef struct dynamic_programming_iterator
{
   DPSUBSOL*             dpsubsol;
   STP_Vectype(int)      stack;              /**< general purpose stack */
   STP_Vectype(TTRIPLET) tripletstack;       /**< special stack */
   STP_Vectype(SOLTRACE) sol_traces;         /**< traces of current sub-solution */
   STP_Bitset            sol_termbits;         /**< marks terminals of sub-solution */
   STP_Vectype(SOLTRACE) valid_traces;       /**< traces of valid extension */
   STP_Bitset            valid_bitset;       /**< marks valid roots */
   SCIP_Real*            nodes_dist;         /**< weight of sub-ST rooted at node */
#ifdef STP_DPTERM_USEDA
   SCIP_Real*            nodes_reddist;      /**< reduced cost weight of sub-ST rooted at node */
#endif
   SCIP_Real*            nodes_ub;           /**< upper bounds for rule-out */
   int*                  nodes_previdx0;     /**< predecessor NOTE: with shift! */
   int*                  nodes_previdx1;     /**< predecessor */
   SCIP_Bool*            nodes_isValidRoot;  /**< is node a valid root? */
   int                   nnodes;             /**< number of nodes */
   int                   sol_nterms;         /**< popcount */
   int                   extterm;             /**< extension terminal */
} DPITER;


/*
 * Local methods
 */


/** prints separator nodes in SCIP_DEBUG mode */
static
void debugPrintSeparator(
   SCIP_Real             maxvaliddist,       /**< maximum value distance */
   const DPITER*         dpiterator          /**< iterator */
)
{
   assert(dpiterator);

#ifdef SCIP_DEBUG
   {
      const SCIP_Real* const nodes_ub = dpiterator->nodes_ub;
      const SCIP_Real* const nodes_dist = dpiterator->nodes_dist;
      const int nnodes = dpiterator->nnodes;

      SCIPdebugMessage("separator nodes: \n");

      for( int i = 0; i < nnodes; i++ )
      {
         if( GT(nodes_ub[i], 0.0) )
         {
            printf("%d (nodes_ub=%f) \n", i, nodes_ub[i]);

         }
         else if( LE(nodes_dist[i], maxvaliddist) )
         {
            printf("%d (dist=%f) \n", i, nodes_dist[i]);
         }
      }
   }
#endif

}



/** assembles final solution nodes */
static
STP_Vectype(int) getSolnodesFinal(
   SCIP*                 scip,               /**< SCIP data structure */
   DPMISC*               dpmisc,             /**< misc */
   DPITER*               dpiterator          /**< iterator */
)
{
   STP_Vectype(int) solnodes = NULL;
   STP_Vectype(int) stack = dpiterator->stack;
   STP_Vectype(SOLTRACE) traces_all = dpmisc->global_traces;

   assert(dpmisc->opt_root);
   StpVecClear(stack);

   SCIPdebugMessage("add solution node %d \n", traces_all[dpmisc->opt_root].root);
   StpVecPushBack(scip, solnodes, traces_all[dpmisc->opt_root].root);

   StpVecPushBack(scip, stack, dpmisc->opt_root);

   if( dpmisc->opt_prev[0] != -1 )
      StpVecPushBack(scip, stack, dpmisc->opt_prev[0]);

   if( dpmisc->opt_prev[1] != -1 )
      StpVecPushBack(scip, stack, dpmisc->opt_prev[1]);

   while( StpVecGetSize(stack) > 0 )
   {
      const int i = stack[StpVecGetSize(stack) - 1];
      const int prev0 = traces_all[i].prevs[0];

      StpVecPopBack(stack);

      if( prev0 != -1 )
      {
         StpVecPushBack(scip, stack, prev0);

         if( traces_all[i].prevs[1] != -1 )
         {
            StpVecPushBack(scip, stack, traces_all[i].prevs[1]);
         }
         else
         {
            SCIPdebugMessage("add solution node %d \n", traces_all[prev0].root);
            StpVecPushBack(scip, solnodes, traces_all[prev0].root);
         }
      }
   }

   dpiterator->stack = stack;

   return solnodes;
}


/** helper */
static inline
SCIP_Bool nodeIsNonSolTerm(
   STP_Bitset            sol_bitset,         /**< bitset */
   const int*            nodes_termId,       /**< ID  */
   int                   node                /**< node to check */
)
{
   return( nodes_termId[node] != -1 && !stpbitset_bitIsTrue(sol_bitset, nodes_termId[node]) );
}


/** helper */
static
SCIP_Bool allExtensionsAreInvalid(
   const GRAPH*          graph,              /**< graph */
   const DPSOLVER*       dpsolver,           /**< solver */
   const DPITER*         dpiterator          /**< iterator */
)
{
   const SCIP_Real* const nodes_ub = dpiterator->nodes_ub;
   STP_Bitset sol_termbits = dpiterator->sol_termbits;
   const int* const terminals = dpsolver->dpgraph->terminals;
   const int nterms = graph->terms;

   for( int i = 0; i < nterms; i++  )
   {
      if( !stpbitset_bitIsTrue(sol_termbits, i)
         && GT(nodes_ub[terminals[i]], 0.0) )
      {
         SCIPdebugMessage("terminal %d has positive UB, all extension are invalid! \n", terminals[i]);
         return TRUE;
      }
   }

   return FALSE;
}


/** gets ordered root indices according to the solution costs */
static
SCIP_RETCODE getOrderedRootIndices(
   SCIP*                 scip,               /**< SCIP data structure */
   DPITER*               dpiterator,         /**< iterator */
   int*                  roots_indices       /**< to initialize */
)
{
   const int nnodes = dpiterator->nnodes;
   SCIP_Real* roots_cost;
   STP_Vectype(SOLTRACE) soltraces = dpiterator->sol_traces;
   const int nsolcands = StpVecGetSize(dpiterator->sol_traces);

   SCIP_CALL( SCIPallocBufferArray(scip, &roots_cost, nnodes) );

   for( int i = 0; i < nsolcands; i++ )
   {
      roots_indices[i] = i;
      roots_cost[i] = -soltraces[i].cost;
   }

   SCIPsortDownRealInt(roots_cost, roots_indices, nsolcands);

   SCIPfreeBufferArray(scip, &roots_cost);

#ifndef NDEBUG
   for( int i = 1; i < nsolcands; i++ )
   {
      const int ind = roots_indices[i];
      const int ind_prev = roots_indices[i - 1];
      assert(LE(soltraces[ind_prev].cost, soltraces[ind].cost));
   }
#endif

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE dpiterInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPITER**              dpiterator          /**< to initialize */
)
{
   DPITER* iter;
   const int nnodes = graph_get_nNodes(graph);

   SCIP_CALL( SCIPallocMemory(scip, dpiterator) );
   iter = *dpiterator;

   iter->stack = NULL;
   iter->tripletstack = NULL;
   iter->sol_traces = NULL;
   iter->sol_termbits = NULL;
   iter->valid_traces = NULL;
   iter->valid_bitset = NULL;
   iter->nnodes = nnodes;
   iter->sol_nterms = -1;
   iter->extterm = -1;

   /* NOTE: could be any positive number */
   StpVecReserve(scip, iter->stack, 2);
   StpVecReserve(scip, iter->valid_traces, 2);

#ifdef STP_DPTERM_USEDA
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_reddist), nnodes) );
#endif
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_dist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_ub), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_previdx0), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_previdx1), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(iter->nodes_isValidRoot), nnodes) );

   return SCIP_OKAY;
}


/** frees */
static
void dpiterFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPITER**              dpiterator          /**< to free */
)
{
   DPITER* iter = *dpiterator;

   assert(!iter->sol_traces && !iter->sol_termbits);
   StpVecFree(scip, iter->valid_traces);
   StpVecFree(scip, iter->tripletstack);
   StpVecFree(scip, iter->stack);

   SCIPfreeMemoryArray(scip, &(iter->nodes_isValidRoot));
   SCIPfreeMemoryArray(scip, &(iter->nodes_previdx1));
   SCIPfreeMemoryArray(scip, &(iter->nodes_previdx0));
   SCIPfreeMemoryArray(scip, &(iter->nodes_ub));
   SCIPfreeMemoryArray(scip, &(iter->nodes_dist));
#ifdef STP_DPTERM_USEDA
   SCIPfreeMemoryArray(scip, &(iter->nodes_reddist));
#endif

   SCIPfreeMemory(scip, dpiterator);
}


/** sets arrays to default values */
static
void dpiterSetDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   DPITER*               dpiterator          /**< to update */
)
{
   SCIP_Real* RESTRICT nodes_dist = dpiterator->nodes_dist;
   SCIP_Real* RESTRICT nodes_ub = dpiterator->nodes_ub;
   int* RESTRICT nodes_pred1 = dpiterator->nodes_previdx0;
   int* RESTRICT nodes_pred2 = dpiterator->nodes_previdx1;
   SCIP_Bool* RESTRICT nodes_isValidRoot = dpiterator->nodes_isValidRoot;
   const int nnodes = dpiterator->nnodes;

#ifdef STP_DPTERM_USEDA
   SCIP_Real* RESTRICT nodes_reddist = dpiterator->nodes_reddist;

   for( int i = 0; i < nnodes; i++ )
      nodes_reddist[i] = 0.0;
#endif

   dpiterator->extterm = -1;


   for( int i = 0; i < nnodes; i++ )
      nodes_dist[i] = FARAWAY;

   for( int i = 0; i < nnodes; i++ )
      nodes_ub[i] = FARAWAY;

   for( int i = 0; i < nnodes; i++ )
      nodes_pred1[i] = -1;

   for( int i = 0; i < nnodes; i++ )
      nodes_pred2[i] = -1;

   BMSclearMemoryArray(nodes_isValidRoot, nnodes);
}


/** gets solution from tree */
static
void dpiterPopSol(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< to update */
)
{
   DPSUBSOL* subsol;

   stpprioqueue_deleteMin((void**) &(dpiterator->sol_termbits), &(dpiterator->sol_nterms), dpsolver->solpqueue);

   if( findSubsol(dpsolver->soltree_root, dpiterator->sol_termbits, &subsol) == 0 )
   {
      dpiterator->sol_traces = subsol->traces;
      SCIPdebugMessage("number of traces: %d \n", StpVecGetSize(dpiterator->sol_traces));

      SCIPrbtreeDelete(&(dpsolver->soltree_root), subsol);

      assert(stpbitset_areEqual(dpiterator->sol_termbits, subsol->bitkey));
      stpbitset_free(scip, &(subsol->bitkey));
   }
   else
   {
      assert(0 && "should never happen");
   }

   assert(stpbitset_getPopcount(dpiterator->sol_termbits) == dpiterator->sol_nterms);
   dpiterator->dpsubsol = subsol;
}


/** updates */
static
void dpiterGetNextSol(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< to update */
)
{
   assert(!dpiterator->sol_termbits);
   assert(!dpiterator->valid_bitset); /* NOTE: should be moved */

   dpiterSetDefault(scip, dpiterator);
   dpiterPopSol(scip, dpsolver, dpiterator);

#ifdef SCIP_DEBUG
      SCIPdebugMessage("processing solution: \n");
      stpbitset_print(dpiterator->sol_termbits);
      SCIPdebugMessage("terminals: \n");
      for( int i = 0; i < dpsolver->dpgraph->nterms; i++ )
      {
         if( stpbitset_bitIsTrue(dpiterator->sol_termbits, i) )
            printf("%d  ", dpsolver->dpgraph->terminals[i]);
      }
      printf(" \n");
#endif
}


/** (Implicitly) constructs sub-Steiner trees */
static
SCIP_RETCODE subtreesBuild(
   SCIP*                 scip,               /**< SCIP data structure */
   DPMISC*               dpmisc,             /**< misc */
   DPITER*               dpiterator          /**< to update */
)
{
   SCIP_Real* RESTRICT nodes_dist = dpiterator->nodes_dist;
   SCIP_Real* RESTRICT nodes_ub = dpiterator->nodes_ub;
   int* RESTRICT nodes_previdx0 = dpiterator->nodes_previdx0;
   int* RESTRICT nodes_previdx1 = dpiterator->nodes_previdx1;
   SCIP_Bool* RESTRICT nodes_isValidRoot = dpiterator->nodes_isValidRoot;
   STP_Vectype(TTRIPLET) tripletstack = dpiterator->tripletstack;
   STP_Vectype(SOLTRACE) soltraces = dpiterator->sol_traces;
   STP_Vectype(SOLTRACE) global_traces = dpmisc->global_traces;
   int* nodes_nvisits;
   int* roots_indices;
   const int nsolcands = StpVecGetSize(soltraces);
   const int nnodes = dpiterator->nnodes;
   int nvalidroots = 0;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &nodes_nvisits, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roots_indices, nnodes) );
   SCIP_CALL( getOrderedRootIndices(scip, dpiterator, roots_indices) );

   SCIPdebugMessage("building Steiner trees from %d candidates \n", nsolcands);

   for( int i = 0; i < nsolcands; i++ )
   {
      const SOLTRACE sol_trace = soltraces[roots_indices[i]];
      const int sol_root = sol_trace.root;

      // todo GE enough?
      if( GT(sol_trace.cost, nodes_dist[sol_root]) )
         continue;

      SCIPdebugMessage("building solution with root %d and prev0=%d prev1=%d \n", sol_root,
            sol_trace.prevs[0], sol_trace.prevs[1]);

#ifdef STP_DPTERM_USEDA
      dpiterator->nodes_reddist[sol_root] = sol_trace.redcost;
#endif
      nodes_isValidRoot[sol_root] = TRUE;
      nodes_dist[sol_root] = sol_trace.cost;
      nodes_ub[sol_root] = 0.0;
      nodes_nvisits[sol_root]++;
      nvalidroots++;
      nodes_previdx0[sol_root] = sol_trace.prevs[0];
      nodes_previdx1[sol_root] = sol_trace.prevs[1];

      assert(sol_trace.prevs[0] != -1 || sol_trace.prevs[1] == -1);

      if( sol_trace.prevs[0] == -1 )
         continue;

      assert(StpVecGetSize(tripletstack) == 0);
      StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, 0.0, sol_trace.prevs[0]}));
      StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, 0.0, sol_trace.prevs[1]}));

      while( StpVecGetSize(tripletstack) > 0 )
      {
         TTRIPLET triplet = tripletstack[StpVecGetSize(tripletstack) - 1];
         SOLTRACE parent_trace = global_traces[triplet.index];
         const int previdx0 = parent_trace.prevs[0];

         StpVecPopBack(tripletstack);

         if( previdx0 != -1 )
         {
            const SCIP_Real bdist_local = triplet.bdist_local;
            const SCIP_Real bdist_global = triplet.bdist_global;
            const int previdx1 = parent_trace.prevs[1];

            assert(previdx0 >= 0);

            /* merged solution? */
            if( previdx1 != -1 )
            {
               assert(previdx0 != -1);

               StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, bdist_global, previdx0}));
               StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, bdist_global, previdx1}));

               SCIPdebugMessage("...merged solution from prev0=%d prev1=%d", previdx0, previdx1);
            }
            else
            {
               const int curr_root = global_traces[previdx0].root;
               const SCIP_Real curr_bdist_local = bdist_local + parent_trace.cost - global_traces[previdx0].cost;
               const SCIP_Real curr_bdist_global = MAX(bdist_global, curr_bdist_local);

               SCIPdebugMessage("at pred. solution with index=%d root=%d \n", previdx0, global_traces[previdx0].root);

               assert(GE(parent_trace.cost - global_traces[previdx0].cost, 0.0));

               if( LT(sol_trace.cost, nodes_dist[curr_root]) )
               {
                  nodes_dist[curr_root] = sol_trace.cost;
#ifdef STP_DPTERM_USEDA
                  dpiterator->nodes_reddist[curr_root] = sol_trace.redcost;
#endif
               }

               nodes_nvisits[curr_root]++;
               nodes_ub[curr_root] = MIN(nodes_ub[curr_root], curr_bdist_global);

               StpVecPushBack(scip, tripletstack, ((TTRIPLET) {curr_bdist_local, curr_bdist_global, previdx0}));
            }
         }
      }
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_nvisits[i] < nvalidroots )
         nodes_ub[i] = 0.0;
   }

   SCIPfreeBufferArray(scip, &roots_indices);
   SCIPfreeBufferArray(scip, &nodes_nvisits);
   dpiterator->tripletstack = tripletstack;

   return SCIP_OKAY;
}


/** extends sub-Steiner trees */
static
SCIP_RETCODE subtreesExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   const CSR* const csr = graph->csr_storage;
   DHEAP* const dheap = dpsolver->dheap;
   const SCIP_Real* const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   int* RESTRICT nodes_previdx0 = dpiterator->nodes_previdx0;
   int* RESTRICT nodes_previdx1 = dpiterator->nodes_previdx1;
   SCIP_Real* RESTRICT nodes_dist = dpiterator->nodes_dist;
   SCIP_Bool* RESTRICT nodes_isValidRoot = dpiterator->nodes_isValidRoot;
   const int nnodes = dpiterator->nnodes;
   STP_Bitset sol_termbits = dpiterator->sol_termbits;
   const int* const nodes_termId = dpsolver->dpgraph->nodes_termId;
   int* terms_adjCount;
   const SCIP_Bool breakEarly = (graph->terms - dpiterator->sol_nterms > 1);
   const int* const grad = graph->grad;
#ifdef STP_DPTERM_USEDA
   SCIP_Real* RESTRICT nodes_reddist = dpiterator->nodes_reddist;
   const SCIP_Real* const csr_redcost = dpsolver->dpredcosts->csr_redcosts;
#endif

   assert(nnodes == graph->knots);
   assert(dpiterator->extterm == -1);
   assert(dpiterator->sol_nterms > 0 && dpiterator->sol_nterms <= graph->terms);

   graph_heap_clean(TRUE, dheap);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &terms_adjCount, graph->terms) );

   for( int i = 0; i < nnodes; i++ )
   {
      if( LT(nodes_dist[i], FARAWAY) )
         graph_heap_correct(i, nodes_dist[i], dheap);
   }

   SCIPdebugMessage("starting Dijkstra with %d roots \n", dheap->size);

   /* run Dijkstra */
   while( dheap->size > 0 )
   {
      int k;
      SCIP_Real k_dist;

      graph_heap_deleteMin(&k, &k_dist, dheap);
      assert(EQ(nodes_dist[k], k_dist));

      SCIPdebugMessage("updated node %d: dist=%f \n", k, k_dist);

      if( nodeIsNonSolTerm(sol_termbits, nodes_termId, k) )
      {
         assert(dpiterator->extterm == -1);
         dpiterator->extterm = k;
         break;
      }
      else
      {
         const int k_start = start_csr[k];
         const int k_end = start_csr[k + 1];
         const SCIP_Bool k_isValid = nodes_isValidRoot[k];

         for( int e = k_start; e != k_end; e++ )
         {
            const int m = head_csr[e];
            const SCIP_Real distnew = k_dist + cost_csr[e];

            if( LT(distnew, nodes_dist[m]) )
            {
               nodes_isValidRoot[m] = k_isValid;
               nodes_dist[m] = distnew;
               nodes_previdx0[m] = k;
               nodes_previdx1[m] = -1;
#ifdef STP_DPTERM_USEDA
               nodes_reddist[m] = nodes_reddist[k] + csr_redcost[e];
#endif

               graph_heap_correct(m, distnew, dheap);
            }
            else if( EQ(distnew, nodes_dist[m]) && !k_isValid )
            {
               nodes_isValidRoot[m] = FALSE;
            }

            if( nodeIsNonSolTerm(sol_termbits, nodes_termId, m) && breakEarly )
            {
               assert(nodes_termId[m] >= 0);
               terms_adjCount[nodes_termId[m]]++;

               /* all neighbors hit? */
               if( terms_adjCount[nodes_termId[m]] == grad[m] )
               {
                  graph_heap_clean(FALSE, dheap);
                  assert(dheap->size == 0);
                  assert(dpiterator->extterm == -1);
                  dpiterator->extterm = m;
                  break;
               }
            }
         }
      }
   }

   SCIPdebugMessage("...ending extension with root %d \n", dpiterator->extterm);

   SCIPfreeBufferArray(scip, &terms_adjCount);

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE dpiterAddNewPrepare(
   SCIP*                 scip,               /**< SCIP data structure */
   DPMISC*               dpmisc,             /**< DP misc data structure */
   DPITER*               dpiterator,         /**< iterator */
   SCIP_Bool*            hasExtension        /**< extensions existing? */
)
{
   STP_Vectype(SOLTRACE) valid_traces = dpiterator->valid_traces;
   const SCIP_Bool* const nodes_isValidRoot = dpiterator->nodes_isValidRoot;
   const SCIP_Real* const nodes_dist = dpiterator->nodes_dist;
   const int* const nodes_previdx0 = dpiterator->nodes_previdx0;
   const int* const nodes_previdx1 = dpiterator->nodes_previdx1;
   const int nnodes = dpiterator->nnodes;
   STP_Bitset valid_bits = stpbitset_new(scip, nnodes);
   int* rootmap;
   const int nall = dpmisc->global_size;

   StpVecClear(valid_traces);

   for( int i = 0; i < nnodes; i++ )
   {
      if( nodes_isValidRoot[i] )
      {
         SCIPdebugMessage("valid trace node %d \n", i);

         assert(GE(nodes_dist[i], 0.0) && LT(nodes_dist[i], FARAWAY));

#ifdef STP_DPTERM_USEDA
         StpVecPushBack(scip, valid_traces,
         ((SOLTRACE) { {nodes_previdx0[i], nodes_previdx1[i]}, nodes_dist[i], dpiterator->nodes_reddist[i], i }) );
#else
         StpVecPushBack(scip, valid_traces,
         ((SOLTRACE) { {nodes_previdx0[i], nodes_previdx1[i]}, nodes_dist[i], i }) );
#endif

         stpbitset_setBitTrue(valid_bits, i);
      }
   }

   /* no valid extensions? */
   if( 0 == StpVecGetSize(valid_traces) )
   {
      SCIPdebugMessage("no valid extensions! \n");
      *hasExtension = FALSE;
      stpbitset_free(scip, &valid_bits);
      return SCIP_OKAY;
   }

   *hasExtension = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &rootmap, nnodes) );
#ifndef NDEBUG
   for( int i = 0; i < nnodes; i++ )
      rootmap[i] = -1;
#endif
   for( int i = 0; i < StpVecGetSize(valid_traces); i++ )
   {
      assert(valid_traces[i].root >= 0);
      rootmap[valid_traces[i].root] = nall + i;
   }

   for( int i = 0; i < StpVecGetSize(valid_traces); i++ )
   {
      if( valid_traces[i].prevs[0] != -1 && valid_traces[i].prevs[1] == -1 )
      {
         assert(rootmap[valid_traces[i].prevs[0]] >= 0);
         valid_traces[i].prevs[0] = rootmap[valid_traces[i].prevs[0]];
      }
   }

   assert(!dpiterator->valid_bitset);
   dpiterator->valid_bitset = valid_bits;
   dpiterator->valid_traces = valid_traces;
   SCIPfreeBuffer(scip, &rootmap);

   return SCIP_OKAY;
}



/** helper */
static
STP_Vectype(SOLTRACE) combineTraces(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Vectype(SOLTRACE) traces1,            /**< ordered traces */
   int                   prevoffset1,
   SOLTRACE*             traces2,            /**< ordered traces */
   int                   ntraces2,           /**< number of traces */
   int                   prevoffset2
)
{
   STP_Vectype(SOLTRACE) combined = NULL;
   int pos1 = 0;
   int pos2 = 0;
   const int ntraces1 = StpVecGetSize(traces1);
   assert(ntraces1 > 0 && ntraces2 > 0);

   while( pos1 < ntraces1 && pos2 < ntraces2 )
   {
      assert(pos1 == ntraces1 - 1 || traces1[pos1].root < traces1[pos1 + 1].root);
      assert(pos2 == ntraces2 - 1 || traces2[pos2].root < traces2[pos2 + 1].root);

      if( traces1[pos1].root == traces2[pos2].root )
      {
#ifdef STP_DPTERM_USEDA
         StpVecPushBack(scip, combined, ((SOLTRACE)
           {  {prevoffset1 + pos1, prevoffset2 + pos2},
              traces1[pos1].cost + traces2[pos2].cost,
              traces1[pos1].redcost + traces2[pos2].redcost,
              traces1[pos1].root
           }  ) );
#else
         StpVecPushBack(scip, combined, ((SOLTRACE)
           {  {prevoffset1 + pos1, prevoffset2 + pos2},
              traces1[pos1].cost + traces2[pos2].cost,
              traces1[pos1].root
           }  ) );
#endif

         pos1++;
         pos2++;
      }
      else if( traces1[pos1].root < traces2[pos2].root )
      {
         pos1++;
      }
      else
      {
         assert(traces1[pos1].root > traces2[pos2].root);
         pos2++;
      }
   }

   assert(combined);

#ifndef NDEBUG
   for( int i = 1; i < StpVecGetSize(combined); i++ )
      assert(combined[i - 1].root <= combined[i].root);
#endif

   return combined;
}


/** updates existing subsol */
static
void subsolUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Vectype(SOLTRACE) combined_traces,     /**< for updating */
   DPITER*               dpiterator,         /**< iterator */
   DPSUBSOL*             dpsubsol            /**< to be updated */
)
{
   STP_Vectype(SOLTRACE) subsol_traces = dpsubsol->traces;
   STP_Vectype(SOLTRACE) updated_traces = NULL;
   int pos1 = 0;
   int pos2 = 0;
   const int nsubsol_traces =  StpVecGetSize(subsol_traces);
   const int ncombined_traces = StpVecGetSize(combined_traces);
   const int nnodes = dpiterator->nnodes;

   assert(subsol_traces && combined_traces);

   while( pos1 < nsubsol_traces || pos2 < ncombined_traces )
   {
      const int root1 = (pos1 == nsubsol_traces) ? nnodes : subsol_traces[pos1].root;
      const int root2 = (pos2 == ncombined_traces) ? nnodes : combined_traces[pos2].root;
      assert(root1 < nnodes || root2 < nnodes);

      if( root1 == root2 )
      {
         if( LE(subsol_traces[pos1].cost, combined_traces[pos2].cost) )
            StpVecPushBack(scip, updated_traces, subsol_traces[pos1]);
         else
            StpVecPushBack(scip, updated_traces, combined_traces[pos2]);
         pos1++;
         pos2++;
      }
      else if( root1 < root2 )
      {
         StpVecPushBack(scip, updated_traces, subsol_traces[pos1]);
         pos1++;
      }
      else
      {
         StpVecPushBack(scip, updated_traces, combined_traces[pos2]);
         pos2++;
      }

   }
   StpVecFree(scip, subsol_traces);
   dpsubsol->traces = updated_traces;
}


/** combines new sub-Steiner trees */
static
SCIP_RETCODE combineWithIntersecting(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   int                   index,              /**< index of intersecting */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   DPSUBSOL* subsol;
   DPMISC* dpmisc = dpsolver->dpmisc;
   STP_Vectype(int) offsets = dpmisc->global_starts;
   STP_Bitset composite_termbits = dpmisc->global_termbits[index];
   SOLTRACE* composite_traces = &(dpmisc->global_traces[offsets[index]]);
   const int composite_ntraces = offsets[index + 1] - offsets[index];
   STP_Bitset combined_termbits = stpbitset_newOr(scip, composite_termbits, dpiterator->sol_termbits);
   int  subsol_pos;

   STP_Vectype(SOLTRACE) combined_traces = combineTraces(scip,
         dpiterator->valid_traces, dpmisc->global_size,
         composite_traces, composite_ntraces, offsets[index]);

   if( (subsol_pos = findSubsol(dpsolver->soltree_root, combined_termbits, &subsol)) == 0 )
   {
      assert(subsol);

      subsolUpdate(scip, combined_traces, dpiterator, subsol);
      StpVecFree(scip, combined_traces);
      stpbitset_free(scip, &combined_termbits);
   }
   else
   {
      DPSUBSOL* newsol;
      const int ncombinedterms = stpbitset_getPopcount(combined_termbits);
      if( 2 * ncombinedterms <= graph->terms )
      {
         SCIP_CALL( stpprioqueue_insert(scip, stpbitset_newCopy(scip, combined_termbits),
               ncombinedterms, dpsolver->solpqueue) );
      }

      SCIP_CALL( dpterms_dpsubsolInit(scip, &newsol) );
      newsol->bitkey = combined_termbits;
      newsol->traces = combined_traces;

      SCIPrbtreeInsert(&(dpsolver->soltree_root), subsol, subsol_pos, newsol);
   }

   return SCIP_OKAY;
}



/** updates best (global) solution */
static
void updateIncumbent(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   DPSUBSOL* dpsubsol;
   DPMISC* dpmisc = dpsolver->dpmisc;
   STP_Bitset sol_termstoggled = stpbitset_newNot(scip, dpiterator->sol_termbits,
         dpsolver->dpgraph->nterms);

   if( findSubsol(dpsolver->soltree_root, sol_termstoggled, &dpsubsol) == 0 )
   {
      STP_Vectype(SOLTRACE) valid_traces = dpiterator->valid_traces;
      STP_Vectype(SOLTRACE) toggled_traces = dpsubsol->traces;
      int pos1 = 0;
      int pos2 = 0;
      const int ntraces1 = StpVecGetSize(valid_traces);
      const int ntraces2 = StpVecGetSize(toggled_traces);
      assert(ntraces1 > 0 && ntraces2 > 0);
      assert(valid_traces && toggled_traces);

#ifdef SCIP_DEBUG
      SCIPdebugMessage("found toggled terminal bits: \n");
      stpbitset_print(dpiterator->sol_termbits);
#endif

      while( pos1 < ntraces1 && pos2 < ntraces2 )
      {
         assert(pos1 == ntraces1 - 1 || valid_traces[pos1].root < valid_traces[pos1 + 1].root);
         assert(pos2 == ntraces2 - 1 || toggled_traces[pos2].root < toggled_traces[pos2 + 1].root);

         if( valid_traces[pos1].root == toggled_traces[pos2].root )
         {
            const SCIP_Real newcost = valid_traces[pos1].cost + toggled_traces[pos2].cost;
            if( LT(newcost, dpmisc->opt_obj) )
            {
               SCIPdebugMessage("updating incumbent obj %f->%f \n", dpmisc->opt_obj, newcost);
               dpmisc->opt_obj = newcost;
               dpmisc->opt_root = dpmisc->global_size + pos1;
               dpmisc->opt_prev[0] = toggled_traces[pos2].prevs[0];
               dpmisc->opt_prev[1] = toggled_traces[pos2].prevs[1];
            }
            pos1++;
            pos2++;
         }
         else if( valid_traces[pos1].root < toggled_traces[pos2].root )
         {
            pos1++;
         }
         else
         {
            assert(valid_traces[pos1].root > toggled_traces[pos2].root);
            pos2++;
         }
      }
   }
   stpbitset_free(scip, &sol_termstoggled);
}


/** finalizes */
static
SCIP_RETCODE subtreesAddNewFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   DPMISC* dpmisc = dpsolver->dpmisc;
   const int nextensions = StpVecGetSize(dpiterator->valid_traces);

   assert(nextensions > 0);
   assert(dpiterator->sol_termbits);
   assert(dpiterator->valid_bitset);

   if( 3 * dpiterator->sol_nterms <= graph->terms )
   {
      const int nsubsets = StpVecGetSize(dpmisc->global_termbits);

      SCIP_CALL( dpterms_streeInsert(scip, stpbitset_newCopy(scip, dpiterator->sol_termbits),
            dpiterator->valid_bitset, nsubsets, dpsolver->dpstree) );
   }
   else
   {
      stpbitset_free(scip, &(dpiterator->valid_bitset));
   }

   for( int i = 0; i < nextensions; i++ )
   {
      StpVecPushBack(scip, dpmisc->global_traces, dpiterator->valid_traces[i]);
   }
   StpVecPushBack(scip, dpmisc->global_termbits, dpiterator->sol_termbits);
   StpVecPushBack(scip, dpmisc->global_termbitscount, dpiterator->sol_nterms);
   dpmisc->global_size += nextensions;
   StpVecPushBack(scip, dpmisc->global_starts, dpmisc->global_size);

   dpiterator->sol_termbits = NULL;
   dpiterator->valid_bitset = NULL;

   return SCIP_OKAY;
}

/** adds new sub-Steiner trees */
static
SCIP_RETCODE subtreesAddNew(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   DPMISC* dpmisc = dpsolver->dpmisc;
   STP_Vectype(int) intersections;
   const int nterms = graph->terms;
   SCIP_Bool hasExtensions;

   SCIP_CALL( dpiterAddNewPrepare(scip, dpmisc, dpiterator, &hasExtensions) );

   if( !hasExtensions )
   {
      return SCIP_OKAY;
   }

   intersections = dpterms_streeCollectIntersects(scip, dpiterator->sol_termbits,
         dpiterator->valid_bitset, dpsolver->dpstree);

   for( int i = 0; i < StpVecGetSize(intersections); i++ )
   {
      const int pos = intersections[i];
      assert(0 <= pos && pos < dpmisc->global_size);
      assert( stpbitset_getPopcount(dpmisc->global_termbits[pos]) == dpmisc->global_termbitscount[pos]);

      if( 2 * dpiterator->sol_nterms + dpmisc->global_termbitscount[pos] > nterms )
         continue;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("intersecting bitset: \n");
      stpbitset_print(dpmisc->global_termbits[pos]);
#endif

      SCIP_CALL( combineWithIntersecting(scip, graph, pos, dpsolver, dpiterator) );
   }

   updateIncumbent(scip, dpsolver, dpiterator);
   SCIP_CALL( subtreesAddNewFinalize(scip, graph, dpsolver, dpiterator) );

   StpVecFree(scip, intersections);
   return SCIP_OKAY;
}


/** remove non-valid sub-Steiner trees */
static
void propagateUBs(
   const GRAPH*          graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   const CSR* const csr = graph->csr_storage;
   DHEAP* const dheap = dpsolver->dheap;
   const SCIP_Real* const cost_csr = csr->cost;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   SCIP_Real* RESTRICT nodes_ub = dpiterator->nodes_ub;
   const int nnodes = dpiterator->nnodes;

   assert(nnodes == graph->knots);

   graph_heap_clean(TRUE, dheap);

   for( int i = 0; i < nnodes; i++ )
   {
      if( GT(nodes_ub[i], 0.0) )
         graph_heap_correct(i, -nodes_ub[i], dheap);
   }

   /* propagate UBs */
   while( dheap->size > 0 )
   {
      int k;
      SCIP_Real k_ub;

      graph_heap_deleteMin(&k, &k_ub, dheap);
      k_ub *= -1.0;
      assert(EQ(nodes_ub[k], k_ub));

      {
         const int k_start = start_csr[k];
         const int k_end = start_csr[k + 1];

         for( int e = k_start; e != k_end; e++ )
         {
            const int m = head_csr[e];
            const SCIP_Real ubnew = k_ub - cost_csr[e];

            if( GT(ubnew, nodes_ub[m]) )
            {
               nodes_ub[m] = ubnew;
               graph_heap_correct(m, -ubnew, dheap);
            }
         }
      }
   }
}


/** remove non-valid sub-Steiner trees */
static
SCIP_RETCODE subtreesRemoveNonValids(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   DPITER*               dpiterator          /**< iterator */
)
{
   const CSR* const csr = graph->csr_storage;
   DHEAP* const dheap = dpsolver->dheap;
   const int* const head_csr = csr->head;
   const int* const start_csr = csr->start;
   const SCIP_Real* const nodes_dist = dpiterator->nodes_dist;
   const SCIP_Real* const nodes_ub = dpiterator->nodes_ub;
   SCIP_Bool* RESTRICT nodes_isValidRoot = dpiterator->nodes_isValidRoot;
   STP_Bitset sol_bitset = dpiterator->sol_termbits;
   const int* const nodes_termId = dpsolver->dpgraph->nodes_termId;
   const int nnodes = dpiterator->nnodes;
   SCIP_Real* nodes_mindist;
   SCIP_Real maxvaliddist = -1.0;
   STP_Vectype(int) stack = dpiterator->stack;
   const int extterm = dpiterator->extterm;
   int termcount = graph->terms - dpiterator->sol_nterms;
#ifdef STP_DPTERM_USEDA
   SCIP_Real* RESTRICT nodes_reddist = dpiterator->nodes_reddist;
   const SCIP_Real cutoffbound = dpsolver->dpredcosts->cutoffbound;
#endif

   assert(nnodes == graph->knots);
   assert(termcount >= 1);
   assert(graph_knot_isInRange(graph, extterm));

   if( allExtensionsAreInvalid(graph, dpsolver, dpiterator) )
   {
      BMSclearMemoryArray(nodes_isValidRoot, nnodes);
      return SCIP_OKAY;
   }

   graph_heap_clean(TRUE, dheap);
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_mindist, nnodes) );
   StpVecClear(stack);

   for( int i = 0; i < nnodes; i++ )
      nodes_mindist[i] = -1.0;

   nodes_mindist[extterm] = nodes_dist[extterm];
   graph_heap_correct(extterm, -nodes_mindist[extterm], dheap);
   SCIPdebugMessage("starting to exclude vertices with separators \n");

   /* propagate UBs */
   while( dheap->size > 0 )
   {
      int heapnode;
      SCIP_Real heapnode_dist;

      graph_heap_deleteMin(&heapnode, &heapnode_dist, dheap);
      heapnode_dist *= -1.0;
      assert(EQ(nodes_mindist[heapnode], heapnode_dist));
      assert(GE(heapnode_dist, 0.0));
      assert(StpVecGetSize(stack) == 0);

      StpVecPushBack(scip, stack, heapnode);
      SCIPdebugMessage("heapnode=%d dist=%f \n", heapnode, heapnode_dist);

      while( StpVecGetSize(stack) > 0 )
      {
         const int m = stack[StpVecGetSize(stack) - 1];
         const int m_start = start_csr[m];
         const int m_end = start_csr[m + 1];

         StpVecPopBack(stack);

         if( nodeIsNonSolTerm(sol_bitset, nodes_termId, m) )
         {
            if( --termcount == 0 )
            {
               maxvaliddist = heapnode_dist;
               graph_heap_clean(FALSE, dheap);
               assert(dheap->size == 0);
               break;
            }
         }

         for( int e = m_start; e != m_end; e++ )
         {
            SCIP_Real dist_new;
            const int q = head_csr[e];

            if( GT(nodes_ub[q], 0.0) )
            {
               continue;
            }

            dist_new = MIN(heapnode_dist, nodes_dist[q]);
            assert(GE(dist_new, 0.0));

            if( GT(dist_new, nodes_mindist[q]) )
            {
               if( LE(heapnode_dist, nodes_dist[q]) )
               {
                  SCIPdebugMessage("update and push-back node %d (dist=%f) \n", q, nodes_dist[q]);
                  StpVecPushBack(scip, stack, q);
               }
               else
               {
                  assert(EQ(dist_new, nodes_dist[q]));
                  graph_heap_correct(q, -dist_new, dheap);
               }

               nodes_mindist[q] = dist_new;
            }
         }
      }
   }

#ifdef STP_DPTERM_USEDA
   maxvaliddist = MIN(maxvaliddist, dpsolver->dpredcosts->upperbound);
#endif

   for( int i = 0; i < nnodes; i++ )
   {
      if( GT(nodes_dist[i], maxvaliddist) )
         nodes_isValidRoot[i] = FALSE;
#ifdef STP_DPTERM_USEDA
      // todo also check dpsolver->dpredcosts->nodes_rootdist[i]
      else if( GT(nodes_reddist[i], cutoffbound) )
      {
         //printf("kill %f > %f \n", nodes_reddist[i], cutoffbound);
         nodes_isValidRoot[i] = FALSE;
      }
#endif
   }

   debugPrintSeparator(maxvaliddist, dpiterator);

   SCIPdebugMessage("maxvaliddist=%f \n", maxvaliddist);

   dpiterator->stack = stack;
   SCIPfreeBufferArray(scip, &nodes_mindist);

   return SCIP_OKAY;
}


/** frees etc. */
static
void dpiterFinalizeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   DPITER*               dpiterator          /**< to update */
)
{
   assert(dpiterator->sol_traces == dpiterator->dpsubsol->traces);

   if( dpiterator->sol_termbits )
   {
      stpbitset_free(scip, &(dpiterator->sol_termbits));
   }

   dpterms_dpsubsolFree(scip, &(dpiterator->dpsubsol));
   dpiterator->sol_traces = NULL;
}


/*
 * Interface methods
 */


/** solves problem */
SCIP_RETCODE dpterms_coreSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph */
   DPSOLVER*             dpsolver,           /**< solver */
   SCIP_Bool*            wasSolved           /**< was problem solved to optimality? */
)
{
   DPITER* dpiterator;
   STP_PQ* const solpqueue = dpsolver->solpqueue;
   DPMISC* const dpmisc = dpsolver->dpmisc;

   assert(scip && graph && dpsolver && wasSolved);
   assert(solpqueue && dpsolver->dpgraph && dpsolver->soltree_root && dpsolver->dpstree && dpmisc);
   assert(!stpprioqueue_isClean(solpqueue));

   *wasSolved = TRUE;

   SCIP_CALL( dpiterInit(scip, graph, &dpiterator) );

   /* DP loop */
   while( !stpprioqueue_isClean(solpqueue) )
   {
      dpiterGetNextSol(scip, dpsolver, dpiterator);

      /*  construct implicit sub-Steiner trees by setting traces */
      SCIP_CALL( subtreesBuild(scip, dpmisc, dpiterator) );

      SCIP_CALL( subtreesExtend(scip, graph, dpsolver, dpiterator) );

      if( graph->terms - dpiterator->sol_nterms > 1 )
      {
         propagateUBs(graph, dpsolver, dpiterator);
         SCIP_CALL( subtreesRemoveNonValids(scip, graph, dpsolver, dpiterator) );
      }

      /* add all valid extensions (including merges) */
      SCIP_CALL( subtreesAddNew(scip, graph, dpsolver, dpiterator) );

      dpiterFinalizeSol(scip, dpiterator);

      if( dpmisc->global_size % 16 == 0 && SCIPisStopped(scip) )
      {
         *wasSolved = FALSE;
         break;
      }
   }

   assert(!dpsolver->solnodes);

   if( *wasSolved )
   {
      SCIPdebugMessage("OBJ=%f \n", dpmisc->opt_obj);
      dpsolver->solnodes = getSolnodesFinal(scip, dpmisc, dpiterator);
   }

   dpiterFree(scip, &dpiterator);

   return SCIP_OKAY;
}
