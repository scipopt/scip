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

/**@file   dpterms_core.c
 * @brief  Core of dynamic programming solver for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Contains core methods.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
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
   STP_Bitset            sol_bitset;         /**< marks terminals of sub-solution */
   SCIP_Real*            nodes_dist;         /**< weight of sub-ST rooted at node */
   SCIP_Real*            nodes_ub;           /**< upper bounds for rule-out */
   int*                  nodes_previdx0;        /**< predecessor NOTE: with shift! */
   int*                  nodes_previdx1;        /**< predecessor */
   SCIP_Bool*            nodes_isValidRoot;  /**< is node a valid root? */
   int                   nnodes;             /**< number of nodes */
   int                   sol_nterms;         /**< popcount */
} DPITER;


/*
 * Local methods
 */


/** gets ordered root indices */
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
   iter->sol_bitset = NULL;
   iter->nnodes = nnodes;
   iter->sol_nterms = -1;

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

   assert(!iter->sol_traces && !iter->sol_bitset);
   StpVecFree(scip, iter->tripletstack);
   StpVecFree(scip, iter->stack);

   SCIPfreeMemoryArray(scip, &(iter->nodes_isValidRoot));
   SCIPfreeMemoryArray(scip, &(iter->nodes_previdx1));
   SCIPfreeMemoryArray(scip, &(iter->nodes_previdx0));
   SCIPfreeMemoryArray(scip, &(iter->nodes_ub));
   SCIPfreeMemoryArray(scip, &(iter->nodes_dist));

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

   dpiterator->sol_bitset = stpprioqueue_deleteMinReturnData(dpsolver->solpqueue);

   if( findSubsol(dpsolver->soltree_root, dpiterator->sol_bitset, &subsol) == 0 )
   {
      STP_Vectype(SOLTRACE) sol_traces = subsol->extensions; // todo needs to be freed
      printf("number of traces: %d \n", StpVecGetSize(sol_traces));

      SCIPrbtreeDelete(&(dpsolver->soltree_root), subsol);

      assert(stpbitset_areEqual(dpiterator->sol_bitset, subsol->bitkey));
      stpbitset_free(scip, &(subsol->bitkey));
   }
   else
   {
      assert(0 && "should never happen");
   }

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
   assert(!dpiterator->sol_bitset);

   dpiterSetDefault(scip, dpiterator);
   dpiterPopSol(scip, dpsolver, dpiterator);

#ifdef SCIP_DEBUG
      SCIPdebugMessage("processing solution: \n");
      stpbitset_print(dpiterator->sol_bitset);
#endif
}



/** (Implicitly) constructs sub-Steiner trees */
static
SCIP_RETCODE dpiterBuildSubTrees(
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
   int* nodes_nvisits;
   int* roots_indices;
   const int nsolcands = StpVecGetSize(soltraces);
   const int nnodes = dpiterator->nnodes;
   int nvalidroots = 0;

   STP_Vectype(SOLTRACE) data = dpmisc->data;
   //STP_Vectype(SOLTRACE) data = NULL;

   // todo at seome test to see that everything works!

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

      nodes_isValidRoot[sol_root] = TRUE;
      nodes_dist[sol_root] = sol_trace.cost;
      nodes_ub[sol_root] = 0.0;
      nodes_nvisits[sol_root]++;
      nvalidroots++;
      nodes_previdx0[sol_root] = sol_trace.prevs[0];
      nodes_previdx1[sol_root] = sol_trace.prevs[1];

      assert(sol_trace.prevs[0] != -1 || sol_trace.prevs[1] == -1);

      if( sol_trace.prevs[0] != -1 )
      {
         assert(StpVecGetSize(tripletstack) == 0);

         StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, 0.0, sol_trace.prevs[0]}));

         if( sol_trace.prevs[1] != -1 )
            StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, 0.0, sol_trace.prevs[1]}));

         while( StpVecGetSize(tripletstack) > 0 )
         {
            TTRIPLET triplet = tripletstack[StpVecGetSize(tripletstack) - 1];
            SOLTRACE parent_trace = data[triplet.index];
            const int previdx0 = parent_trace.prevs[0];
            const int previdx1 = parent_trace.prevs[1];
            const SCIP_Real bdist_local = triplet.bdist_local;
            const SCIP_Real bdist_global = triplet.bdist_global;

            StpVecPopBack(tripletstack);

            assert(previdx0 >= 0);

            /* merged solution? */
            if( previdx1 != -1 )
            {
               assert(previdx0 != -1);

               StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, bdist_global, previdx0}));
               StpVecPushBack(scip, tripletstack, ((TTRIPLET) {0.0, bdist_global, previdx1}));
            }
            else if( previdx0 != -1 )
            {
               const int curr_root = data[previdx0].root;
               const SCIP_Real curr_bdist_local = bdist_local + parent_trace.cost - data[previdx0].cost;
               const SCIP_Real curr_bdist_global = MAX(bdist_global, curr_bdist_local);

               assert(GT(parent_trace.cost - data[previdx0].cost, 0.0));

               nodes_dist[curr_root] = MIN(nodes_dist[curr_root], sol_trace.cost);
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


/** frees etc. */
static
void dpiterFinalizeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   DPITER*               dpiterator          /**< to update */
)
{
   stpbitset_free(scip, &(dpiterator->sol_bitset));
   dpterms_dpsubsolFree(scip, &(dpiterator->dpsubsol));
}


/*
 * Interface methods
 */


/** solves problem */
SCIP_RETCODE dpterms_coreSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   DPITER* dpiterator;
   STP_PQ* const solpqueue = dpsolver->solpqueue;
   DPGRAPH* const dpgraph = dpsolver->dpgraph;
   DPSTREE* const dpstree = dpsolver->dpstree;
   DPMISC* const dpmisc = dpsolver->dpmisc;

   assert(scip && graph && dpsolver);
   assert(solpqueue && dpgraph && dpsolver->soltree_root && dpstree && dpmisc);
   assert(!stpprioqueue_isClean(solpqueue));

   SCIP_CALL( dpiterInit(scip, graph, &dpiterator) );

   /* DP loop */
   while( !stpprioqueue_isClean(solpqueue) )
   {
      dpiterGetNextSol(scip, dpsolver, dpiterator);

      /*  construct implicit sub-Steiner trees by setting traces */
      SCIP_CALL( dpiterBuildSubTrees(scip, dpmisc, dpiterator) );


      // build extended Steiner trees

      // remove non-valid

      // merge subtrees

      dpiterFinalizeSol(scip, dpiterator);
   }


   assert(dpsolver->soltree_root == NULL);
   dpiterFree(scip, &dpiterator);


   assert(0);

#ifdef XXXXX

   // test only!
   {
      STP_Vectype(int) intersect;
      const int nterms = graph->terms;
      const int nnodes = graph->knots;
      STP_Bitset termsmark = stpbitset_new(scip, nterms);
      STP_Bitset rootsmark = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(termsmark, 1);
      stpbitset_setBitTrue(rootsmark, 2);
      SCIP_CALL( dpterms_streeInsert(scip, termsmark, rootsmark, 0, dpsolver->dpstree) );

      termsmark = stpbitset_new(scip, nterms);
      rootsmark = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(termsmark, 0);
      stpbitset_setBitTrue(rootsmark, 0);
      SCIP_CALL( dpterms_streeInsert(scip, termsmark, rootsmark, 1, dpsolver->dpstree) );

      termsmark = stpbitset_new(scip, nterms);
      rootsmark = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(termsmark, 2);
      stpbitset_setBitTrue(rootsmark, 0);
      SCIP_CALL( dpterms_streeInsert(scip, termsmark, rootsmark, 10, dpsolver->dpstree) );


      termsmark = stpbitset_new(scip, nterms);
      rootsmark = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(termsmark, 1);
      stpbitset_setBitTrue(rootsmark, 0);
     // SCIP_CALL( dpterms_streeInsert(scip, termsmark, rootsmark, 1, dpsolver->dpstree) );

     intersect = dpterms_streeCollectIntersects(scip, termsmark, rootsmark, dpsolver->dpstree);

      printf("size=%d \n", StpVecGetSize(intersect));

      for( int i = 0; i < StpVecGetSize(intersect); i++ )
      {
         printf("intersect=%d \n", intersect[i]);
      }



   }

   // todo delete, just for testing
   {
      DPSUBSOL* test;

      FOR_EACH_NODE(DPSUBSOL*, node, soltree_root,
      {
              assert(node);
              stpbitset_print(node->bitkey);
              printf("popcount=%d \n", stpbitset_getPopcount(node->bitkey));
      })

      SCIP_CALL( dpsubsolInit(scip, &test) );
      test->bitkey = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(test->bitkey, 3);

      if( findSubsol(soltree_root, test->bitkey, &soltree_parent) == 0 )
      {
         printf("found \n");
         stpbitset_print(soltree_parent->bitkey);
      }
      else
      {
         printf("not found \n");
      }

      dpsubsolFree(scip, &test);
   }
#endif



   return SCIP_OKAY;
}
