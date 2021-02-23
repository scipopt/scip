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
 * Local methods
 */


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
   STP_PQ* const solpqueue = dpsolver->solpqueue;
   DPGRAPH* const dpgraph = dpsolver->dpgraph;
   DPSUBSOL* soltree_root = dpsolver->soltree_root;
   DPSTREE* const dpstree = dpsolver->dpstree;
   DPMISC* const dpmisc = dpsolver->dpmisc;

   assert(scip && graph && dpsolver);
   assert(solpqueue && dpgraph && soltree_root && soltree_root && dpstree && dpmisc);
   assert(!stpprioqueue_isClean(solpqueue));

   /* DP loop */
   while( !stpprioqueue_isClean(solpqueue) )
   {
      // extra method getNextSol
      DPSUBSOL* subsol;
      STP_Bitset sol_bitset = stpprioqueue_deleteMinReturnData(solpqueue);

      if( findSubsol(soltree_root, sol_bitset, &subsol) == 0 )
      {
         STP_Vectype(SOLTRACE) sol_traces = subsol->extensions; // todo needs to be freed
         printf("number of traces: %d \n", StpVecGetSize(sol_traces));

         SCIPrbtreeDelete(&soltree_root, subsol);

         stpbitset_free(scip, &(subsol->bitkey));
      }
      else
      {
         assert(0 && "should never happen");
      }

    //

    //


#ifdef SCIP_DEBUG
      SCIPdebugMessage("processing solution: \n");
      stpbitset_print(sol_bitset);
#endif

      stpbitset_free(scip, &sol_bitset);
   }


   assert(soltree_root == NULL);

   dpsolver->soltree_root = soltree_root;
   assert(0);


   // test only!
   {
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

      STP_Vectype(int) intersect = dpterms_streeCollectIntersects(scip, termsmark, rootsmark, dpsolver->dpstree);

      printf("size=%d \n", StpVecGetSize(intersect));

      for( int i = 0; i < StpVecGetSize(intersect); i++ )
      {
         printf("intersect=%d \n", intersect[i]);
      }



   }

#ifdef XXXXX
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
