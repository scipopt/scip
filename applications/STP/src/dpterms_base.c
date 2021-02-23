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

/**@file   dpterms_base.c
 * @brief  Dynamic programming solver for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

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



/** initializes */
static
SCIP_RETCODE dpsubsolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSUBSOL**            subsol              /**< solution */
)
{
   DPSUBSOL* sub;
   SCIP_CALL( SCIPallocBlockMemory(scip, subsol) );
   sub = *subsol;
   sub->bitkey = NULL;
   sub->extensions = NULL;

   return SCIP_OKAY;
}


/** frees */
static
void dpsubsolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSUBSOL**            subsol              /**< solution */
)
{
   DPSUBSOL* sub = *subsol;

   if( sub->bitkey )
      stpbitset_free(scip, &(sub->bitkey));

   if( sub->extensions )
   {
      assert(0); // todo! free every single element

   }

   SCIPfreeBlockMemory(scip, subsol);
}


/** initializes data */
static
SCIP_RETCODE dpsolverInitData(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   DPSUBSOL* soltree_root = NULL;
   DPSUBSOL* soltree_parent;
   const int nterms = 3;
   const int nnodes = 4;

   SCIP_CALL( dpterms_streeInit(scip, nterms, nnodes, &(dpsolver->dpstree)) );

   // test only!
   {
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


      assert(0);

   }

   for( int i = 0; i < nterms; i++ )
   {
      DPSUBSOL* singleton;
      const int term = i; // todo
      int pos;

      SCIP_CALL( dpsubsolInit(scip, &singleton) );
      singleton->bitkey = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(singleton->bitkey, term);

      pos = findSubsol(soltree_root, singleton->bitkey, &soltree_parent);
      assert(pos != 0); /* not found */

      SCIPrbtreeInsert(&soltree_root, soltree_parent, pos, singleton);
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


   // todo: remove all

   assert(0);

   dpsolver->soltree_root = soltree_root;

   return SCIP_OKAY;
}


/** frees data*/
static
void dpsolverFreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   assert(dpsolver->soltree_root == NULL);

   dpterms_streeFree(scip, &(dpsolver->dpstree));


   // todo should not be necessary, all should be freed anyway
#ifdef SCIP_DISABLED
   DPSUBSOL* soltree_root = dpsolver->soltree_root;

   FOR_EACH_NODE(DPSUBSOL*, node, soltree_root,
        {
           assert(node);
           //delete!
        })
#endif
}


/** solve problem */
static
SCIP_RETCODE dpsolverSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph of sub-problem */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   // compress...

   SCIP_CALL( dpterms_coreSolve(scip, g, dpsolver) );

   return SCIP_OKAY;
}


/** gets optimal solution */
static
SCIP_RETCODE dpsolverGetSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER*             dpsolver,           /**< the solver */
   int*                  solution            /**< to store solution */
)
{

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE dpsolverInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPSOLVER**            dpsolver            /**< solver */
)
{
   SCIP_CALL( SCIPallocMemory(scip, dpsolver) );

   SCIP_CALL( dpsolverInitData(scip, graph, *dpsolver) );

   return SCIP_OKAY;
}


/** frees */
static
void dpsolverFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSOLVER**            dpsolver            /**< solver */
)
{
   dpsolverFreeData(scip, *dpsolver);

   SCIPfreeMemory(scip, dpsolver);
}



/*
 * Interface methods
 */


/** solves problem given by graph */
SCIP_RETCODE dpterms_solve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   int*                  solution            /**< was sub-problem solved to optimality? */
)
{
   DPSOLVER* dpsolver;

   assert(scip && graph && solution);

   SCIP_CALL( dpsolverInit(scip, graph, &dpsolver) );

   SCIP_CALL( dpsolverSolve(scip, graph, dpsolver) );
   SCIP_CALL( dpsolverGetSolution(scip, dpsolver, solution) );

   dpsolverFree(scip, &dpsolver);

   return SCIP_OKAY;
}


/** solves sub-problem */
SCIP_Bool dpterms_isPromising(
   const GRAPH*          graph               /**< graph */
)
{
   assert(graph);

   // todo
   return FALSE;
}
