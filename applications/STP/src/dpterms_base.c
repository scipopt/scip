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
SCIP_RETCODE dpgraphInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPGRAPH**             dpgraph             /**< DP graph */
)
{
   DPGRAPH* dpg;
   int* terminals;
   int* nodes_termId;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);
   const int nterms = graph_get_nTerms(graph);
   int termscount = 0;

   assert(nterms >= 2);

   SCIP_CALL( SCIPallocMemory(scip, dpgraph) );
   dpg = *dpgraph;
   dpg->nnodes = nnodes;
   dpg->nterms = nterms;
   dpg->nedges = nedges;

   SCIP_CALL( SCIPallocMemoryArray(scip, &terminals, nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodes_termId, nnodes) );

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(graph->term[i]) )
      {
         assert(termscount < nterms);
         nodes_termId[i] = termscount;
         terminals[termscount++] = i;
      }
      else
      {
         nodes_termId[i] = -1;
      }
   }
   assert(termscount == nterms);

   dpg->terminals = terminals;
   dpg->nodes_termId = nodes_termId;

#ifndef NDEBUG
   for( int i = 0; i < nterms; i++ )
   {
      const int term = terminals[i];
      assert(graph_knot_isInRange(graph, term));
      assert(nodes_termId[term] == i);
   }
#endif

   return SCIP_OKAY;
}


/** frees */
static
void dpgraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPGRAPH**             dpgraph             /**< DP graph */
)
{
   DPGRAPH* dpg = *dpgraph;

   assert(dpg);
   assert(dpg->terminals && dpg->nodes_termId);

   SCIPfreeMemoryArray(scip, &(dpg->nodes_termId));
   SCIPfreeMemoryArray(scip, &(dpg->terminals));

   SCIPfreeMemory(scip, dpgraph);
}


/** initializes */
static
SCIP_RETCODE dpmiscInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPMISC**              dpmisc              /**< to initialize */
)
{
   DPMISC* misc;
   SCIP_CALL( SCIPallocMemory(scip, dpmisc) );
   misc = *dpmisc;

   misc->min = FARAWAY;
   misc->min_x = -1;
   misc->min_prev[0] = -1;
   misc->min_prev[1] = -1;
   misc->bits = NULL;
   misc->bits_count = NULL;
   misc->data = NULL;
   misc->offsets = NULL;
   misc->total_size = 0;

   StpVecPushBack(scip, misc->offsets, 0);

   return SCIP_OKAY;
}


/** frees */
static
void dpmiscFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPMISC**              dpmisc              /**< to free */
)
{
   DPMISC* misc = *dpmisc;
   int todo; // also free members of Vectors
   // ALSO RENAME MEMBERS!!!
   // ALSO RENAME MEMBERS!!!
   // ALSO RENAME MEMBERS!!!


   assert(misc);

   if( misc->bits_count )
   {
      StpVecFree(scip, misc->bits_count);
   }

   if( misc->bits )
   {
      StpVecFree(scip, misc->bits);
   }

   if( misc->offsets )
   {
      StpVecFree(scip, misc->offsets);
   }

   if( misc->data )
   {
      StpVecFree(scip, misc->data);
   }

   SCIPfreeMemory(scip, dpmisc);
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
   const int nterms = graph->terms;
   const int nnodes = graph->knots;
   const int* terminals;
   assert(nterms >= 2);

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &(dpsolver->dheap) ));
   SCIP_CALL( stpprioqueue_create(scip, nnodes, &(dpsolver->solpqueue)) );
   SCIP_CALL( dpterms_streeInit(scip, nterms, nnodes, &(dpsolver->dpstree)) );
   SCIP_CALL( dpmiscInit(scip, graph,  &(dpsolver->dpmisc)) );
   SCIP_CALL( dpgraphInit(scip, graph, &(dpsolver->dpgraph)) );
   terminals = dpsolver->dpgraph->terminals;

   for( int i = 0; i < nterms; i++ )
   {
      DPSUBSOL* singleton_sol;
      const int term = terminals[i];
      int pos;
      SOLTRACE trace = { .prevs = {-1,-1},
                         .cost = 0.0,
                         .root = term};

      assert(Is_term(graph->term[term]));
      SCIPdebugMessage("add term %d (index=%d) \n", term, i);

      SCIP_CALL( dpterms_dpsubsolInit(scip, &singleton_sol) );
      singleton_sol->bitkey = stpbitset_new(scip, nnodes);
      stpbitset_setBitTrue(singleton_sol->bitkey, i);

      SCIP_CALL( stpprioqueue_insert(scip, ((void*) stpbitset_newCopy(scip, singleton_sol->bitkey)),
            1, dpsolver->solpqueue) );

      assert(NULL == singleton_sol->extensions);
      StpVecPushBack(scip, singleton_sol->extensions, trace);

      pos = findSubsol(soltree_root, singleton_sol->bitkey, &soltree_parent);
      assert(pos != 0); /* not found */

      SCIPrbtreeInsert(&soltree_root, soltree_parent, pos, singleton_sol);
   }

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
   assert(dpsolver->dpgraph && dpsolver->dpstree && dpsolver->solpqueue);
   assert(dpsolver->dheap);
   assert(stpprioqueue_isClean(dpsolver->solpqueue));

   if( dpsolver->soltree_root )
   {
      FOR_EACH_NODE(DPSUBSOL*, node, dpsolver->soltree_root,
      {
         assert(node);
         SCIPrbtreeDelete(&(dpsolver->soltree_root), node);
         dpterms_dpsubsolFree(scip, &node);
      })

      assert(!dpsolver->soltree_root);
   }

   dpgraphFree(scip, &(dpsolver->dpgraph));
   dpmiscFree(scip, &(dpsolver->dpmisc));
   dpterms_streeFree(scip, &(dpsolver->dpstree));
   stpprioqueue_free(scip, &(dpsolver->solpqueue));
   graph_heap_free(scip, TRUE, TRUE, &(dpsolver->dheap));
}


/** solve problem */
static
SCIP_RETCODE dpsolverSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph of sub-problem */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   // todo compress?


   SCIP_CALL( graph_init_csr(scip, g) );

   SCIP_CALL( dpterms_coreSolve(scip, g, dpsolver) );

   graph_free_csr(scip, g);

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

   // todo we might want to pack the graph here...

   SCIP_CALL( dpsolverInit(scip, graph, &dpsolver) );

   SCIP_CALL( dpsolverSolve(scip, graph, dpsolver) );
   SCIP_CALL( dpsolverGetSolution(scip, dpsolver, solution) );

   dpsolverFree(scip, &dpsolver);

  // assert(0);


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
