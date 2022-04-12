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
#include "solstp.h"
#ifdef STP_DPTERM_USEDA
#include "dualascent.h"
#include "heur_ascendprune.h"
#include "heur_local.h"
#endif


// todo more or less random values, tune them!
#define PROMISING_FULL_MAXNTERMS              20
#define PROMISING_FULL_MAXNTERMS_LARGE        15
#define PROMISING_PARTLY_MAXDENSITY          2.2
#define PROMISING_PARTLY_SMALL_MAXNTERMS      45
#define PROMISING_PARTLY_SMALL_MINNEDGES    1100
#define PROMISING_PARTLY_MEDIUM_MAXNTERMS     55
#define PROMISING_PARTLY_MEDIUM_MINNEDGES   1500
#define PROMISING_PARTLY_LARGE_MAXNTERMS      65
#define PROMISING_PARTLY_LARGE_MINNEDGES    3000
#define PROMISING_FULL_LARGE_MINNEDGES    100000
#define PROMISING_FULL_MAXAVGDEG             5.0


/*
 * Local methods
 */

#ifdef STP_DPTERM_USEDA

/** initializes */
static
SCIP_RETCODE dpredcostsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< original graph */
   DPREDCOST**           dpredcosts          /**< DP graph */
)
{
   SCIP_Real* redcosts_tmp;
   int* soledges_tmp;
   DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = FALSE, .root = graph->source,
         .is_pseudoroot = FALSE, .damaxdeviation = -1.0 };
   DPREDCOST* dprc;
   SCIP_Real dualobjval = -1.0;
   SCIP_Real primalobjval;
   SCIP_Bool success;
   int* pathedge;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   SCIP_CALL( SCIPallocMemory(scip, dpredcosts) );
   dprc = *dpredcosts;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(redcosts_tmp), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(soledges_tmp), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dprc->csr_redcosts), nedges) );

   SCIP_CALL( dualascent_exec(scip, graph, NULL, &daparams, redcosts_tmp, &dualobjval) );

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, graph, redcosts_tmp, soledges_tmp, graph->source, &success, FALSE));
   assert(success);
   SCIP_CALL(SCIPStpHeurLocalRun(scip, graph, soledges_tmp));
   assert(solstp_isValid(scip, graph, soledges_tmp));

   primalobjval = solstp_getObj(graph, soledges_tmp, 0.0);
   assert(GE(primalobjval, dualobjval));
   dprc->cutoffbound = primalobjval - dualobjval;
   dprc->upperbound = primalobjval;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(dprc->nodes_rootdist), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes + 1) );
   graph_mark(graph);
   graph_path_execX(scip, graph, graph->source, redcosts_tmp, dprc->nodes_rootdist, pathedge);
   SCIPfreeBufferArray(scip, &pathedge);

   printf("dual=%f primal =%f \n", dualobjval, primalobjval);
   printf("cutoffbound=%f \n", dprc->cutoffbound);

   // todo extra method
   {
      const CSR* const csr = graph->csr_storage;
      const int* const edgeid_csr = csr->edge_id;
      const int* const start_csr = csr->start;

      assert(csr && edgeid_csr && start_csr);

      for( int k = 0; k < nnodes; k++ )
      {
         for( int j = start_csr[k]; j != start_csr[k + 1]; j++  )
         {
            const int edge = edgeid_csr[j];
            const SCIP_Real rc = MIN(redcosts_tmp[edge], redcosts_tmp[flipedge(edge)]);

            assert(EQ(graph->cost[edge], csr->cost[j]));
            assert(GE(rc, 0.0));
            assert(0 <= j && j < nedges);

            dprc->csr_redcosts[j] = rc;
         }
      }
   }

   SCIPfreeMemoryArray(scip, &soledges_tmp);
   SCIPfreeMemoryArray(scip, &redcosts_tmp);

   return SCIP_OKAY;
}


/** frees */
static
void dpredcostsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPREDCOST**           dpredcosts          /**< DP graph */
)
{
   DPREDCOST* dprc = *dpredcosts;

   assert(dprc);
   assert(dprc->csr_redcosts);

   SCIPfreeMemoryArray(scip, &(dprc->nodes_rootdist));
   SCIPfreeMemoryArray(scip, &(dprc->csr_redcosts));

   SCIPfreeMemory(scip, dpredcosts);
}
#endif


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

   misc->opt_obj = FARAWAY;
   misc->opt_root = -1;
   misc->opt_prev[0] = -1;
   misc->opt_prev[1] = -1;
   misc->global_termbits = NULL;
   misc->global_termbitscount = NULL;
   misc->global_traces = NULL;
   misc->global_starts = NULL;
   misc->global_size = 0;

   StpVecPushBack(scip, misc->global_starts, 0);

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
   assert(misc);

   if( misc->global_termbitscount )
   {
      StpVecFree(scip, misc->global_termbitscount);
   }

   if( misc->global_termbits )
   {
      const int size = StpVecGetSize(misc->global_termbits);
      for( int i = 0; i < size; i++ )
         stpbitset_free(scip, &(misc->global_termbits[i]));

      StpVecFree(scip, misc->global_termbits);
   }

   if( misc->global_starts )
   {
      StpVecFree(scip, misc->global_starts);
   }

   if( misc->global_traces )
   {
      StpVecFree(scip, misc->global_traces);
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

#ifdef STP_DPTERM_USEDA
   SCIP_CALL( dpredcostsInit(scip, graph, &(dpsolver->dpredcosts)) );
#endif

   for( int i = 0; i < nterms; i++ )
   {
      DPSUBSOL* singleton_sol;
      const int term = terminals[i];
      int pos;
#ifdef STP_DPTERM_USEDA
      SOLTRACE trace = { .prevs = {-1,-1},
                         .cost = 0.0,
                         .redcost = 0.0,
                         .root = term};
#else
      SOLTRACE trace = { .prevs = {-1,-1},
                         .cost = 0.0,
                         .root = term};
#endif


      assert(Is_term(graph->term[term]));
      SCIPdebugMessage("add term %d (index=%d) \n", term, i);

      SCIP_CALL( dpterms_dpsubsolInit(scip, &singleton_sol) );
      singleton_sol->bitkey = stpbitset_new(scip, nterms);
      stpbitset_setBitTrue(singleton_sol->bitkey, i);

      SCIP_CALL( stpprioqueue_insert(scip, ((void*) stpbitset_newCopy(scip, singleton_sol->bitkey)),
            1, dpsolver->solpqueue) );

      assert(NULL == singleton_sol->traces);
      StpVecPushBack(scip, singleton_sol->traces, trace);

      pos = findSubsol(soltree_root, singleton_sol->bitkey, &soltree_parent);
      assert(pos != 0); /* not found */

      SCIPrbtreeInsert(&soltree_root, soltree_parent, pos, singleton_sol);
   }

   dpsolver->soltree_root = soltree_root;
   dpsolver->solnodes = NULL;

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
   assert(SCIPisStopped(scip) || stpprioqueue_isClean(dpsolver->solpqueue));

   StpVecFree(scip, dpsolver->solnodes);

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

#ifdef STP_DPTERM_USEDA
   dpredcostsFree(scip, &(dpsolver->dpredcosts));
#endif

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
   DPSOLVER*             dpsolver,           /**< solver */
   SCIP_Bool*            wasSolved           /**< was problem solved to optimality? */
)
{
   SCIP_CALL( dpterms_coreSolve(scip, g, dpsolver, wasSolved) );

   return SCIP_OKAY;
}


/** gets optimal solution */
static
SCIP_RETCODE dpsolverGetSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPSOLVER*             dpsolver,           /**< the solver */
   int*                  solution            /**< to store solution */
)
{
   STP_Bool* connected;
   STP_Vectype(int) solnodes = dpsolver->solnodes;

   assert(solnodes);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &connected, graph->knots) );

   for( int i = 0; i < StpVecGetSize(solnodes); i++ )
   {
      const int node = solnodes[i];
      assert(graph_knot_isInRange(graph, node));
      assert(!connected[node]);

      connected[node] = TRUE;
   }

   SCIP_CALL( solstp_pruneFromNodes(scip, graph, solution, connected) );

   SCIPfreeBufferArray(scip, &connected);

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
   int*                  solution,           /**< optimal solution (out) */
   SCIP_Bool*            wasSolved           /**< was problem solved to optimality? */
)
{
   DPSOLVER* dpsolver;

   assert(scip && graph && solution);

   SCIP_CALL( graph_init_csrWithEdgeId(scip, graph) );
   SCIP_CALL( dpsolverInit(scip, graph, &dpsolver) );

   SCIP_CALL( dpsolverSolve(scip, graph, dpsolver, wasSolved) );

   if( *wasSolved )
   {
      SCIP_CALL( dpsolverGetSolution(scip, graph, dpsolver, solution) );
      assert(solstp_isValid(scip, graph, solution));
   }

   dpsolverFree(scip, &dpsolver);
   graph_free_csr(scip, graph);

   return SCIP_OKAY;
}

/** is DP at least partly promising? */
SCIP_Bool dpterms_isPromisingPartly(
   const GRAPH*          graph               /**< graph */
)
{
   /* NOTE: we count the undirected edges here */
   const int nedges = (graph->edges / 2);
   SCIP_Real density;

   assert(graph);

   if( dpterms_isPromisingFully(graph) )
      return TRUE;

   density = nedges / graph->knots;

   if( GT(density, PROMISING_PARTLY_MAXDENSITY) )
      return FALSE;

   if( graph->terms <= PROMISING_PARTLY_SMALL_MAXNTERMS && nedges >= PROMISING_PARTLY_SMALL_MINNEDGES )
      return TRUE;

   if( graph->terms <= PROMISING_PARTLY_MEDIUM_MAXNTERMS && nedges >= PROMISING_PARTLY_MEDIUM_MINNEDGES )
      return TRUE;

   if( graph->terms < PROMISING_PARTLY_LARGE_MAXNTERMS && nedges >= PROMISING_PARTLY_LARGE_MINNEDGES )
      return TRUE;

   // todo just a test, remove!
   if( graph->terms == 73 && nedges >= 3500 && nedges <= 4000  )
       return TRUE;

   return FALSE;
}


/** is DP embarrassingly promising? */
SCIP_Bool dpterms_isPromisingEmbarrassingly(
   const GRAPH*          graph               /**< graph */
)
{
   return ( graph->terms <= 3 );
}


/** is DP fully promising? */
SCIP_Bool dpterms_isPromisingFully(
   const GRAPH*          graph               /**< graph */
)
{
   assert(graph);

   if( graph->terms <= PROMISING_FULL_MAXNTERMS_LARGE )
      return TRUE;

   if( graph->terms <= PROMISING_FULL_MAXNTERMS )
   {
      int nedges;
      int nnodes;
      graph_get_nVET(graph, &nnodes, &nedges, NULL);

      if( nedges < PROMISING_FULL_LARGE_MINNEDGES && LT((SCIP_Real) nedges / (SCIP_Real) nnodes, PROMISING_FULL_MAXAVGDEG) )
         return TRUE;
   }

   return FALSE;
}
