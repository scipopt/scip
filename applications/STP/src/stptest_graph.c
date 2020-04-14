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

/**@file   stptest_graph.c
 * @brief  tests for Steiner tree problem methods
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problems.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "scip/scip.h"
#include "stptest.h"
#include "graph.h"
#include "portab.h"


/** are the CSRs equal? */
static
SCIP_Bool csrdepoCSRsAreEqual(
   const CSR*            csr1,                /**< csr */
   const CSR*            csr2                 /**< csr */
   )
{
   if( csr1->nedges_max != csr2->nedges_max )
   {
      SCIPdebugMessage("wrong edge count \n");
      return FALSE;
   }

   if( csr1->nnodes != csr2->nnodes )
   {
      SCIPdebugMessage("wrong node count \n");
      return FALSE;
   }

   for( int i = 0; i <= csr1->nnodes; ++i )
   {
      if( csr1->start[i] != csr2->start[i] )
      {
         SCIPdebugMessage("wrong start array \n");
         return FALSE;
      }
   }

   for( int i = 0; i < csr1->nedges_max; ++i )
   {
      if( !EQ(csr1->cost[i], csr2->cost[i]) )
      {
         SCIPdebugMessage("wrong cost array \n");
         return FALSE;
      }

      if( csr1->head[i] != csr2->head[i] )
      {
         SCIPdebugMessage("wrong head array \n");
         return FALSE;
      }
   }

   return TRUE;
}


/** simple pseudo random fill of CSRs */
static
void csrdepoFillRandom(
   int                   seed,               /**< seed */
   CSR*                  csrd                /**< csr */
   )
{
   int* start_csr;
   int* head_csr;
   SCIP_Real* cost_csr;
   const int nnodes = csrd->nnodes;
   const int nedges = csrd->nedges_max;

   assert(nnodes >= 1 && nedges >= 1);
   assert(seed >= 0);

   start_csr = csrd->start;
   head_csr = csrd->head;
   cost_csr = csrd->cost;

   start_csr[0] = 0;

   for( int i = 1; i <= nnodes; i++ )
   {
      start_csr[i] = start_csr[i - 1] + (i + seed) % 8;
   }

   for( int i = 0; i < nedges; i++ )
   {
      head_csr[i] = (i + seed) % 7;
      cost_csr[i] = (SCIP_Real)((i + seed) % 12) + 0.5 ;
   }
}


/** first test: Simple additions and deletions */
static
SCIP_RETCODE csrdepoTest1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   const int nnodes1 = 5;
   const int nedges1 = 12;
   const int nnodes2 = 9;
   const int nedges2 = 24;

   CSRDEPO* depo;
   CSR csr;

   SCIP_CALL( graph_csrdepo_init(scip, 8, 1000, &depo) );

   /* add first */
   graph_csrdepo_addEmptyTop(depo, nnodes1, nedges1);
   graph_csrdepo_getEmptyTop(depo, &csr);

   if( csr.nedges_max != nedges1 || csr.nnodes != nnodes1 )
   {
      SCIPdebugMessage("wrong number of nodes/edges (%d, %d) stored! \n", csr.nnodes, csr.nedges_max);
      return SCIP_ERROR;
   }

   csrdepoFillRandom(5522, &csr);
   graph_csrdepo_emptyTopSetMarked(depo);

   /* add second */
   graph_csrdepo_addEmptyTop(depo, nnodes2, nedges2);
   graph_csrdepo_getEmptyTop(depo, &csr);

   if( csr.nedges_max != nedges2 || csr.nnodes != nnodes2 )
   {
      SCIPdebugMessage("wrong number of nodes/edges (%d, %d) stored! \n", csr.nnodes, csr.nedges_max);
      return SCIP_ERROR;
   }

   /* remove both */
   graph_csrdepo_removeTop(depo);
   graph_csrdepo_removeTop(depo);

   if( !graph_csrdepo_isEmpty(depo) )
   {
      SCIPdebugMessage("depo is not empty! \n");
      return SCIP_ERROR;
   }

   graph_csrdepo_free(scip, &depo);

   return SCIP_OKAY;
}


/** second test: Make sure that the CSRs are stored correctly */
static
SCIP_RETCODE csrdepoTest2(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   const int nnodes0 = 5;
   const int nedges0 = 12;
   const int nnodes1 = 10;
   const int nedges1 = 24;
   const int nnodes2 = 15;
   const int nedges2 = 28;

   CSRDEPO* depo;
   CSR csr_in;
   CSR csr_out;
   CSR* csr0;
   CSR* csr1;
   CSR* csr2;

   SCIP_CALL( graph_csrdepo_init(scip, 8, 1000, &depo) );

   /* add first */
   graph_csrdepo_addEmptyTop(depo, nnodes0, nedges0);
   graph_csrdepo_getEmptyTop(depo, &csr_in);

   SCIP_CALL( graph_csr_alloc(scip, nnodes0, nedges0, &csr0) );
   csrdepoFillRandom(55, &csr_in);
   csrdepoFillRandom(55, csr0);

   graph_csrdepo_emptyTopSetMarked(depo);

   /* dummy add */
   graph_csrdepo_addEmptyTop(depo, nnodes1, nedges1);
   graph_csrdepo_getEmptyTop(depo, &csr_in);
   csrdepoFillRandom(200, &csr_in);

   graph_csrdepo_removeTop(depo);

   /* add second */
   graph_csrdepo_addEmptyTop(depo, nnodes1, nedges1);
   graph_csrdepo_getEmptyTop(depo, &csr_in);

   SCIP_CALL( graph_csr_alloc(scip, nnodes1, nedges1, &csr1) );
   csrdepoFillRandom(551, &csr_in);
   csrdepoFillRandom(551, csr1);

   graph_csrdepo_emptyTopSetMarked(depo);

   /* add third */
   graph_csrdepo_addEmptyTop(depo, nnodes2, nedges2);
   graph_csrdepo_getEmptyTop(depo, &csr_in);

   SCIP_CALL( graph_csr_alloc(scip, nnodes2, nedges2, &csr2) );
   csrdepoFillRandom(44, &csr_in);
   csrdepoFillRandom(44, csr2);

   graph_csrdepo_emptyTopSetMarked(depo);

   /* now check: */

   graph_csrdepo_getCSR(depo, 0, &csr_out);

   if( !csrdepoCSRsAreEqual(&csr_out, csr0) )
   {
      SCIPdebugMessage("CSRs 0 not equal! \n");
      return SCIP_ERROR;
   }

   graph_csrdepo_getCSR(depo, 2, &csr_out);

   if( !csrdepoCSRsAreEqual(&csr_out, csr2) )
   {
      SCIPdebugMessage("CSRs 2 not equal! \n");
      return SCIP_ERROR;
   }

   graph_csrdepo_removeTop(depo);

   graph_csrdepo_getCSR(depo, 1, &csr_out);

   if( !csrdepoCSRsAreEqual(&csr_out, csr1) )
   {
      SCIPdebugMessage("CSRs 1 not equal! \n");
      return SCIP_ERROR;
   }

   /* remove all */
   graph_csrdepo_removeTop(depo);
   graph_csrdepo_removeTop(depo);

   assert(graph_csrdepo_isEmpty(depo));

   /* clean-up*/
   graph_csrdepo_free(scip, &depo);
   graph_csr_free(scip, &csr2);
   graph_csr_free(scip, &csr1);
   graph_csr_free(scip, &csr0);

   return SCIP_OKAY;
}

/** frees, etc. */
void stptest_graphTearDown(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
)
{
   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);
}


/** sets up graph */
SCIP_RETCODE stptest_graphSetUp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   graph_mark(graph);

   return SCIP_OKAY;
}


/** sets up graph for (undirected) PC */
SCIP_RETCODE stptest_graphSetUpPcOrg(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   int*                  nnodes_new,         /**< to store new number of nodes (if != NULL)  */
   int*                  nedges_new          /**< to store new number of edge (if != NULL) */
   )
{
   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, nnodes_new, nedges_new) );

   graph_pc_2org(scip, graph);

   return SCIP_OKAY;
}


/** sets up graph for RMW */
SCIP_RETCODE stptest_graphSetUpRmwOrg(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   int*                  nnodes_new,         /**< to store new number of nodes (if != NULL)  */
   int*                  nedges_new          /**< to store new number of edge (if != NULL) */
   )
{
   SCIP_CALL( stptest_graphSetUpRmwExtended(scip, graph, nnodes_new, nedges_new) );

   graph_pc_2org(scip, graph);

   return SCIP_OKAY;
}


/** sets up graph for RMW */
SCIP_RETCODE stptest_graphSetUpRmwExtended(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   int*                  nnodes_new,         /**< to store new number of nodes (if != NULL)  */
   int*                  nedges_new          /**< to store new number of edge (if != NULL) */
   )
{
   graph->stp_type = STP_RMWCSP;

   SCIP_CALL( graph_transRmw(scip, graph) );

   stptest_graphSetUp(scip, graph);

   if( nnodes_new )
      *nnodes_new = graph->knots;

   if( nedges_new )
      *nedges_new = graph->edges;

   return SCIP_OKAY;
}


/** sets up graph for (undirected) PC */
SCIP_RETCODE stptest_graphSetUpPcExtended(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   int*                  nnodes_new,         /**< to store new number of nodes (if != NULL)  */
   int*                  nedges_new          /**< to store new number of edge (if != NULL) */
   )
{
   graph->stp_type = STP_PCSPG;

   SCIP_CALL( graph_transPc(scip, graph) );

   stptest_graphSetUp(scip, graph);

   if( nnodes_new )
      *nnodes_new = graph->knots;

   if( nedges_new )
      *nedges_new = graph->edges;

   return SCIP_OKAY;
}


/** tests CSR depository */
SCIP_RETCODE stptest_csrdepo(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( csrdepoTest1(scip) );
   SCIP_CALL( csrdepoTest2(scip) );

   printf("csrdepo test: all ok \n");

   return SCIP_OKAY;
}
