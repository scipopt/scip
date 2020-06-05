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

/**@file   stptest_reduceutils.c
 * @brief  tests for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements unit tests for Steiner problem reduction utility methods.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "graph.h"
#include "reduce.h"
#include "portab.h"


/** builds up 3x3 MST */
static
SCIP_RETCODE dcmstTest1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DCMST* dcmst;
   CSR* csr_base;
   CSR* csr_extended;
   SCIP_Real adjcosts[] = { -1.0, -1.0 };

   SCIP_CALL( reduce_dcmstInit(scip, 3, &dcmst) );
   SCIP_CALL( graph_csr_alloc(scip, 2, 2, &csr_base) );
   SCIP_CALL( graph_csr_alloc(scip, 3, 4, &csr_extended) );

   csr_base->nnodes = 1;
   csr_base->nedges_max = 0;

   reduce_dcmstGet1NodeMst(scip, csr_base);

   csr_extended->nnodes = 2;
   csr_extended->nedges_max = 2;

   adjcosts[0] = 2.5;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 2.5) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   csr_base->nnodes = 2;
   csr_base->nedges_max = 2;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 3;
   csr_extended->nedges_max = 4;

   adjcosts[0] = 1.0;
   adjcosts[1] = 1.2;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 2.2) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   adjcosts[0] = 3.0;
   adjcosts[1] = 3.2;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 5.5) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   reduce_dcmstFree(scip, &dcmst);
   graph_csr_free(scip, &csr_base);
   graph_csr_free(scip, &csr_extended);

   return SCIP_OKAY;
}


/** builds up 4x4 MST */
static
SCIP_RETCODE dcmstTest2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DCMST* dcmst;
   CSR* csr_base;
   CSR* csr_extended;
   SCIP_Real adjcosts[] = { -1.0, -1.0, -1.0 };

   SCIP_CALL( reduce_dcmstInit(scip, 4, &dcmst) );
   SCIP_CALL( graph_csr_alloc(scip, 3, 4, &csr_base) );
   SCIP_CALL( graph_csr_alloc(scip, 4, 6, &csr_extended) );

   csr_base->nnodes = 2;
   csr_base->nedges_max = 2;

   reduce_dcmstGet2NodeMst(scip, 1.0, csr_base);

   csr_extended->nnodes = 3;
   csr_extended->nedges_max = 4;

   adjcosts[0] = 1.5;
   adjcosts[1] = 0.7;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 1.7) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   csr_base->nnodes = 3;
   csr_base->nedges_max = 4;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 4;
   csr_extended->nedges_max = 6;

   adjcosts[0] = 0.5;
   adjcosts[1] = 2.7;
   adjcosts[2] = 2.7;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 2.2) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }


   reduce_dcmstFree(scip, &dcmst);
   graph_csr_free(scip, &csr_base);
   graph_csr_free(scip, &csr_extended);

   return SCIP_OKAY;
}



/** builds up a 4x4 MST */
static
SCIP_RETCODE dcmstTest2b(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DCMST* dcmst;
   CSR* csr_base;
   CSR* csr_extended;
   SCIP_Real adjcosts[] = { -1.0, -1.0, -1.0 };

   SCIP_CALL( reduce_dcmstInit(scip, 4, &dcmst) );
   SCIP_CALL( graph_csr_alloc(scip, 3, 4, &csr_base) );
   SCIP_CALL( graph_csr_alloc(scip, 4, 6, &csr_extended) );

   csr_base->nnodes = 2;
   csr_base->nedges_max = 2;

   reduce_dcmstGet2NodeMst(scip, 395.0, csr_base);

   csr_extended->nnodes = 3;
   csr_extended->nedges_max = 4;

   adjcosts[0] = 200.0;
   adjcosts[1] = 302.0;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 502.0) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   csr_base->nnodes = 3;
   csr_base->nedges_max = 4;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 4;
   csr_extended->nedges_max = 6;

   adjcosts[0] = 197.0;
   adjcosts[1] = FARAWAY;
   adjcosts[2] = 197.0;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 696.0) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }


   reduce_dcmstFree(scip, &dcmst);
   graph_csr_free(scip, &csr_base);
   graph_csr_free(scip, &csr_extended);

   return SCIP_OKAY;
}


/** builds up 6x6 MST */
static
SCIP_RETCODE dcmstTest3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DCMST* dcmst;
   CSR* csr_base;
   CSR* csr_extended;
   SCIP_Real adjcosts[] = { -1.0, -1.0, -1.0, -1.0, -1.0 };

   SCIP_CALL( reduce_dcmstInit(scip, 6, &dcmst) );
   SCIP_CALL( graph_csr_alloc(scip, 5, 8, &csr_base) );
   SCIP_CALL( graph_csr_alloc(scip, 6, 10, &csr_extended) );

   csr_base->nnodes = 2;
   csr_base->nedges_max = 2;

   reduce_dcmstGet2NodeMst(scip, 3.1, csr_base);

   /* build on 3x3 */
   csr_extended->nnodes = 3;
   csr_extended->nedges_max = 4;

   adjcosts[0] = 1.2;
   adjcosts[1] = 7.3;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 4.3) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   /* build on 4x4 */
   csr_base->nnodes = 3;
   csr_base->nedges_max = 4;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 4;
   csr_extended->nedges_max = 6;

   adjcosts[0] = 1.6;
   adjcosts[1] = 0.1;
   adjcosts[2] = 0.3;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 1.6) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   /* build on 5x5 */
   csr_base->nnodes = 4;
   csr_base->nedges_max = 6;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 5;
   csr_extended->nedges_max = 8;

   adjcosts[0] = 0.2;
   adjcosts[1] = 0.01;
   adjcosts[2] = 2.6;
   adjcosts[3] = 4.1;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 0.61) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   /* build on 6x6 */
   csr_base->nnodes = 5;
   csr_base->nedges_max = 8;

   graph_csr_copy(csr_extended, csr_base);

   csr_extended->nnodes = 6;
   csr_extended->nedges_max = 10;

   adjcosts[0] = 0.3;
   adjcosts[1] = 0.2;
   adjcosts[2] = 0.6;
   adjcosts[3] = 4.1;
   adjcosts[4] = 5.1;

   reduce_dcmstAddNode(scip, csr_base, adjcosts, dcmst, csr_extended);

   if( !EQ(reduce_dcmstGetWeight(scip, csr_extended), 0.81) )
   {
      SCIPdebugMessage("wrong MST weight \n");

      return SCIP_ERROR;
   }

   reduce_dcmstFree(scip, &dcmst);
   graph_csr_free(scip, &csr_base);
   graph_csr_free(scip, &csr_extended);

   return SCIP_OKAY;
}


/** STAR of degree 3 should be correctly computed */
static
SCIP_RETCODE testStar3IsCorrect(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   STAR* star;
   GRAPH* graph;
   int nnodes = 4;
   int nedges = 6;
   int nstaredges = -1;
   const int* staredges;

   SCIP_CALL( reduce_starInit(scip, 3, &star) );
   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   reduce_starReset(graph, 0, star);
   staredges = reduce_starGetNext(star, &nstaredges);

   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 4);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);
   STPTEST_ASSERT(reduce_starAllAreChecked(star));

   stptest_graphTearDown(scip, graph);
   reduce_starFree(scip, &star);

   return SCIP_OKAY;
}


/** STAR of degree 4 should be correctly computed */
static
SCIP_RETCODE testStar4IsCorrect(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   STAR* star;
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 8;
   int nstaredges = -1;
   const int* staredges;

   SCIP_CALL( reduce_starInit(scip, 5, &star) );
   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 6


   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   reduce_starReset(graph, 2, star);
   reduce_starReset(graph, 0, star);

   /* full star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);
   STPTEST_ASSERT(staredges[3] == 0);

   /* 1st degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);

   /* 2nd degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 3rd degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 4th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 4);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);

   STPTEST_ASSERT(reduce_starAllAreChecked(star));

   stptest_graphTearDown(scip, graph);
   reduce_starFree(scip, &star);

   return SCIP_OKAY;
}


/** STAR of degree 5 should be correctly computed */
static
SCIP_RETCODE testStar5IsCorrect(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   STAR* star;
   GRAPH* graph;
   int nnodes = 6;
   int nedges = 10;
   int nstaredges = -1;
   const int* staredges;

   SCIP_CALL( reduce_starInit(scip, 5, &star) );
   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 6
   graph_edge_addBi(scip, graph, 0, 5, 1.0); // 8

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   reduce_starReset(graph, 0, star);

   /* full star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 5);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 4);
   STPTEST_ASSERT(staredges[3] == 2);
   STPTEST_ASSERT(staredges[4] == 0);

   /* 1st degree 4 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 4);
   STPTEST_ASSERT(staredges[3] == 2);

   /* 2nd degree 4 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 4);
   STPTEST_ASSERT(staredges[3] == 0);

   /* 3rd degree 4 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 2);
   STPTEST_ASSERT(staredges[3] == 0);

   /* 4th degree 4 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);
   STPTEST_ASSERT(staredges[3] == 0);

   /* 5th degree 4 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 4);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);
   STPTEST_ASSERT(staredges[3] == 0);


   /* 1st degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 4);

   /* 2nd degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 2);

   /* 3rd degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 6);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 4th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);

   /* 5th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 6th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 8);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 7th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 2);

   /* 8th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 4);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 9th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 6);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);

   /* 10th degree 3 star */
   staredges = reduce_starGetNext(star, &nstaredges);
   STPTEST_ASSERT(nstaredges == 3);
   STPTEST_ASSERT(staredges[0] == 4);
   STPTEST_ASSERT(staredges[1] == 2);
   STPTEST_ASSERT(staredges[2] == 0);

   STPTEST_ASSERT(reduce_starAllAreChecked(star));

   stptest_graphTearDown(scip, graph);
   reduce_starFree(scip, &star);

   return SCIP_OKAY;
}


/** tests that BLC works for a degree 3 star */
static
SCIP_RETCODE testBLCworksFor3Star(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   BLCTREE* blctree;
   GRAPH* graph;
   int nnodes = 4;
   int nedges = 10;
   int mstedges[3];
   SCIP_Real bottlenecks[3];

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 1, 2, 2.0); // 6
   graph_edge_addBi(scip, graph, 2, 3, 3.0); // 8

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_blctreeInit(scip, graph, &blctree) );

   reduce_blctreeGetMstEdges(graph, blctree, mstedges);
   reduce_blctreeGetMstBottlenecks(graph, blctree, bottlenecks);

   for( int i = 0; i < 3; i++ )
   {
      const int edge = mstedges[i];
      assert(0 <= edge && edge < nedges);

      if( (edge / 2) * 2 == 0 )
      {
         STPTEST_ASSERT(EQ(bottlenecks[i], 2.0));
      }
      else if( (edge / 2) * 2 == 2 )
      {
         STPTEST_ASSERT(EQ(bottlenecks[i], 2.0));
      }
      else
      {
         STPTEST_ASSERT((edge / 2) * 2 == 4);
         STPTEST_ASSERT(EQ(bottlenecks[i], 3.0));
      }
   }

   stptest_graphTearDown(scip, graph);
   reduce_blctreeFree(scip, &blctree);

   return SCIP_OKAY;
}


/** tests that BLC works for a tree with five nodes */
static
SCIP_RETCODE testBLCworksFor5Tree(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   BLCTREE* blctree;
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 16;
   int mstedges[4];
   SCIP_Real bottlenecks[4];

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;
   graph_knot_chg(graph, 3, STP_TERM);

   /* the MST edges */
   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.1); // 2
   graph_edge_addBi(scip, graph, 2, 3, 1.2); // 4
   graph_edge_addBi(scip, graph, 2, 4, 1.3); // 6

   /* additional edges */
   graph_edge_addBi(scip, graph, 0, 4, 1.4);
   graph_edge_addBi(scip, graph, 1, 2, 1.5);
   graph_edge_addBi(scip, graph, 3, 4, 1.6);
   graph_edge_addBi(scip, graph, 1, 3, 1.7);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_blctreeInit(scip, graph, &blctree) );

   reduce_blctreeGetMstEdges(graph, blctree, mstedges);
   reduce_blctreeGetMstBottlenecks(graph, blctree, bottlenecks);

   for( int i = 0; i < 4; i++ )
   {
      const int edge = mstedges[i];
      assert(0 <= edge && edge < nedges);

      if( (edge / 2) * 2 == 0 )
      {
         STPTEST_ASSERT(EQ(bottlenecks[i], 1.5));
      }
      else if( (edge / 2) * 2 == 2 )
      {
         STPTEST_ASSERT(EQ(bottlenecks[i], 1.4));
      }
      else if( (edge / 2) * 2 == 4 )
      {
         STPTEST_ASSERT(EQ(bottlenecks[i], 1.6));
      }
      else
      {
         STPTEST_ASSERT((edge / 2) * 2 == 6);
         STPTEST_ASSERT(EQ(bottlenecks[i], 1.4));
      }
   }

   stptest_graphTearDown(scip, graph);
   reduce_blctreeFree(scip, &blctree);

   return SCIP_OKAY;
}


/** tests DCMST */
SCIP_RETCODE stptest_dcmst(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( dcmstTest1(scip) );
   SCIP_CALL( dcmstTest2(scip) );
   SCIP_CALL( dcmstTest2b(scip) );
   SCIP_CALL( dcmstTest3(scip) );

   printf("dcmst test: all ok \n");

   return SCIP_OKAY;
}


/** tests STAR methods */
SCIP_RETCODE stptest_reduceStar(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testStar3IsCorrect(scip) );
   SCIP_CALL( testStar4IsCorrect(scip) );
   SCIP_CALL( testStar5IsCorrect(scip) );

   printf("reduce star test: all ok \n");

   return SCIP_OKAY;
}


/** tests bottleneck tree methods */
SCIP_RETCODE stptest_reduceBLCtree(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testBLCworksFor3Star(scip) );
   SCIP_CALL( testBLCworksFor5Tree(scip) );

   printf("reduce BLC tree test: all ok \n");

   return SCIP_OKAY;
}
