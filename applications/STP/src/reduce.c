/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce.c
 * @brief  Reduction tests for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 *
 * This file includes several packages of reduction techniques for different Steiner problem variants.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h)           */
#define REDUCE_C
#define STP_RED_SDSPBOUND    200         /**< visited edges bound for SDSP test  */
#define STP_RED_SDSPBOUND2   500         /**< visited edges bound for SDSP test  */
#define STP_RED_BD3BOUND     400         /**< visited edges bound for BD3 test  */
#define STP_RED_EXTENSIVE FALSE
#define STP_RED_MWTERMBOUND 400
#define STP_RED_MAXNROUNDS 15

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "heur_tm.h"
#include "misc_stp.h"
#include "scip/scip.h"
#include "probdata_stp.h"
#include "prop_stp.h"

/** iterate NV and SL test while at least minelims many contractions are being performed */
static
SCIP_RETCODE nvreduce_sl(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 vnoi,
   SCIP_Real*            nodearrreal,
   SCIP_Real*            fixed,
   int*                  edgearrint,
   int*                  heap,
   int*                  state,
   int*                  vbase,
   int*                  neighb,
   int*                  distnode,
   int*                  solnode,
   STP_Bool*             visited,
   int*                  nelims,
   int                   minelims
   )
{
   int elims;
   int nvelims;
   int slelims;
   int degelims;
   int totalelims;

   assert(g != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(vnoi != NULL);
   assert(nodearrreal != NULL);
   assert(visited != NULL);
   assert(minelims >= 0);

   *nelims = 0;
   totalelims = 0;

   do
   {
      elims = 0;
      degelims = 0;

      /* NV-reduction */
      SCIP_CALL( reduce_nvAdv(scip, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, neighb, distnode, solnode, &nvelims) );
      elims += nvelims;

      SCIPdebugMessage("NV-reduction (in NVSL): %d \n", nvelims);

      /* SL-reduction */
      SCIP_CALL( reduce_sl(scip, g, vnoi, fixed, heap, state, vbase, neighb, visited, solnode, &slelims) );
      elims += slelims;

      SCIPdebugMessage("SL-reduction (in NVSL): %d \n", slelims);

      /* trivial reductions */
      if( elims > 0 )
      {
         if( g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG )
            SCIP_CALL( reduce_simple_pc(scip, g, fixed, &degelims, solnode, FALSE) );
         else
            SCIP_CALL( reduce_simple(scip, g, fixed, solnode, &degelims, NULL) );
      }
      else
      {
         degelims = 0;
      }

      elims += degelims;

      SCIPdebugMessage("Degree Test-reduction (in NVSL): %d \n", degelims);

      totalelims += elims;
   }while( elims > minelims );

   *nelims = totalelims;

   assert(graph_valid(g));

   return SCIP_OKAY;
}

/* remove unconnected vertices, overwrites g->mark */
SCIP_RETCODE level0(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   int k;
   int nnodes;

   assert(scip != NULL);
   assert(g != NULL);

   nnodes = g->knots;

   for( k = nnodes - 1; k >= 0 ; k-- )
      g->mark[k] = FALSE;

   SCIP_CALL( graph_trail_arr(scip, g, g->source) );

   for( k = nnodes - 1; k >= 0 ; k-- )
   {
      if( !g->mark[k] && (g->grad[k] > 0) )
      {
         assert(!Is_term(g->term[k]));

         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
      }
   }
   return SCIP_OKAY;
}

/* remove unconnected vertices, keep g->mark */
SCIP_RETCODE level0save(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   int* savemark;
   const int nnodes = g->knots;

   assert(scip != NULL);
   assert(g != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &savemark, nnodes) );
   BMScopyMemoryArray(savemark, g->mark, nnodes);

   for( int k = nnodes - 1; k >= 0; k-- )
      g->mark[k] = FALSE;

   SCIP_CALL( graph_trail_arr(scip, g, g->source) );

   for( int k = nnodes - 1; k >= 0 ; k-- )
   {
      if( !g->mark[k] && (g->grad[k] > 0) )
      {
         assert(!Is_term(g->term[k]));

         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
      }
   }

   BMScopyMemoryArray(g->mark, savemark, nnodes);

   SCIPfreeBufferArray(scip, &savemark);

   return SCIP_OKAY;
}

/** basic reduction package for the STP */
SCIP_RETCODE reduceStp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             dualascent,         /**< perform dualascent reductions? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int*                  edgestate,          /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   PATH* vnoi;
   PATH* path;
   GRAPH* g;
   GNODE** gnodearr;
   SCIP_Real*  nodearrreal;
   SCIP_Real*  edgearrreal;
   SCIP_Real*  edgearrreal2;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     i;
   int     nnodes;
   int     nedges;
   int     nterms;
   int     reductbound;
   STP_Bool* nodearrchar;

   SCIP_Bool    bred = FALSE;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(*graph != NULL);
   assert(minelims >= 0);

   g = *graph;

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, 0.03) )
      bred = TRUE;

   if( dualascent )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
      }
   }
   else
   {
      gnodearr = NULL;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   if( bred || dualascent )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, nedges) );
   }
   else
   {
      edgearrreal2 = NULL;
   }

   reductbound = MAX(nedges / 1000, minelims);

   /* reduction loop */
   SCIP_CALL( redLoopStp(scip, g, vnoi, path, gnodearr, nodearrreal, edgearrreal, edgearrreal2, heap, state,
         vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, fixed, -1.0, dualascent, bred, nodereplacing, reductbound, edgestate, userec) );

   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n", *fixed);

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &edgearrreal2);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &edgearrreal);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &edgearrint);

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   return SCIP_OKAY;
}

/** basic reduction package for the (R)PCSTP */
static
SCIP_RETCODE reducePc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   STP_Bool              dualascent,         /**< perform dual ascent reductions? */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   PATH* vnoi;
   PATH* path;
   GRAPH* g = *graph;
   GNODE** gnodearr;
   SCIP_Real* exedgearrreal;
   SCIP_Real* nodearrreal;
   SCIP_Real* exedgearrreal2;
   SCIP_Real timelimit;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     i;
   int     nnodes;
   int     nterms;
   int     nedges;
   int     extnedges;
   int     reductbound;
   STP_Bool*   nodearrchar;
   SCIP_Bool    bred = FALSE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;

   /* for PCSPG more memory is necessary */
   if( g->stp_type == STP_RPCSPG || !dualascent )
      extnedges = nedges;
   else
      extnedges = nedges + 2 * (g->terms - 1);

   /* get timelimit parameter*/
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exedgearrreal, extnedges ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, 0.03) )
      bred = TRUE;

   if( bred || dualascent )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &exedgearrreal2, extnedges) );
   }
   else
   {
      exedgearrreal2 = NULL;
   }

   if( dualascent )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, extnedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
      }
   }
   else
   {
      gnodearr = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   }

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = MAX(nnodes / 1000, minelims);

   /* reduction loop */
   SCIP_CALL( redLoopPc(scip, g, vnoi, path, gnodearr, nodearrreal, exedgearrreal, exedgearrreal2, heap, state,
         vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, fixed, dualascent, bred, reductbound, userec) );

   /* free memory */

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }
   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArrayNull(scip, &exedgearrreal2);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &exedgearrreal);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);

   return SCIP_OKAY;
}

/** reduction package for the MWCSP */
static
SCIP_RETCODE reduceMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   STP_Bool              advanced,           /**< perform advanced reductions? */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   GRAPH* g = *graph;
   PATH* vnoi;
   PATH* path;
   GNODE** gnodearr;
   SCIP_Real* nodearrreal;
   SCIP_Real* edgearrreal;
   SCIP_Real* edgearrreal2;
   int* state;
   int* vbase;
   int* edgearrint;
   int* nodearrint;
   int* nodearrint2;
   int* nodearrint3;
   int i;
   int nterms;
   int nnodes;
   int nedges;
   int redbound;
   int extnedges;
   STP_Bool* nodearrchar;
   STP_Bool bred = FALSE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   nnodes = g->knots;
   nedges = g->edges;
   nterms = g->terms;
   redbound = MAX(nnodes / 1000, minelims);

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, 0.1) )
      bred = TRUE;

   if( advanced )
   {
      extnedges = nedges + 2 * (g->terms - 1);

      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, extnedges) );
   }
   else
   {
      extnedges = nedges;
      edgearrint = NULL;
      gnodearr = NULL;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint3, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes + 1) );

   if( bred || advanced )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, extnedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, extnedges) );
   }
   else
   {
      nodearrreal = NULL;
      edgearrreal = NULL;
      edgearrreal2 = NULL;
   }

   /* reduction loop */
   SCIP_CALL( redLoopMw(scip, g, vnoi, path, gnodearr, nodearrreal, edgearrreal, edgearrreal2, state,
         vbase, nodearrint, edgearrint, nodearrint2, nodearrint3, NULL, nodearrchar, fixed, advanced, bred, advanced, redbound, userec) );

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &edgearrreal2);
   SCIPfreeBufferArrayNull(scip, &edgearrreal);
   SCIPfreeBufferArrayNull(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint3);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArrayNull(scip, &edgearrint);

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   return SCIP_OKAY;
}

/** basic reduction package for the HCDSTP */
static
SCIP_RETCODE reduceHc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   GRAPH* g = *graph;
   PATH* vnoi;
   SCIP_Real*  cost;
   SCIP_Real*  radius;
   SCIP_Real*  costrev;
   SCIP_Real timelimit;
   SCIP_Real upperbound;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    pathedge;
   int     nnodes;
   int     nedges;
   int     redbound;
#if 0
   int     danelims;
#endif
   int     degnelims;
   int     brednelims;
   int     hbrednelims;
   int     hcrnelims;
   int     hcrcnelims;
   STP_Bool*   nodearrchar;
#if 0
   DOES NOT WORK for HC!
      STP_Bool    da = !TRUE;
#endif
   STP_Bool    bred = TRUE;
   STP_Bool    hbred = TRUE;
   STP_Bool    rbred = TRUE;
   STP_Bool    rcbred = TRUE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nnodes = g->knots;
   nedges = g->edges;
   degnelims = 0;
   upperbound = -1.0;
   redbound = MAX(g->knots / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &radius, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );

   SCIP_CALL( reduce_simple_hc(scip, g, fixed, &degnelims) );

   while( (bred || hbred || rbred || rcbred) && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      upperbound = -1.0;

      if( rbred )
      {
         SCIP_CALL( reduce_boundHopR(scip, g, vnoi, cost, costrev, radius, heap, state, vbase, &hcrnelims, pathedge) );
         if( hcrnelims <= redbound )
            rbred = FALSE;
      }

      if( rcbred )
      {
         SCIP_CALL( reduce_boundHopRc(scip, g, vnoi, cost, costrev, radius, -1.0, heap, state, vbase, &hcrcnelims, pathedge, FALSE) );
         if( hcrcnelims <= redbound )
            rcbred = FALSE;
      }

      if( bred )
      {
         SCIP_CALL( reduce_bound(scip, g, vnoi, cost, NULL, radius, costrev, fixed, &upperbound, heap, state, vbase, &brednelims) );
         if( brednelims <= redbound )
            bred = FALSE;
      }

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( hbred )
      {
         SCIP_CALL( reduce_boundHop(scip, g, vnoi, cost, radius, costrev, heap, state, vbase, &hbrednelims) );
         if( hbrednelims <= redbound )
            hbred = FALSE;
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &radius);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);

   return SCIP_OKAY;
}

/** basic reduction package for the SAP */
static
SCIP_RETCODE reduceSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH*   vnoi;
   PATH*   path;
   SCIP_Real ub = FARAWAY;
   SCIP_Real timelimit;
   SCIP_Real*  nodearrreal;
   SCIP_Real*  edgearrreal;
   SCIP_Real*  edgearrreal2;
   GRAPH*  g = *graph;
   GNODE** gnodearr;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     e;
   int     i;
   int     nnodes;
   int     nedges;
   int     nterms;
   int     danelims;
   int     sdnelims;
   int     rptnelims;
   int     degtnelims;
   int     redbound;
   STP_Bool    da = TRUE;
   STP_Bool    sd = !TRUE;
   STP_Bool* nodearrchar;
   STP_Bool    rpt = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1) );

   nnodes = g->knots;
   nedges = g->edges;
   nterms = g->terms;

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, nedges) );

   /* @todo change .stp file format for SAP! */
   for( e = 0; e < g->edges; e++ )
      if( SCIPisEQ(scip, g->cost[e], 20000.0) )
         g->cost[e] = FARAWAY;

   SCIP_CALL( reduce_simple_sap(scip, g, fixed, &degtnelims) );

   /* main loop */
   while( (sd || rpt || da) && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( sd )
      {
         SCIP_CALL( reduce_sdspSap(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdnelims, 300) );
         if( sdnelims <= redbound )
            sd = FALSE;
      }

      if( rpt )
      {
         SCIP_CALL( reduce_rpt(scip, g, fixed, &rptnelims) );
         if( rptnelims <= redbound )
            rpt = FALSE;
      }

      SCIP_CALL( reduce_simple_sap(scip, g, fixed, &degtnelims) );

      if( da )
      {
         SCIP_CALL( reduce_da(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, &ub, fixed, edgearrint, vbase, state, heap, nodearrint,
               nodearrint2, nodearrchar, &danelims, 0, randnumgen, TRUE, NULL, FALSE) );

         if( danelims <= 2 * redbound )
            da = FALSE;
      }
   }

   SCIP_CALL( reduce_simple_sap(scip, g, fixed, &degtnelims) );

   SCIPfreeBufferArray(scip, &edgearrreal2);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &edgearrreal);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &edgearrint);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBlockMemory(scip, &gnodearr[i]);
   SCIPfreeBufferArray(scip, &gnodearr);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}


static
SCIP_RETCODE reduceNw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH*   vnoi;
   SCIP_Real*  nodearrreal;
   SCIP_Real*  edgearrreal;
   SCIP_Real*  edgearrreal2;
   SCIP_Real   ub = FARAWAY;
   SCIP_Real   timelimit;
   GRAPH*  g = *graph;
   GNODE** gnodearr;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     i;
   int     nnodes;
   int     nedges;
   int     nterms;
   int     danelims;
   int     redbound;

   STP_Bool*   nodearrchar;
   STP_Bool    da = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1) );

   nnodes = g->knots;
   nedges = g->edges;
   nterms = g->terms;

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
   for( i = 0; i < nterms - 1; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, nedges) );

   while( (da) && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      SCIP_CALL( reduce_da(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, &ub, fixed, edgearrint, vbase, state, heap, nodearrint,
            nodearrint2, nodearrchar, &danelims, 0, randnumgen, TRUE, NULL, FALSE) );

      if( danelims <= 2 * redbound )
         da = FALSE;
   }

   SCIPfreeBufferArray(scip, &edgearrreal2);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &edgearrreal);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &edgearrint);

   for( i = nterms - 2; i >= 0; i-- )
      SCIPfreeBlockMemory(scip, &gnodearr[i]);
   SCIPfreeBufferArray(scip, &gnodearr);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** MWCS loop */
SCIP_RETCODE redLoopMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   SCIP_Real*            edgearrreal,        /**< edges-sized array  */
   SCIP_Real*            edgearrreal2,       /**< edges-sized array  */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< voronoi base array  */
   int*                  nodearrint,         /**< nodes-sized array  */
   int*                  edgearrint,         /**< edges-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   int*                  nodearrint3,        /**< nodes-sized array  */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   STP_Bool              advanced,           /**< do advanced reduction? */
   STP_Bool              bred,               /**< do bound-based reduction? */
   STP_Bool              tryrmw,             /**< try to convert problem to RMWCSP? Only possible if advanced = TRUE and userec = TRUE */
   int                   redbound,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   SCIP_Real timelimit;
   int daelims;
   int anselims;
   int nnpelims;
   int degelims;
   int npvelims;
   int bredelims;
   int ansadelims;
   int ansad2elims;
   int chain2elims;

   STP_Bool da = advanced;
   STP_Bool ans = TRUE;
   STP_Bool nnp = TRUE;
   STP_Bool npv = TRUE;
   STP_Bool rerun = TRUE;
   STP_Bool ansad = TRUE;
   STP_Bool ansad2 = TRUE;
   STP_Bool chain2 = TRUE;
   STP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real prizesum;


   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(advanced || !tryrmw);

   tryrmw = tryrmw && userec;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1) );

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   graph_pc_2org(g);

   degelims = 0;

   SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
   assert(graph_pc_term2edgeConsistent(g));

   prizesum = graph_pc_getPosPrizeSum(scip, g);

   while( rerun && !SCIPisStopped(scip) )
   {
      daelims = 0;
      anselims = 0;
      nnpelims = 0;
      degelims = 0;
      npvelims = 0;
      bredelims = 0;
      ansadelims = 0;
      ansad2elims = 0;
      chain2elims = 0;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( ans || extensive )
      {
         reduce_ans(scip, g, nodearrint2, &anselims);

         if( anselims <= redbound )
            ans = FALSE;

         SCIPdebugMessage("ans deleted: %d \n", anselims);
      }

      if( ansad || extensive )
      {
         reduce_ansAdv(scip, g, nodearrint2, &ansadelims, FALSE);

         if( ansadelims <= redbound )
            ansad = FALSE;

         SCIPdebugMessage("ans advanced deleted: %d \n", ansadelims);
      }

      if( ans || ansad || nnp || npv || extensive )
         SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

      if( (da || (advanced && extensive)) )
      {
         SCIP_CALL( reduce_daPcMw(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, vbase, nodearrint, edgearrint,
               state, nodearrchar, &daelims, TRUE, FALSE, FALSE, FALSE, userec, randnumgen, prizesum) );

         if( daelims <= 2 * redbound )
            da = FALSE;
         else
            SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

         SCIPdebugMessage("Dual-Ascent Elims: %d \n", daelims);
      }

      if( nnp )
      {
         reduce_nnp(scip, g, nodearrint2, &nnpelims);

         if( nnpelims <= redbound )
            nnp = FALSE;

         SCIPdebugMessage("nnp deleted: %d \n", nnpelims);
      }

      if( nnp || extensive )
      {
         SCIP_CALL(reduce_chain2(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &chain2elims, 500));

         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 delete: %d \n", chain2elims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( npv || extensive )
      {
         SCIP_CALL(reduce_npv(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &npvelims, 400));

         if( npvelims <= redbound )
            npv = FALSE;

         SCIPdebugMessage("npv delete: %d \n", npvelims);
      }

      if( chain2 || extensive )
      {
         SCIP_CALL(reduce_chain2(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &chain2elims, 300));

         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 delete: %d \n", chain2elims);
      }

      SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

      if( ansad2 || extensive )
      {
         reduce_ansAdv2(scip, g, nodearrint2, &ansad2elims);

         if( ansad2elims <= redbound )
            ansad2 = FALSE;
         else
            SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &ansad2elims) );

         SCIPdebugMessage("ans advanced 2 deleted: %d (da? %d ) \n", ansad2elims, da);
      }

      if( bred )
      {
         SCIP_CALL( reduce_boundMw(scip, g, vnoi, path, edgearrreal, nodearrreal, edgearrreal2, fixed, nodearrint, state, vbase, NULL, &bredelims) );

         if( bredelims <= redbound )
            bred = FALSE;


         SCIPdebugMessage("reduce_bound: %d \n", bredelims);
      }

      if( anselims + nnpelims + chain2elims + bredelims + npvelims + ansadelims + ansad2elims + daelims <= redbound )
         rerun = FALSE;

      if( !rerun && advanced && g->terms > 2 )
      {
         int cnsadvelims = 0;

         SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
#if 0
         SCIP_CALL( reduce_cnsAdv(scip, g, nodearrint2, &cnsadvelims) );
#endif

         SCIP_CALL( reduce_daPcMw(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, vbase, nodearrint, edgearrint,
               state, nodearrchar, &daelims, TRUE, (g->terms > STP_RED_MWTERMBOUND), FALSE, tryrmw, userec, randnumgen, prizesum) );

         userec = FALSE;

         if( cnsadvelims + daelims >= redbound || (extensive && (cnsadvelims + daelims > 0))  )
         {
            ans = TRUE;
            nnp = TRUE;
            npv = TRUE;
            ansad = TRUE;
            ansad2 = TRUE;
            chain2 = TRUE;
            rerun = TRUE;
            advanced = FALSE;

            SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
            SCIPdebugMessage("Restarting reduction loop! (%d eliminations) \n\n ", cnsadvelims + daelims);
            if( extensive )
               advanced = TRUE;
         }
      }
   }

   SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

#if 0
   if( tryrmw )
   {
      SCIP_CALL( reduce_daPcMw(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, vbase, nodearrint, edgearrint,
            state, nodearrchar, &daelims, TRUE, FALSE, FALSE, TRUE, userec, randnumgen, prizesum) );

      SCIP_CALL( reduce_npv(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &npvelims, 400) );
      reduce_ans(scip, g, nodearrint2, &anselims);
      reduce_ansAdv(scip, g, nodearrint2, &ansadelims, FALSE);
   }
#endif

   /* go back to the extended graph */
   graph_pc_2trans(g);

   SCIP_CALL( level0(scip, g) );

   if( tryrmw )
      SCIP_CALL( graph_pc_mw2rmw(scip, g, prizesum) );

   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** (R)PC loop */
SCIP_RETCODE redLoopPc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   SCIP_Real*            exedgearrreal,      /**< edges-sized array  */
   SCIP_Real*            exedgearrreal2,     /**< edges-sized array  */
   int*                  heap,               /**< shortest path array  */
   int*                  state,              /**< voronoi base array  */
   int*                  vbase,              /**< nodes-sized array  */
   int*                  nodearrint,         /**< edges-sized array  */
   int*                  edgearrint,         /**< nodes-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   int*                  solnode,            /**< solution nodes array (or NULL) */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   SCIP_Bool             dualascent,         /**< do dual-ascent reduction? */
   SCIP_Bool             bred,               /**< do bound-based reduction? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   SCIP_Real ub;
   SCIP_Real fix;
   SCIP_Real timelimit;
   SCIP_Real rootprize = 0.0;
   const SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   int nelims;
   int danelims;
   int sdnelims;
   int sdcnelims;
   int bd3nelims;
   int degnelims;
   int nvslnelims;
   int brednelims;
   SCIP_Bool da = dualascent;
   SCIP_Bool sd = TRUE;
   SCIP_Bool sdc = TRUE;
   SCIP_Bool bd3 = TRUE;
   SCIP_Bool nvsl = TRUE;
   SCIP_Bool rerun = TRUE;
   SCIP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_Bool advancedrun = dualascent;
   SCIP_Real prizesum;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1) );

   if( rpc )
   {
      rootprize = g->prize[g->source];
      g->prize[g->source] = FARAWAY;
   }

   ub = -1.0;
   fix = 0.0;

   graph_pc_2org(g);

   assert(graph_pc_term2edgeConsistent(g));

   SCIP_CALL( reduce_simple_pc(scip, g, &fix, &degnelims, solnode, FALSE) );

   assert(graph_pc_term2edgeConsistent(g));

   prizesum = graph_pc_getPosPrizeSum(scip, g);

   /* get timelimit parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   for( int rounds = 0; rounds < STP_RED_MAXNROUNDS && !SCIPisStopped(scip) && rerun; rounds++ )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      nelims = 0;
      danelims = 0;
      sdnelims = 0;
      sdcnelims = 0;
      bd3nelims = 0;
      nvslnelims = 0;
      degnelims = 0;
      brednelims = 0;

      if( sd || extensive )
      {
         SCIP_CALL( reduce_sdPc(scip, g, vnoi, heap, state, vbase, nodearrint, nodearrint2, &sdnelims) );

         if( sdnelims <= reductbound )
            sd = FALSE;

         SCIPdebugMessage("SDpc: %d \n", sdnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sdc || extensive )
      {
         SCIP_CALL( reduce_sdsp(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims,
               ((rounds > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND), NULL) );

         if( sdcnelims <= reductbound )
            sdc = FALSE;

         SCIPdebugMessage("SDsp: %d \n", sdcnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      SCIP_CALL( reduce_simple_pc(scip, g, &fix, &nelims, solnode, FALSE) );

      degnelims += nelims;

      if( bd3 && dualascent )
      {
         SCIP_CALL( reduce_bd34(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &bd3nelims, STP_RED_BD3BOUND) );
         if( bd3nelims <= reductbound )
         {
            bd3 = FALSE;
         }
         else if( !rpc )
         {
            SCIP_CALL( reduce_sdPc(scip, g, vnoi, heap, state, vbase, nodearrint, nodearrint2, &sdnelims) );
            SCIP_CALL( reduce_sdsp(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims, ((rounds > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND), NULL) );
         }

         SCIPdebugMessage("bd3: %d \n", bd3nelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( nvsl || extensive )
      {
         SCIP_CALL( nvreduce_sl(scip, g, vnoi, nodearrreal, &fix, edgearrint, heap, state, vbase, nodearrint, nodearrint2, solnode, nodearrchar, &nvslnelims, reductbound) );

         if( nvslnelims <= 0.5 * reductbound )
            nvsl = FALSE;

         SCIPdebugMessage("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bred )
      {
         ub = -1.0;
         SCIP_CALL( reduce_bound(scip, g, vnoi, exedgearrreal, g->prize, nodearrreal, exedgearrreal2, &fix, &ub, heap, state, vbase, &brednelims) );
         if( brednelims <= reductbound )
            bred = FALSE;

         SCIPdebugMessage("bndelims %d \n", brednelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      ub = -1.0;

      assert(graph_pc_term2edgeConsistent(g));

      if( da || (dualascent && extensive) )
      {
         if( rpc )
            SCIP_CALL( reduce_da(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, &ub, &fix, edgearrint, vbase, state, heap,
                  nodearrint, nodearrint2, nodearrchar, &danelims, 0, randnumgen, TRUE, NULL, FALSE) );
         else
            SCIP_CALL( reduce_daPcMw(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, vbase, heap, edgearrint,
                  state, nodearrchar, &danelims, TRUE, FALSE, FALSE, FALSE, userec, randnumgen, prizesum) );

         if( danelims <= reductbound )
            da = FALSE;

         SCIPdebugMessage("da: %d \n", danelims);
      }

      SCIP_CALL( reduce_simple_pc(scip, g, &fix, &degnelims, solnode, TRUE) );

      if( ub >= 0 )
      {
         SCIP_CALL( reduce_bound(scip, g, vnoi, exedgearrreal, g->prize, nodearrreal, exedgearrreal2, &fix, &ub, heap, state, vbase, &brednelims) );
         if( brednelims <= reductbound )
            bred = FALSE;

         SCIPdebugMessage("bndelims %d \n", brednelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( degnelims + sdnelims + sdcnelims + bd3nelims + danelims + brednelims + nvslnelims <= reductbound )
         rerun = FALSE;

      if( !rerun && advancedrun && g->terms > 2 )
      {
         rerun = TRUE;
         danelims = 0;
         degnelims = 0;
         advancedrun = FALSE;
         if( rpc )
         {
            SCIP_CALL( reduce_da(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, &ub, &fix, edgearrint, vbase, state, heap,
                  nodearrint, nodearrint2, nodearrchar, &danelims, 0, randnumgen, TRUE, NULL, FALSE) );
         }
         else
         {
            SCIP_CALL( reduce_daPcMw(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, vbase, heap, edgearrint,
                  state, nodearrchar, &danelims, TRUE, TRUE, FALSE, FALSE, userec, randnumgen, prizesum) );
         }
         SCIP_CALL( reduce_simple_pc(scip, g, &fix, &degnelims, solnode, TRUE) );
         if( danelims + degnelims > reductbound || (extensive && (danelims + degnelims > 0)) )
         {
            da = dualascent;
            sd = TRUE;
            sdc = TRUE;
            bd3 = TRUE;
            nvsl = TRUE;
            if( extensive )
               advancedrun = TRUE;
         }
         else
         {
            rerun = FALSE;
         }
      }
   }
   SCIP_CALL( reduce_simple_pc(scip, g, &fix, &degnelims, solnode, TRUE) );

   assert(graph_pc_term2edgeConsistent(g));


   if( rpc )
      g->prize[g->source] = rootprize;

   graph_pc_2trans(g);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   *fixed += fix;

   SCIPdebugMessage("Reduction Level PC 1: Fixed Cost = %.12e\n", *fixed);
   return SCIP_OKAY;
}

/** STP loop */
SCIP_RETCODE redLoopStp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   SCIP_Real*            edgearrreal,        /**< edges-sized array  */
   SCIP_Real*            edgearrreal2,       /**< edges-sized array  */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< Voronoi base array  */
   int*                  nodearrint,         /**< edges-sized array  */
   int*                  edgearrint,         /**< nodes-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   int*                  solnode,            /**< solution nodes array (or NULL) */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   SCIP_Real             upperbound,         /**< upper bound */
   SCIP_Bool             da,                 /**< do dual-ascent reduction? */
   SCIP_Bool             bred,               /**< do bound-based reduction? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   int*                  edgestate,          /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   SCIP_Real    ub;
   SCIP_Real    fix;
   SCIP_Real    timelimit;
   SCIP_Bool    le = TRUE;
   SCIP_Bool    sd = TRUE;
   SCIP_Bool    sdc = TRUE;
   SCIP_Bool    bd3 = nodereplacing;
   SCIP_Bool    nvsl = nodereplacing;
   SCIP_Bool    rerun = TRUE;
   const SCIP_Bool extensive = STP_RED_EXTENSIVE;
   int i = 0;
   int rounds = 0;
   SCIP_RANDNUMGEN* randnumgen;

   assert(graph_valid(g));

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1) );

   ub = upperbound;
   fix = 0.0;

   SCIP_CALL( reduce_contractZeroEdges(scip, g) );

   SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &i, edgestate) );

   /* get timelimit parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   while( rerun && !SCIPisStopped(scip) )
   {
      int danelims = 0;
      int lenelims = 0;
      int sdnelims = 0;
      int sdcnelims = 0;
      int bd3nelims = 0;
      int nvslnelims = 0;
      int brednelims = 0;
      int degtnelims = 0;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( le || extensive )
      {
         SCIP_CALL( reduce_ledge(scip, g, vnoi, heap, state, vbase, &lenelims, edgestate) );

         if( lenelims <= reductbound )
            le = FALSE;
         else
            SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &degtnelims, edgestate) );

         SCIPdebugMessage("le: %d \n", lenelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd || extensive  )
      {
         SCIP_CALL( reduce_sd(scip, g, vnoi, edgearrreal, nodearrreal, heap, state, vbase, nodearrint, nodearrint2, edgearrint, &sdnelims, nodereplacing, edgestate) );

         if( sdnelims <= reductbound )
            sd = FALSE;

         SCIPdebugMessage("sd: %d, \n", sdnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sdc || extensive  )
      {
         SCIP_CALL( reduce_sdsp(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims,
               ((rounds > 0) ? (STP_RED_SDSPBOUND2 / 2) : (STP_RED_SDSPBOUND / 2)), edgestate) );

         if( sdcnelims <= reductbound )
            sdc = FALSE;

         SCIPdebugMessage("sdsp: %d \n", sdcnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd || sdc )
         SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &degtnelims, edgestate) );

      if( bd3 || extensive )
      {
         SCIP_CALL( reduce_bd34(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &bd3nelims, STP_RED_BD3BOUND) );
         if( bd3nelims <= reductbound )
            bd3 = FALSE;
         else
         {
            SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &degtnelims, edgestate) );
         }

         SCIPdebugMessage("bd3: %d \n", bd3nelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( nvsl || extensive  )
      {
         SCIP_CALL( nvreduce_sl(scip, g, vnoi, nodearrreal, &fix, edgearrint, heap, state, vbase, nodearrint, NULL, solnode, nodearrchar, &nvslnelims, reductbound) );

         if( nvslnelims <= reductbound )
            nvsl = FALSE;

         SCIPdebugMessage("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      ub = -1.0;

      if( da )
      {
         SCIP_CALL( reduce_da(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, &ub, &fix, edgearrint, vbase, state, heap,
               nodearrint, nodearrint2, nodearrchar, &danelims, rounds, randnumgen, nodereplacing, edgestate, userec) );

         if( danelims <= 3 * reductbound )
            da = FALSE;

         SCIPdebugMessage("da: %d \n", danelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bred && nodereplacing )
      {
         SCIP_CALL( reduce_bound(scip, g, vnoi, edgearrreal, NULL, nodearrreal, edgearrreal2, &fix, &ub, heap, state, vbase, &brednelims) );

         SCIP_CALL( level0(scip, g) );

         if( brednelims <= 2 * reductbound )
            bred = FALSE;

         SCIPdebugMessage("bnd: %d \n\n", brednelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }
      SCIP_CALL( level0(scip, g) );
      SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &degtnelims, edgestate) );

      if( (danelims + sdnelims + bd3nelims + nvslnelims + lenelims + brednelims + sdcnelims) <= 2 * reductbound  )
         rerun = FALSE;

      if( extensive && (danelims + sdnelims + bd3nelims + nvslnelims + lenelims + brednelims + sdcnelims) > 0 )
         rerun = TRUE;

      rounds++;
   }

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   *fixed += fix;

   return SCIP_OKAY;
}


/** reduces the graph */
SCIP_RETCODE reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph structure */
   SCIP_Real*            offset,             /**< pointer to store offset generated by reductions */
   int                   level,              /**< reduction level 0: none, 1: basic, 2: advanced */
   int                   minelims,           /**< minimal amount of reductions to reiterate reduction methods */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   int stp_type;

   assert((*graph)      != NULL);
   assert((*graph)->fixedges == NULL);
   assert(level  >= 0 && level <= 2);
   assert(minelims >= 0);
   assert((*graph)->layers == 1);

   *offset = 0.0;

   stp_type = (*graph)->stp_type;

   /* initialize ancestor list for each edge */
   SCIP_CALL( graph_init_history(scip, (*graph)) );

   /* initialize shortest path algorithms */
   SCIP_CALL( graph_path_init(scip, (*graph)) );

   SCIP_CALL( level0(scip, (*graph)) );

   /* if no reduction methods available, return */
   if( (*graph)->stp_type == STP_DCSTP || (*graph)->stp_type == STP_RMWCSP )
   {
      graph_path_exit(scip, (*graph));
      return SCIP_OKAY;
   }

   if( level == 1 )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reducePc(scip, (graph), offset, minelims, FALSE, FALSE) );
      }
      else if( stp_type == STP_MWCSP )
      {
         SCIP_CALL( reduceMw(scip, (graph), offset, minelims, FALSE, FALSE) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduceHc(scip, (graph), offset, minelims) );
      }
      else if( stp_type == STP_SAP )
      {
         SCIP_CALL( reduceSap(scip, (graph), offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduceNw(scip, (graph), offset, minelims) );
      }
      else
      {
         SCIP_CALL( reduceStp(scip, (graph), offset, minelims, FALSE, TRUE, NULL, FALSE) );
      }
   }
   else if( level == 2 )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reducePc(scip, (graph), offset, minelims, TRUE, userec) );
      }
      else if( stp_type == STP_MWCSP )
      {
         SCIP_CALL( reduceMw(scip, (graph), offset, minelims, TRUE, userec) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduceHc(scip, (graph), offset, minelims) );
      }
      else if( stp_type == STP_SAP )
      {
         SCIP_CALL( reduceSap(scip, (graph), offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduceNw(scip, (graph), offset, minelims) );
      }
      else
      {
         SCIP_CALL( reduceStp(scip, (graph), offset, minelims, TRUE, TRUE, NULL, userec) );
      }
   }
   SCIPdebugMessage("offset : %f \n", *offset);

   assert(graph_valid(*graph));

   graph_path_exit(scip, (*graph));

   return SCIP_OKAY;
}
