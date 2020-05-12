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

/**@file   reduce_base.c
 * @brief  Reduction tests for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes several packages of reduction techniques for different Steiner problem variants.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h)           */
#define REDUCE_C
#define STP_RED_SDSPBOUND    200         /**< visited edges bound for SDSP test  */
#define STP_RED_SDSPBOUND2   1000        /**< visited edges bound for SDSP test  */
#define STP_RED_BD3BOUND     500         /**< visited edges bound for BD3 test  */
#define STP_RED_EXTENSIVE FALSE
#define STP_RED_MWTERMBOUND 400
#define STP_RED_MAXNROUNDS 15
#define STP_RED_EXFACTOR   2
#define STP_RED_GLBFACTOR  2
#define STP_RED_EDGELIMIT 100000
#define STP_BND_THRESHOLD 0.03

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "misc_stp.h"
#include "scip/scip.h"
#include "probdata_stp.h"
#include "prop_stp.h"


SCIP_Bool show;

enum PC_REDTYPE {pc_sdc, pc_sdstar, pc_sdw1, pc_sdw2, pc_bd3};


static
int getWorkLimits_mw(
   const GRAPH* g,
   int roundnumber,
   enum PC_REDTYPE redtype
)
{
   // todo: fill for chain2 and npv2 perhaps
   return 0;
}


static
int getWorkLimits_pc(
   const GRAPH* g,
   int roundnumber,
   enum PC_REDTYPE redtype
)
{
   const int nedges = g->edges;
   const int sqrtnedges = (int) sqrt(nedges);
   int limit = 0;
   int todo; // pc_bd3 far too aggressive! for handbi04 takes a lot of time...try to go down to STP_RED_SDSPBOUND2 / 2? or / 3 ?

   assert(roundnumber >= 0);
   assert(nedges >= 0);
   assert(sqrtnedges >= 0 && sqrtnedges <= nedges);

   switch (redtype)
   {
   case pc_sdc:
      limit = (roundnumber > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND;
      break;
   case pc_sdstar:
      limit = (roundnumber > 0) ? STP_RED_SDSPBOUND2 : 2 * STP_RED_SDSPBOUND;
      break;
   case pc_sdw1:
      limit = (roundnumber > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND;
      break;
   case pc_sdw2:
      limit = (roundnumber > 0) ? STP_RED_SDSPBOUND2 : 0;
      break;
   case pc_bd3:
      limit = (roundnumber > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND / 2;
      break;
   default:
      assert(0);

   }

   if( nedges >= STP_RED_EDGELIMIT && roundnumber == 0 )
      limit = (int) MAX(limit, limit * sqrtnedges / 5000.0);
   else
      limit = (int) MAX(limit, limit * sqrtnedges / 150.0);

 //  printf("limit %d \n", limit);
 //  printf("avg degree %f \n", graph_get_avgDeg(g));

   assert(limit >= 0);

   return limit;
}


/** print reduction information */
static
void reduceStatsPrint(
   SCIP_Bool             print,
   const char*           method,
   int                   nelims
)
{
   assert(nelims >= 0);

#ifdef STP_PRINT_STATS
    if( print )
       printf("%s: %d \n", method, nelims);
#endif

}


/** iterate NV and SL test while at least minelims many contractions are being performed */
static
SCIP_RETCODE nvreduce_sl(
   SCIP*                 scip,
   const int*            edgestate,          /**< for propagation or NULL */
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
      SCIP_CALL( reduce_nvAdv(scip, edgestate, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, neighb, distnode, solnode, &nvelims) );
      elims += nvelims;

      SCIPdebugMessage("NV-reduction (in NVSL): %d \n", nvelims);

      /* SL-reduction */
      SCIP_CALL( reduce_sl(scip, edgestate, g, vnoi, fixed, heap, state, vbase, neighb, visited, solnode, &slelims) );
      elims += slelims;

      SCIPdebugMessage("SL-reduction (in NVSL): %d \n", slelims);

      /* trivial reductions */
      if( elims > 0 )
      {
         if( g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG )
            SCIP_CALL( reduce_simple_pc(scip, edgestate, g, fixed, &degelims, NULL, solnode) );
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

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}

static
SCIP_RETCODE execPc_SD(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nodesid,            /**< array */
   int*                  nodesorg,           /**< array */
   int*                  nelims,             /**< pointer to store number of eliminated edges */
   int                   redbound,           /**< reduction bound */
   SCIP_Bool             verbose,            /**< be verbose? */
   SCIP_Bool*            rerun               /**< use again? */
)
{
   SCIP_CALL( reduce_sdPc(scip, g, vnoi, heap, state, vbase, nodesid, nodesorg, nelims) );

   if( verbose )
      printf("pc_SD eliminations: %d \n", *nelims);

   if( *nelims <= redbound )
      *rerun = FALSE;

   return SCIP_OKAY;
}

#if 0
static
SCIP_RETCODE execPc_SDSP(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit,
   int*                  edgestate,
   int                   redbound,          /**< reduction bound */
   SCIP_Bool             verbose,           /**< be verbose? */
   SCIP_Bool*            rerun              /**< use again? */
)
{
   SCIP_CALL( reduce_sdsp(scip, g, pathtail, pathhead, heap, statetail, statehead, memlbltail, memlblhead, nelims,
         limit, edgestate) );

   if( verbose )
      printf("pc_SDSP eliminations: %d \n", *nelims);

   if( *nelims <= redbound )
      *rerun = FALSE;

   return SCIP_OKAY;
}
#endif

static
SCIP_RETCODE execPc_BDk(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   PATH*                 pathtail,           /**< array for internal use */
   PATH*                 pathhead,           /**< array for internal use */
   int*                  heap,               /**< array for internal use */
   int*                  statetail,          /**< array for internal use */
   int*                  statehead,          /**< array for internal use */
   int*                  memlbltail,         /**< array for internal use */
   int*                  memlblhead,         /**< array for internal use */
   int*                  nelims,             /**< point to return number of eliminations */
   int                   limit,              /**< limit for edges to consider for each vertex */
   SCIP_Real*            offset,             /**< offset */
   int                   redbound,           /**< reduction bound */
   SCIP_Bool             verbose,            /**< be verbose? */
   SCIP_Bool*            rerun               /**< use again? */
)
{
   SCIP_CALL( reduce_bd34(scip, g, pathtail, pathhead, heap, statetail, statehead, memlbltail,
         memlblhead, nelims, limit, offset) );

   if( verbose )
      printf("pc_BDk eliminations: %d \n", *nelims);

   if( *nelims <= redbound )
      *rerun = FALSE;

   return SCIP_OKAY;
}

static
SCIP_RETCODE execPc_NVSL(
   SCIP*                 scip,
   const int*            edgestate,          /**< for propagation or NULL */
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
   int                   redbound,           /**< reduction bound */
   SCIP_Bool             verbose,            /**< be verbose? */
   SCIP_Bool*            rerun               /**< use again? */
)
{
   SCIP_CALL( nvreduce_sl(scip, edgestate, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, neighb,
         distnode, solnode, visited, nelims, redbound) );

   if( verbose )
      printf("pc_NVSL eliminations: %d \n", *nelims);

   if( *nelims <= redbound / 2 )
      *rerun = FALSE;

   return SCIP_OKAY;
}

static
SCIP_RETCODE execPc_BND(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            radius,             /**< radius array                       */
   SCIP_Real*            offset,             /**< pointer to the offset              */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nelims,             /**< pointer to store number of eliminated edges */
   int                   redbound,           /**< reduction bound */
   SCIP_Bool             verbose,            /**< be verbose? */
   SCIP_Bool*            rerun               /**< use again? */
)
{
   SCIP_Real ub = -1.0;

   SCIP_CALL( reduce_bound(scip, graph, vnoi, radius, offset, &ub, heap, state, vbase, nelims) );

   if( verbose )
      printf("pc_BND eliminations: %d \n", *nelims);

   if( *nelims <= redbound )
      *rerun = FALSE;

   return SCIP_OKAY;
}


/** inner MWCS loop */
static
SCIP_RETCODE redLoopMw_inner(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< voronoi base array  */
   int*                  nodearrint,         /**< nodes-sized array  */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   STP_Bool              dualascent,         /**< do dual-ascent reduction? */
   STP_Bool              bred,               /**< do bound-based reduction? */
   int                   redbound,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Real             prizesum            /**< prize sum */
   )
{
   SCIP_Real timelimit;
   STP_Bool da = dualascent;
   STP_Bool ans = TRUE;
   STP_Bool nnp = TRUE;
   STP_Bool npv = TRUE;
   STP_Bool rerun = TRUE;
   STP_Bool ansad = TRUE;
   STP_Bool ansad2 = TRUE;
   STP_Bool chain2 = TRUE;
   STP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_RANDNUMGEN* randnumgen;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   for( int rounds = 0; rounds < STP_RED_MAXNROUNDS && !SCIPisStopped(scip) && rerun; rounds++ )
   {
      int daelims = 0;
      int anselims = 0;
      int nnpelims = 0;
      int npvelims = 0;
      int bredelims = 0;
      int ansadelims = 0;
      int ansad2elims = 0;
      int chain2elims = 0;
      int degelims = 0;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( ans || extensive )
      {
         SCIP_CALL( reduce_ans(scip, g, &anselims) );
         if( anselims <= redbound )
            ans = FALSE;
      }

      if( ansad || extensive )
      {
         SCIP_CALL( reduce_ansAdv(scip, g, &ansadelims, FALSE) );
         if( ansadelims <= redbound )
            ansad = FALSE;

         SCIPdebugMessage("ans advanced deleted: %d \n", ansadelims);
      }

      if( ans || ansad || nnp || npv || extensive )
         SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

      if( (da || (dualascent && extensive)) )
      {
         const RPDA paramsda = { .prevrounds = 0, .useRec = userec, .useExtRed = FALSE, .nodereplacing = FALSE,
               .pcmw_solbasedda = TRUE, .pcmw_useMultRoots = FALSE, .pcmw_markroots = FALSE, .pcmw_fastDa = (rounds == 0) };

         if( graph_pc_isRootedPcMw(g) )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, NULL, fixed, vbase, state, nodearrint, nodearrchar, &daelims, randnumgen) );
         }
         else
         {
            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, vnoi, gnodearr, nodearrreal, vbase, nodearrint, state, nodearrchar, &daelims, randnumgen, prizesum) );
         }

         if( daelims <= 2 * redbound )
            da = FALSE;
         else
            SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

         SCIPdebugMessage("Dual-ascent eliminations: %d \n", daelims);
      }

      if( nnp )
      {
         SCIP_CALL( reduce_nnp(scip, g, &nnpelims) );
         if( nnpelims <= redbound )
            nnp = FALSE;

         SCIPdebugMessage("nnp eliminations: %d \n", nnpelims);
      }

      if( nnp || extensive )
      {
         SCIP_CALL(reduce_chain2(scip, g, vnoi, state, vbase, nodearrint, &chain2elims, 500));
         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 eliminations: %d \n", chain2elims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( npv || extensive )
      {
         SCIP_CALL(reduce_npv(scip, g, vnoi, state, vbase, nodearrint, &npvelims, 400));
         if( npvelims <= redbound )
            npv = FALSE;

         SCIPdebugMessage("npv eliminations: %d \n", npvelims);
      }

      if( chain2 || extensive )
      {
         SCIP_CALL(reduce_chain2(scip, g, vnoi, state, vbase, nodearrint, &chain2elims, 300));
         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 eliminations: %d \n", chain2elims);
      }

      SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

      if( ansad2 || extensive )
      {
         SCIP_CALL( reduce_ansAdv2(scip, g, &ansad2elims) );
         if( ansad2elims <= redbound )
            ansad2 = FALSE;
         else
            SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &ansad2elims) );

         SCIPdebugMessage("ans advanced 2 eliminations: %d (da? %d ) \n", ansad2elims, da);
      }

      if( bred )
      {
         SCIP_CALL( reduce_boundMw(scip, g, vnoi, fixed, nodearrint, state, vbase, NULL, &bredelims) );
         if( bredelims <= redbound )
            bred = FALSE;

         SCIPdebugMessage("reduce_bound eliminations: %d \n", bredelims);
      }

      if( anselims + nnpelims + chain2elims + bredelims + npvelims + ansadelims + ansad2elims + daelims <= redbound )
         rerun = FALSE;

#ifdef SCIP_DEBUG
      if(  advanced && tryrmw )
      {
         int nreal, mreal;
         graph_get_nVET(g, &nreal, &mreal, NULL);

         printf("round %d  nreal=%d mreal=%d\n", rounds, nreal, mreal);
      }
#endif
   }

   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/* remove unconnected vertices */
SCIP_RETCODE reduceLevel0(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nnodes = g->knots;
   SCIP_Bool* nodevisited;

   assert(scip && g);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );

   SCIP_CALL( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = nnodes - 1; k >= 0 ; k-- )
   {
      if( !nodevisited[k] && (g->grad[k] > 0) )
      {
         if( Is_term(g->term[k]) )
         {
            assert(graph_pc_isPc(g));
            assert(graph_pc_termIsNonLeafTerm(g, k));

            continue;
         }

         graph_knot_del(scip, g, k, TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}


/* remove unconnected vertices, including pseudo terminals, and checks for feasibility (adapts g->mark) */
SCIP_RETCODE reduceLevel0RpcRmwInfeas(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offsetp,            /**< pointer to offset */
   SCIP_Bool*            infeas              /**< is problem infeasible? */
)
{
   int* stackarr;
   SCIP_Bool* gmark;
   int stacksize;
   const int nnodes = g->knots;
   const SCIP_Bool isExtended = g->extended;

   assert(scip != NULL);
   assert(g != NULL);
   assert(offsetp != NULL);
   assert(graph_pc_isRootedPcMw(g));

   if( !isExtended )
      graph_pc_2trans(scip, g);

   SCIP_CALL( SCIPallocBufferArray(scip, &gmark, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stackarr, nnodes) );

   BMSclearMemoryArray(gmark, nnodes);

   *infeas = FALSE;
   stacksize = 0;
   stackarr[stacksize++] = g->source;
   assert(!gmark[g->source]);
   gmark[g->source] = TRUE;

   /* DFS loop */
   while( stacksize != 0 )
   {
      const int node = stackarr[--stacksize];

      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( !gmark[head] )
         {
            /* don't mark pseudo-terminals from the root */
            if( node == g->source && Is_term(g->term[head]) && !graph_pc_knotIsFixedTerm(g, head) )
            {
               assert(g->cost[flipedge(a)] == FARAWAY && g->grad[head] == 2);
               continue;
            }

            gmark[head] = TRUE;
            stackarr[stacksize++] = head;
         }
      }
   }

   SCIPfreeBufferArray(scip, &stackarr);

   for( int k = 0; k < nnodes; k++ )
   {
      if( !gmark[k] && Is_term(g->term[k]) )
      {
         assert(k != g->source);
         assert(graph_pc_knotIsFixedTerm(g, k) || (g->grad[k] > 0));

         if( graph_pc_knotIsFixedTerm(g, k) )
         {
            *infeas = TRUE;
            SCIPfreeBufferArray(scip, &gmark);

            if( !isExtended )
               graph_pc_2org(scip, g);

            return SCIP_OKAY;
         }
         else if( graph_pc_termIsNonLeafTerm(g, k) )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage("delete: \n");
            graph_knot_printInfo(g, k);
#endif

            /* offset is not necessary, because graph is in extended mode */
            graph_knot_del(scip, g, k, TRUE);
         }
         else
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage("delete: \n");
            graph_knot_printInfo(g, k);
#endif

            assert(g->term2edge[k] >= 0);
            assert(!gmark[graph_pc_getTwinTerm(g, k)]);

            (void) graph_pc_deleteTerm(scip, g, k, offsetp);
         }
      }
   }

   for( int k = 0; k < nnodes; k++ )
   {
      if( !gmark[k] && (g->grad[k] > 0) )
      {
         assert(!graph_pc_knotIsFixedTerm(g, k));
         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
         g->mark[k] = FALSE;
      }
   }

   SCIPfreeBufferArray(scip, &gmark);

   if( !isExtended )
      graph_pc_2org(scip, g);

   return SCIP_OKAY;
}

/* remove unconnected vertices, including pseudo terminals, adapts g->mark */
SCIP_RETCODE reduceLevel0RpcRmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offsetp             /**< pointer to offset */
)
{
   SCIP_Bool infeas;

   SCIP_CALL( reduceLevel0RpcRmwInfeas(scip, g, offsetp, &infeas) );

   if( infeas )
   {
      printf("level0RpcRmw detected infeasibility \n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/* remove unconnected vertices and checks whether problem is infeasible  */
SCIP_RETCODE reduceLevel0infeas(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool*            infeas              /**< is problem infeasible? */
)
{
   const int nnodes = g->knots;
   SCIP_Bool* nodevisited;

   assert(scip && g && infeas);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );

   *infeas = FALSE;

   SCIP_CALL( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
      if( !nodevisited[k] && Is_term(g->term[k]) )
      {
         assert(k != g->source);
         *infeas = TRUE;
         break;
      }

   for( int k = 0; k < nnodes; k++ )
   {
      if( !nodevisited[k] && (g->grad[k] > 0) )
      {
         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}

/** basic reduction package for the STP */
SCIP_RETCODE reduceStp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             dualascent,         /**< perform dual-ascent reductions? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   PATH* vnoi;
   PATH* path;
   SCIP_Real*  nodearrreal;
   SCIP_Real*  edgearrreal;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     nnodes;
   int     nedges;
   int     nterms;
   int     reductbound;
   STP_Bool* nodearrchar;
   SCIP_Bool bred = FALSE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, 0.03) )
      bred = TRUE;

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

   reductbound = MAX(nedges / 1000, minelims);

   /* reduction loop */
   SCIP_CALL( redLoopStp(scip, g, vnoi, path, nodearrreal, edgearrreal, heap, state,
         vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, fixed, -1.0, dualascent, bred, nodereplacing, reductbound, userec, (dualascent && userec)) );

   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n", *fixed);

   /* free memory */
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

   return SCIP_OKAY;
}

/** basic reduction package for the (R)PCSTP */
SCIP_RETCODE reducePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             advanced,           /**< perform advanced (e.g. dual ascent) reductions? */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             nodereplacing       /**< should node replacement (by edges) be performed? */
   )
{
   PATH* vnoi;
   PATH* path;
   GNODE** gnodearr;
   SCIP_Real* nodearrreal;
   SCIP_Real timelimit;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     nnodes;
   int     nterms;
   int     nedges;
   int     reductbound;
   STP_Bool* nodearrchar;
   SCIP_Bool bred = FALSE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;

   /* get timelimit parameter*/
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, STP_BND_THRESHOLD) )
   {
      bred = TRUE;

      if( SCIPisGT(scip, (double) (nterms + graph_pc_nNonLeafTerms(g)) / (double) nnodes, STP_BND_THRESHOLD ) )
         bred = FALSE;
   }

   if( advanced )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );

      for( int i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &gnodearr[i]) ); /*lint !e866*/
      }
   }
   else
   {
      gnodearr = NULL;
   }

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = MAX(nnodes / 1000, minelims);

   /* reduction loop */
   SCIP_CALL( redLoopPc(scip, edgestate, g, vnoi, path, gnodearr, nodearrreal, heap, state,
         vbase, nodearrint, edgearrint, nodearrint2, NULL, nodearrchar, fixed, advanced, bred, userec && advanced, reductbound, userec, nodereplacing) );

   /* free memory */

   if( gnodearr != NULL )
   {
      for( int i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }
   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);

   return SCIP_OKAY;
}

/** reduction package for the MWCSP */
SCIP_RETCODE reduceMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             advanced,           /**< perform advanced reductions? */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   PATH* vnoi;
   GNODE** gnodearr;
   SCIP_Real* nodearrreal;
   int* state;
   int* vbase;
   int* nodearrint;
   int i;
   int nterms;
   int nnodes;
   int redbound;
   STP_Bool* nodearrchar;
   STP_Bool bred = FALSE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   nnodes = g->knots;
   nterms = g->terms;
   redbound = MAX(nnodes / 500, minelims);

   if( SCIPisLE(scip, (double) nterms / (double) nnodes, 0.1) )
      bred = TRUE;

   if( advanced )
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

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );

   if( bred || advanced )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 2) );
   }
   else
   {
      nodearrreal = NULL;
   }

   /* reduction loop */
   SCIP_CALL( redLoopMw(scip, g, vnoi, gnodearr, nodearrreal, state,
         vbase, nodearrint, NULL, nodearrchar, fixed, advanced, bred, advanced, redbound, userec) );

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint);

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBlockMemory(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   return SCIP_OKAY;
}

/** basic reduction package for the HCDSTP */
SCIP_RETCODE reduceHc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
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
         SCIP_CALL( reduce_bound(scip, g, vnoi, radius, fixed, &upperbound, heap, state, vbase, &brednelims) );
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
SCIP_RETCODE reduceSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH*   vnoi;
   PATH*   path;
   SCIP_Real ub = FARAWAY;
   SCIP_Real timelimit;
   SCIP_Real*  nodearrreal;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    nodearrint2;
   int     e;
   int     nnodes;
   int     degtnelims;
   int     redbound;
   STP_Bool    da = TRUE;
   STP_Bool    sd = !TRUE;
   STP_Bool* nodearrchar;
   STP_Bool    rpt = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   nnodes = g->knots;

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   /* @todo change .stp file format for SAP! */
   for( e = 0; e < g->edges; e++ )
      if( SCIPisEQ(scip, g->cost[e], 20000.0) )
         g->cost[e] = FARAWAY;

   SCIP_CALL( reduce_simple_sap(scip, g, fixed, &degtnelims) );

   /* main loop */
   while( (sd || rpt || da) && !SCIPisStopped(scip) )
   {
      int danelims = 0;
      int sdnelims = 0;
      int rptnelims = 0;

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
         const RPDA paramsda = { .prevrounds = 0, .useRec = FALSE, .useExtRed = FALSE, .nodereplacing = FALSE};

         ub = -1.0;

         SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, &ub, fixed, vbase, state,
               heap, nodearrchar, &danelims, randnumgen) );

         if( danelims <= 2 * redbound )
            da = FALSE;
      }
   }

   SCIP_CALL( reduce_simple_sap(scip, g, fixed, &degtnelims) );

   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}


/** reduce node-weighted Steiner tree problem */
SCIP_RETCODE reduceNw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH *vnoi;
   SCIP_Real *nodearrreal;
   SCIP_Real timelimit;
   int *heap;
   int *state;
   int *vbase;
   int *nodearrint2;
   int nnodes;
   int redbound;

   STP_Bool*   nodearrchar;
   STP_Bool    da = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   nnodes = g->knots;

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   while( (da) && !SCIPisStopped(scip) )
   {
      int danelims = 0;
      const RPDA paramsda = { .prevrounds = 0, .useRec = FALSE, .useExtRed = FALSE, .nodereplacing = FALSE};

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, NULL, fixed, vbase, state, heap,
            nodearrchar, &danelims, randnumgen) );

      if( danelims <= 2 * redbound )
         da = FALSE;
   }

   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** MWCS loop */
SCIP_RETCODE redLoopMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< voronoi base array  */
   int*                  nodearrint,         /**< nodes-sized array  */
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
   SCIP_Real da = advanced;
   SCIP_Real timelimit;
   int degelims;
   SCIP_Bool fullrerun = TRUE;
   STP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real prizesum;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(advanced || !tryrmw);

   tryrmw = tryrmw && userec;
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   graph_pc_2org(scip, g);
   degelims = 0;
   SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
   assert(graph_pc_term2edgeIsConsistent(scip, g));

   prizesum = graph_pc_getPosPrizeSum(scip, g);

   /* todo currently run a most two times */
   for( int rounds = 0; rounds < 3 && !SCIPisStopped(scip) && fullrerun; rounds++ )
   {
      const RPDA paramsda = { .prevrounds = 0, .useRec = userec, .useExtRed = FALSE, .nodereplacing = FALSE,
            .pcmw_solbasedda = TRUE, .pcmw_useMultRoots =  (g->terms > STP_RED_MWTERMBOUND),
            .pcmw_markroots = tryrmw, .pcmw_fastDa = FALSE };
      int daelims = 0;

      fullrerun = FALSE;

      SCIP_CALL( redLoopMw_inner(scip, g, vnoi, gnodearr, nodearrreal, state,
            vbase, nodearrint, solnode, nodearrchar, fixed, da, bred, redbound, userec, prizesum) );

      if( advanced && g->terms > 2 )
      {
         SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

         if( graph_pc_isRootedPcMw(g) )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, NULL, fixed, vbase, state, nodearrint, nodearrchar, &daelims, randnumgen) );
         }
         else
         {
            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, vnoi, gnodearr, nodearrreal, vbase, nodearrint, state, nodearrchar, &daelims, randnumgen, prizesum) );
         }

         da = FALSE;
         userec = FALSE;
         advanced = FALSE;
         fullrerun = (daelims >= redbound || (extensive && (daelims > 0)));

         if( !graph_pc_isRootedPcMw(g) && tryrmw && g->terms > 2 )
         {
#ifdef WITH_UG
             SCIP_CALL( graph_transPcmw2rooted(scip, g, prizesum, FALSE) );
#else
             SCIP_CALL( graph_transPcmw2rooted(scip, g, prizesum, TRUE) );
#endif
             if( graph_pc_isRootedPcMw(g) )
             {
                if( fullrerun )
                   da = TRUE;
             }
         }

         if( fullrerun )
         {
#ifndef WITH_UG
            printf("Restarting reduction loop after %d rounds! (%d DA eliminations) \n\n ", rounds, daelims);
#endif
            if( extensive )
               advanced = TRUE;
         }
      }

      SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
   }

   /* go back to the extended graph */
   graph_pc_2trans(scip, g);

   SCIP_CALL( reduceLevel0(scip, g) );

   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** (R)PC loop */
SCIP_RETCODE redLoopPc(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   GNODE**               gnodearr,           /**< nodes-sized array  */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  heap,               /**< shortest path array  */
   int*                  state,              /**< voronoi base array  */
   int*                  vbase,              /**< nodes-sized array  */
   int*                  nodearrint,         /**< node-sized array  */
   int*                  edgearrint,         /**< edge-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   int*                  solnode,            /**< solution nodes array (or NULL) */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   SCIP_Bool             dualascent,         /**< do dual-ascent reduction? */
   SCIP_Bool             bred,               /**< do bound-based reduction? */
   SCIP_Bool             tryrpc,             /**< try to transform to rpc? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             nodereplacing       /**< should node replacement (by edges) be performed? */
   )
{
   DHEAP* dheap;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real fix = 0.0;
   SCIP_Real timelimit;
   SCIP_Bool rerun = TRUE;
   const SCIP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_Bool advancedrun = dualascent;
   SCIP_Real prizesum;
   const SCIP_Bool verbose = show && dualascent && userec && nodereplacing;
   int nelims;
   int degnelims;
   const int reductbound_global = reductbound * STP_RED_GLBFACTOR;

   /* already finished? */
   if( g->grad[g->source] == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( graph_heap_create(scip, g->knots, NULL, NULL, &dheap) );
   graph_pc_2org(scip, g);

   assert(graph_pc_term2edgeIsConsistent(scip, g));

   SCIP_CALL( graph_pc_presolInit(scip, g) );

   SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &degnelims, NULL, solnode) );
   if( verbose ) printf("initial degnelims: %d \n", degnelims);

   prizesum = graph_pc_getPosPrizeSum(scip, g);
   assert(prizesum < FARAWAY);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* main reduction loop */
   for( int outterrounds = 0; outterrounds < 3 && rerun; outterrounds++ )
   {
      int ninnerelims = 0;
      SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
      SCIP_Bool da = dualascent;
      SCIP_Bool sd = TRUE;
      SCIP_Bool sdc = TRUE;
      SCIP_Bool sdw = TRUE;
      SCIP_Bool sdstar = TRUE;
      SCIP_Bool bd3 = nodereplacing;
      SCIP_Bool nvsl = TRUE;

      /* inner loop todo outsource! */
      for( int rounds = 0; rounds < STP_RED_MAXNROUNDS && !SCIPisStopped(scip) && rerun; rounds++ )
      {
         int danelims = 0;
         int sdnelims = 0;
         int sdcnelims = 0;
         int bd3nelims = 0;
         int nvslnelims = 0;
         int sdwnelims = 0;
         int sdstarnelims = 0;
         int brednelims = 0;
         degnelims = 0;
         nelims = 0;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

         if( sdstar || extensive )
         {
            int sdstarpcnelims = 0;
            SCIP_CALL( reduce_sdStar(scip, getWorkLimits_pc(g, rounds, pc_sdstar), NULL, g, nodearrreal, nodearrint, nodearrint2, nodearrchar, dheap, &sdstarnelims));

            if( verbose ) printf("sdstarnelims %d \n", sdstarnelims);

            SCIP_CALL( reduce_sdStarPc2(scip, getWorkLimits_pc(g, rounds, pc_sdstar), NULL, g, nodearrreal, nodearrint, nodearrint2, nodearrchar, dheap, &sdstarpcnelims));

            if( verbose )  printf("sdstarpcnelims %d \n", sdstarpcnelims);

            sdstarnelims += sdstarpcnelims;

            if( sdstarnelims <= reductbound )
               sdstar = FALSE;
         }

         SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

         if( sd || extensive )
         {
            SCIP_CALL( execPc_SD(scip, g, vnoi, heap, state, vbase, nodearrint, nodearrint2, &sdnelims,
                  reductbound, verbose, &sd) );
         }

         if( sd )
         SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

         if( sdw || extensive )
         {
            int sdwnelims2 = 0;

            SCIP_CALL( reduce_sdWalkTriangle(scip, getWorkLimits_pc(g, rounds, pc_sdw1), NULL, g, nodearrint, nodearrreal, vbase, nodearrchar, dheap, &sdwnelims));
            SCIP_CALL( reduce_sdWalkExt(scip, getWorkLimits_pc(g, rounds, pc_sdw2), NULL, g, nodearrreal, heap, state, vbase, nodearrchar, &sdwnelims2) );

            if( verbose )
               printf("SDw: %d, SDwExt: %d\n", sdwnelims, sdwnelims2);

            sdwnelims += sdwnelims2;

            if( sdwnelims <= reductbound )
               sdw = FALSE;
         }

         if( sdc || extensive )
         {
#if 0
            SCIP_CALL( execPc_SDSP(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims,
               getWorkLimits_pc(g, rounds, pc_sdc), NULL, reductbound, verbose, &sdc) );
#endif
         }
         if( SCIPgetTotalTime(scip) > timelimit )
             break;

         if( sdw )
            SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

         if( bd3 && dualascent )
         {
            SCIP_CALL( execPc_BDk(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2,
                  &bd3nelims, getWorkLimits_pc(g, rounds, pc_bd3), &fix, reductbound, verbose, &bd3) );
         }

         if( nvsl || extensive )
         {
            SCIP_CALL( execPc_NVSL(scip, edgestate, g, vnoi, nodearrreal, &fix, edgearrint, heap, state, vbase,
                  nodearrint, nodearrint2, solnode, nodearrchar, &nvslnelims, reductbound, verbose, &nvsl) );
         }


         if( bred )
         {
            SCIP_CALL( execPc_BND(scip, g, vnoi, nodearrreal, &fix, heap, state, vbase, &brednelims, reductbound, verbose, &bred) );
         }

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

         if( da || (dualascent && extensive) )
         {
            const RPDA paramsda = { .prevrounds = 0, .useRec = userec, .useExtRed = FALSE, .nodereplacing = nodereplacing,
                  .pcmw_solbasedda = TRUE, .pcmw_useMultRoots = FALSE, .pcmw_markroots = FALSE, .pcmw_fastDa = (rounds == 0) };

            SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

            if( rpc )
            {
               SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, NULL, &fix, vbase, state, heap,
                     nodearrchar, &danelims, randnumgen) );
            }
            else
            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, vnoi, gnodearr, nodearrreal, vbase, heap,
                  state, nodearrchar, &danelims, randnumgen, prizesum) );

            if( danelims <= reductbound )
               da = FALSE;

            if( verbose ) printf("daX: %d \n", danelims);
         }

         SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

         nelims = degnelims + sdnelims + sdcnelims + bd3nelims + danelims + brednelims + nvslnelims + sdwnelims + sdstarnelims;
         ninnerelims += nelims;

         if( nelims <= reductbound_global )
            rerun = FALSE;
      } /* inner reduction loop */

      if( advancedrun && g->terms > 2 )
      {
         const RPDA paramsda = { .prevrounds = 0, .useRec = userec, .useExtRed = TRUE, .nodereplacing = nodereplacing,
                                 .pcmw_solbasedda = TRUE, .pcmw_useMultRoots = TRUE, .pcmw_markroots = TRUE, .pcmw_fastDa = FALSE };
         int danelims = 0;
         degnelims = 0;
         nelims = 0;
         advancedrun = FALSE;

         if( rpc )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, NULL, &fix, vbase, state, heap,
                  nodearrchar, &danelims, randnumgen) );
         }
         else
         {
            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, vnoi, gnodearr, nodearrreal, vbase, heap,
                  state, nodearrchar, &danelims, randnumgen, prizesum) );
         }

         SCIP_CALL( reduce_simple_pc(scip, edgestate, g, &fix, &nelims, &degnelims, solnode) );

         ninnerelims += danelims + degnelims;

         if( ninnerelims > reductbound_global )
         {
            if( danelims > reductbound_global )
               advancedrun = TRUE;

            rerun = TRUE;
         }
      }

      if( !rpc && tryrpc && g->terms > 2 )
      {
         assert(graph_pc_term2edgeIsConsistent(scip, g));

#ifdef WITH_UG
         SCIP_CALL(graph_transPcmw2rooted(scip, g, prizesum, FALSE));
#else
         SCIP_CALL(graph_transPcmw2rooted(scip, g, prizesum, TRUE));
#endif
         rpc = (g->stp_type == STP_RPCSPG);

         if( rpc )
         {
            SCIP_CALL(reduceLevel0RpcRmw(scip, g, &fix));
            rerun = TRUE;
            advancedrun = dualascent;
            outterrounds = 0;
         }
      }
   } // outter loop; todo refactor

   if( dualascent && tryrpc)
   {
      //reduce_simple_aritculations(scip, g, NULL, NULL);
      SCIP_CALL( reduce_deleteConflictEdges(scip, g) );
   }

   assert(graph_pc_term2edgeIsConsistent(scip, g));
   graph_pc_2trans(scip, g);
   graph_pc_presolExit(scip, g);

   graph_heap_free(scip, TRUE, TRUE, &dheap);
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
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   SCIP_Real*            edgearrreal,        /**< edges-sized array  */
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
   SCIP_Bool             dualascent,         /**< do dual-ascent reduction? */
   SCIP_Bool             boundreduce,        /**< do bound-based reduction? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             fullreduce          /**< use full reductions? (including extended techniques) */
   )
{
   SCIP_Real    ub;
   SCIP_Real    fix;
   SCIP_Real    timelimit;
   SCIP_Bool    le = TRUE;
   SCIP_Bool    sd = TRUE;
   SCIP_Bool    da = dualascent;
   SCIP_Bool    sdc = TRUE;
   SCIP_Bool    bd3 = nodereplacing;
   SCIP_Bool    bred = boundreduce;
   SCIP_Bool    nvsl = nodereplacing;
   SCIP_Bool    rerun = TRUE;

   const SCIP_Bool extensive = STP_RED_EXTENSIVE;
   int i = 0;

   SCIP_RANDNUMGEN* randnumgen;

   assert(reductbound > 0);
   assert(graph_valid(scip, g));

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   ub = upperbound;
   fix = 0.0;

   SCIP_CALL( reduce_contract0Edges(scip, g, TRUE) );
   SCIP_CALL( reduce_simple(scip, g, &fix, solnode, &i, NULL) );

   /* get timelimit parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   do
   {
      int inner_rounds = 0;
      int inner_restarts = 0;

      /* inner reduction loop */
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

         // for debugging of extended reductions todo deleteme
#if 0
         if( fullreduce && da)
         {
           const RPDA paramsda = { .prevrounds = inner_rounds, .useRec = userec, .useExtRed = TRUE, .nodereplacing = nodereplacing};
           int extendedelims = 0;
           SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, &ub, &fix, vbase, state, heap, nodearrchar, &extendedelims, randnumgen) );
           printf("debug extendedelims=%d \n", extendedelims);
         }
#endif

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

#if 1
#if 0
         SCIP_CALL( reduce_sdsp(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims,
                           ((inner_rounds > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND), NULL));



         SCIP_CALL( reduce_sd(scip, g, vnoi, edgearrreal, nodearrreal, heap, state, vbase, nodearrint, nodearrint2, edgearrint, &sdnelims,
                          nodereplacing, NULL));

         SCIP_CALL(reduce_ledge(scip, g, vnoi, heap, state, vbase, &lenelims, NULL));
         printf("lenelims=%d, sdnelims=%d, sdcnelims=%d \n", lenelims, sdnelims, sdcnelims);

#endif

         int nelimsx = 0;
         SCIP_CALL( reduce_sdStarBiased(scip, ((inner_rounds > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND)
               , NULL, g, &nelimsx));

#if 0
         FILE *fp;
                       fp = fopen("/nfs/optimi/kombadon/bzfrehfe/projects/scip/applications/STP/redsNew2.txt", "a+");
                       fprintf(fp, "%s %d \n", SCIPgetProbName(scip), nelimsx);
                       fclose(fp);
                       exit(1);
#endif
#endif


         if( le || extensive )
         {
            SCIP_CALL(reduce_ledge(scip, g, vnoi, heap, state, vbase, &lenelims, NULL));

            if( lenelims <= reductbound )
               le = FALSE;
            else
               SCIP_CALL(reduce_simple(scip, g, &fix, solnode, &degtnelims, NULL));

            reduceStatsPrint(fullreduce, "le", lenelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         if( sd || extensive )
         {
            SCIP_CALL( reduce_sd(scip, g, vnoi, edgearrreal, nodearrreal, heap, state, vbase, nodearrint, nodearrint2, edgearrint, &sdnelims,
                  nodereplacing, NULL));

            if( sdnelims <= reductbound )
               sd = FALSE;

            reduceStatsPrint(fullreduce, "sd", sdnelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         if( sdc || extensive )
         {
            SCIP_CALL( reduce_sdsp(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims,
                  ((inner_rounds > 0) ? STP_RED_SDSPBOUND2 : STP_RED_SDSPBOUND), NULL));

            if( sdcnelims <= reductbound )
               sdc = FALSE;

            reduceStatsPrint(fullreduce, "sdsp", sdcnelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         if( sd || sdc )
            SCIP_CALL(reduce_simple(scip, g, &fix, solnode, &degtnelims, NULL));

         if( bd3 || extensive )
         {
            SCIP_CALL(reduce_bd34(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &bd3nelims, STP_RED_BD3BOUND, &fix));
            if( bd3nelims <= reductbound )
               bd3 = FALSE;
            else
               SCIP_CALL(reduce_simple(scip, g, &fix, solnode, &degtnelims, NULL));

            reduceStatsPrint(fullreduce, "bd3", bd3nelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         if( nvsl || extensive )
         {
            SCIP_CALL( nvreduce_sl(scip, NULL, g, vnoi, nodearrreal, &fix, edgearrint, heap, state, vbase, nodearrint, NULL,
                  solnode, nodearrchar, &nvslnelims, reductbound));

            if( nvslnelims <= reductbound )
               nvsl = FALSE;

            reduceStatsPrint(fullreduce, "nvsl", nvslnelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         ub = -1.0;

         if( da )
         {
            const RPDA paramsda = { .prevrounds = inner_rounds, .useRec = userec, .useExtRed = FALSE, .nodereplacing = nodereplacing};
            SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, &ub, &fix, vbase,
                  state, heap, nodearrchar, &danelims, randnumgen));

            if( danelims <= STP_RED_EXFACTOR * reductbound )
               da = FALSE;

            reduceStatsPrint(fullreduce, "da", danelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         if( bred && nodereplacing )
         {
            SCIP_CALL(reduce_bound(scip, g, vnoi, nodearrreal, &fix, &ub, heap, state, vbase, &brednelims));

            SCIP_CALL(reduceLevel0(scip, g));

            if( brednelims <= reductbound )
               bred = FALSE;

            reduceStatsPrint(fullreduce, "bnd", brednelims);

            if( SCIPgetTotalTime(scip) > timelimit )
               break;
         }

         SCIP_CALL(reduceLevel0(scip, g));
         SCIP_CALL(reduce_simple(scip, g, &fix, solnode, &degtnelims, NULL));

         if( (danelims + sdnelims + bd3nelims + nvslnelims + lenelims + brednelims + sdcnelims) <= 2 * reductbound )
         {
            // at least one successful round and full reduce and no inner_restarts yet?
            if( inner_rounds > 0 && fullreduce && inner_restarts == 0 )
            {
               inner_restarts++;
               le = TRUE;
               sd = TRUE;
               sdc = TRUE;
               da = TRUE;
               bd3 = nodereplacing;
               nvsl = nodereplacing;

#ifdef STP_PRINT_STATS
               printf("RESTART reductions (restart %d) \n", inner_restarts);
#endif
            }
            else
               rerun = FALSE;
         }

         if( extensive && (danelims + sdnelims + bd3nelims + nvslnelims + lenelims + brednelims + sdcnelims) > 0 )
            rerun = TRUE;

         inner_rounds++;
      } /* inner reduction loop */

      if( fullreduce && !SCIPisStopped(scip) )
      {
         int extendedelims = 0;
         const RPDA paramsda = { .prevrounds = inner_rounds, .useRec = userec, .useExtRed = TRUE, .nodereplacing = nodereplacing};

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

         assert(!rerun);

         SCIP_CALL( reduce_da(scip, g, &paramsda, vnoi, nodearrreal, &ub, &fix, vbase, state, heap,
                     nodearrchar, &extendedelims, randnumgen) );

         reduceStatsPrint(fullreduce, "ext", extendedelims);

         SCIP_CALL(reduce_simple(scip, g, &fix, solnode, &extendedelims, NULL));

         if( nodereplacing )
         {
            int conflictnelims = 0;

            SCIP_CALL( reduce_fixedConflicts(scip, NULL, g, &conflictnelims) );

            reduceStatsPrint(fullreduce, "fixedconflict", conflictnelims);
            extendedelims += conflictnelims;
         }

         if( extendedelims > STP_RED_EXFACTOR * reductbound )
         {
            le = TRUE;
            sd = TRUE;
            sdc = TRUE;
            da = TRUE;
            bd3 = nodereplacing;
            nvsl = nodereplacing;
            rerun = TRUE;
         }
      }
   }
   while( rerun && !SCIPisStopped(scip) ); /* extensive reduction loop */

   assert(graph_valid_ancestors(scip, g));

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   *fixed += fix;

   return SCIP_OKAY;
}


/** reduces the graph */
SCIP_RETCODE reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph structure */
   SCIP_Real*            offset,             /**< pointer to store offset generated by reductions */
   int                   reductionlevel,     /**< reduction level 0: none, 1: basic, 2: advanced */
   int                   minelims,           /**< minimal amount of reductions to reiterate reduction methods */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   int stp_type;

   assert(graph != NULL);
   assert(graph_get_fixedges(graph) == NULL);
   assert(reductionlevel == STP_REDUCTION_NONE || reductionlevel == STP_REDUCTION_BASIC || reductionlevel == STP_REDUCTION_ADVANCED );
   assert(minelims >= 0);

   *offset = 0.0;
   show = FALSE;
   stp_type = graph->stp_type;

   /* initialize ancestor list for each edge */
   SCIP_CALL( graph_init_history(scip, graph) );

   /* initialize shortest path algorithms */
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL( reduceLevel0(scip, graph) );

   /* if no reduction methods available, return */
   if( graph->stp_type == STP_DCSTP || graph->stp_type == STP_RMWCSP || graph->stp_type == STP_NWPTSPG || graph->stp_type == STP_BRMWCSP )
   {
      graph_path_exit(scip, graph);
      return SCIP_OKAY;
   }

   if( reductionlevel == STP_REDUCTION_BASIC )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reducePc(scip, NULL, graph, offset, minelims, FALSE, FALSE, TRUE) );
      }
      else if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
      {
         SCIP_CALL( reduceMw(scip, graph, offset, minelims, FALSE, FALSE) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduceHc(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_SAP )
      {
         SCIP_CALL( reduceSap(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduceNw(scip, graph, offset, minelims) );
      }
      else
      {
         assert(graph_typeIsSpgLike(graph));

         SCIP_CALL( reduceStp(scip, graph, offset, minelims, FALSE, TRUE, FALSE) );
      }
   }
   else if( reductionlevel == STP_REDUCTION_ADVANCED )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reducePc(scip, NULL, graph, offset, minelims, TRUE, userec, TRUE) );
      }
      else if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
      {
         SCIP_CALL( reduceMw(scip, graph, offset, minelims, TRUE, userec) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduceHc(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_SAP )
      {
         SCIP_CALL( reduceSap(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduceNw(scip, graph, offset, minelims) );
      }
      else
      {
         assert(graph_typeIsSpgLike(graph));

         SCIP_CALL( reduceStp(scip, graph, offset, minelims, TRUE, TRUE, userec) );
      }
   }
   SCIPdebugMessage("offset : %f \n", *offset);

   SCIP_CALL( reduceLevel0(scip, graph) );
   show = FALSE;

   assert(graph_valid(scip, graph));

   graph_path_exit(scip, graph);

   return SCIP_OKAY;
}
