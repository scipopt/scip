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
#define STP_RED_EXTENSIVE FALSE          /**< just for testing */
#define STP_REDBOUND_SDSP    200         /**< visited edges bound for SDSP test  */
#define STP_REDBOUND_SDSP2   800         /**< visited edges bound for SDSP test  */
#define STP_REDBOUND_PCBDK   80          /**< visited edges bound for PC BDK test  */
#define STP_REDBOUND_PCBDK2  400         /**< visited edges bound for PC BDK test  */

#define STP_REDBOUND_BDK     400         /**< visited edges bound for BDK test  */
#define STP_REDBOUND_BDK2    600         /**< visited edges bound for BDK test  */
#define STP_REDBOUND_SDSTAR     400         /**< visited edges bound for SD Star test  */
#define STP_REDBOUND_SDSTAR2    800         /**< visited edges bound for SD Star test  */
#define STP_RED_MWTERMBOUND 400
#define STP_RED_EXPENSIVEFACTOR   2
#define STP_RED_GLBFACTOR         3
#define STP_RED_EDGELIMIT       100000
#define STP_RED_EDGELIMIT_HUGE 1000000
#define STP_RED_MAXNINNROUNDS    15
#define STP_RED_MAXNOUTROUNDS_PC  4
#define STP_RED_MAXNOUTROUNDS_MW  4

#define STP_BND_THRESHOLD 0.03
#define USE_FULLSEPA

#include "mincut.h"
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
#include "portab.h"
#include "dpterms.h"
#include "relax_stpenum.h"



enum PC_REDTYPE {pc_sdc, pc_sdstar, pc_sdstar2, pc_sdw1, pc_sdw2, pc_bd3};
enum STP_REDTYPE {stp_bdk, stp_sdstar, stp_sdstarbot};


/** returns limit parameter for SPG method */
static
int getWorkLimitsStp(
   const GRAPH* g,
   int roundnumber,
   SCIP_Bool fullreduce,
   enum STP_REDTYPE redtype
)
{
   int limit = -1;

   assert(roundnumber >= 0);
   switch (redtype)
   {
      case stp_bdk:
         limit = (roundnumber == 0) ? STP_REDBOUND_BDK : STP_REDBOUND_BDK2;
         break;
      case stp_sdstarbot:
         limit = (roundnumber == 0) ? STP_REDBOUND_SDSTAR : STP_REDBOUND_SDSTAR2;
         break;
      case stp_sdstar:
         limit = STP_REDBOUND_SDSTAR;
         break;
      default:
         assert(0);
   }

   if( fullreduce )
   {
      if( stp_bdk == redtype )
         limit /= 2;
   }
   else
   {
      limit /= 3;
   }

   assert(limit >= 0);

   return limit;
}


/** returns limit parameters for PCSTP method */
static
int getWorkLimitsPc(
   const GRAPH* g,
   int roundnumber,
   SCIP_Bool beFast,
   enum PC_REDTYPE redtype
)
{
   const int nedges = g->edges;
   const int sqrtnedges = (int) sqrt(nedges);
   int limit = 0;

   assert(roundnumber >= 0);
   assert(nedges >= 0);
   assert(sqrtnedges >= 0 && sqrtnedges <= nedges);

   // todo try mabe 100 500?
   // or somethin similar

   switch (redtype)
   {
      case pc_sdc:
         limit = (roundnumber > 0) ? STP_REDBOUND_SDSP2 : STP_REDBOUND_SDSP;
         break;
      case pc_sdstar:
         limit = (roundnumber > 0) ? STP_REDBOUND_SDSP2 : STP_REDBOUND_SDSP;
         break;
      case pc_sdstar2:
         limit = (roundnumber > 0) ? STP_REDBOUND_SDSP :  STP_REDBOUND_SDSP / 2;
         break;
      case pc_sdw1:
         limit = (roundnumber > 0) ? STP_REDBOUND_SDSP2 : STP_REDBOUND_SDSP / 2;
         break;
      case pc_sdw2:
         limit = (roundnumber > 0) ? STP_REDBOUND_SDSP : 0;
         break;
      case pc_bd3:
         limit = (roundnumber > 0) ? STP_REDBOUND_PCBDK2 : STP_REDBOUND_PCBDK;
         break;
      default:
         assert(0);
   }

   // todo more or less random number, tune them at least a bit
   if( beFast )
   {
      limit /= 2;
   }
   else
   {
      limit = (int) MAX(limit, limit * sqrtnedges / 200.0);
   }

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
SCIP_RETCODE execNvSl(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph */
   PATH*                 vnoi,               /**< Voronoi */
   SCIP_Real*            nodearrreal,        /**< array */
   SCIP_Real*            fixed,              /**< offset */
   int*                  edgearrint,         /**< array */
   int*                  vbase,              /**< array */
   int*                  neighb,             /**< array */
   int*                  distnode,           /**< array */
   int*                  solnode,            /**< array */
   STP_Bool*             visited,            /**< array */
   int*                  nelims,             /**< number of eliminations */
   int                   minelims            /**< minimum number of eliminations */
   )
{
   int elims;
   int nvelims;
   int slelims;
   int degelims;
   int totalelims;

   assert(g != NULL);
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
      SCIP_CALL( reduce_nvAdv(scip, edgestate, g, vnoi, nodearrreal, fixed, edgearrint, vbase, distnode, solnode, &nvelims) );
      elims += nvelims;

      SCIPdebugMessage("NV-reduction (in NVSL): %d \n", nvelims);

      /* SL-reduction */
      SCIP_CALL( reduce_sl(scip, edgestate, g, vnoi, fixed, vbase, neighb, visited, solnode, &slelims) );
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
   // todo test properly
#ifdef SCIP_DISABLED
   SCIP_CALL( reduce_sdPc(scip, g, vnoi, heap, state, vbase, nodesid, nodesorg, nelims) );

   if( verbose )
      printf("pc_SD eliminations: %d \n", *nelims);

   if( *nelims <= redbound )
      *rerun = FALSE;
#else
   *rerun = FALSE;
#endif

   return SCIP_OKAY;
}


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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             usestrongreds,      /**< allow strong reductions? */
   GRAPH*                g,                  /**< graph */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            nodearrreal,        /**< array */
   SCIP_Real*            fixed,              /**< offset */
   int*                  edgearrint,         /**< array */
   int*                  vbase,              /**< array */
   int*                  neighb,             /**< array */
   int*                  distnode,           /**< array */
   int*                  solnode,            /**< array */
   STP_Bool*             visited,            /**< array */
   int*                  nelims,             /**< number of eliminations */
   int                   redbound,           /**< reduction bound */
   SCIP_Bool             verbose,            /**< be verbose? */
   SCIP_Bool*            rerun               /**< use again? */
)
{
   // todo propagate usestrongreds and use it!
   SCIP_CALL( execNvSl(scip, NULL, g, vnoi, nodearrreal, fixed, edgearrint, vbase, neighb,
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



/** returns pointer */
static
SCIP_Real* redbaseGetOffsetPointer(
   REDBASE*              redbase             /**< base */
   )
{
   assert(redbase);
   assert(redbase->redsol);

   return reduce_solGetOffsetPointer(redbase->redsol);
}


/** returns pointer
 *  todo: remove once redsol is properly tested */
static
int* redbaseGetSolnode(
   REDSOLLOCAL*          redsollocal,        /**< data structure for retaining primal solution */
   REDBASE*              redbase             /**< base */
   )
{
   REDSOL* redsol;
   assert(redbase && redsollocal);

   redsol = redbase->redsol;
   assert(redsol);

   if( reduce_solUsesNodesol(redsol) )
   {
      assert(!redbase->solnode);
      assert(reduce_sollocalUsesNodesol(redsollocal));

      return reduce_sollocalGetSolnode(redsollocal);
   }

   return redbase->solnode;
}

/** inner STP reduction loop */
static
SCIP_RETCODE redLoopInnerStp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< generator */
   GRAPH*                g,                  /**< graph data structure */
   REDSOLLOCAL*          redsollocal,        /**< data structure for retaining primal solution */
   REDBASE*              redbase,            /**< parameters */
   SCIP_Bool*            wasDecomposed       /**< pointer to mark whether to exit early */
)
{
   const RPARAMS* const redparameters = redbase->redparameters;
   BIDECPARAMS* const bidecompparams = redbase->bidecompparams;
   SCIP_Real timelimit;
   SCIP_Bool rerun = TRUE;
   int inner_rounds = 0;
   int inner_restarts = 0;
   const int reductbound = redparameters->reductbound;
   const SCIP_Bool fullreduce = redparameters->fullreduce;
   const SCIP_Bool nodereplacing = redparameters->nodereplacing;
   const SCIP_Bool extensive = STP_RED_EXTENSIVE;
   const SCIP_Bool usestrongreds = redparameters->usestrongreds;
   SCIP_Bool sd = TRUE;
   SCIP_Bool sdbiased = redparameters->usestrongreds;
   SCIP_Bool sdc = FALSE;
   SCIP_Bool pathrep = usestrongreds && redparameters->dualascent;
   SCIP_Bool sdstar = TRUE;
   SCIP_Bool da = redparameters->dualascent;
   SCIP_Bool bdk = redparameters->nodereplacing;
   SCIP_Bool bred = redparameters->boundreduce;
   SCIP_Bool nvsl = redparameters->nodereplacing;

   *wasDecomposed = FALSE;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* inner reduction loop */
   while( rerun && !SCIPisStopped(scip) )
   {
      SCIP_Bool skiptests = FALSE;
      int danelims = 0;
      int sdnelims = 0;
      int sdcnelims = 0;
      int sdstarnelims = 0;
      int bdknelims = 0;
      int nvslnelims = 0;
      int brednelims = 0;
      int degtnelims = 0;
      int sdbiasnelims = 0;
      int pathrepnelims = 0;

      if( bidecompparams && bidecompparams->newLevelStarted )
      {
         skiptests = TRUE;
         bidecompparams->newLevelStarted = FALSE;
      }

      assert(*wasDecomposed == FALSE);

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( sd && !skiptests )
      {
         SCIP_CALL( reduce_sd(scip, g, redbase, &sdnelims) );

         if( sdnelims <= reductbound && !extensive )
            sd = FALSE;

         reduceStatsPrint(fullreduce, "sd", sdnelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd && !skiptests )
         SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &degtnelims, NULL));

      if( sdstar && !skiptests )
      {
         SCIP_CALL( reduce_sdStarBiased(scip, getWorkLimitsStp(g, inner_rounds, fullreduce, stp_sdstar),
               usestrongreds, g, &sdstarnelims) );

         if( sdstarnelims <= reductbound && !extensive )
            sdstar = FALSE;

         reduceStatsPrint(fullreduce, "sdstar", sdstarnelims);
      }

      if( sdc && !skiptests )
      {
         SCIP_CALL( reduce_sdsp(scip, g, redbase->vnoi,
               redbase->heap, redbase->state, redbase->vbase, redbase->nodearrint,
               redbase->nodearrint2, &sdcnelims,
               ((inner_rounds > 0) ? STP_REDBOUND_SDSP2 : STP_REDBOUND_SDSP), usestrongreds));

         if( sdcnelims <= reductbound && !extensive )
            sdc = FALSE;

         reduceStatsPrint(fullreduce, "sdsp", sdcnelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( (sdc || sdstar) && !skiptests )
         SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &degtnelims, NULL));

      if( pathrep )
      {
         SCIP_CALL( reduce_pathreplace(scip, g, &pathrepnelims) );
         if( pathrepnelims <= reductbound )
            pathrep = FALSE;
//printf("pathrep=%d \n", pathrepnelims);

         reduceStatsPrint(fullreduce, "pathrep", pathrepnelims);
      }

      if( bdk && !skiptests )
      {
         SCIP_CALL( reduce_bdk(scip, getWorkLimitsStp(g, inner_rounds, fullreduce, stp_bdk), g, &bdknelims) );

         if( bdknelims <= STP_RED_EXPENSIVEFACTOR * reductbound && !extensive )
            bdk = FALSE;
         else
            SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &degtnelims, NULL));

         reduceStatsPrint(fullreduce, "bdk", bdknelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sdbiased && !skiptests )
      {
         SCIP_CALL( reduce_impliedProfitBased(scip, getWorkLimitsStp(g, inner_rounds, fullreduce, stp_sdstarbot), g,
               redbaseGetSolnode(redsollocal, redbase), redbaseGetOffsetPointer(redbase), &sdbiasnelims) );

         if( sdbiasnelims <= reductbound && !extensive  )
            sdbiased = FALSE;

         reduceStatsPrint(fullreduce, "sdbiasnelims", sdbiasnelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( nvsl )
      {
         SCIP_CALL( execNvSl(scip, NULL, g, redbase->vnoi,
               redbase->nodearrreal, redbaseGetOffsetPointer(redbase), redbase->edgearrint,
               redbase->vbase, redbase->nodearrint, NULL,
               redbaseGetSolnode(redsollocal, redbase),
               redbase->nodearrchar, &nvslnelims, reductbound));

         if( nvslnelims <= reductbound && !extensive )
            nvsl = FALSE;

         reduceStatsPrint(fullreduce, "nvsl", nvslnelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bidecompparams && bidecompparams->depth < bidecompparams->maxdepth )
      {
         SCIPdebugMessage("go with depth %d \n", bidecompparams->depth);
         SCIP_CALL( reduce_bidecomposition(scip, g, redbase, redbaseGetSolnode(redsollocal, redbase), wasDecomposed) );
         SCIPdebugMessage("wasDecomposed=%d \n", *wasDecomposed);

         if( *wasDecomposed )
            break;
      }

      if( (inner_rounds > 0) && bred && nodereplacing )
      {
         SCIP_Real ub;
         reduce_sollocalSetOffset(*redbaseGetOffsetPointer(redbase), redsollocal);
         ub = reduce_sollocalGetUpperBound(redsollocal);
         SCIP_CALL(reduce_bound(scip, g,
               redbase->vnoi, redbase->nodearrreal, redbaseGetOffsetPointer(redbase), &ub,
               redbase->heap, redbase->state, redbase->vbase, &brednelims));

         if( brednelims <= reductbound )
            bred = FALSE;

         reduceStatsPrint(fullreduce, "bnd", brednelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( da )
      {
         const RPDA paramsda = {
           .damode = (!redparameters->userec && inner_rounds > 0) ? STP_DAMODE_MEDIUM : STP_DAMODE_FAST,
           .useRec = FALSE,
           .useSlackPrune = FALSE, .extredMode = extred_none, .nodereplacing = nodereplacing };

         SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal,
             redbaseGetOffsetPointer(redbase), &danelims, randnumgen));

         if( danelims <= STP_RED_EXPENSIVEFACTOR * reductbound )
            da = FALSE;

         reduceStatsPrint(fullreduce, "da", danelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      SCIP_CALL(reduce_unconnected(scip, g));
      SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &degtnelims, NULL));

      /* too few eliminations? */
      if( (sdbiasnelims + danelims + sdnelims + bdknelims + nvslnelims + pathrepnelims +
           brednelims + sdcnelims + sdstarnelims) <= STP_RED_GLBFACTOR * reductbound )
      {
         // at least one successful round and full reduce and no inner_restarts yet?
         if( inner_rounds > 0 && fullreduce && inner_restarts == 0 )
         {
            inner_restarts++;
            sd = TRUE;
            sdstar = TRUE;
            sdbiased = TRUE;
            sdc = TRUE;
            nvsl = nodereplacing;

            assert(extensive || sdcnelims == 0);

#ifdef STP_PRINT_STATS
            printf("RESTART inner reductions (restart number %d) \n", inner_restarts);
#endif
         }
         else
         {
            rerun = FALSE;
         }
      }

      if( extensive && (sdbiasnelims + danelims + sdnelims + bdknelims + nvslnelims + brednelims + sdcnelims + sdstarnelims) > 0 )
         rerun = TRUE;

      inner_rounds++;
   } /* inner reduction loop */

   return SCIP_OKAY;
}

/** inner MWCS loop
 *  todo use REDBASE here instead of all the parameters and also allocate the memory locally! */
static
SCIP_RETCODE redLoopInnerMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDSOLLOCAL*          redsollocal,        /**< data structure for retaining primal solution */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< voronoi base array  */
   int*                  nodearrint,         /**< nodes-sized array  */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   STP_Bool              dualascent,         /**< do dual-ascent reduction? */
   STP_Bool              bred,               /**< do bound-based reduction? */
   int                   redbound,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Real             prizesum            /**< prize sum */
   )
{
   int* solnode = reduce_sollocalGetSolnode(redsollocal);
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

   /* inner loop */
   for( int rounds = 0; rounds < STP_RED_MAXNINNROUNDS && !SCIPisStopped(scip) && rerun; rounds++ )
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
         const RPDA paramsda = {
            .damode = STP_DAMODE_FAST , .useRec = userec, .extredMode = extred_none,
            .nodereplacing = FALSE, .useSlackPrune = FALSE,
            .pcmw_solbasedda = (rounds > 0), .pcmw_useMultRoots = FALSE, .pcmw_markroots = FALSE, .pcmw_fastDa = (rounds == 0) };

         if( graph_pc_isRootedPcMw(g) )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &daelims, randnumgen) );
         }
         else
         {
            reduce_sollocalSetOffset(*fixed, redsollocal);
            SCIP_CALL( reduce_sollocalRebuildTry(scip, g, redsollocal) );

            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, redsollocal,
                  vnoi, nodearrreal, vbase, nodearrint, state, nodearrchar, &daelims, randnumgen, prizesum) );
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

         printf("round %d (of %d) nreal=%d mreal=%d\n", rounds, STP_RED_MAXNINNROUNDS, nreal, mreal);
      }
#endif
   } /* inner loop */

   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}


/** inner (R)PC loop
 * todo use REDBASE here instead of all the parameters and also allocate the memory locally! */
static
SCIP_RETCODE redLoopInnerPc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDSOLLOCAL*          redsollocal,        /**< solution store */
   DHEAP*                dheap,              /**< heap data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  heap,               /**< shortest path array  */
   int*                  state,              /**< voronoi base array  */
   int*                  vbase,              /**< nodes-sized array  */
   int*                  nodearrint,         /**< node-sized array  */
   int*                  edgearrint,         /**< edge-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Real*            fixed,              /**< offset */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             prizesum,           /**< prize sum */
   SCIP_Bool             dualascent,         /**< do dual-ascent reduction? */
   SCIP_Bool             bred,               /**< do bound-based reduction? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   SCIP_Bool             usestrongreds,      /**< allow strong reductions? */
   SCIP_Bool             beFast,             /**< fast mode? */
   int*                  ninnerelims         /**< number of eliminations */
   )
{
   int* solnode = reduce_sollocalGetSolnode(redsollocal);
   SCIP_Real timelimit;
   const int reductbound_global = reductbound * STP_RED_GLBFACTOR;
   SCIP_Bool dapaths = usestrongreds;
   SCIP_Bool da = dualascent;
   SCIP_Bool sd = usestrongreds;
   SCIP_Bool sdw = usestrongreds;
   SCIP_Bool sdstar = TRUE;
   SCIP_Bool bd3 = nodereplacing;
   SCIP_Bool nvsl = TRUE;
   SCIP_Bool pathrep = usestrongreds;
   SCIP_Bool rerun = TRUE;
   const SCIP_Bool verbose = FALSE && dualascent && userec && nodereplacing;
   const STP_Bool extensive = STP_RED_EXTENSIVE;

   *ninnerelims = 0;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* inner loop! */
   for( int rounds = 0; rounds < STP_RED_MAXNINNROUNDS && !SCIPisStopped(scip) && rerun; rounds++ )
   {
      int dapathelims = 0;
      int danelims = 0;
      int sdnelims = 0;
      int sdcnelims = 0;
      int bd3nelims = 0;
      int nvslnelims = 0;
      int sdwnelims = 0;
      int sdstarnelims = 0;
      int brednelims = 0;
      int degnelims = 0;
      int nelims = 0;
      int pathrepnelims = 0;


      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( sdstar || extensive )
      {
         int sdstarpcnelims = 0;

         SCIP_CALL( reduce_sdStarBiased(scip, getWorkLimitsPc(g, rounds, beFast, pc_sdstar), usestrongreds, g, &sdstarnelims));
         if( verbose ) printf("sdstarnelims %d \n", sdstarnelims);

         if( sdstarnelims <= reductbound )
         {
            sdstar = FALSE;
         }
         else
         {
            SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &nelims, &degnelims, solnode) );
            SCIP_CALL( reduce_sdStarPc2(scip, getWorkLimitsPc(g, rounds, beFast, pc_sdstar2), usestrongreds, g, nodearrreal, nodearrint, nodearrint2, nodearrchar, dheap, &sdstarpcnelims));
            if( verbose ) printf("sdstarpcnelims %d \n", sdstarpcnelims);

            sdstarnelims += sdstarpcnelims;

         }
      }

      if( (graph_pc_isRootedPcMw(g) && (dapaths || extensive)) )
      {
         reduce_dapaths(scip, g, fixed, &dapathelims);

         if( dapathelims <= reductbound )
            dapaths = FALSE;
      }

      if( sd || extensive )
      {
         SCIP_CALL( execPc_SD(scip, g, vnoi, heap, state, vbase, nodearrint, nodearrint2, &sdnelims,
               reductbound, verbose, &sd) );
      }

      SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &nelims, &degnelims, solnode) );

      if( sdw || extensive )
      {
         SCIP_CALL( reduce_sdWalkTriangle(scip, getWorkLimitsPc(g, rounds, beFast, pc_sdw1), usestrongreds, g, nodearrint, nodearrreal, vbase, nodearrchar, dheap, &sdwnelims));

         if( verbose ) printf("SDw:  %d\n", sdwnelims);

         if( sdwnelims <= reductbound )
         {
            sdw = FALSE;
         }
         else
         {
            int sdwnelims2 = 0;
            SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &nelims, &degnelims, solnode) );
            SCIP_CALL( reduce_sdWalkExt(scip, getWorkLimitsPc(g, rounds, beFast, pc_sdw2), usestrongreds, g, nodearrreal, heap, state, vbase, nodearrchar, &sdwnelims2) );
            sdwnelims += sdwnelims2;

            if( verbose ) printf("SDwExt: %d\n", sdwnelims2);
         }
      }

      if( SCIPgetTotalTime(scip) > timelimit )
          break;

      if( bd3 || extensive )
      {
         SCIP_CALL( execPc_BDk(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2,
               &bd3nelims, getWorkLimitsPc(g, rounds, beFast, pc_bd3), fixed, reductbound, verbose, &bd3) );
      }

      if( nvsl || extensive )
      {
         SCIP_CALL( execPc_NVSL(scip, usestrongreds, g, vnoi, nodearrreal, fixed, edgearrint, vbase,
               nodearrint, nodearrint2, solnode, nodearrchar, &nvslnelims, reductbound, verbose, &nvsl) );
      }

      if( bred || extensive)
      {
         SCIP_CALL( execPc_BND(scip, g, vnoi, nodearrreal, fixed, heap, state, vbase, &brednelims, reductbound, verbose, &bred) );
      }

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( da || (dualascent && extensive) )
      {
         RPDA paramsda = { .damode = STP_DAMODE_FAST, .useRec = userec, .extredMode = extred_none, .nodereplacing = nodereplacing, .useSlackPrune = FALSE,
             .pcmw_solbasedda = !beFast,
             .pcmw_useMultRoots = FALSE, .pcmw_markroots = FALSE,
             .pcmw_fastDa = beFast };

         if( g->edges > STP_RED_EDGELIMIT_HUGE )
         {
            paramsda.pcmw_solbasedda = TRUE;
            paramsda.pcmw_fastDa = (beFast && rounds == 0);
         }

         SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &nelims, &degnelims, solnode) );

         if( g->stp_type == STP_RPCSPG )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &danelims, randnumgen) );
         }
         else
         {
            reduce_sollocalSetOffset(*fixed, redsollocal);
            SCIP_CALL( reduce_sollocalRebuildTry(scip, g, redsollocal) );

            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, redsollocal, vnoi,  nodearrreal, vbase, heap,
               state, nodearrchar, &danelims, randnumgen, prizesum) );
         }

         if( danelims <= reductbound )
            da = FALSE;

         if( verbose ) printf("daX: %d \n", danelims);
      }

      SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &nelims, &degnelims, solnode) );

      if( pathrep )
      {
         SCIP_CALL( reduce_pathreplace(scip, g, &pathrepnelims) );

         if( pathrepnelims <= reductbound )
            pathrep = FALSE;
         if( verbose ) printf("pathrep: %d\n", pathrepnelims);
      }

      nelims = degnelims + sdnelims + sdcnelims + bd3nelims + danelims + brednelims + nvslnelims
            + sdwnelims + sdstarnelims + dapathelims + pathrepnelims;
      *ninnerelims += nelims;

      if( nelims <= reductbound_global )
         rerun = FALSE;
   } /* inner reduction loop */

   return SCIP_OKAY;
}




/*
 * Interface methods
 */



/** gets reduction bound */
int reduce_getMinNreductions(
   const GRAPH*          g,                  /**< graph data structure */
   int                   lowerbound          /**< lower bound on number of reductions (>= 1) */
   )
{
   int min;
   const int nedges = graph_get_nEdges(g);
   const int nnodes = graph_get_nNodes(g);

   assert(lowerbound >= 1);

   if( graph_typeIsSpgLike(g) )
   {
      min = MAX(nedges / 1000, lowerbound);
   }
   else if( graph_pc_isPc(g) )
   {
      min = MAX(nnodes / 500, lowerbound);
   }
   else if( graph_pc_isMw(g) )
   {
      min = MAX(nnodes / 500, lowerbound);
   }
   else
   {
      min = MAX(nnodes / 1000, lowerbound);
   }

   return min;
}


/** initializes */
SCIP_RETCODE reduce_baseInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDBASE**             redbase             /**< base */
   )
{
   REDBASE* red;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, redbase) );
   red = *redbase;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->edgearrint), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->nodearrchar), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->heap), nnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->state), 4 * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->nodearrreal), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->vbase), 4 * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->nodearrint), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->nodearrint2), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(red->vnoi), 4 * nnodes) );
   red->path = NULL;

   red->redparameters = NULL;
   red->solnode = NULL;
   red->redsol = NULL;

   return SCIP_OKAY;
}


/** frees */
void reduce_baseFree(
   SCIP*                 scip,               /**< SCIP data structure */
   REDBASE**             redbase             /**< base */
   )
{
   REDBASE* red;
   assert(scip && redbase);
   assert(*redbase);

   red = *redbase;

   SCIPfreeMemoryArray(scip, &(red->vnoi));
   SCIPfreeMemoryArray(scip, &(red->nodearrint2));
   SCIPfreeMemoryArray(scip, &(red->nodearrint));
   SCIPfreeMemoryArray(scip, &(red->vbase));
   SCIPfreeMemoryArray(scip, &(red->nodearrreal));
   SCIPfreeMemoryArray(scip, &(red->state));
   SCIPfreeMemoryArray(scip, &(red->heap));
   SCIPfreeMemoryArray(scip, &(red->nodearrchar));
   SCIPfreeMemoryArray(scip, &(red->edgearrint));

   SCIPfreeMemory(scip, redbase);
}


/** basic reduction package for the STP */
SCIP_RETCODE reduce_stp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDSOL*               redsol,             /**< primal solution container */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             dualascent,         /**< perform dual-ascent reductions? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             usestrongreds       /**< allow strong reductions?
                                                  NOTE: needed for propagation, because arcs might have been fixed to 0 */
   )
{
   PATH* vnoi;
   PATH* path;
   SCIP_Real* nodearrreal;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   STP_Bool* nodearrchar;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms = graph_get_nTerms(g);
   const int reductbound = reduce_getMinNreductions(g, minelims);
   SCIP_Bool bred = FALSE;

   assert(scip && redsol);

   if( SCIPisLE(scip, (SCIP_Real) nterms / (SCIP_Real) nnodes, STP_BND_THRESHOLD) )
      bred = TRUE;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   path = NULL;

   /* reduction loop */
   {
      RPARAMS parameters = { .dualascent = dualascent, .boundreduce = bred, .nodereplacing = nodereplacing,
                             .reductbound_min = minelims,
                             .reductbound = reductbound, .userec = userec, .fullreduce = (dualascent && userec),
                             .usestrongreds = usestrongreds };
      BIDECPARAMS decparameters = { .depth = 0, .maxdepth = 3, .newLevelStarted = FALSE };
      REDBASE redbase = { .redparameters = &parameters, .bidecompparams = &decparameters,
                          .solnode = NULL, .redsol = redsol,
                          .vnoi = vnoi, .path = path, .heap = heap,
                          .nodearrreal = nodearrreal,
                          .state = state, .vbase = vbase, .nodearrint = nodearrint,
                          .edgearrint = edgearrint, .nodearrint2 = nodearrint2, .nodearrchar = nodearrchar };

      SCIP_CALL( reduce_redLoopStp(scip, g, &redbase) );
   }

   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n", reduce_solGetOffset(redsol));

   /* free memory */
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &edgearrint);

   return SCIP_OKAY;
}

/** basic reduction package for the (R)PCSTP */
SCIP_RETCODE reduce_pc(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOL*               redsol,             /**< solution storage */
   GRAPH*                g,                  /**< graph data structure */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             advanced,           /**< perform advanced (e.g. dual ascent) reductions? */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   SCIP_Bool             usestrongreds       /**< allow strong reductions?
                                                  NOTE: needed for propagation, because arcs might have been fixed to 0 */
   )
{
   PATH* vnoi;
   PATH* path;
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

   if( SCIPisLE(scip, (SCIP_Real) nterms / (SCIP_Real) nnodes, STP_BND_THRESHOLD) )
   {
      bred = TRUE;

      if( SCIPisGT(scip, (SCIP_Real) (nterms + graph_pc_nNonLeafTerms(g)) / (SCIP_Real) nnodes, STP_BND_THRESHOLD ) )
         bred = FALSE;
   }

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = reduce_getMinNreductions(g, minelims);

   /* reduction loop */
   SCIP_CALL( reduce_redLoopPc(scip, redsol, g, vnoi, path, nodearrreal, heap, state,
         vbase, nodearrint, edgearrint, nodearrint2, nodearrchar,
         advanced, bred, userec && advanced, reductbound, userec, nodereplacing, usestrongreds) );

   /* free memory */
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
SCIP_RETCODE reduce_mw(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOL*               redsol,             /**< solution storage */
   GRAPH*                g,                  /**< graph data structure */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             advanced,           /**< perform advanced reductions? */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             usestrongreds       /**< allow strong reductions?
                                                  NOTE: needed for propagation, because arcs might have been fixed to 0 */
   )
{
   PATH* vnoi;
   SCIP_Real* nodearrreal;
   int* state;
   int* vbase;
   int* nodearrint;
   const int nnodes = graph_get_nNodes(g);
   const int nterms = graph_get_nTerms(g);
   const int redbound = reduce_getMinNreductions(g, minelims);
   STP_Bool* nodearrchar;
   STP_Bool bred = FALSE;

   assert(scip != NULL);
   assert(redsol != NULL);

   if( SCIPisLE(scip, (SCIP_Real) nterms / (SCIP_Real) nnodes, 0.1) )
      bred = TRUE;

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
   SCIP_CALL( reduce_redLoopMw(scip, redsol, g, vnoi, nodearrreal, state,
         vbase, nodearrint, nodearrchar, advanced, bred, advanced, redbound, userec, usestrongreds) );

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint);

   return SCIP_OKAY;
}

/** basic reduction package for the HCDSTP */
SCIP_RETCODE reduce_hc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   REDSOLLOCAL* redsollocal;
   SCIP_RANDNUMGEN* randnumgen;
   PATH* vnoi;
   SCIP_Real*  cost;
   SCIP_Real*  radius;
   SCIP_Real*  costrev;
   SCIP_Real timelimit;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    pathedge;
   STP_Bool* nodearrchar;
   int     redbound;
   int     degnelims;
   STP_Bool bred = FALSE; // todo reduction method is not correct! ignores direction
   STP_Bool hbred = TRUE;
   STP_Bool rbred = TRUE;
   STP_Bool rcbred = TRUE;
   STP_Bool da = TRUE;
   STP_Bool dahop = TRUE;
   int nrounds = 0;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   degnelims = 0;
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
   SCIP_CALL( reduce_unconnectedForDirected(scip, g) );

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( reduce_sollocalInit(scip, g, &redsollocal) );

   while( (da || bred || hbred || rbred || rcbred || dahop ) && !SCIPisStopped(scip) )
   {
      int ndahopelims = 0;
      int danelims = 0;
      int brednelims = 0;
      int hbrednelims = 0;
      int hcrnelims = 0;
      int hcrcnelims = 0;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;


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

      SCIP_CALL( reduce_unconnectedForDirected(scip, g) );

      if( bred )
      {
         SCIP_Real upperbound = -1.0;

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

      if( da )
      {
         const RPDA paramsda = { .damode = STP_DAMODE_FAST , .useSlackPrune = FALSE, .useRec = FALSE, .extredMode = extred_none, .nodereplacing = FALSE};

         SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &danelims, randnumgen) );
         if( danelims <= redbound )
            da = FALSE;
      }

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( dahop )
      {
         SCIP_CALL( reduce_boundHopDa(scip, g, &ndahopelims, randnumgen) );

         if( ndahopelims <= redbound )
            dahop = FALSE;
      }


      if( da && nrounds == 0 )
      {
         rbred = TRUE;
         hbred = TRUE;
         rcbred = TRUE;
      }

      nrounds++;
   }

   reduce_sollocalFree(scip, &redsollocal);
   SCIPfreeRandom(scip, &randnumgen);

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
SCIP_RETCODE reduce_sap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool             dualascent,         /**< perform dual-ascent reductions? */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH*   vnoi;
   PATH*   path;
   SCIP_Real timelimit;
   SCIP_Real*  nodearrreal;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    nodearrint2;
   int     degtnelims;
   int     redbound;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   STP_Bool da = dualascent; // && (g->stp_type != STP_NWPTSPG);
   STP_Bool sd = TRUE;
   STP_Bool rpt = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   assert(!graph_typeIsUndirected(g));

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   if( g->stp_type == STP_SAP )
   {
      /* todo change .stp file format for SAP! */
      for( int e = 0; e < nedges; e++ )
      {
         if( EQ(g->cost[e], 20000.0) )
            g->cost[e] = FARAWAY;
      }
   }

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
         const RPDA paramsda = { .damode = STP_DAMODE_FAST, .useSlackPrune = FALSE, .useRec = FALSE, .extredMode = extred_none, .nodereplacing = FALSE};
         SCIP_CALL( reduce_da(scip, g, &paramsda, NULL, fixed, &danelims, randnumgen) );

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

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}


/** reduce node-weighted Steiner tree problem */
SCIP_RETCODE reduce_nw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH *vnoi;
   SCIP_Real *nodearrreal;
   SCIP_Real timelimit;
   int *vbase;
   int *nodearrint2;
   int nnodes;
   int redbound;

   STP_Bool    da = TRUE;
   SCIP_RANDNUMGEN* randnumgen;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   nnodes = g->knots;

   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   while( (da) && !SCIPisStopped(scip) )
   {
      int danelims = 0;
      const RPDA paramsda = { .damode = STP_DAMODE_FAST , .useSlackPrune = FALSE, .useRec = FALSE, .extredMode = extred_none, .nodereplacing = FALSE};

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      SCIP_CALL( reduce_da(scip, g, &paramsda, NULL, fixed, &danelims, randnumgen) );

      if( danelims <= 2 * redbound )
         da = FALSE;
   }

   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrreal);

   /* free random number generator */
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}



/** reduce degree constrained Steiner tree problem */
SCIP_RETCODE reduce_dc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   REDSOLLOCAL* redsollocal;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real timelimit;
   const int nnodes = graph_get_nNodes(g);
   const int redbound = MAX(nnodes / 1000, minelims);
   STP_Bool da = TRUE;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( reduce_sollocalInit(scip, g, &redsollocal) );

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   reduce_simple_dc(scip, g);

   while( (da) && !SCIPisStopped(scip) )
   {
      int danelims = 0;
      const RPDA paramsda = { .damode = STP_DAMODE_FAST , .useSlackPrune = FALSE, .useRec = FALSE, .extredMode = extred_none, .nodereplacing = FALSE};

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &danelims, randnumgen) );

      if( danelims <= 2 * redbound )
         da = FALSE;
   }

   reduce_sollocalFree(scip, &redsollocal);
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

/** MWCS loop */
SCIP_RETCODE reduce_redLoopMw(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOL*               redsol,             /**< solution contained */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  state,              /**< shortest path array  */
   int*                  vbase,              /**< voronoi base array  */
   int*                  nodearrint,         /**< nodes-sized array  */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   STP_Bool              advanced,           /**< do advanced reduction? */
   STP_Bool              bred,               /**< do bound-based reduction? */
   STP_Bool              tryrmw,             /**< try to convert problem to RMWCSP? Only possible if advanced = TRUE and userec = TRUE */
   int                   redbound,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             usestrongreds       /**< allow strong reductions? */

   )
{
   SCIP_Real* fixed;
   int* solnode;
   REDSOLLOCAL* redsollocal;
   SCIP_Real da = advanced;
   SCIP_Real timelimit;
   int degelims;
   SCIP_Bool fullrerun = TRUE;
   STP_Bool extensive = STP_RED_EXTENSIVE;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real prizesum;

   assert(scip && g);
   assert(advanced || !tryrmw);

   tryrmw = tryrmw && userec;

   SCIP_CALL( reduce_solInitLocal(scip, g, redsol, &redsollocal));
   fixed = reduce_solGetOffsetPointer(redsol);
   solnode = reduce_sollocalGetSolnode(redsollocal);

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
      const RPDA paramsda = { .damode = STP_DAMODE_FAST , .useRec = userec, .extredMode = extred_none, .nodereplacing = FALSE, .useSlackPrune = TRUE,
            .pcmw_solbasedda = TRUE, .pcmw_useMultRoots =  (g->terms > STP_RED_MWTERMBOUND),
            .pcmw_markroots = tryrmw, .pcmw_fastDa = FALSE };
      int daelims = 0;

      fullrerun = FALSE;

      SCIP_CALL( redLoopInnerMw(scip, g, redsollocal, vnoi, nodearrreal, state,
            vbase, nodearrint, nodearrchar, fixed, da, bred, redbound, userec, prizesum) );

      if( advanced && g->terms > 2 )
      {
         SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );

         if( graph_pc_isRootedPcMw(g) )
         {
            SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &daelims, randnumgen) );
         }
         else
         {
            reduce_sollocalSetOffset(*fixed, redsollocal);
            SCIP_CALL( reduce_sollocalRebuildTry(scip, g, redsollocal) );

            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, redsollocal, vnoi, nodearrreal, vbase,
                  nodearrint, state, nodearrchar, &daelims, randnumgen, prizesum) );
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
                {
                   advanced = TRUE;
                   da = TRUE;
                }
             }
         }

         if( fullrerun )
         {
            SCIPdebugMessage("Restarting reduction loop after %d rounds! (%d DA eliminations) \n\n ", rounds, daelims);

            if( extensive )
               advanced = TRUE;
         }
      }

      SCIP_CALL( reduce_simple_mw(scip, g, solnode, fixed, &degelims) );
   }

   /* go back to the extended graph */
   graph_pc_2trans(scip, g);

   SCIP_CALL( reduce_unconnected(scip, g) );

   SCIPfreeRandom(scip, &randnumgen);

   reduce_solFinalizeLocal(scip, g, redsol);

   return SCIP_OKAY;
}

/** (R)PC loop */
SCIP_RETCODE reduce_redLoopPc(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOL*               redsol,             /**< solution contained */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   PATH*                 path,               /**< path data structure */
   SCIP_Real*            nodearrreal,        /**< nodes-sized array  */
   int*                  heap,               /**< shortest path array  */
   int*                  state,              /**< voronoi base array  */
   int*                  vbase,              /**< nodes-sized array  */
   int*                  nodearrint,         /**< node-sized array  */
   int*                  edgearrint,         /**< edge-sized array  */
   int*                  nodearrint2,        /**< nodes-sized array  */
   STP_Bool*             nodearrchar,        /**< nodes-sized array  */
   SCIP_Bool             dualascent,         /**< do dual-ascent reduction? */
   SCIP_Bool             bred,               /**< do bound-based reduction? */
   SCIP_Bool             tryrpc,             /**< try to transform to rpc? */
   int                   reductbound,        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec,             /**< use recombination heuristic? */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   SCIP_Bool             usestrongreds       /**< allow strong reductions? */
   )
{
   SCIP_Real* fixed;
   int* solnode;
   REDSOLLOCAL* redsollocal;
   DHEAP* dheap;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real timelimit;
   SCIP_Bool rerun = TRUE;
   SCIP_Bool advancedrun = dualascent;
   SCIP_Real prizesum;
   // todo remove
   const SCIP_Bool verbose = FALSE && dualascent && userec && nodereplacing;
   int degnelims;
   int nadvruns = 0;
   SCIP_Bool isSpg = FALSE;
   SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   const SCIP_Bool trySpgTrans = userec; // todo extra parameter
   const int reductbound_global = reductbound * STP_RED_GLBFACTOR;

   SCIP_CALL( reduce_solInitLocal(scip, g, redsol, &redsollocal));
   fixed = reduce_solGetOffsetPointer(redsol);
   solnode = reduce_sollocalGetSolnode(redsollocal);

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( graph_heap_create(scip, g->knots, NULL, NULL, &dheap) );
   graph_pc_2org(scip, g);

   assert(graph_pc_term2edgeIsConsistent(scip, g));

   SCIP_CALL( graph_pc_presolInit(scip, g) );

   SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &degnelims, NULL, solnode) );
   if( verbose ) printf("initial degnelims: %d \n", degnelims);

   prizesum = graph_pc_getPosPrizeSum(scip, g);
   assert(prizesum < FARAWAY);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* main reduction loop */
   for( int outterrounds = 0; outterrounds < STP_RED_MAXNOUTROUNDS_PC && rerun; outterrounds++ )
   {
      int nouterelims = 0;
      const SCIP_Bool beFast = (outterrounds == 0 && advancedrun);
      rerun = FALSE;

      SCIP_CALL( redLoopInnerPc(scip, g, redsollocal, dheap, vnoi, path, nodearrreal, heap, state,
             vbase, nodearrint, edgearrint, nodearrint2, nodearrchar, fixed, randnumgen, prizesum,
             dualascent, bred, reductbound, userec, nodereplacing, usestrongreds, beFast, &nouterelims) );

      if( advancedrun && g->terms > 2 )
      {
         const RPDA paramsda = {
               .damode = STP_DAMODE_MEDIUM,
               .useRec = userec,
               .extredMode = (nadvruns == 0 ? extred_fast : extred_full),
               .nodereplacing = nodereplacing,
               .useSlackPrune = FALSE,
               .pcmw_solbasedda = TRUE, .pcmw_useMultRoots = TRUE, .pcmw_markroots = TRUE, .pcmw_fastDa = FALSE };

         int danelims = 0;
         int implnelims = 0;

         degnelims = 0;
         advancedrun = FALSE;
         nadvruns++;

         if( rpc )
         {
            if( trySpgTrans && graph_pc_nNonLeafTerms(g) == 0 && graph_pc_nProperPotentialTerms(g) == 0 )
            {
               isSpg = TRUE;
               break;
            }

            SCIP_CALL( reduce_impliedProfitBasedRpc(scip, g, redsollocal, fixed, &implnelims) );
            if( verbose ) printf("implnelims: %d \n", implnelims);

            SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, fixed, &danelims, randnumgen) );
         }
         else
         {
            reduce_sollocalSetOffset(*fixed, redsollocal);
            SCIP_CALL( reduce_sollocalRebuildTry(scip, g, redsollocal) );
            SCIP_CALL( reduce_daPcMw(scip, g, &paramsda, redsollocal, vnoi, nodearrreal, vbase, heap,
                  state, nodearrchar, &danelims, randnumgen, prizesum) );
         }

         SCIP_CALL( reduce_simple_pc(scip, NULL, g, fixed, &degnelims, NULL, solnode) );
         nouterelims += danelims + degnelims + implnelims;

         if( nouterelims + implnelims > reductbound_global )
         {
            if( danelims + implnelims > reductbound_global )
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
            SCIP_CALL(reduce_unconnectedRpcRmw(scip, g, fixed));
            rerun = TRUE;
            advancedrun = dualascent;
            outterrounds = 1;

            if( trySpgTrans && graph_pc_nNonLeafTerms(g) == 0 && graph_pc_nProperPotentialTerms(g) == 0 )
            {
               isSpg = TRUE;
               break;
            }
         }
      }
   }

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

   reduce_solFinalizeLocal(scip, g, redsol);

   if( isSpg )
   {
#ifndef WITH_UG
      printf("changing problem to SPG \n");
      graph_printInfo(g);
#endif
      SCIP_CALL( graph_transRpc2SpgTrivial(scip, g) );
   }

   SCIPdebugMessage("Reduction Level PC 1: Fixed Cost = %.12e\n", *fixed);
   return SCIP_OKAY;
}

/** STP loop */
SCIP_RETCODE reduce_redLoopStp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase             /**< parameters */
)
{
   REDSOLLOCAL* redsollocal;
   SCIP_RANDNUMGEN* randnumgen;
   const RPARAMS* redparameters = redbase->redparameters;
   SCIP_Real timelimit;
   const int reductbound = redparameters->reductbound;
   const SCIP_Bool fullreduce = redparameters->fullreduce;
   const SCIP_Bool nodereplacing = redparameters->nodereplacing;
   int dummy = 0;
   int nruns = 0;
   SCIP_Bool rerun = TRUE;
#ifdef USE_FULLSEPA
   SCIP_Bool runFullSepa = TRUE;
#endif

   assert(scip && g && redbase);
   assert(reductbound > 0);
   assert(graph_valid(scip, g));
   assert(graph_typeIsSpgLike(g));

   SCIP_CALL( reduce_solInitLocal(scip, g, redbase->redsol, &redsollocal));
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   SCIP_CALL( reduce_contract0Edges(scip, g, redbaseGetSolnode(redsollocal, redbase), TRUE) );
   SCIP_CALL( reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &dummy, NULL) );

   /* reduction loop */
   do
   {
      SCIP_Bool wasDecomposed;
      rerun = FALSE;

      SCIP_CALL( redLoopInnerStp(scip, randnumgen, g, redsollocal, redbase, &wasDecomposed) );

      if( wasDecomposed )
      {
         SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &dummy, NULL));
         break;
      }

      if( fullreduce && !SCIPisStopped(scip) )
      {
         int extendedelims = 0;
         const RPDA paramsda = {
                         .damode = (nruns == 0 ? STP_DAMODE_FAST : STP_DAMODE_MEDIUM),
                         .useRec = redparameters->userec,
                         .extredMode = (nruns == 0 ? extred_fast : extred_full),
                         .useSlackPrune = FALSE,
                         .nodereplacing = nodereplacing};

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

         if( redparameters->userec && (SCIPStpEnumRelaxIsPromising(g) || dpterms_isPromisingFully(g)) )
            break;

         assert(!rerun);

         SCIP_CALL( reduce_da(scip, g, &paramsda, redsollocal, redbaseGetOffsetPointer(redbase), &extendedelims, randnumgen) );
         reduceStatsPrint(fullreduce, "ext", extendedelims);

         SCIP_CALL(reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &extendedelims, NULL));

         if( nodereplacing )
         {
            int conflictnelims = 0;

            SCIP_CALL( reduce_fixedConflicts(scip, NULL, g, &conflictnelims) );
            reduceStatsPrint(fullreduce, "fixedconflict", conflictnelims);
            extendedelims += conflictnelims;
         }

         if( extendedelims > STP_RED_EXPENSIVEFACTOR * reductbound )
         {
            rerun = TRUE;
         }

#ifdef USE_FULLSEPA
         // todo tune!
         if( !rerun && runFullSepa )
         {
            int sepanelims = 0;

            if( dpterms_isPromisingFully(g) )
               break;

            SCIP_CALL( reduce_termsepaFull(scip, g, redbaseGetSolnode(redsollocal, redbase), redbase, &sepanelims) );

            if( sepanelims > STP_RED_EXPENSIVEFACTOR * reductbound )
               SCIP_CALL( redLoopInnerStp(scip, randnumgen, g, redsollocal, redbase, &wasDecomposed) );
            else if( sepanelims > reductbound )
               SCIP_CALL( reduce_simple(scip, g, redbaseGetOffsetPointer(redbase), redbaseGetSolnode(redsollocal, redbase), &dummy, NULL));

            runFullSepa = FALSE;
         }
#endif
      }

      nruns++;
   }
   while( rerun && !SCIPisStopped(scip) ); /* reduction loop */

   assert(graph_valid_ancestors(scip, g));
   SCIPfreeRandom(scip, &randnumgen);

   reduce_solFinalizeLocal(scip, g, redbase->redsol);

   return SCIP_OKAY;
}


/** reduces the graph */
SCIP_RETCODE reduce_exec(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph structure */
   REDSOL*               redsol,             /**< primal solution container */
   int                   reductionlevel,     /**< reduction level 0: none, 1: basic, 2: advanced */
   int                   minelims,           /**< minimal amount of reductions to reiterate reduction methods */
   SCIP_Bool             userec              /**< use recombination heuristic? */
   )
{
   int stp_type;
   SCIP_Real* offset = reduce_solGetOffsetPointer(redsol);
   SCIP_Bool graphWasInitialized = FALSE;

   assert(scip && offset);
   assert(graph_getFixedges(graph) == NULL);
   assert(reductionlevel == STP_REDUCTION_NONE || reductionlevel == STP_REDUCTION_BASIC || reductionlevel == STP_REDUCTION_ADVANCED );
   assert(minelims >= 0);
   assert(EQ(*offset, 0.0));

   stp_type = graph->stp_type;

   if( !graph_isSetUp(graph) )
   {
      SCIP_CALL( graph_initHistory(scip, graph) );
   }

   if( !graph_path_exists(graph) )
   {
      SCIP_CALL( graph_path_init(scip, graph) );
      graphWasInitialized = TRUE;
   }

   /* start with trivial preprocessing step */
   SCIP_CALL( reduce_unconnected(scip, graph) );

   /* if no reduction methods available for given problem, return */
   if( graph->stp_type == STP_BRMWCSP )
   {
      graph_path_exit(scip, graph);
      return SCIP_OKAY;
   }

   if( reductionlevel == STP_REDUCTION_BASIC )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reduce_pc(scip, redsol, graph, minelims, FALSE, FALSE, TRUE, TRUE) );
      }
      else if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
      {
         SCIP_CALL( reduce_mw(scip, redsol, graph, minelims, FALSE, FALSE, TRUE) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduce_hc(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_SAP || stp_type == STP_NWPTSPG )
      {
         SCIP_CALL( reduce_sap(scip, graph, FALSE, offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduce_nw(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_DCSTP )
      {
         SCIP_CALL( reduce_dc(scip, graph, offset, minelims) );
      }
      else
      {
         assert(graph_typeIsSpgLike(graph));

         SCIP_CALL( reduce_stp(scip, graph, redsol, minelims, FALSE, TRUE, FALSE, TRUE) );
      }
   }
   else if( reductionlevel == STP_REDUCTION_ADVANCED )
   {
      if( stp_type == STP_PCSPG || stp_type == STP_RPCSPG )
      {
         SCIP_CALL( reduce_pc(scip, redsol, graph, minelims, TRUE, userec, TRUE, TRUE) );

         if( graph->stp_type == STP_SPG )
         {
            reduce_solReInitLocal(graph, redsol);
            SCIP_CALL( reduce_stp(scip, graph, redsol, minelims, TRUE, TRUE, userec, TRUE) );
         }
      }
      else if( stp_type == STP_MWCSP || stp_type == STP_RMWCSP )
      {
         SCIP_CALL( reduce_mw(scip, redsol, graph, minelims, TRUE, userec, TRUE) );
      }
      else if( stp_type == STP_DHCSTP )
      {
         SCIP_CALL( reduce_hc(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_SAP || stp_type == STP_NWPTSPG )
      {
         SCIP_CALL( reduce_sap(scip, graph, TRUE, offset, minelims) );
      }
      else if( stp_type == STP_NWSPG )
      {
         SCIP_CALL( reduce_nw(scip, graph, offset, minelims) );
      }
      else if( stp_type == STP_DCSTP )
      {
         SCIP_CALL( reduce_dc(scip, graph, offset, minelims) );
      }
      else
      {
         assert(graph_typeIsSpgLike(graph));

         SCIP_CALL( reduce_stp(scip, graph, redsol, minelims, TRUE, TRUE, userec, TRUE) );
      }
   }

   SCIP_CALL( reduce_unconnected(scip, graph) );

   /* NOTE: ugly, but necessary to allow for clean execution of decomposition routine during solving process  */
   if( graph_typeIsSpgLike(graph) && userec )
   {
      int nartelims = 0;
      SCIP_CALL( reduce_articulations(scip, graph, offset, &nartelims) );
   }

   SCIPdebugMessage("offset : %f \n", *offset);

   assert(graph_valid(scip, graph));

   if( graphWasInitialized )
      graph_path_exit(scip, graph);

   return SCIP_OKAY;
}
