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

/**@file   heur_rec.c
 * @brief  Primal recombination heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a recombination heuristic for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" (2017) by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_rec.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_rec.h"
#include "heur_prune.h"
#include "heur_slackprune.h"
#include "heur_local.h"
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"
#include "solstp.h"
#include "math.h"

#define HEUR_NAME             "rec"
#define HEUR_DESC             "recombination heuristic for Steiner problems"
#define HEUR_DISPCHAR         'R'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_MAXFREQREC     FALSE         /**< executions of the rec heuristic at maximum frequency?                             */
#define DEFAULT_MAXNSOLS       15            /**< default maximum number of (good) solutions be regarded in the subproblem          */
#define DEFAULT_NUSEDSOLS      4             /**< number of solutions that will be taken into account                               */
#define DEFAULT_RANDSEED       1984          /**< random seed                                                                       */
#define DEFAULT_NTMRUNS        100           /**< number of runs in TM heuristic                                                    */
#define DEFAULT_NWAITINGSOLS   4             /**< max number of new solutions to be available before executing the heuristic again  */

#define BOUND_MAXNTERMINALS 1000
#define BOUND_MAXNEDGES     20000
#define RUNS_RESTRICTED     3
#define RUNS_NORMAL         10
#define CYCLES_PER_RUN      3
#define REC_MAX_FAILS       4
#define REC_MIN_NSOLS       4
#define REC_ADDEDGE_FACTOR 1.0

#define COST_MAX_POLY_x0 500
#define COST_MIN_POLY_x0 100


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   maxnsols;           /**< maximum number of (good) solutions be regarded in the subproblem  */
   SCIP_Longint          ncalls;             /**< number of calls                                                   */
   SCIP_Longint          nlastsols;          /**< number of solutions during the last run                           */
   int                   ntmruns;            /**< number of runs in TM heuristic                                    */
   int                   nusedsols;          /**< number of solutions that will be taken into account               */
   int                   nselectedsols;      /**< number of solutions actually selected                             */
   int                   nwaitingsols;       /**< number of new solutions before executing the heuristic again      */
   int                   nfailures;          /**< number of failures since last successful call                     */
   unsigned int          randseed;           /**< seed value for random number generator                            */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                                           */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};


/*
 * Local methods
 */

/** information method for a parameter change of random seed */
static
SCIP_DECL_PARAMCHGD(paramChgdRandomseed)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int newrandseed;

   newrandseed = SCIPparamGetInt(param);

   heurdata = (SCIP_HEURDATA*)SCIPparamGetData(param);
   assert(heurdata != NULL);

   heurdata->randseed = (int) newrandseed;

   return SCIP_OKAY;
}


/** edge cost multiplier */
static
int edgecostmultiplier(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Real             avg                 /**< number of solutions containing this edge */
   )
{
   SCIP_Real normed;
   int maxpolyx0;
   int minpolyx0;

   int upper;
   int lower;
   int factor;

   /* if STP, then avg >= 2.0 */
   assert(SCIPisGE(scip, avg, 1.0));
   assert(heurdata->nusedsols  >= 2);

   maxpolyx0 = COST_MAX_POLY_x0 * (heurdata->nusedsols - 1);
   minpolyx0 = COST_MIN_POLY_x0 * (heurdata->nusedsols - 1);

   /* norm to [0,1] (STP) or [-1,1] and evaluate polynomials */
   normed = (avg - 2.0) / ((double) heurdata->nusedsols - 1.0);
   upper = (int) (maxpolyx0 - 2 * maxpolyx0 * normed + maxpolyx0 * (normed * normed));
   lower = (int) (minpolyx0 - 2 * minpolyx0 * normed + minpolyx0 * (normed * normed));

   assert(SCIPisGE(scip, normed, -1.0) && SCIPisLE(scip, normed, 1.0));

   lower = MAX(0, lower);
   upper = MAX(0, upper);

   factor = SCIPrandomGetInt(heurdata->randnumgen, lower, upper);
   factor++;

   assert(factor <= COST_MAX_POLY_x0 * 3 + 1 || avg <= 2.0);
   assert(factor >= 1);

   return factor;
}


static
void marksolverts(
   GRAPH* g,
   IDX* curr,
   int* unodemap,
   STP_Bool* stvertex
   )
{
   while( curr != NULL )
   {
      int i = curr->index;

      stvertex[unodemap[ g->orghead[i] ]] = TRUE;
      stvertex[unodemap[ g->orgtail[i] ]] = TRUE;

      curr = curr->parent;
   }
}


static
void setHeurdataNUsedSols(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   int                   nsols,              /**< number of solutions */
   SCIP_Bool             restrictheur        /**< restrict heuristic? */
)
{
   int randn;

   if( restrictheur )
      randn = SCIPrandomGetInt(heurdata->randnumgen, 0, 3);
   else
      randn = SCIPrandomGetInt(heurdata->randnumgen, 0, 6);

   if( (randn <= 2) || (nsols < 3) )
      heurdata->nusedsols = 2;
   else if( (randn <= 4) || (nsols < 4) )
      heurdata->nusedsols = 3;
   else
      heurdata->nusedsols = 4;
}


/** compute (heuristic) solution on reduced (and merged) problem */
static
SCIP_RETCODE retransformReducedProbSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data */
   const GRAPH*          solgraph,           /**< graph data */
   IDX**                 ancestors,          /**< ancestors */
   const int*            edgeancestor,       /**< ancestor to edge */
   const int*            soledges,           /**< old solution edges */
   int*                  newsoledges,        /**< new solution edges */
   STP_Bool*             stnodes             /**< new solution nodes */
)
{
   STP_Bool* solnodes;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   for( int i = 0; i < nedges; i++ )
       newsoledges[i] = UNKNOWN;

    for( int i = 0; i < nnodes; i++ )
       stnodes[i] = FALSE;

    if( pcmw )
    {
       SCIP_CALL( SCIPallocBufferArray(scip, &solnodes, solgraph->orgknots) );

       for( int i = 0; i < solgraph->orgknots; i++ )
          solnodes[i] = FALSE;
    }

   if( solgraph->terms > 1 )
   {
      assert(soledges != NULL);
      for( int e = 0; e < solgraph->edges; e++ )
      {
         if( soledges[e] == CONNECT )
         {
            /* iterate through list of ancestors */
            if( graph->stp_type != STP_DCSTP )
            {
               IDX* curr = ancestors[e];

               while( curr != NULL )
               {
                  const int i = edgeancestor[curr->index];

                  stnodes[graph->head[i]] = TRUE;
                  stnodes[graph->tail[i]] = TRUE;

                  if( pcmw )
                  {
                     assert(solgraph->orghead[curr->index] < solgraph->orgknots);
                     assert(solgraph->orgtail[curr->index] < solgraph->orgknots);

                     solnodes[solgraph->orghead[curr->index]] = TRUE;
                     solnodes[solgraph->orgtail[curr->index]] = TRUE;
                  }

                  curr = curr->parent;
               }
            }
            else
            {
               IDX* curr = ancestors[e];
               while( curr != NULL )
               {
                  int i = edgeancestor[curr->index];
                  newsoledges[i] = CONNECT;

                  curr = curr->parent;
               }
            }
         }
      }
   } /* if( solgraph->terms > 1 ) */

   /* retransform edges fixed during graph reduction */
   if( graph->stp_type != STP_DCSTP )
   {
      IDX* curr = graph_get_fixedges(solgraph);

      while( curr != NULL )
      {
         const int i = edgeancestor[curr->index];

         stnodes[graph->head[i]] = TRUE;
         stnodes[graph->tail[i]] = TRUE;
         if( pcmw )
         {
            assert(solgraph->orghead[curr->index] < solgraph->orgknots);
            assert(solgraph->orgtail[curr->index] < solgraph->orgknots);

            solnodes[solgraph->orghead[curr->index]] = TRUE;
            solnodes[solgraph->orgtail[curr->index]] = TRUE;
         }

         curr = curr->parent;
      }
   }
   else
   {
      IDX* curr = graph_get_fixedges(solgraph);

      while( curr != NULL )
      {
         const int i = edgeancestor[curr->index];
         newsoledges[i] = CONNECT;
      }
   }

   if( pcmw )
   {
      STP_Bool* pcancestoredges;
      SCIP_CALL(SCIPallocBufferArray(scip, &pcancestoredges, solgraph->orgedges));

      for( int k = 0; k < solgraph->orgedges; k++ )
         pcancestoredges[k] = FALSE;

      SCIP_CALL(
            solstp_markPcancestors(scip, solgraph->pcancestors, solgraph->orgtail, solgraph->orghead, solgraph->orgknots, solnodes, pcancestoredges, NULL, NULL, NULL ));

      for( int k = 0; k < solgraph->orgedges; k++ )
         if( pcancestoredges[k] )
         {
            const int i = edgeancestor[k];
            stnodes[graph->tail[i]] = TRUE;
            stnodes[graph->head[i]] = TRUE;
         }

      SCIPfreeBufferArray(scip, &pcancestoredges);
      SCIPfreeBufferArrayNull(scip, &solnodes);
   }

   return SCIP_OKAY;
}


/** compute (heuristic) solution on reduced (and merged) problem */
static
SCIP_RETCODE computeReducedProbSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   const GRAPH*          graph,              /**< graph data */
   GRAPH*                solgraph,           /**< graph data */
   IDX**                 ancestors,          /**< ancestors */
   const int*            edgeweight,         /**< for each edge: number of solution that contain this edge */
   const int*            edgeancestor,       /**< ancestor to edge edge */
   SCIP_VAR**            vars,               /**< variables or NULL */
   SCIP_Bool             usestppool,         /**< using STP pool? */
   int*                  soledges            /**< solution edges */
)
{
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real* prize = NULL;
   SCIP_Real* nodepriority;
   const int nsolnodes = solgraph->knots;
   const int nsoledges = solgraph->edges;
   const int probtype = graph->stp_type;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(soledges != NULL);
   assert(graph_valid(scip, solgraph));

   SCIPdebugMessage("REC: graph not completely reduced, nodes: %d, edges: %d, terminals: %d \n", solgraph->knots, nsoledges, solgraph->terms);

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nsoledges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nsoledges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodepriority, nsolnodes) );

   for( int i = 0; i < nsolnodes; i++ )
      nodepriority[i] = 0.0;

   BMScopyMemoryArray(cost, solgraph->cost, nsoledges);

   for( int e = 0; e < nsoledges; e++ )
   {
      IDX* curr = ancestors[e];
      SCIP_Real avg = 0.0;
      int i = 0;
      SCIP_Bool fixed = FALSE;

      while( curr != NULL )
      {
         i++;
         avg += edgeweight[curr->index];
         if( !usestppool && SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5
               && SCIPvarGetUbGlobal(vars[flipedge(edgeancestor[curr->index])] ) < 0.5 )
            fixed = TRUE;

         curr = curr->parent;
      }

      if( i > 0 )
      {
         avg /= (double) i;
         assert(avg >= 1.0);
      }

      /* is an ancestor edge fixed? */
      if( fixed )
      {
         cost[e] = MAX(cost[e], BLOCKED);

         nodepriority[solgraph->head[e]] /= 2.0;
         nodepriority[solgraph->tail[e]] /= 2.0;
      }

      nodepriority[solgraph->head[e]] += avg - 1.0;
      nodepriority[solgraph->tail[e]] += avg - 1.0;

      cost[e] *= edgecostmultiplier(scip, heurdata, avg);
   }

   /* adapted prizes */
   if( pcmw )
   {
      assert(solgraph->extended);

      SCIP_CALL( SCIPallocBufferArray(scip, &prize, nsolnodes) );

      for( int k = 0; k < nsolnodes; k++ )
      {
         if( Is_pseudoTerm(solgraph->term[k]) )
         {
            const int term = graph_pc_getTwinTerm(solgraph, k);
            const int root2termedge = graph_pc_getRoot2PtermEdge(solgraph, term);

            assert(term >= 0 && root2termedge >= 0);
            assert(solgraph->source != k);
            assert(Is_term(solgraph->term[term]));
            assert(solgraph->cost[root2termedge] > 0.0);

            prize[k] = cost[root2termedge];
         }
         else
         {
            prize[k] = solgraph->prize[k];
         }
      }
   }

   for( int e = 0; e < nsoledges; e++ )
   {
      costrev[e] = cost[flipedge(e)];
      soledges[e] = UNKNOWN;
   }

   SCIP_CALL( graph_path_init(scip, solgraph) );

   // todo: run prune heuristic with changed weights!
   {
      SCIP_Bool success;
      SCIP_Real hopfactor = 0.1;

      SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fromheurdata, solgraph, NULL, prize, soledges, heurdata->ntmruns,
         solgraph->source, cost, costrev, &hopfactor, nodepriority, &success) );

      assert(SCIPisStopped(scip) || success);
      assert(SCIPisStopped(scip) || solstp_isValid(scip, solgraph, soledges));
   }

   SCIPfreeBufferArrayNull(scip, &prize);

   /* run local heuristic with original costs! */
   if( !SCIPisStopped(scip) && probtype != STP_DHCSTP && probtype != STP_DCSTP
         && probtype != STP_SAP && probtype != STP_NWSPG && probtype != STP_NWPTSPG )
   {
      SCIP_CALL( SCIPStpHeurLocalRun(scip, solgraph, soledges) );
      assert(solstp_isValid(scip, solgraph, soledges));
   }

   graph_path_exit(scip, solgraph);

   SCIPfreeBufferArray(scip, &nodepriority);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/** initialize solution vector and objective of incumbent */
static
void initializeIncumbent(
   SCIP*                 scip,               /**< SCIP data structure */
   const STPSOLPOOL*     pool,               /**< solution pool or NULL */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nsols,              /**< number of solutions in pool (SCIP or STP) */
   int                   newsolindex,        /**< index of new solution */
   double*               incumentobj,        /**< pointer to incumbent value */
   int*                  incumbentedges      /**< edges of solution to be used as recombination root */
)
{
   const STP_Bool usestppool = (pool != NULL);
   const int nedges = graph->edges;

   if( usestppool )
   {
      int i;
      STPSOL** poolsols = pool->sols;

      assert(pool->nedges == graph->edges);

      for( i = 0; i < nsols; i++ )
         if( newsolindex == poolsols[i]->index )
            break;
      assert(i < nsols);

      BMScopyMemoryArray(incumbentedges, poolsols[i]->soledges, nedges);
   }
   else
   {
      SCIP_SOL* newsol;
      SCIP_SOL** sols = SCIPgetSols(scip);
      int i;

      assert(sols != NULL);

      for( i = 0; i < nsols; i++ )
         if( newsolindex == SCIPsolGetIndex(sols[i]) )
            break;
      assert(i < nsols);

      newsol = sols[i];

      for( int e = 0; e < nedges; e++ )
         if( SCIPisEQ(scip, SCIPgetSolVal(scip, newsol, vars[e]), 1.0) )
            incumbentedges[e] = CONNECT;
         else
            incumbentedges[e] = UNKNOWN;
   }

   /* get objective value of incumbent */
   *incumentobj = solstp_getObjBounded(graph, incumbentedges, 0.0, nedges);
}


#ifndef NDEBUG
/** any edge of given solution deleted? */
static
SCIP_Bool poolSolIsUnreduced(
   const GRAPH*          graph,              /**< graph structure */
   const int             solindex,           /**< the index */
   const STPSOLPOOL*     pool                /**< the pool */
)
{
   const int nedges = graph->edges;
   STPSOL* const poolsol = pool->sols[solindex];

   assert(pool);
   assert(solindex >= 0 && solindex < pool->size);

   for( int i = 0; i < nedges; i += 2 )
   {
      if( graph->oeat[i] != EAT_FREE )
         continue;

      if( (poolsol->soledges[i] == CONNECT || poolsol->soledges[i + 1] == CONNECT) )
         return FALSE;
   }

   return TRUE;
}
#endif


/** store solution in (SCIP or STP) pool */
static
SCIP_RETCODE poolAddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL*           pool,               /**< STP solution pool or NULL */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   const SCIP_Real       incumentobj,        /**< objective of solution to be added */
   const int*            incumbentedges,     /**< edge array of solution to be added */
   int                   nedges,             /**< nedges */
   int*                  newsolindex,        /**< index of new solution */
   SCIP_Bool*            soladded            /**< solution added? */
)
{
   SCIPdebugMessage("REC: better solution found ...  ");

   if( pool != NULL )
   {
      const int maxindex = pool->maxindex;
      SCIP_CALL( solpool_addSol(scip, incumentobj, incumbentedges, pool, soladded) );

      if( *soladded )
      {
         assert(pool->maxindex == maxindex + 1);
         *newsolindex = maxindex + 1;
         SCIPdebugMessage("and added! \n");
      }
      else
      {
         *newsolindex = -1;
      }
   }
   else
   {
      SCIP_Real* nval = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );

      for( int e = 0; e < nedges; e++ )
      {
         if( incumbentedges[e] == CONNECT )
            nval[e] = 1.0;
         else
            nval[e] = 0.0;
      }

      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, heur, soladded) );

      if( *soladded )
      {
         SCIP_SOL** sols = SCIPgetSols(scip);
         const int nsols = SCIPgetNSols(scip);
         int idx = -1;

         SCIPdebugMessage("and added! \n");

         assert(nsols > 0);

         for( int i = 1; i < nsols; i++ )
            if( SCIPsolGetIndex(sols[i]) > idx )
               idx = SCIPsolGetIndex(sols[i]);

         assert(idx >= 0);
         *newsolindex = idx;
      }
      SCIPfreeBufferArray(scip, &nval);

   }
   return SCIP_OKAY;
}


/** select solutions to be merged */
static
SCIP_RETCODE solgraphSelectSolsDiff(
   SCIP*                 scip,               /**< SCIP data structure */
   const STPSOLPOOL*     pool,               /**< solution pool or NULL */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   const int*            incumbentedges,     /**< edges of solution to be used as recombination root */
   int*                  incumbentindex,     /**< index of ancestor of incumbent solution */
   int*                  selection,          /**< selected solutions */
   SCIP_Bool*            success             /**< could at least two solutions be selected? */
   )
{
   SCIP_SOL** sols = NULL;
   STPSOL** poolsols = NULL;
   int* perm;
   int* soledgestmp;
   STP_Bool* solselected;
   STP_Bool* soledges;
   int pos;
   int nsols;
   int maxnsols = heurdata->maxnsols;
   int nselectedsols = 0;
   const int nedges = graph->edges;
   const int nusedsols = heurdata->nusedsols;
   const SCIP_Bool usestppool = (pool != NULL);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(selection != NULL);
   assert(graph != NULL);

   if( !usestppool )
   {
      assert(vars != NULL);

      /* get solution data */
      sols = SCIPgetSols(scip);
      nsols = SCIPgetNSols(scip);
   }
   else
   {
      poolsols = pool->sols;
      nsols = pool->size;
      assert(nsols > 1);
   }

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledges, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledgestmp, nedges / 2) );

   for( int i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      solselected[i] = FALSE;
   }

   if( usestppool )
   {
      for( pos = 0; pos < nsols; pos++ )
         if( *incumbentindex == poolsols[pos]->index )
            break;
   }
   else
   {
      for( pos = 0; pos < nsols; pos++ )
         if( *incumbentindex == SCIPsolGetIndex(sols[pos]) )
            break;
   }
   assert(pos < nsols);
   solselected[pos] = TRUE;
   selection[nselectedsols++] = pos;

   // mark current solution
   for( int e = 0; e < nedges; e += 2 )
      soledges[e / 2] = (incumbentedges[e] == CONNECT
                      || incumbentedges[e + 1] == CONNECT);

   maxnsols = MIN(nsols, maxnsols);

   if( !usestppool )
   {
      const int sqrtnallsols = (int) sqrt((double) nsols);

      if( sqrtnallsols >= REC_MIN_NSOLS && sqrtnallsols < maxnsols )
         maxnsols = sqrtnallsols;
   }

   SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, maxnsols);

   // check solutions to be merged
   for( int i = 0; i < nsols; i++ )
   {
      if( solselected[perm[i]] == FALSE )
      {
         int eqnedges = 0;
         int diffnedges = 0;
         const int k = perm[i];
         SCIP_SOL* solk = NULL;
         const int* solKedges = NULL;

         if( usestppool )
            solKedges = poolsols[k]->soledges;
         else
            solk = sols[k];

         for( int e = 0; e < nedges; e += 2 )
         {
            SCIP_Bool hit;

            if( usestppool )
               hit = (solKedges[e] == CONNECT || solKedges[e + 1] == CONNECT);
            else
               hit = SCIPisEQ(scip, SCIPgetSolVal(scip, solk, vars[e]), 1.0)
                  || SCIPisEQ(scip, SCIPgetSolVal(scip, solk, vars[e + 1]), 1.0);

            if( hit )
            {
               if( soledges[e / 2] == FALSE )
                  soledgestmp[diffnedges++] = e / 2;
               else
               {
                  const int tail = graph->tail[e];
                  const int head = graph->head[e];

                  /* no edge between root and dummy terminal? */
                  if( !pcmw || !(Is_term(graph->term[tail]) && Is_term(graph->term[head])) )
                     eqnedges++;
               }
            }
         }

         /* enough similarities and differences with new solution? */
         if( diffnedges > 5 && eqnedges > 0  )
         {
            SCIPdebugMessage("REC: different edges: %d same edges: %d \n", diffnedges, eqnedges);
            selection[nselectedsols++] = k;
            solselected[k] = TRUE;
            *success = TRUE;

            for( int j = 0; j < diffnedges; j++ )
               soledges[soledgestmp[j]] = TRUE;

            if( nselectedsols >= nusedsols )
               break;
         }
      }
   }

   assert(nselectedsols <= nusedsols);

   heurdata->nselectedsols = nselectedsols;

   SCIPdebugMessage("REC: selected %d sols \n", nselectedsols);

   /* free memory */
   SCIPfreeBufferArray(scip, &soledgestmp);
   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &solselected);

   return SCIP_OKAY;
}

/** select solutions to be merged */
static
SCIP_RETCODE solgraphSelectSols(
   SCIP*                 scip,               /**< SCIP data structure */
   const STPSOLPOOL*     pool,               /**< solution pool or NULL */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   int*                  incumbentindex,     /**< index of ancestor of incumbent solution */
   int*                  selection,          /**< selected solutions */
   SCIP_Bool             randomize           /**< select solutions randomly? */
   )
{
   SCIP_SOL** sols = NULL;                          /* array of all solutions found so far         */
   STPSOL** poolsols = NULL;
   int* perm;
   STP_Bool* solselected;
   int i;
   int nsols;                                /* number of all available solutions        */
   int maxnsols;
   int nusedsols;                            /* number of solutions to use in rec           */
   int nselectedsols;

   const SCIP_Bool usestppool = (pool != NULL);
   assert(selection != NULL);

   if( usestppool )
   {
      poolsols = pool->sols;
      nsols = pool->size;
      assert(nsols > 1);
   }
   else
   {
      /* get solution data */
      sols = SCIPgetSols(scip);
      nsols = SCIPgetNSols(scip);
   }

   maxnsols = heurdata->maxnsols;
   nusedsols = heurdata->nusedsols;
   assert(nusedsols <= nsols);
   nselectedsols = 0;

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      solselected[i] = FALSE;
   }

   if( usestppool )
   {
      for( i = 0; i < nsols; i++ )
         if( *incumbentindex == poolsols[i]->index )
            break;
   }
   else
   {
      for( i = 0; i < nsols; i++ )
         if( *incumbentindex == SCIPsolGetIndex(sols[i]) )
            break;
   }

   assert(i < nsols);
   solselected[i] = TRUE;
   selection[nselectedsols++] = i;

   if( !randomize )
   {
      const int end = SCIPrandomGetInt(heurdata->randnumgen, 1, nusedsols - 1);
      int shift = SCIPrandomGetInt(heurdata->randnumgen, end, 2 * nusedsols - 1);

      if( shift > nsols )
         shift = nsols;

      SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, shift);

      for( i = 0; i < end; i++ )
      {
         if( solselected[perm[i]] == FALSE )
         {
            selection[nselectedsols++] = perm[i];
            solselected[perm[i]] = TRUE;
         }
      }
   }

   maxnsols = MIN(nsols, maxnsols);

   if( !usestppool )
   {
      const int sqrtnallsols = (int) sqrt(nsols);

      if( sqrtnallsols >= REC_MIN_NSOLS && sqrtnallsols < maxnsols )
         maxnsols = sqrtnallsols;
   }

   SCIPdebugMessage("maxnsols in REC %d \n", maxnsols);

   if( nselectedsols < nusedsols )
   {
      SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, maxnsols);
      for( i = 0; i < maxnsols; i++ )
      {
         if( solselected[perm[i]] == FALSE )
         {
            selection[nselectedsols++] = perm[i];
            if( nselectedsols >= nusedsols )
               break;
         }
      }
   }

   assert(nselectedsols <= nusedsols);

   SCIPdebugMessage("REC: selected %d sols \n", nselectedsols);

   heurdata->nselectedsols = nselectedsols;
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &solselected);

   return SCIP_OKAY;
}


/** adapt merged solution graph for PC/MW */
static
void solgraphAdaptForPcMw(
   const GRAPH*          graph,              /**< graph structure */
   STP_Bool*             solnode,            /**< solution nodes array */
   int*                  nsolnodes,          /**< number of nodes pointer */
   STP_Bool*             soledge,            /**< solution edges array  */
   int*                  nsoledges           /**< number of edges pointer */
   )
{
   const int rpcmw = graph_pc_isRootedPcMw(graph);
   const int oldroot = graph->source;
   const int nnodes = graph->knots;

   assert(solnode && soledge && nsolnodes && nsoledges);
   assert(graph->extended);
   assert(solnode[oldroot]);

   /* remove unnecessary dummy edges and add necessary ones */
   for( int k = 0; k < nnodes; k++ )
   {
      /* dummy terminal other than the root? */
      if( solnode[k] && graph_pc_knotIsDummyTerm(graph, k) && k != oldroot )
      {
         const int edgeroot = graph_pc_getRoot2PtermEdge(graph, k);
         const int edge2pterm = graph->term2edge[k];
         const int pterm = graph->head[edge2pterm];

         assert(graph->grad[k] == 2 && edge2pterm >= 0 && Is_pseudoTerm(graph->term[pterm]));

         if( !solnode[pterm] )
         {
            assert(soledge[edgeroot / 2]);
            assert(!soledge[edge2pterm / 2]);

            soledge[edgeroot / 2] = FALSE;
            solnode[k] = FALSE;
            (*nsoledges)--;
            (*nsolnodes)--;
         }
         else if( !soledge[edgeroot / 2] )
         {
            soledge[edgeroot / 2] = TRUE;
            (*nsoledges)++;
         }
      }

      /* pseudo terminal? */
      if( solnode[k] && Is_pseudoTerm(graph->term[k]) )
      {
         const int edge2term = graph->term2edge[k];
         const int term = graph->head[edge2term];

         if( !rpcmw )
         {
            const int edgeroot = graph_pc_getRoot2PtermEdge(graph, k);

            if( !soledge[edgeroot / 2] )
            {
               soledge[edgeroot / 2] = TRUE;
               (*nsoledges)++;
            }
         }

         if( !solnode[term] )
         {
            const int edgeroot2term = graph_pc_getRoot2PtermEdge(graph, term);

            solnode[term] = TRUE;
            (*nsolnodes)++;

            assert(!soledge[edge2term / 2]);
            assert(!soledge[edgeroot2term / 2]);

            soledge[edge2term / 2] = TRUE;
            soledge[edgeroot2term / 2] = TRUE;
            (*nsoledges) += 2;
         }
      }
   }

   assert(*nsolnodes >= 1);
   assert(*nsoledges >= 0);
}


/** adapt merged solution graph for DCSTP */
static
void solgraphAdaptForDCSTP(
   const GRAPH*          graph,              /**< graph structure */
   const STP_Bool*       solnode,            /**< solution nodes array */
   STP_Bool*             soledge,            /**< solution edges array  */
   int*                  nsoledges           /**< number of edges pointer */
   )
{
   const int nnodes = graph->knots;


   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) )
      {
         assert(solnode[k]);
         for( int i = graph->outbeg[k]; i != EAT_LAST; i = graph->oeat[i] )
         {
            if( solnode[graph->head[i]] && !soledge[i / 2] )
            {
               soledge[i / 2] = TRUE;
               (*nsoledges)++;
            }
         }
      }
   }
}



/** add additional edges to solution graph */
static
void solgraphAddEdges(
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   const STP_Bool*       solnode,            /**< solution nodes array */
   STP_Bool*             soledge,            /**< solution edges array  */
   int*                  nsoledges           /**< number of edges pointer */
   )
{
   const int nedges = graph->edges;

   int naddedges = 0;
   const int offset = SCIPrandomGetInt(heurdata->randnumgen, 0, nedges - 1);

   for( int e = 0; e < nedges && naddedges <= (int)(REC_ADDEDGE_FACTOR * (*nsoledges)); e += 2 )
   {
      const int newedge = (e + offset) % nedges;
      int tail;
      int head;

      if( soledge[newedge / 2] )
         continue;

      // todo if fixed to zero, continue?
      if( graph->oeat[newedge] == EAT_FREE )
         continue;

      tail = graph->tail[newedge];
      head = graph->head[newedge];

      if( solnode[tail] && solnode[head] )
      {
         soledge[newedge / 2] = TRUE;
         naddedges++;
      }
   }
   SCIPdebugMessage("additional edges %d (orig: %d )\n", naddedges, *nsoledges);

   (*nsoledges) += naddedges;
}


/** merge selected solutions to a new graph */
static
SCIP_RETCODE buildsolgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const STPSOLPOOL*     pool,               /**< solution pool or NULL */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   GRAPH**               solgraph,           /**< pointer to store new graph */
   const int*            incumbentedges,     /**< edges of solution to be used as recombination root */
   int*                  incumbentindex,     /**< index of ancestor of incumbent solution */
   int**                 edgeancestor,       /**< ancestor to edge edge */
   int**                 edgeweight,         /**< for each edge: number of solution that contain this edge */
   SCIP_Bool*            success,            /**< new graph constructed? */
   SCIP_Bool             randomize,          /**< select solution randomly? */
   SCIP_Bool             addedges            /**< add additional edges between solution vertices? */
   )
{
   GRAPH* newgraph = NULL;
   SCIP_SOL** sols = NULL;
   STPSOL** poolsols = NULL;
   SCIP_VAR** vars = NULL;
   int* solselection;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   const SCIP_Bool usestppool = (pool != NULL);

   assert(scip != NULL);
   assert(graph != NULL);

   if( !usestppool )
   {
      sols = SCIPgetSols(scip);
      vars = SCIPprobdataGetEdgeVars(scip);
      assert(vars != NULL);
      assert(sols != NULL);
   }
   else
   {
      poolsols = pool->sols;
   }

   *success = TRUE;
   *edgeweight = NULL;
   *edgeancestor = NULL;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselection, heurdata->nusedsols) );

   /* select solutions to be merged */
   if( pcmw || graph->stp_type == STP_DCSTP )
      SCIP_CALL( solgraphSelectSolsDiff(scip, pool, graph, heurdata, vars, incumbentedges, incumbentindex, solselection, success) );
   else
      SCIP_CALL( solgraphSelectSols(scip, pool, heurdata, incumbentindex, solselection, randomize) );

   if( *success )
   {
      int* dnodemap;
      STP_Bool* solnode;               /* marks nodes contained in at least one solution */
      STP_Bool* soledge;               /* marks edges contained in at least one solution */
      int j;
      int nsoledges = 0;
      int nsolnodes = 0;
      const int nedges = graph->edges;
      const int nnodes = graph->knots;
      const int selectedsols = heurdata->nselectedsols;

      SCIPdebugMessage("REC: successfully selected new solutions \n");

      assert(selectedsols > 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &dnodemap, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &soledge, nedges / 2) );

      for( int i = 0; i < nnodes; i++ )
      {
         solnode[i] = FALSE;
         dnodemap[i] = UNKNOWN;
      }

      /* count and mark selected nodes and edges */
      for( int i = 0; i < nedges; i += 2 )
      {
         const int ihalf = i / 2;
         soledge[ihalf] = FALSE;

         if( graph->oeat[i] == EAT_FREE )
            continue;

         for( j = 0; j < selectedsols; j++ )
         {
            SCIP_Bool hit;

            if( j == 0 )
            {
               hit = (incumbentedges[i] == CONNECT
                   || incumbentedges[i + 1] == CONNECT);
            }
            else if( usestppool )
               hit = (poolsols[solselection[j]]->soledges[i] == CONNECT
                   || poolsols[solselection[j]]->soledges[i + 1] == CONNECT);
            else
               hit = (SCIPisEQ(scip, SCIPgetSolVal(scip, sols[solselection[j]], vars[i]), 1.0)
                   || SCIPisEQ(scip, SCIPgetSolVal(scip, sols[solselection[j]], vars[i + 1]), 1.0));

            if( hit )
            {
               nsoledges++;
               soledge[ihalf] = TRUE;
               if( !solnode[graph->tail[i]] )
               {
                  solnode[graph->tail[i]] = TRUE;
                  nsolnodes++;
               }
               if( !solnode[graph->head[i]] )
               {
                  solnode[graph->head[i]] = TRUE;
                  nsolnodes++;
               }
               break;
            }
         }
      }

      /* need to adapt solution graph for PC/MW? */
      if( pcmw )
      {
         solgraphAdaptForPcMw(graph, solnode, &nsolnodes, soledge, &nsoledges);
      }

      /* add additional edges? */
      if( addedges )
      {
         solgraphAddEdges(heurdata, graph, solnode, soledge, &nsoledges);
      }

      if( graph->stp_type == STP_GSTP )
      {
         solgraphAdaptForDCSTP(graph, solnode, soledge, &nsoledges);
      }

      /* initialize new graph */
      SCIP_CALL( graph_init(scip, &newgraph, nsolnodes, 2 * nsoledges, 1) );

      if( graph->stp_type == STP_RSMT || graph->stp_type == STP_OARSMT || graph->stp_type == STP_GSTP )
         newgraph->stp_type = STP_SPG;
      else
         newgraph->stp_type = graph->stp_type;

      if( pcmw )
      {
         SCIP_CALL( graph_pc_initSubgraph(scip, newgraph) );
      }

      newgraph->hoplimit = graph->hoplimit;
      j = 0;
      for( int i = 0; i < nnodes; i++ )
      {
         if( solnode[i] )
         {
            if( pcmw )
            {
               assert(graph->extended);

               newgraph->prize[j] = graph->prize[i];
#ifndef NDEBUG
               if( graph_pc_knotIsDummyTerm(graph, i) && i != graph->source )
                  assert(0.0 == newgraph->prize[j]);
#endif
            }

            graph_knot_add(newgraph, graph->term[i]);
            dnodemap[i] = j++;
         }
      }

      if( pcmw )
      {
         newgraph->extended = TRUE;
         newgraph->norgmodelknots = newgraph->knots - newgraph->terms;

         if( graph_pc_isRootedPcMw(graph) )
            newgraph->norgmodelknots += graph_pc_nFixedTerms(graph);
      }

      SCIPdebugMessage("REC: sol graph with nodes: %d, edges: %d, terminals: %d  \n", nsolnodes, 2 * nsoledges, newgraph->terms);

      /* set root */
      newgraph->source = dnodemap[graph->source];

      assert(newgraph->source >= 0);
      assert(!graph_pc_isRootedPcMw(graph) || newgraph->prize[newgraph->source] == FARAWAY);

      /* copy max degrees*/
      if( graph->stp_type == STP_DCSTP )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(newgraph->maxdeg), nsolnodes) );
         for( int i = 0; i < nnodes; i++ )
            if( solnode[i] )
               newgraph->maxdeg[dnodemap[i]] = graph->maxdeg[i];
      }
      /* allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, edgeancestor, 2 * nsoledges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, edgeweight, 2 * nsoledges) );

      for( int i = 0; i < 2 * nsoledges; i++ )
         (*edgeweight)[i] = 1;

      /* store original ID of each new edge (i.e. edge in the merged graph) */
      j = 0;
      assert(selectedsols == heurdata->nselectedsols);
      for( int i = 0; i < nedges; i += 2 )
      {
         if( soledge[i / 2] )
         {
            (*edgeancestor)[j++] = i;
            (*edgeancestor)[j++] = i + 1;

            graph_edge_addSubgraph(scip, graph, dnodemap, i, newgraph);

            /* (*edgeweight)[e]: number of solutions containing edge e */
            for( int k = 0; k < selectedsols; k++ )
            {
               SCIP_Bool hit;

               if( k == 0 )
               {
                  hit = (incumbentedges[i] == CONNECT
                      || incumbentedges[i + 1] == CONNECT);
               }
               else if( usestppool )
                  hit = (poolsols[solselection[k]]->soledges[i] == CONNECT
                      || poolsols[solselection[k]]->soledges[i + 1] == CONNECT);
               else
                  hit = (SCIPisEQ(scip, SCIPgetSolVal(scip, sols[solselection[k]], vars[i]), 1.0)
                      || SCIPisEQ(scip, SCIPgetSolVal(scip, sols[solselection[k]], vars[i + 1]), 1.0));

               /* is edge i or reversed edge in current solution? */
               if( hit )
               {
                  (*edgeweight)[j - 2]++;
                  (*edgeweight)[j - 1]++;
               }
            }
         }
      }

      assert(j == 2 * nsoledges);

      SCIP_CALL( graph_pc_finalizeSubgraph(scip, newgraph) );

      SCIPfreeBufferArray(scip, &soledge);
      SCIPfreeBufferArray(scip, &dnodemap);
      SCIPfreeBufferArray(scip, &solnode);

      if( !graph_valid(scip, newgraph) )
      {
         assert(usestppool);

         for( int k = 1; k < selectedsols; k++ )
         {
            assert(!poolSolIsUnreduced(graph, solselection[k], pool));
         }

         graph_free(scip, &newgraph, TRUE);

         *success = FALSE;
      }

      assert(!pcmw || graph_pc_nFixedTerms(graph) == graph_pc_nFixedTerms(newgraph));
   }

   SCIPfreeBufferArray(scip, &solselection);

   assert(*success || newgraph == NULL);

   *solgraph = newgraph;

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */


/** runs STP recombination heuristic */
SCIP_RETCODE SCIPStpHeurRecRun(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL*           pool,               /**< STP solution pool or NULL */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data or NULL */
   const GRAPH*          graph,              /**< graph data */
   SCIP_VAR**            vars,               /**< variables or NULL */
   int*                  newsolindex,        /**< index of new solution */
   int                   runs,               /**< number of runs */
   int                   nsols,              /**< number of solutions in pool (SCIP or STP) */
   SCIP_Bool             restrictheur,       /**< use restricted version of heuristic? */
   SCIP_Bool*            solfound            /**< new solution found? */
)
{
   int* newsoledges;
   int* incumbentedges;
   SCIP_Real incumentobj;

   const STP_Bool usestppool = (pool != NULL);
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const int probtype = graph->stp_type;
   STP_Bool* stnodes;

   assert(runs >= 0);
   assert(graph != NULL);
   assert(scip != NULL);
   assert(*newsolindex >= 0);

   SCIPdebugMessage("REC: start \n");

   *solfound = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &newsoledges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &incumbentedges, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, nnodes) );

   /* initialize solution vector (edges) and objective of incumbent */
   initializeIncumbent(scip, pool, graph, vars, nsols, *newsolindex, &incumentobj, incumbentedges);
   SCIPdebugMessage("REC: incumbent obj: %f \n", incumentobj);

   if( heur == NULL )
   {
      heur = SCIPfindHeur(scip, "rec");
      assert(heur != NULL);
      heurdata = SCIPheurGetData(heur);
      assert(heurdata != NULL);
   }

   /* main loop (for recombination) */
   for( int v = 0, failcount = 0; v < CYCLES_PER_RUN * runs && !SCIPisStopped(scip); v++ )
   {
      GRAPH* solgraph = NULL;
      int* edgeweight;
      int* edgeancestor;
      SCIP_Bool success;
      const SCIP_Bool randomize = (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 1);
      const SCIP_Bool addmoreedges = (SCIPrandomGetInt(heurdata->randnumgen, 0, 2) == 1);

      /* compute number of solutions to merge (and save in heurdata->nsols) */
      setHeurdataNUsedSols(heurdata, nsols, restrictheur);
      SCIPdebugMessage("REC: merge %d solutions \n", heurdata->nusedsols);

      /* build a new graph, consisting of several solutions */
      SCIP_CALL( buildsolgraph(scip, pool, heurdata, graph, &solgraph, incumbentedges, newsolindex,
            &edgeancestor, &edgeweight, &success, randomize, addmoreedges) );

      /* valid graph built? */
      if( success )
      {
         IDX** ancestors;
         GRAPH* psolgraph;
         int* soledges = NULL;
         SCIP_Real pobj;
         SCIP_Real offsetdummy;

         SCIPdebugMessage("REC: solution successfully built \n");
         assert(graph_valid(scip, solgraph));

         /* reduce new graph */
         if( probtype == STP_DHCSTP || probtype == STP_DCSTP || probtype == STP_NWSPG || probtype == STP_SAP )
            SCIP_CALL( reduce(scip, solgraph, &offsetdummy, 0, 5, FALSE) );
         else
         {
            SCIP_CALL( reduce(scip, solgraph, &offsetdummy, 2, 5, FALSE) );
         }

         SCIP_CALL( graph_pack(scip, solgraph, &psolgraph, &offsetdummy, FALSE) );

         solgraph = psolgraph;
         ancestors = solgraph->ancestors;

         /* if graph reduction solved the whole problem, solgraph has only one terminal */
         if( solgraph->terms > 1 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &soledges, solgraph->edges) );

            SCIP_CALL( computeReducedProbSolution(scip, heurdata, graph, solgraph, ancestors, edgeweight,
                  edgeancestor, vars, usestppool, soledges) );
         }

         /*  retransform solution (soledges or single vertex) found above */
         SCIP_CALL( retransformReducedProbSolution(scip, graph, solgraph, ancestors, edgeancestor, soledges, newsoledges, stnodes) );

         SCIPfreeBufferArrayNull(scip, &soledges);
         SCIPfreeMemoryArray(scip, &edgeancestor);

         graph_free(scip, &solgraph, TRUE);

         /* prune solution (in the original graph) */
         if( probtype == STP_DCSTP )
            SCIP_CALL( SCIPStpHeurTMBuildTreeDc(scip, graph, newsoledges, stnodes) );
         else
            SCIP_CALL( solstp_prune(scip, graph, newsoledges, stnodes) );

         assert(solstp_isValid(scip, graph, newsoledges) || SCIPisStopped(scip));
         pobj = solstp_getObjBounded(graph, newsoledges, 0.0, nedges);

         SCIPdebugMessage("REC: new obj: %f \n", pobj);

         /* new solution better than incumbent? */
         if( !SCIPisStopped(scip) && SCIPisLT(scip, pobj, incumentobj) )
         {
            SCIPdebugMessage("...is better! \n");
            incumentobj = pobj;
            BMScopyMemoryArray(incumbentedges, newsoledges, nedges);

            if( (!restrictheur && failcount >= 1) || !(*solfound) )
               failcount--;

            *solfound = TRUE;
         }
         else
         {
            failcount++;
         }
      }
      else // no new graph could be built
      {
         assert(solgraph == NULL);
         failcount++;
         SCIPdebugMessage("REC: no new graph could be build \n");
      }

      SCIPfreeMemoryArray(scip, &edgeweight);

      /* too many fails? */
      if( failcount >= REC_MAX_FAILS )
      {
         SCIPdebugMessage("REC: too many fails, exit! \n");
         break;
      }
   } /* main recombination loop */


   /* improving solution found? */
   if( *solfound )
   {
      SCIP_CALL( poolAddSol(scip, pool, heur, incumentobj, incumbentedges, nedges, newsolindex, solfound) );
   }
   else
   {
      *newsolindex = -1;
   }

   SCIPdebugMessage("incumbentobj=%f newsolobj=%f \n", solstp_getObjBounded(graph, incumbentedges, 0.0, nedges), solstp_getObjBounded(graph, newsoledges, 0.0, nedges));
   SCIPdebugMessage("incumentobj=%f \n", incumentobj);

   SCIPfreeBufferArray(scip, &stnodes);
   SCIPfreeBufferArray(scip, &incumbentedges);
   SCIPfreeBufferArray(scip, &newsoledges);

   return SCIP_OKAY;
}

/** heuristic to exclude vertices or edges from a given solution (and inserting other edges) to improve objective */
SCIP_RETCODE SCIPStpHeurRecExclude(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   const int*            result,             /**< edge solution array (UNKNOWN/CONNECT) */
   int*                  newresult,          /**< new edge solution array (UNKNOWN/CONNECT) */
   int*                  dnodemap,           /**< node array for internal use */
   STP_Bool*             stvertex,           /**< node array for internally marking solution vertices */
   SCIP_Bool*            success             /**< solution improved? */
   )
{
   GRAPH* newgraph;
   int* unodemap;
   STP_Bool* solnodes;
   SCIP_Real dummy;
   int j;
   int nsolterms;
   int nsoledges;
   int nsolnodes;
   const int root = graph->source;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(result != NULL);
   assert(stvertex != NULL);
   assert(success != NULL);
   assert(dnodemap != NULL);
   assert(newresult != NULL);

   *success = FALSE;
   assert(graph->stp_type == STP_MWCSP);
   assert(graph_valid(scip, graph));

   /* killed solution edge? */
   for( int e = 0; e < nedges; e++ )
      if( result[e] == CONNECT && graph->oeat[e] == EAT_FREE )
         return SCIP_OKAY;

   /*** 1. step: for solution S and original graph (V,E) initialize new graph (V[S], (V[S] X V[S]) \cup E) ***/

   BMSclearMemoryArray(stvertex, nnodes);

   nsolnodes = 1;
   nsolterms = 0;
   stvertex[root] = 1;

   /* mark nodes in solution */
   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         int tail = graph->tail[e];
         int head = graph->head[e];

         if( tail == root )
         {
            /* there might be only one node */
            if( Is_pseudoTerm(graph->term[head]) )
            {
               stvertex[head] = 1;
               nsolterms++;
               nsolnodes++;
            }

            continue;
         }

         assert(!stvertex[head] );
         stvertex[head] = 1;

         if( Is_pseudoTerm(graph->term[head]) )
            nsolterms++;
         nsolnodes++;
      }
   }

   assert(nsolterms > 0);

   nsoledges = 0;

   /* count edges of new graph */
   for( int i = 0; i < nedges; i += 2 )
      if( stvertex[graph->tail[i]] && stvertex[graph->head[i]] && graph->oeat[i] != EAT_FREE  )
         nsoledges++;

   nsoledges *= 2;

   /* create new graph */
   SCIP_CALL( graph_init(scip, &newgraph, nsolnodes, nsoledges, 1) );

   newgraph->stp_type = graph->stp_type;
   newgraph->extended = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &unodemap, nsolnodes) );
   SCIP_CALL( graph_pc_initSubgraph(scip, newgraph) );

   j = 0;
   for( int i = 0; i < nnodes; i++ )
   {
      if( stvertex[i] )
      {
         if( (!Is_term(graph->term[i])) )
            newgraph->prize[j] = graph->prize[i];
         else
         {
            assert(SCIPisZero(scip, graph->prize[i]));

            newgraph->prize[j] = 0.0;
         }

         graph_knot_add(newgraph, graph->term[i]);
         unodemap[j] = i;
         dnodemap[i] = j++;
      }
      else
      {
         dnodemap[i] = -1;
      }
   }

   assert(j == nsolnodes);

   /* set root */
   newgraph->source = dnodemap[root];
   assert(newgraph->source >= 0 && newgraph->source < nsolnodes);

   /* add edges */
   for( int i = 0; i < nedges; i += 2 )
      if( stvertex[graph->tail[i]] && stvertex[graph->head[i]] && graph->oeat[i] != EAT_FREE )
         graph_edge_addSubgraph(scip, graph, dnodemap, i, newgraph);

   assert(newgraph->edges == nsoledges);

   SCIP_CALL( graph_pc_finalizeSubgraph(scip, newgraph) );


   /*** step 2: presolve ***/

   newgraph->norgmodelknots = nsolnodes;

   dummy = 0.0;
   SCIP_CALL( reduce(scip, newgraph, &dummy, 1, 5, FALSE) );


   /*** step 3: compute solution on new graph ***/

   SCIP_CALL( graph_path_init(scip, newgraph) );

   /* compute Steiner tree to obtain upper bound */
   SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fromheurdata,
      newgraph, NULL, NULL, newresult, MIN(50, nsolterms), newgraph->source, newgraph->cost, newgraph->cost, &dummy, NULL, success) );

   graph_path_exit(scip, newgraph);

   assert(*success);
   assert(solstp_isValid(scip, newgraph, newresult));


   /*** step 4: retransform solution to original graph ***/


   BMSclearMemoryArray(stvertex, nnodes);

   for( int e = 0; e < nsoledges; e++ )
      if( newresult[e] == CONNECT )
         marksolverts(newgraph, newgraph->ancestors[e], unodemap, stvertex);

   /* retransform edges fixed during graph reduction */
   marksolverts(newgraph, graph_get_fixedges(newgraph), unodemap, stvertex);

   for( int k = 0; k < nsolnodes; k++ )
      if( stvertex[unodemap[k]] )
         marksolverts(newgraph, newgraph->pcancestors[k], unodemap, stvertex);

   SCIP_CALL(SCIPallocBufferArray(scip, &solnodes, nsolnodes));

   for( int k = 0; k < nsolnodes; k++ )
      solnodes[k] = FALSE;

   for( int k = 0; k < nnodes; k++ )
      if( stvertex[k] )
      {
         assert(dnodemap[k] < nsolnodes);
         solnodes[dnodemap[k]] = TRUE;
      }

   SCIP_CALL( solstp_markPcancestors(scip, newgraph->pcancestors, newgraph->orgtail, newgraph->orghead, nsolnodes, solnodes, NULL, NULL, NULL, NULL ) );

   for( int k = 0; k < nsolnodes; k++ )
      if( solnodes[k] )
         stvertex[unodemap[k]] = TRUE;

   SCIPfreeBufferArray(scip, &solnodes);

   assert(graph_pc_isPcMw(graph));

   SCIP_CALL( solstp_prune(scip, graph, newresult, stvertex) );

   /* solution better than original one?  */

   if( SCIPisLT(scip, solstp_getObjBounded(graph, newresult, 0.0, nedges),
         solstp_getObjBounded(graph, result, 0.0, nedges)) )
   {
      *success = TRUE;
      SCIPdebugMessage("success %f < %f \n", solstp_getObjBounded(graph, newresult, 0.0, nedges), solstp_getObjBounded(graph, result, 0.0, nedges));
   }
   else
   {
      SCIPdebugMessage("no improvements %f >= %f \n", solstp_getObjBounded(graph, newresult, 0.0, nedges), solstp_getObjBounded(graph, result, 0.0, nedges));
      *success = FALSE;
   }

   assert(solstp_isValid(scip, graph, newresult));

   if( !solstp_isValid(scip, graph, newresult) )
      *success = FALSE;

   SCIPfreeBufferArray(scip, &unodemap);
   graph_free(scip, &newgraph, TRUE);

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRec)
{
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRec)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurRec(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRec)
{  /*lint --e{715}*/

   assert(heur != NULL);
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->nselectedsols = 0;
   heurdata->nusedsols = DEFAULT_NUSEDSOLS;
   heurdata->randseed = DEFAULT_RANDSEED;
   heurdata->nselectedsols = 0;
   heurdata->ncalls = 0;
   heurdata->nlastsols = 0;
   heurdata->lastsolindex = -1;
   heurdata->bestsolindex = -1;
   heurdata->nfailures = 0;


   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRec)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL** sols;
   GRAPH* graph;
   SCIP_Longint nallsols;
   SCIP_Bool pcmw;
   SCIP_Bool solfound;
   SCIP_Bool restrictheur;

   int runs;
   int nsols;
   int solindex;
   int probtype;
   int newsolindex;
   int nreadysols;
   int bestsolindex;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get graph */
   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   probtype = graph->stp_type;
   *result = SCIP_DIDNOTRUN;

   if( probtype == STP_BRMWCSP )
      return SCIP_OKAY;

   SCIPdebugMessage("REC: checking ... \n");

   pcmw = graph_pc_isPcMw(graph);
   nallsols = SCIPgetNSolsFound(scip);
   nreadysols = SCIPgetNSols(scip);

   /* only call heuristic if sufficiently many solutions are available */
   if( nreadysols < DEFAULT_NUSEDSOLS )
      return SCIP_OKAY;

   if( probtype == STP_DCSTP )
      return SCIP_OKAY;

   /* suspend heuristic? */
   if( pcmw || probtype == STP_DHCSTP || probtype == STP_DCSTP )
   {
      int i;
      if( heurdata->ncalls == 0 )
         i = 0;
      else if( heurdata->maxfreq )
         i = 1;
      else if( probtype == STP_DCSTP )
         i = MIN(2 * heurdata->nwaitingsols, 2 * heurdata->nfailures);
      else
         i = MIN(heurdata->nwaitingsols, heurdata->nfailures);

      if( nallsols <= heurdata->nlastsols + i )
         return SCIP_OKAY;
   }
   else
   {
      int i;
      if( heurdata->maxfreq )
         i = 1;
      else
         i = MIN(heurdata->nwaitingsols, heurdata->nfailures);

      if( nallsols <= heurdata->nlastsols + i && heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         return SCIP_OKAY;
   }

   /* get edge variables */
   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);
   assert(vars[0] != NULL);

   heurdata->ncalls++;

   restrictheur = (graph->terms > BOUND_MAXNTERMINALS && graph->edges > BOUND_MAXNEDGES);

   if( restrictheur )
      runs = RUNS_RESTRICTED;
   else
      runs = RUNS_NORMAL;

   if( runs > nreadysols )
      runs = nreadysols;

   assert(runs > 0);

   sols = SCIPgetSols(scip);
   assert(sols != NULL);

   bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* first run or exotic variant? */
   if( heurdata->lastsolindex == -1 || probtype == STP_DHCSTP || probtype == STP_DCSTP )
      newsolindex = bestsolindex;
   else
   {
      /* get best new solution */

      SCIP_Real minobj = FARAWAY;
      const int lastindex = heurdata->lastsolindex;

      newsolindex = -1;
      assert(lastindex >= 0);

      for( int i = 0; i < nreadysols; i++ )
      {
         const int currindex = SCIPsolGetIndex(sols[i]);
         if( currindex > lastindex && SCIPisLT(scip, SCIPgetSolOrigObj(scip, sols[i]), minobj) )
         {
            newsolindex = currindex;
            minobj = SCIPgetSolOrigObj(scip, sols[i]);
            SCIPdebugMessage("REC: better start solution, obj: %f \n", minobj);
         }
      }
      if( newsolindex == -1 )
      {
         SCIPdebugMessage("REC: random start solution\n");
         newsolindex = SCIPsolGetIndex(sols[SCIPrandomGetInt(heurdata->randnumgen, 0, nreadysols)]);
      }
   }

   /* run the actual heuristic */
   SCIP_CALL( SCIPStpHeurRecRun(scip, NULL, heur, heurdata, graph, vars, &newsolindex, runs, nreadysols, restrictheur, &solfound) );

   /* save latest solution index */
   solindex = 0;
   nsols = SCIPgetNSols(scip);
   assert(nsols > 0);
   sols = SCIPgetSols(scip);

   for( int i = 1; i < nsols; i++ )
      if( SCIPsolGetIndex(sols[i]) > SCIPsolGetIndex(sols[solindex]) )
         solindex = i;

   if( solfound )
      *result = SCIP_FOUNDSOL;

   if( SCIPsolGetIndex(SCIPgetBestSol(scip)) == bestsolindex )
      heurdata->nfailures++;
   else
      heurdata->nfailures = 0;

   heurdata->lastsolindex = SCIPsolGetIndex(sols[solindex]);
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));
   heurdata->nlastsols = SCIPgetNSolsFound(scip);

   return SCIP_OKAY;
}


/** creates the rec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurRec(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rec primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRec, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRec) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRec) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRec) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRec) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRec) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRec) );

   /* add rec primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nwaitingsols",
         "number of solution findings to be in abeyance",
         &heurdata->nwaitingsols, FALSE, DEFAULT_NWAITINGSOLS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/randseed",
         "random seed for heuristic",
         NULL, FALSE, DEFAULT_RANDSEED, 1, INT_MAX, paramChgdRandomseed, (SCIP_PARAMDATA*)heurdata) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnsols",
         "max size of solution pool for heuristic",
         &heurdata->maxnsols, FALSE, DEFAULT_MAXNSOLS, 5, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/ntmruns",
         "number of runs in TM",
         &heurdata->ntmruns, FALSE, DEFAULT_NTMRUNS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_MAXFREQREC, NULL, NULL) );

   heurdata->nusedsols = DEFAULT_NUSEDSOLS;
   heurdata->randseed = DEFAULT_RANDSEED;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, heurdata->randseed, TRUE) );

   return SCIP_OKAY;
}
