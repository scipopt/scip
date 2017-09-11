/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_rec.c
 * @brief  Primal recombination heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a recombination heuristic for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" (2016) by
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
#include "grph.h"
#include "heur_tm.h"
#include "cons_stp.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"
#include "math.h"

#define HEUR_NAME             "rec"
#define HEUR_DESC             "recombination heuristic for Steiner problems"
#define HEUR_DISPCHAR         'R'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_MAXFREQREC     FALSE         /**< executions of the rec heuristic at maximum frequency?                             */
#define DEFAULT_MAXNSOLS       12            /**< default maximum number of (good) solutions be regarded in the subproblem                  */
#define DEFAULT_NUSEDSOLS      4             /**< number of solutions that will be taken into account                               */
#define DEFAULT_RANDSEED       177           /**< random seed                                                                       */
#define DEFAULT_NTMRUNS        75            /**< number of runs in TM heuristic                                                    */
#define DEFAULT_NWAITINGSOLS   4             /**< max number of new solutions to be available before executing the heuristic again  */

#define BOUND_MAXNTERMINALS 1000
#define BOUND_MAXNEDGES     20000
#define RUNS_RESTRICTED     3
#define RUNS_NORMAL         10

#ifndef DEBUG
#define DEBUG
#endif

#ifdef WITH_UG
extern
int getUgRank(void);
#endif

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
   int                   randseed;           /**< seed value for random number generator                            */
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
SCIP_Real costMultiplier(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   SCIP_Real             avg                 /**< number of solutions containing this edge */
   )
{
   int factor = 1;
   int nusedsols = heurdata->nusedsols;
   SCIP_Real mult = 1;

   assert(SCIPisGE(scip, avg, 1.0));

   if( nusedsols <= 3 )
   {
      if( SCIPisLT(scip, avg, 1.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 1000, 1400);
      }
      else if( SCIPisLT(scip, avg, 2.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 200, 1000);
      }
      else if( nusedsols == 3 && SCIPisLT(scip, avg, 3.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 40, 100);
      }
   }
   else
   {
      if( SCIPisLT(scip, avg, 1.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 1400, 1800);
      }
      else if( SCIPisLT(scip, avg, 2.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 400, 1400);
      }
      else if( SCIPisLT(scip, avg, 3.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 150, 250);
      }
      else if( SCIPisLT(scip, avg, 4.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 60, 90);
      }
      else if( nusedsols >= 5 && SCIPisLT(scip, avg, 5.6) )
      {
         factor = SCIPrandomGetInt(heurdata->randnumgen, 20, 40);
      }

   }

   mult =  (double) factor * (1.0 / avg);


   return mult;
}

/** select solutions to be merged */
static
SCIP_RETCODE selectdiffsols(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_SOL**            newsol,             /**< new solution */
   int*                  selection,          /**< selected solutions */
   SCIP_Bool*            success             /**< could at least two solutions be selected? */
   )
{
   int* solidx;
   SCIP_SOL** sols;                          /* array of all solutions found so far         */
   SCIP_Real varsolval;
   SCIP_Real varrevsolval;
   int i;
   int e;
   int k;
   int head;
   int tail;
   int eqnedges;
   int diffnedges;
   int maxnsols;
   int nselectedsols;
   int nsols;                                /* number of all solutions found so far        */
   int nedges;
   int nusedsols;                            /* number of solutions to use in rec     */
   int* perm;
   int* solselected;
   STP_Bool* soledges;
   assert(selection != NULL);
   assert(graph != NULL);
   assert(vars != NULL);

   /* get solution data */
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   maxnsols = heurdata->maxnsols;
   nusedsols = heurdata->nusedsols;
   nedges = graph->edges;
   assert(nusedsols <= nsols);
   nselectedsols = 0;

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solidx, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledges, nedges / 2) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      solidx[i] = SCIPsolGetIndex(sols[i]);
      solselected[i] = FALSE;
   }

   if( *newsol == NULL )
   {
      SCIPsortIntInt(solidx, perm, nsols);
      i = nsols - 1;

      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         *newsol = sols[perm[i]];
      }
      else
      {
         i = SCIPrandomGetInt(heurdata->randnumgen, 0, nsols - 1);
      }
      *newsol = sols[perm[i]];
      solselected[perm[i]] = TRUE;
      selection[nselectedsols++] = perm[i];

      for( i = 0; i < nsols; i++ )
         perm[i] = i;
   }
   else
   {
      for( i = 0; i < nsols; i++ )
	 if( *newsol == sols[i] )
            break;
      assert(i < nsols);
      solselected[i] = TRUE;
      selection[nselectedsols++] = i;
   }

   for( e = 0; e < nedges; e += 2 )
   {
      varsolval = SCIPgetSolVal(scip, sols[selection[0]], vars[e]);
      varrevsolval = SCIPgetSolVal(scip, sols[selection[0]], vars[e + 1]);

      if( SCIPisEQ(scip, varsolval, 1.0) || SCIPisEQ(scip, varrevsolval, 1.0) )
      {
         soledges[e / 2] = TRUE;
      }
      else
      {
         soledges[e / 2] = FALSE;
      }
   }
   maxnsols = MIN(nsols, maxnsols);

   SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, maxnsols);
   for( i = 0; i < maxnsols; i++ )
   {
      if( solselected[perm[i]] == FALSE )
      {
         k = perm[i];
         eqnedges = 0;
         diffnedges = 0;
         for( e = 0; e < nedges; e += 2 )
         {
            varsolval = SCIPgetSolVal(scip, sols[k], vars[e]);
            varrevsolval = SCIPgetSolVal(scip, sols[k], vars[e + 1]);

            if( SCIPisEQ(scip, varsolval, 1.0) || SCIPisEQ(scip, varrevsolval, 1.0)  )
            {
               head = graph->head[e];
               tail = graph->tail[e];
               if( soledges[e / 2] == FALSE )
               {
                  soledges[e / 2] = TRUE;
                  if( !(Is_term(graph->term[tail]) && Is_term(graph->term[head])) )
                     diffnedges++;
               }
               else
               {
                  if( !(Is_term(graph->term[tail]) && Is_term(graph->term[head])) )
                     eqnedges++;
               }
            }
         }

         /* enough similarities and differences with new solution? */
         if( diffnedges > 3 && eqnedges > 0  )
         {
            selection[nselectedsols++] = k;
            solselected[k] = TRUE;
            *success = TRUE;
            if( nselectedsols >= nusedsols )
               break;
         }
      }
   }

   assert(nselectedsols <= nusedsols);

   heurdata->nselectedsols = nselectedsols;

   /* free memory */
   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &solidx);
   SCIPfreeBufferArray(scip, &solselected);

   return SCIP_OKAY;
}

/** select solutions to be merged */
static
SCIP_RETCODE selectsols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   SCIP_SOL**            newsol,             /**< new solution */
   int*                  selection,          /**< selected solutions */
   SCIP_Bool             randomize           /**< select solutions randomly? */
   )
{
   int* solidx;
   SCIP_SOL** sols;                          /* array of all solutions found so far         */
   SCIP_Longint nallsols;
   int i;
   int end;
   int maxnsols;
   int nselectedsols;
   int shift;
   int nsols;                                /* number of all available solutions        */
   int nusedsols;                            /* number of solutions to use in rec           */
   int sqrtnallsols;
   int* perm;
   int* solselected;

   assert(selection != NULL);

   /* get solution data */
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   nallsols = SCIPgetNSolsFound(scip);
   maxnsols = heurdata->maxnsols;
   nusedsols = heurdata->nusedsols;
   assert(nusedsols <= nsols);
   nselectedsols = 0;

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solidx, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      solidx[i] = SCIPsolGetIndex(sols[i]);
      solselected[i] = FALSE;
   }

   if( *newsol == NULL )
   {
      SCIPsortIntInt(solidx, perm, nsols);
      i = nsols - 1;

      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         int j = i - 1;

         /* get best new solution */
         while( j >= 0 && heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[j]]) )
         {
            if( SCIPisLT(scip, SCIPgetSolOrigObj(scip, sols[perm[j]]), SCIPgetSolOrigObj(scip, sols[perm[i]])) )
               i = j;
            j--;
         }

         *newsol = sols[perm[i]];
      }
      else
      {
         i = SCIPrandomGetInt(heurdata->randnumgen, 0, nsols - 1);
      }

      *newsol = sols[perm[i]];
      solselected[perm[i]] = TRUE;
      selection[nselectedsols++] = perm[i];

      for( i = 0; i < nsols; i++ )
         perm[i] = i;
   }
   else
   {
      for( i = 0; i < nsols; i++ )
         if( *newsol == sols[i] )
            break;

      assert(i < nsols);

      if( i >= nsols )
         i = 0;

      solselected[i] = TRUE;
      selection[nselectedsols++] = i;
   }

   if( !randomize )
   {
      end = (int) (SCIPrandomGetReal(heurdata->randnumgen, 1.0, (SCIP_Real) nusedsols - 1));

      shift = SCIPrandomGetInt(heurdata->randnumgen, end, 2 * nusedsols - 1);

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
   sqrtnallsols = (int) sqrt(nallsols / 4);

   if( sqrtnallsols > maxnsols && sqrtnallsols < nsols )
      maxnsols = sqrtnallsols;

   SCIPdebugMessage("maxnsols in rec %d \n", maxnsols);

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

   heurdata->nselectedsols = nselectedsols;
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &solidx);
   SCIPfreeBufferArray(scip, &solselected);

   return SCIP_OKAY;
}

/** merge selected solutions to a new graph */
static
SCIP_RETCODE buildsolgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   SCIP_SOL**            newsol,             /**< new solution */
   GRAPH*                graph,              /**< graph structure */
   GRAPH**               solgraph,           /**< pointer to store new graph */
   int**                 edgeancestor,       /**< ancestor to edge edge */
   int**                 edgeweight,         /**< fore each edge: number of solution that contain this edge */
   SCIP_Bool*            success,            /**< new graph constructed? */
   SCIP_Bool             randomize           /**< select solution randomly? */
   )
{
   GRAPH* newgraph;
   SCIP_SOL**  sols;
   SCIP_VAR** vars;
   SCIP_Real varsolval;
   SCIP_Real varrevsolval;

   int    i;
   int    j;
   int    k;
   int    ihalf;
   int    nedges;
   int    nnodes;
   int    nsoledges;
   int    nsolnodes;
   int*   dnodemap;
   int*   solselection;          /* pool of solutions rec will use */
   SCIP_Bool pcmw;
   STP_Bool*  solnode;               /* marks nodes contained in at least one solution */
   STP_Bool*  soledge;               /* marks edges contained in at least one solution */

   assert(scip != NULL);
   assert(graph != NULL);

   sols = SCIPgetSols(scip);
   nedges = graph->edges;
   nnodes = graph->knots;
   nsoledges = 0;
   nsolnodes = 0;
   *success = TRUE;
   *edgeweight = NULL;
   *edgeancestor = NULL;
   newgraph = NULL;
   pcmw = (graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RPCSPG);
   assert(sols != NULL);

   vars = SCIPprobdataGetEdgeVars(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselection, heurdata->nusedsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dnodemap, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledge, nedges / 2) );

   for( i = 0; i < nnodes; i++ )
   {
      solnode[i] = FALSE;
      dnodemap[i] = UNKNOWN;
   }

   /* select solutions to be merged */
   if( pcmw || graph->stp_type == STP_DCSTP )
      SCIP_CALL( selectdiffsols(scip, graph, heurdata, vars, newsol, solselection, success) );
   else
      SCIP_CALL( selectsols(scip, heurdata, newsol, solselection, randomize) );

   if( *success )
   {
      assert(heurdata->nselectedsols > 0);

      /* count and mark selected nodes and edges */
      for( i = 0; i < nedges; i += 2 )
      {
         ihalf = i / 2;
         soledge[ihalf] = FALSE;
         for( j = 0; j < heurdata->nselectedsols; j++ )
         {
            varsolval = SCIPgetSolVal(scip, sols[solselection[j]], vars[i]);
            varrevsolval = SCIPgetSolVal(scip, sols[solselection[j]], vars[i + 1]);

            if( SCIPisEQ(scip, varsolval, 1.0) || SCIPisEQ(scip, varrevsolval, 1.0) )
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
      if( pcmw ) /* todo this probably won't work for RMWCSP */
      {
         int e2;
         int head;
         int oldroot = graph->source[0];
         for( i = graph->outbeg[oldroot]; i != EAT_LAST; i = graph->oeat[i] )
         {
            if( Is_gterm(graph->term[graph->head[i]]) )
            {
               ihalf = i / 2;
               head = graph->head[i];
               if( soledge[ihalf] == FALSE )
               {
                  nsoledges++;
                  soledge[ihalf] = TRUE;
                  if( !solnode[head] && SCIPisEQ(scip, graph->cost[flipedge(i)], FARAWAY) )
                  {
                     solnode[head] = TRUE;
                     nsolnodes++;
                  }
                  assert(solnode[graph->head[i]]);
               }

               if( Is_pterm(graph->term[head]) )
               {
                  for( e2 = graph->outbeg[head]; e2 != EAT_LAST; e2 = graph->oeat[e2] )
                     if( Is_term(graph->term[graph->head[e2]]) &&  graph->head[e2] != oldroot )
                        break;

                  assert(e2 != EAT_LAST);

                  if( soledge[e2 / 2] == FALSE )
                  {
                     nsoledges++;
                     soledge[e2 / 2] = TRUE;
                  }
               }
               else
               {
                  assert(Is_term(graph->term[head]));
                  for( e2 = graph->outbeg[head]; e2 != EAT_LAST; e2 = graph->oeat[e2] )
                      if( Is_pterm(graph->term[graph->head[e2]]) && graph->head[e2] != oldroot )
                         break;

                   assert(e2 != EAT_LAST);

                   if( soledge[e2 / 2] == FALSE )
                   {
                      nsoledges++;
                      soledge[e2 / 2] = TRUE;
                   }
               }
            }
         }
      }

      if( graph->stp_type == STP_GSTP )
      {
         for( k = 0; k < nnodes; k++ )
         {
            if( Is_term(graph->term[k]) )
            {
               assert(solnode[k]);
               for( i = graph->outbeg[k]; i != EAT_LAST; i = graph->oeat[i] )
               {
                  if( solnode[graph->head[i]] && !soledge[i / 2] )
                  {
                     soledge[i / 2] = TRUE;
                     nsoledges++;
                  }
               }
            }
         }
      }
      /* initialize new graph */
      SCIP_CALL( graph_init(scip, &newgraph, nsolnodes, 2 * nsoledges, 1, 0) );

      if( graph->stp_type == STP_RSMT || graph->stp_type == STP_OARSMT || graph->stp_type == STP_GSTP )
         newgraph->stp_type = STP_SPG;
      else
         newgraph->stp_type = graph->stp_type;

      if( pcmw )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(newgraph->prize), nsolnodes) );
      }

      newgraph->hoplimit = graph->hoplimit;
      j = 0;
      for( i = 0; i < nnodes; i++ )
      {
         if( solnode[i] )
         {
            if( pcmw )
            {
               if( (!Is_term(graph->term[i])) )
                  newgraph->prize[j] = graph->prize[i];
               else
                  newgraph->prize[j] = 0.0;
            }

            graph_knot_add(newgraph, graph->term[i]);
            dnodemap[i] = j++;
         }
      }

      if( pcmw )
         newgraph->norgmodelknots = newgraph->knots - newgraph->terms;

      /* set root */
      newgraph->source[0] = dnodemap[graph->source[0]];
      if( newgraph->stp_type == STP_RPCSPG )
         newgraph->prize[newgraph->source[0]] = FARAWAY;

      assert(newgraph->source[0] >= 0);

      /* copy max degrees*/
      if( graph->stp_type == STP_DCSTP )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(newgraph->maxdeg), nsolnodes) );
         for( i = 0; i < nnodes; i++ )
            if( solnode[i] )
               newgraph->maxdeg[dnodemap[i]] = graph->maxdeg[i];
      }
      /* allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, edgeancestor, 2 * nsoledges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, edgeweight, 2 * nsoledges) );

      for( i = 0; i < 2 * nsoledges; i++ )
         (*edgeweight)[i] = 1;

      /* store original ID of each new edge (i.e. edge in the merged graph) */
      j = 0;
      for( i = 0; i < nedges; i += 2 )
      {
         if( soledge[i / 2] )
         {
            (*edgeancestor)[j++] = i;
            (*edgeancestor)[j++] = i + 1;
            graph_edge_add(scip, newgraph, dnodemap[graph->tail[i]], dnodemap[graph->head[i]], graph->cost[i], graph->cost[i + 1]);

            /* (*edgeweight)[e]: number of solutions containing edge e */
            for( k = 0; k < heurdata->nselectedsols; k++ )
            {
               varsolval = SCIPgetSolVal(scip, sols[solselection[k]], vars[i]);
               varrevsolval = SCIPgetSolVal(scip, sols[solselection[k]], vars[i + 1]);
               if( SCIPisEQ(scip, varsolval, 1.0) ||  SCIPisEQ(scip, varrevsolval, 1.0) )
               {
                  (*edgeweight)[j - 2]++;
                  (*edgeweight)[j - 1]++;
               }
            }
         }
      }

      assert(j == 2 * nsoledges);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &soledge);
   SCIPfreeBufferArray(scip, &dnodemap);
   SCIPfreeBufferArray(scip, &solnode);
   SCIPfreeBufferArray(scip, &solselection);
   *solgraph = newgraph;
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

   /* free random number generator */
   SCIPrandomFree(&heurdata->randnumgen);

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

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->nselectedsols = 0;
   heurdata->ncalls = 0;
   heurdata->ntmruns = 100;
   heurdata->nlastsols = 0;
   heurdata->lastsolindex = -1;
   heurdata->bestsolindex = -1;
   heurdata->nfailures = 0;
   heurdata->nusedsols = DEFAULT_NUSEDSOLS;
   heurdata->randseed = DEFAULT_RANDSEED;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#endif

   /* create random number generator */
   SCIP_CALL( SCIPrandomCreate(&heurdata->randnumgen, SCIPblkmem(scip),
         SCIPinitializeRandomSeed(scip, heurdata->randseed)) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_HEURDATA* tmheurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL** sols;
   GRAPH* graph;
   GRAPH* solgraph;
   GRAPH* psolgraph;
   SCIP_SOL* sol;
   SCIP_SOL*  newsol;
   SCIP_SOL*  bestsol;
   SCIP_Real* nval;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real* nodepriority = NULL;
   SCIP_Real pobj;
   SCIP_Real avg;
   SCIP_Real maxcost = 0.0;
   SCIP_Real hopfactor = 0.1;
   SCIP_Longint nallsols;
   SCIP_Bool pcmw;
   SCIP_Bool bisuccess;
   SCIP_Bool success;
   SCIP_Bool fixed;
   SCIP_Bool randomize;
   SCIP_Bool solfound;
   SCIP_Bool restrictheur;
   IDX* curr;
   IDX** ancestors;
   int i;
   int e;
   int v;
   int runs;
   int nsols;
   int randn;
   int count;
   int nedges;
   int nnodes;
   int nrounds;
   int solindex;
   int probtype;
   int roundstop;
   int nsoledges;
   int nreadysols;
   int best_start;
   int bestsolindex;
   int* results;
   int* edgeweight;
   int* orgresults;
   int* edgeancestor;
   STP_Bool* stnodes;

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

   if( probtype == STP_RMWCSP )
      return SCIP_OKAY;

   pcmw = (probtype == STP_PCSPG || probtype == STP_MWCSP || probtype == STP_RPCSPG);
   nallsols = SCIPgetNSolsFound(scip);
   nreadysols = SCIPgetNSols(scip);

   /* only call heuristic if sufficiently many solutions are available */
   if( nreadysols < DEFAULT_NUSEDSOLS )
      return SCIP_OKAY;

   /* suspend heuristic? */
   if( pcmw || probtype == STP_DHCSTP || probtype == STP_DCSTP )
   {
      if( heurdata->ncalls == 0 )
         i = 0;
      else if( heurdata->maxfreq )
         i = 1;
      else if( probtype == STP_RPCSPG || probtype == STP_DCSTP )
         i = MIN(2 * heurdata->nwaitingsols, 2 * heurdata->nfailures);
      else
         i = MIN(heurdata->nwaitingsols, heurdata->nfailures);

      if( nallsols <= heurdata->nlastsols + i )
         return SCIP_OKAY;
   }
   else
   {
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

   nedges = graph->edges;
   nnodes = graph->knots;
   bestsol = SCIPgetBestSol(scip);
   results = NULL;

   heurdata->ncalls++;

   restrictheur = (graph->terms > BOUND_MAXNTERMINALS && nedges > BOUND_MAXNEDGES);

   if( restrictheur)
      runs = RUNS_RESTRICTED;
   else
      runs = RUNS_NORMAL;

   if( runs > nreadysols )
      runs = nreadysols;

   assert(runs > 0);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgresults, nedges) );

   if( probtype == STP_MWCSP || probtype == STP_DHCSTP || probtype == STP_DCSTP )
   {
      newsol = (SCIPgetSols(scip))[0];
   }
   /* first run? */
   else if( heurdata->lastsolindex == -1 )
   {
      newsol = (SCIPgetSols(scip))[SCIPrandomGetInt(heurdata->randnumgen, 0, heurdata->nusedsols - 1)];
   }
   else
   {
      newsol = NULL;
   }

   bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   /* save latest solution index */
   solindex = 0;
   nsols = SCIPgetNSols(scip);
   assert(nsols > 0);
   sols = SCIPgetSols(scip);

   for( i = 1; i < nsols; i++ )
      if( SCIPsolGetIndex(sols[i]) > SCIPsolGetIndex(sols[solindex]) )
         solindex = i;

   if( SCIPsolGetIndex(SCIPgetBestSol(scip)) == bestsolindex )
      heurdata->nfailures++;
   else
      heurdata->nfailures = 0;

   heurdata->lastsolindex = SCIPsolGetIndex(sols[solindex]);
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));
   heurdata->nlastsols = SCIPgetNSolsFound(scip);
   SCIPfreeBufferArray(scip, &orgresults);
   SCIPfreeBufferArray(scip, &nval);

   return SCIP_OKAY;
}


static
inline void markSolVerts(
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
SCIP_Bool isInPool(
   const int*            soledges,           /**< edge array of solution to be checked */
   const STPSOLPOOL*     pool                /**< the pool */
)
{
   STPSOL** poolsols = pool->sols;
   const int poolsize = pool->size;
   const int nedges = pool->nedges;

   for( int i = 0; i < poolsize; i++ )
   {
      int j;
      const int* pooledges = poolsols[i]->soledges;
      assert(pooledges != NULL);

      for( j = 0; j < nedges; j++ )
         if( pooledges[j] != soledges[j] )
            break;

      /* pooledges == soledges? */
      if( j == nedges )
         return TRUE;
   }
   return FALSE;
}

/*
 * primal heuristic specific interface methods
 */

/** initializes STPSOL pool */
SCIP_RETCODE SCIPStpHeurRecInitPool(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool,               /**< the pool */
   const int             maxsize             /**< capacity of pool */
   )
{
   STPSOLPOOL* dpool;

   assert(pool != NULL);
   assert(maxsize > 0);

   SCIP_CALL( SCIPallocBuffer(scip, pool) );

   dpool = *pool;

   SCIP_CALL( SCIPallocBufferArray(scip, &(dpool->sols), maxsize) );

   for( int i = 0; i < maxsize; i++ )
      dpool->sols[i] = NULL;

   dpool->maxsize = maxsize;
   dpool->size = 0;

   return SCIP_OKAY;
}


/** frees STPSOL pool */
void SCIPStpHeurRecFreePool(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool                /**< the pool */
   )
{
   STPSOLPOOL* dpool = *pool;
   const int poolsize = dpool->size;

   assert(pool != NULL);
   assert(dpool != NULL);
   assert(poolsize == dpool->maxsize || dpool->sols[poolsize] == NULL);

   for( int i = poolsize - 1; i >= 0; i-- )
   {
      STPSOL* sol = dpool->sols[i];

      assert(sol != NULL);

      SCIPfreeBufferArray(scip, &(sol->soledges));
      SCIPfreeBuffer(scip, &sol);
   }

   SCIPfreeBuffer(scip, pool);
}


/** tries to add STPSOL to pool */
SCIP_RETCODE SCIPStpHeurRecAddToPool(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real       obj,                /**< objective of solution to be added */
   const int*            soledges,           /**< edge array of solution to be added */
   STPSOLPOOL*           pool,               /**< the pool */
   SCIP_Bool*            success             /**< has solution been added? */
   )
{
   STPSOL** poolsols = pool->sols;
   STPSOL* sol;
   int i;
   int poolsize = pool->size;
   const int nedges = pool->nedges;
   const int poolmaxsize = pool->maxsize;

   assert(scip != NULL);
   assert(pool != NULL);
   assert(poolsols != NULL);
   assert(poolsize >= 0);
   assert(poolmaxsize >= 0);
   assert(poolsize <= poolmaxsize);

   *success = FALSE;

   /* is solution in pool? */
   if( !isInPool(soledges, pool) )
      return SCIP_OKAY;

   /* enlarge pool if possible */
   if( poolsize < poolmaxsize )
   {
      sol = pool->sols[poolsize];
      SCIP_CALL( SCIPallocBuffer(scip, &sol) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(sol->soledges), nedges) );

      poolsize++;
      pool->size++;
   }
   /* pool is full; new solution worse than worst solution in pool? */
   else if( SCIPisGT(scip, obj, poolsols[poolsize - 1]->obj) )
   {
      return SCIP_OKAY;
   }

   /* overwrite last element of pool (either empty or inferior to current solution) */
   sol = poolsols[poolsize - 1];
   assert(sol != NULL);
   sol->obj = obj;
   BMScopyMemoryArray(sol->soledges, soledges, nedges);

   /* shift solution up */
   for( i = poolsize - 1; i >= 1; i-- )
   {
      if( SCIPisGT(scip, obj, poolsols[i - 1]->obj) )
         break;

      poolsols[i] = poolsols[i - 1];
   }

   poolsols[i] = sol;
   *success = TRUE;

   return SCIP_OKAY;
}
#if 0
/** runs STP recombination heuristic */
SCIP_RETCODE SCIPStpHeurRecRun(


)
{
   count = 0;
   nsols = nreadysols;
   solfound = FALSE;

   for( v = 0; v < 2 * runs && !SCIPisStopped(scip); v++ )
   {
      if( SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 1 )
         randomize = TRUE;
      else
         randomize = FALSE;

      /* one cycle finished? */
      if( count++ >= runs )
      {
         if( solfound )
         {
            count = 0;
            solfound = FALSE;
         }
         else
         {
            break;
         }
      }

      if( restrictheur )
         randn = SCIPrandomGetInt(heurdata->randnumgen, 0, 3);
      else
         randn = SCIPrandomGetInt(heurdata->randnumgen, 0, 5);

      if( (randn <= 2) || (nsols < 3) )
         heurdata->nusedsols = 2;
      else if( (randn <= 4) || (nsols < 4) )
         heurdata->nusedsols = 3;
      else
         heurdata->nusedsols = 4;

      /* build up a new graph, consisting of several solutions */
      SCIP_CALL( buildsolgraph(scip, heurdata, &newsol, graph, &solgraph, &edgeancestor, &edgeweight, &success, randomize) );

      /* valid graph built? */
      if( success )
      {
         assert(newsol != NULL);

         /* get TM heuristic data */
         assert(SCIPfindHeur(scip, "TM") != NULL);
         tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

         /* presolve new graph */
         assert(graph_valid(solgraph));

         /* reduce new graph */
         if( probtype == STP_RPCSPG || probtype == STP_DHCSTP || probtype == STP_DCSTP || probtype == STP_MWCSP
             || probtype == STP_NWSPG || probtype == STP_SAP )
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 0, 5) );
         else
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 2, 5) );

         SCIP_CALL( graph_pack(scip, solgraph, &psolgraph, FALSE) );

         solgraph = psolgraph;
         ancestors = solgraph->ancestors;
         nsoledges = solgraph->edges;

         /* if graph reduction solved the whole problem, solgraph has only one node */
         if( solgraph->terms > 1 )
         {
            /* edge multiplier */
            SCIP_Real mult;

            /* allocate memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &results, nsoledges) );
            SCIP_CALL( SCIPallocBufferArray(scip, &cost, nsoledges) );
            SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nsoledges) );
            SCIP_CALL( SCIPallocBufferArray(scip, &nodepriority, solgraph->knots) );

            for( i = 0; i < solgraph->knots; i++ )
               nodepriority[i] = 0.0;

            /* copy edge costs */
            BMScopyMemoryArray(cost, solgraph->cost, nsoledges);

            maxcost = 0.0;
            for( e = 0; e < nsoledges; e++ )
            {
               curr = ancestors[e];
               avg = 0.0;
               i = 0;
               fixed = FALSE;
               if( curr != NULL )
               {
                  while( curr != NULL )
                  {
                     i++;
                     avg += edgeweight[curr->index];
                     if( SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
                        fixed = TRUE;

                     curr = curr->parent;
                  }
                  avg = (double) avg / (double) i;
                  assert(avg >= 1);
               }
               /* is an ancestor edge fixed? */
               if( fixed )
               {
                  cost[e] = BLOCKED;
               }
               else
               {
                  nodepriority[solgraph->head[e]] += avg - 1.0;
                  nodepriority[solgraph->tail[e]] += avg - 1.0;
                  mult = costMultiplier(scip, heurdata, avg);
                  cost[e] = cost[e] * mult;
               }

               if( probtype == STP_DHCSTP && SCIPisLT(scip, cost[e], BLOCKED) && SCIPisGT(scip, cost[e], maxcost) )
                  maxcost = cost[e];
            }

            for( e = 0; e < nsoledges; e++ )
            {
               costrev[e] = cost[flipedge(e)];
               results[e] = UNKNOWN;
            }

            /* initialize shortest path algorithm */
            SCIP_CALL( graph_path_init(scip, solgraph) );

            /* run TM heuristic */
            SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, solgraph, NULL, &best_start, results, heurdata->ntmruns,
                  solgraph->source[0], cost, costrev, &hopfactor, nodepriority, maxcost, &success, FALSE) );

            if( !success )
            {
               SCIPdebugMessage("failed to build tree in rec heuristic \n");

               graph_path_exit(scip, solgraph);
               SCIPfreeBufferArrayNull(scip, &nodepriority);
               SCIPfreeBufferArrayNull(scip, &costrev);
               SCIPfreeBufferArrayNull(scip, &cost);
               SCIPfreeBufferArrayNull(scip, &results);
               SCIPfreeMemoryArray(scip, &edgeweight);
               SCIPfreeMemoryArray(scip, &edgeancestor);
               graph_free(scip, solgraph, TRUE);
               continue;
            }

            assert(graph_valid(solgraph));
            assert(graph_sol_valid(scip, solgraph, results));

            /* run local heuristic */
            if( probtype != STP_DHCSTP && probtype != STP_DCSTP && probtype != STP_MWCSP
                  && probtype != STP_SAP && probtype != STP_NWSPG )
               SCIP_CALL( SCIPStpHeurLocalRun(scip, solgraph, cost, costrev, results) );

            if( probtype == STP_SPG )
               assert(graph_sol_valid(scip, solgraph, results));
         }

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, nnodes) );

         bisuccess = FALSE;

#if 0 /* too expensive */
         roundstop = SCIPrandomGetInt(heurdata->randnumgen, 1, 2);
#endif

         roundstop = 1;

         for( nrounds = 0; nrounds < roundstop; nrounds++ )
         {
            for( i = 0; i < nedges; i++ )
               orgresults[i] = UNKNOWN;

            if( nrounds == 1 )
            {
               assert(results != NULL);
               SCIP_CALL( SCIPStpHeurSlackPruneRun(scip, NULL, solgraph, results, &success, FALSE) );

               assert(graph_sol_valid(scip, solgraph, results));

            }
            if( (probtype != STP_SPG || nrounds == roundstop - 1) && (solgraph->terms > 1) )
            {
               graph_path_exit(scip, solgraph);
            }

            for( i = 0; i < nnodes; i++ )
               stnodes[i] = FALSE;

            /* retransform solution found by heuristic */
            if( solgraph->terms > 1 )
            {
               assert(results != NULL);
               for( e = 0; e < nsoledges; e++ )
               {
                  if( results[e] == CONNECT )
                  {
                     /* iterate through list of ancestors */
                     if( probtype != STP_DCSTP )
                     {
                        curr = ancestors[e];

                        while( curr != NULL )
                        {
                           i = edgeancestor[curr->index];

                           stnodes[graph->head[i]] = TRUE;
                           stnodes[graph->tail[i]] = TRUE;

                           curr = curr->parent;
                        }
                     }
                     else
                     {
                        curr = ancestors[e];
                        while( curr != NULL )
                        {
                           i = edgeancestor[curr->index];
                           orgresults[i] = CONNECT;

                           curr = curr->parent;
                        }
                     }
                  }
               }
            }

            /* retransform edges fixed during graph reduction */
            if( probtype != STP_DCSTP )
            {
               curr = solgraph->fixedges;

               while( curr != NULL )
               {
                  i = edgeancestor[curr->index];

                  stnodes[graph->head[i]] = TRUE;
                  stnodes[graph->tail[i]] = TRUE;

                  curr = curr->parent;
               }
            }
            else
            {
               curr = solgraph->fixedges;
               while( curr != NULL )
               {
                  i = edgeancestor[curr->index];
                  orgresults[i] = CONNECT;
               }
            }

            if( pcmw )
            {
               for( i = 0; i < solgraph->knots; i++ )
               {
                  if( stnodes[i] == TRUE )
                  {
                     curr = solgraph->pcancestors[i];
                     while (curr != NULL)
                     {
                        if( stnodes[graph->tail[curr->index]] == FALSE )
                           stnodes[graph->tail[curr->index]] = TRUE;

                        if( stnodes[graph->head[curr->index]] == FALSE )
                           stnodes[graph->head[curr->index]] = TRUE;

                        curr = curr->parent;
                     }
                  }
               }
            }

            /* prune solution (in the original graph) */
            if( pcmw )
               SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, orgresults, stnodes) );
            else if( probtype == STP_DCSTP )
               SCIP_CALL( SCIPStpHeurTMBuildTreeDc(scip, graph, orgresults, stnodes) );
            else
               SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, orgresults, stnodes) );

            pobj = 0.0;

            for( e = 0; e < nedges; e++ )
            {
               if( orgresults[e] == CONNECT )
               {
                  nval[e] = 1.0;
                  pobj += graph->cost[e];
               }
               else
               {
                  nval[e] = 0.0;
               }
            }

            if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, newsol) - SCIPprobdataGetOffset(scip), pobj) )
            {
               SCIPdebugMessage("better solution found ...      ");
               sol = NULL;
               SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

               if( success )
               {
                  SCIPdebugMessage("and added! \n");
                  *result = SCIP_FOUNDSOL;
                  solfound = TRUE;
                  nsols = SCIPgetNSols(scip);

                  assert(nsols > 0);

                  sols = SCIPgetSols(scip);
                  solindex = 0;

                  for( i = 1; i < nsols; i++ )
                     if( SCIPsolGetIndex(sols[i]) > SCIPsolGetIndex(sols[solindex]) )
                        solindex = i;

                  newsol = sols[solindex];

                  assert(graph_sol_valid(scip, graph, orgresults));

                  if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, bestsol) - SCIPprobdataGetOffset(scip), pobj) )
                     heurdata->nfailures = 0;
               }
            }

            bisuccess = success || bisuccess;

            if( solgraph->terms <= 1 || probtype != STP_SPG )
            {
               break;
            }
         }

         if( bisuccess )
            SCIPdebugMessage("solution found by rec! \n");

         SCIPfreeBufferArray(scip, &stnodes);

         /* free memory */
         SCIPfreeBufferArrayNull(scip, &nodepriority);
         SCIPfreeBufferArrayNull(scip, &costrev);
         SCIPfreeBufferArrayNull(scip, &cost);
         SCIPfreeBufferArrayNull(scip, &results);
         SCIPfreeMemoryArray(scip, &edgeweight);
         SCIPfreeMemoryArray(scip, &edgeancestor);
         graph_free(scip, solgraph, TRUE);
      }
   }

   return SCIP_OKAY;
}
#endif

/** heuristic to exclude vertices or edges from a given solution (and inserting other edges) to improve objective (also prunes original solution) */
SCIP_RETCODE SCIPStpHeurRecExclude(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   int*                  result,             /**< edge solution array (UNKNOWN/CONNECT) */
   int*                  result2,            /**< second edge solution array or NULL */
   int*                  newresult,          /**< new edge solution array (UNKNOWN/CONNECT) */
   int*                  dnodemap,           /**< node array for internal use */
   STP_Bool*             stvertex,           /**< node array for internally marking solution vertices */
   SCIP_Bool*            success             /**< solution improved? */
   )
{
   SCIP_HEURDATA* tmheurdata;
   GRAPH* newgraph;
   SCIP_Real dummy;
   int*   unodemap;
   int    j;
   int    root;
   int    nedges;
   int    nnodes;
   int    nsolterms;
   int    nsoledges;
   int    nsolnodes;
   int    best_start;
   SCIP_Bool pcmw;
   SCIP_Bool mergesols;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(result != NULL);
   assert(stvertex != NULL);
   assert(success != NULL);
   assert(dnodemap != NULL);
   assert(newresult != NULL);
   assert(graph->stp_type != STP_DHCSTP);

   pcmw = (graph->stp_type == STP_PCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP || graph->stp_type == STP_RPCSPG);
   nedges = graph->edges;
   nnodes = graph->knots;
   *success = TRUE;

   // todo
   assert(graph->stp_type == STP_MWCSP);

   /* killed solution edge? */
   for( int e = 0; e < nedges; e++ )
      if( result[e] == CONNECT && graph->oeat[e] == EAT_FREE )
         return SCIP_OKAY;

   BMSclearMemoryArray(stvertex, nnodes);

   SCIP_Real cc = graph_computeSolVal(graph->cost, result, 0.0, nedges);

   int *x1;
   int *x2;
   SCIP_CALL( SCIPallocBufferArray(scip, &x1, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &x2, nnodes) );
   BMSclearMemoryArray(x1, nnodes);
   BMSclearMemoryArray(x2, nnodes);
   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         x1[graph->tail[e]] = TRUE;
         x1[graph->head[e]] = TRUE;
      }
   }

   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         stvertex[graph->tail[e]] = TRUE;
         stvertex[graph->head[e]] = TRUE;
      }
      result[e] = UNKNOWN;
   }

   if( pcmw )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, result, stvertex) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, result, stvertex) );

#if 1
   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         x2[graph->tail[e]] = TRUE;
         x2[graph->head[e]] = TRUE;
      }
   }

   for( int k = 0; k < nnodes; k++ )
      if(x1[k] != x2[k])
      {
         printf("FAIL for %d \n", k);
         return SCIP_ERROR;
      }

   if( !SCIPisEQ(scip, cc, graph_computeSolVal(graph->cost, result, 0.0, nedges)) )
   {
      printf("fail2 %f %f \n", cc, graph_computeSolVal(graph->cost, result, 0.0, nedges));
      return SCIP_ERROR;
   }

   SCIPfreeBufferArray(scip, &x1);
   SCIPfreeBufferArray(scip, &x2);
#endif


   /*** 1. step: for solution S and original graph (V,E) initialize new graph (V[S], (V[S] X V[S]) \cup E) ***/

   BMSclearMemoryArray(stvertex, nnodes);

   root = graph->source[0];
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
            if( Is_pterm(graph->term[head]) )
            {
               stvertex[head] = 1;
               nsolterms++;
               nsolnodes++;
            }

            continue;
         }

         //todo
         if( stvertex[head] )
         {
            printf("ohoh %d \n", 0);
            return SCIP_ERROR;
         }

         stvertex[head] = 1;

         if( Is_pterm(graph->term[head]) )
            nsolterms++;
         nsolnodes++;
      }
   }

   mergesols = FALSE;

   /* if there is second solution, check whether it can be merged with first one */
   if( result2 != NULL )
   {
      int ed = 0;
      for( ; ed < nedges; ed++ )
         if( result2[ed] == CONNECT )
         {
            int e;
            int k = graph->head[ed];

            if( Is_term(graph->term[k]) )
               continue;

            if( stvertex[k] )
               break;

            for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               if( stvertex[graph->head[e]] )
                  break;

            if( e != EAT_LAST )
               break;
         }

      if( ed != nedges )
         mergesols = TRUE;
   }

   /* mark nodes in second solution */
   if( mergesols )
   {
      for( int e = 0; e < nedges; e++ )
      {
         if( result2[e] == CONNECT )
         {
            int tail = graph->tail[e];
            int head = graph->head[e];

            if( stvertex[head] )
            {
               stvertex[head]++;
               continue;
            }

            if( tail == root )
            {
               /* there might be only one node */
               if( Is_pterm(graph->term[head]) )
               {
                  stvertex[head] = 1;
                  nsolterms++;
                  nsolnodes++;
               }

               continue;
            }

            stvertex[head] = 1;

            if( Is_pterm(graph->term[head]) )
               nsolterms++;
            nsolnodes++;
         }
      }
   }
   else if( result2 != NULL )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }


//todo
   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         int tail = graph->tail[e];
         int head = graph->head[e];

         if( tail == root )
         {
            continue;
         }


         if(!stvertex[head] )
                        {
                           printf("FAILLL HEAD %d \n", head);
                           return SCIP_ERROR;

                        }

         if( !stvertex[tail] )
         {
            printf("FAILLL TAIL %d \n", tail);
            printf(" %d \n", 0);
            return SCIP_ERROR;

         }


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
   SCIP_CALL( graph_init(scip, &newgraph, nsolnodes, nsoledges, 1, 0) );

   SCIP_CALL( SCIPallocBufferArray(scip, &unodemap, nsolnodes) );

   if( graph->stp_type == STP_RSMT || graph->stp_type == STP_OARSMT || graph->stp_type == STP_GSTP )
      newgraph->stp_type = STP_SPG;
   else
      newgraph->stp_type = graph->stp_type;

   if( pcmw )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(newgraph->prize), nsolnodes) );
   }

   j = 0;
   for( int i = 0; i < nnodes; i++ )
   {
      if( stvertex[i] )
      {
         if( pcmw )
         {
            if( (!Is_term(graph->term[i])) )
               newgraph->prize[j] = graph->prize[i];
            else
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
   newgraph->source[0] = dnodemap[root];
   if( newgraph->stp_type == STP_RPCSPG )
      newgraph->prize[newgraph->source[0]] = FARAWAY;

   assert(newgraph->source[0] >= 0);

   /* add edges */
   for( int i = 0; i < nedges; i += 2 )
      if( stvertex[graph->tail[i]] && stvertex[graph->head[i]] && graph->oeat[i] != EAT_FREE )
         graph_edge_add(scip, newgraph, dnodemap[graph->tail[i]], dnodemap[graph->head[i]], graph->cost[i], graph->cost[i + 1]);

   assert(newgraph->edges == nsoledges);

#if 0
   if( pcmw && result2 != NULL )
   {
      const int newroot = newgraph->source[0];
      srand(nsolnodes);

      for( int i = 0; i < nsolnodes; i++ )
      {
         if( stvertex[unodemap[i]] > 1 && !Is_term(newgraph->term[i]) )
         {
            if( SCIPisPositive(scip, newgraph->prize[i]) )
            {
               int e;
               int term;
               const SCIP_Real factor = rand() % 4 + 4;

               assert(Is_pterm(newgraph->term[i]));

               newgraph->prize[i] *= factor;

               for( e = newgraph->outbeg[i]; e != EAT_LAST; e = newgraph->oeat[e] )
                  if( newgraph->head[e] != newroot && Is_term(newgraph->term[newgraph->head[e]]) )
                     break;

               assert(e != EAT_LAST);

               // todo
               if(e == EAT_LAST)
                  return SCIP_ERROR;

               term = newgraph->head[e];
               for( e = newgraph->inpbeg[term]; e != EAT_LAST; e = newgraph->ieat[e] )
                  if( newgraph->tail[e] == newroot )
                     break;

               assert(e != EAT_LAST);

               // todo
               if(e == EAT_LAST)
                  return SCIP_ERROR;

               newgraph->cost[e] *= factor;
            }
            else if( SCIPisNegative(scip, newgraph->prize[i]) )
            {
               const SCIP_Real factor = 1.0 / (rand() % 4 + 4);

               newgraph->prize[i] *= factor;

               for( int e = newgraph->inpbeg[i]; e != EAT_LAST; e = newgraph->ieat[e] )
                  newgraph->cost[e] *= factor;
            }
         }
      }
   }
#endif

#ifdef DEBUG
   if( j != nsolnodes )
   {
      printf("miscal %d != %d \n", j, nsolnodes);
      return SCIP_ERROR;
   }

   for( int k = 0; k < newgraph->knots; k++ )
   {
      if( Is_pterm(newgraph->term[k]) )
      {
         int e = newgraph->outbeg[k];
         for( ; e != EAT_LAST; e = newgraph->oeat[e] )
         {
            int head = newgraph->head[e];
            if( newgraph->source[0] != head && Is_term(newgraph->term[head]) )
               break;

         }
         if( e == EAT_LAST )
         {
            printf("1graph construction fail in heur_rec \n");
            return SCIP_ERROR;
         }
      }
      if( Is_term(newgraph->term[k]) )
      {
         int e = newgraph->outbeg[k];
         for( ; e != EAT_LAST; e = newgraph->oeat[e] )
         {
            int head = newgraph->head[e];
            if( newgraph->source[0] != head && Is_pterm(newgraph->term[head]) )
               break;

         }
         if( e == EAT_LAST )
         {
            printf("2graph construction fail in heur_rec \n");
            return SCIP_ERROR;
         }
      }
   }

   for( int e = newgraph->outbeg[newgraph->source[0]] ; e != EAT_LAST; e = newgraph->oeat[e] )
   {
      if( Is_term(newgraph->term[newgraph->head[e]]) &&  SCIPisZero(scip, newgraph->cost[e] ) )
      {
         printf("TERM FAIL with %d \n", newgraph->head[e]);
         return SCIP_ERROR;
      }
      if( Is_pterm(newgraph->term[newgraph->head[e]]) &&  !SCIPisZero(scip, newgraph->cost[e] ) )
      {
         printf("PTERM FAIL with %d \n", newgraph->head[e]);
         return SCIP_ERROR;
      }
      if( newgraph->term[newgraph->head[e]] == -1 )
      {
         printf("NEWTERM FAIL with %d \n", newgraph->head[e]);
         return SCIP_ERROR;
      }
   }

   if( !graph_valid(newgraph) )
   {
      printf("GRAPH NOT VALID %d \n", 0);
      return SCIP_ERROR;
   }
#endif


   /*** step 2: presolve ***/

   newgraph->norgmodelknots = nsolnodes;

   dummy = 0.0;
   SCIP_CALL( reduce(scip, &newgraph, &dummy, 1, 5) );


   /*** step 3: compute solution on new graph ***/


   /* get TM heuristic data */
   assert(SCIPfindHeur(scip, "TM") != NULL);
   tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   graph_path_init(scip, newgraph);

   /* compute Steiner tree to obtain upper bound */
   best_start = newgraph->source[0];
   SCIP_CALL( SCIPStpHeurTMRun(scip, tmheurdata, newgraph, NULL, &best_start, newresult, MIN(50, nsolterms), newgraph->source[0], newgraph->cost, newgraph->cost, &dummy, NULL, 0.0, success, FALSE) );

   graph_path_exit(scip, newgraph);

   assert(*success);
   assert(graph_sol_valid(scip, newgraph, newresult));

#ifdef DEBUG
   if( !graph_sol_valid(scip, newgraph, newresult) )
   {
      printf("FAIL %d \n", 0);
      return SCIP_ERROR;
   }
#endif

   /*** step 4: retransform solution to original graph ***/


   BMSclearMemoryArray(stvertex, nnodes);

   for( int e = 0; e < nsoledges; e++ )
      if( newresult[e] == CONNECT )
         markSolVerts(newgraph, newgraph->ancestors[e], unodemap, stvertex);

   /* retransform edges fixed during graph reduction */
   markSolVerts(newgraph, newgraph->fixedges, unodemap, stvertex);

   if( pcmw )
      for( int k = 0; k < nsolnodes; k++ )
         if( stvertex[unodemap[k]] )
            markSolVerts(newgraph, newgraph->pcancestors[k], unodemap, stvertex);

   for( int e = 0; e < nedges; e++ )
      newresult[e] = UNKNOWN;

   if( pcmw )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, newresult, stvertex) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, graph, graph->cost, 0, newresult, stvertex) );

   /* solution better than original one?  */

   if( SCIPisLT(scip, graph_computeSolVal(graph->cost, newresult, 0.0, nedges),
         graph_computeSolVal(graph->cost, result, 0.0, nedges)) )
      *success = TRUE;
   else
      *success = FALSE;

   if( !graph_sol_valid(scip, graph, newresult) )
      *success = FALSE;

   SCIPfreeBufferArray(scip, &unodemap);
   graph_free(scip, newgraph, TRUE);

#ifdef DEBUG
   if( !graph_sol_valid(scip, graph, newresult) )
   {
      printf("invalid sol in REC \n");
      return SCIP_ERROR;
   }
#endif
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

   return SCIP_OKAY;
}
