/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_rec.c
 * @brief  primal recombination heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a recombination heuristic for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
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
#include "heur_local.h"
#include "grph.h"
#include "heur_tm.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"

#define HEUR_NAME             "rec"
#define HEUR_DESC             "LNS heuristic fixing all variables corresponding to edges used in at least one of several selected solutions"
#define HEUR_DISPCHAR         'R'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_MAXNSOLS       50            /**< maximum number of (good) solutions be regarded in the subproblem                  */
#define DEFAULT_NUSEDSOLS      4             /**< number of solutions that will be taken into account                               */
#define DEFAULT_RANDSEED       0             /**< random seed                                                                       */
#define DEFAULT_NTMRUNS        50            /**< number of runs in TM heuristic                                                    */
#define DEFAULT_NWAITINGSOLS   2             /**< max number of new solutions to be available before executing the heuristic again  */


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
   unsigned int          randseed;           /**< seed value for random number generator                            */
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

   heurdata->randseed = (unsigned int)newrandseed;

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
         factor = SCIPgetRandomInt(1000, 1400, &(heurdata->randseed));
         mult =  (double) factor * (1.0 / avg);
      }
      else if( SCIPisLT(scip, avg, 2.6) )
      {
         factor = SCIPgetRandomInt(200, 1000, &(heurdata->randseed));
         mult =  (double) factor * (1.0 / avg);
      }
      else if( nusedsols == 3 && SCIPisLT(scip, avg, 3.6) )
      {
         factor = SCIPgetRandomInt(40, 100, &(heurdata->randseed));
         mult =  (double) factor * (1.0 / avg);
      }
   }
   else
   {
      if( SCIPisLT(scip, avg, 1.6) )
      {
         factor = SCIPgetRandomInt(1400, 1800, &(heurdata->randseed));
      }
      else if( SCIPisLT(scip, avg, 2.6) )
      {
         factor = SCIPgetRandomInt(400, 1400, &(heurdata->randseed));
      }
      else if( SCIPisLT(scip, avg, 3.6) )
      {
         factor = SCIPgetRandomInt(150, 250, &(heurdata->randseed));
      }
      else if( SCIPisLT(scip, avg, 4.6) )
      {
         factor = SCIPgetRandomInt(60, 90, &(heurdata->randseed));
      }
      else if( nusedsols >= 5 && SCIPisLT(scip, avg, 5.6) )
      {
         factor = SCIPgetRandomInt(20, 40, &(heurdata->randseed));
      }
      mult =  (double) factor * (1.0 / avg);
   }

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
   SCIP_Real* soltimes;
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
   char* soledges;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &soltimes, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledges, nedges / 2) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      soltimes[i] = SCIPgetSolTime(scip, sols[i]);
      solselected[i] = FALSE;
   }

   if( *newsol == NULL )
   {
      SCIPsortRealInt(soltimes, perm, nsols);
      i = nsols - 1;
      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         *newsol = sols[perm[i]];
      }
      else
      {
	 i = SCIPgetRandomInt(0, nsols - 1, &(heurdata->randseed));
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

   SCIPpermuteIntArray(perm, 0, maxnsols, &(heurdata->randseed));
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
   SCIPfreeBufferArray(scip, &soledges);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &soltimes);
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
   SCIP_Real* soltimes;
   SCIP_SOL** sols;                          /* array of all solutions found so far         */
   int i;
   int end;
   int maxnsols;
   int nselectedsols;
   int shift;
   int nsols;                                /* number of all solutions found so far        */
   int nusedsols;                            /* number of solutions to use in rec           */
   int* perm;
   int* solselected;

   assert(selection != NULL);

   /* get solution data */
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   maxnsols = heurdata->maxnsols;
   nusedsols = heurdata->nusedsols;
   assert(nusedsols <= nsols);
   nselectedsols = 0;

   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soltimes, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      soltimes[i] = SCIPgetSolTime(scip, sols[i]);
      solselected[i] = FALSE;
   }

   if( *newsol == NULL )
   {
      SCIPsortRealInt(soltimes, perm, nsols);
      i = nsols - 1;

      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         *newsol = sols[perm[i]];
      }
      else
      {
	 i = SCIPgetRandomInt(0, nsols - 1, &(heurdata->randseed));
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
      if( nsols >= nsols )
	 i = 0;
      solselected[i] = TRUE;
      selection[nselectedsols++] = i;
   }

   if( !randomize )
   {
      end = (int) ((SCIPgetRandomReal(1.0, (double) nusedsols - 1, &(heurdata->randseed))) );

      shift = SCIPgetRandomInt(end, 2 * nusedsols - 1, &(heurdata->randseed));

      if( shift > nsols )
         shift = nsols;
      SCIPpermuteIntArray(perm, 0, shift, &(heurdata->randseed));

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
   if( nselectedsols < nusedsols )
   {
      SCIPpermuteIntArray(perm, 0, maxnsols, &(heurdata->randseed));
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
   SCIPfreeBufferArray(scip, &soltimes);
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
   int    nedges;
   int    nnodes;
   int    nsoledges;
   int    nsolnodes;
   int*   dnodemap;
   int*   solselection;          /* pool of solutions rec will use */
   char*  solnode;               /* marks nodes contained in at least one solution */
   char*  soledge;               /* marks edges contained in at least one solution */

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

   assert(sols != NULL);

   vars = SCIPprobdataGetEdgeVars(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solselection, heurdata->nusedsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solnode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dnodemap, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soledge, nedges / 2) );

   for( i = 0; i < nedges / 2; i++ )
      soledge[i] = FALSE;
   for( i = 0; i < nnodes; i++ )
   {
      solnode[i] = FALSE;
      dnodemap[i] = UNKNOWN;
   }

   /* select solutions to be merged */
   if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING
      || graph->stp_type == STP_DEG_CONS )
      SCIP_CALL( selectdiffsols(scip, graph, heurdata, vars, newsol, solselection, success) );
   else
      SCIP_CALL( selectsols(scip, heurdata, newsol, solselection, randomize) );

   if( *success )
   {
      assert(heurdata->nselectedsols > 0);

      /* count and mark selected nodes and edges */
      for( i = 0; i < nedges; i += 2 )
      {
         for( j = 0; j < heurdata->nselectedsols; j++ )
         {
            varsolval = SCIPgetSolVal(scip, sols[solselection[j]], vars[i]);
            varrevsolval = SCIPgetSolVal(scip, sols[solselection[j]], vars[i + 1]);

            if( SCIPisEQ(scip, varsolval, 1.0) || SCIPisEQ(scip, varrevsolval, 1.0) )
            {
               nsoledges++;
               soledge[i / 2] = TRUE;
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
      if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
      {
         for( i = graph->outbeg[graph->source[0]]; i != EAT_LAST; i = graph->oeat[i] )
         {
            if( soledge[i / 2] == FALSE && Is_term(graph->term[graph->head[i]]) )
            {
               nsoledges++;
               soledge[i / 2] = TRUE;
               assert(solnode[graph->head[i]]);
            }
         }
      }

      if( graph->stp_type == GSTP )
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

      if( graph->stp_type == STP_GRID || graph->stp_type == STP_OBSTACLES_GRID )
         newgraph->stp_type = STP_UNDIRECTED;
      else
         newgraph->stp_type = graph->stp_type;

      newgraph->hoplimit = graph->hoplimit;
      j = 0;
      for( i = 0; i < nnodes; i++ )
      {
         if( solnode[i] )
         {
            dnodemap[i] = j++;
            if( Is_term(graph->term[i]) )
               graph_knot_add(newgraph, 0);
            else
               graph_knot_add(newgraph, -1);
         }
      }
      /* set root */
      newgraph->source[0] = dnodemap[graph->source[0]];
      assert(newgraph->source[0] >= 0);

      /* copy max degrees*/
      if( graph->stp_type == STP_DEG_CONS )
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
      for( i = 0; i < 2 * nsoledges; i++ )
         assert((*edgeweight)[i] >= 1);
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


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRec)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRec(scip) );

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

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#else
   heurdata->randseed = 0;
#endif

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
   SCIP_Real* nval;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real* nodepriority = NULL;
   SCIP_Real pobj;
   SCIP_Real avg;
   SCIP_Real maxcost = 0.0;
   SCIP_Real hopfactor = 0.1;
   SCIP_Longint nsols;
   SCIP_Bool success;
   SCIP_Bool fixed;
   SCIP_Bool randomize;
   SCIP_Bool modcost;
   SCIP_Bool solfound;
   IDX* curr;
   IDX** ancestors;
   int i;
   int e;
   int v;
   int head;
   int tail;
   int runs;
   int count;
   int solindex;
   int nedges;
   int nnodes;
   int probtype;
   int nsoledges;
   int best_start;
   int lastsolindex;
   int* perm;
   int* results;
   int* edgeweight;
   int* orgresults;
   int* edgeancestor;
   char* stnodes;

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
   nsols = SCIPgetNSolsFound(scip);

   /* only call heuristic, if sufficiently many solutions are available */
   if( SCIPgetNSols(scip) < DEFAULT_NUSEDSOLS + 1 )
      return SCIP_OKAY;

   /* suspend heuristic? */
   if( probtype == STP_PRIZE_COLLECTING || probtype == STP_MAX_NODE_WEIGHT || probtype == STP_ROOTED_PRIZE_COLLECTING
      || probtype == STP_HOP_CONS || probtype == STP_DEG_CONS )
   {
      if( heurdata->ncalls == 0 )
	 i = 0;
      else if( probtype == STP_ROOTED_PRIZE_COLLECTING || probtype == STP_DEG_CONS )
         i = MAX(heurdata->nwaitingsols, 2 * heurdata->nfailures);
      else
	 i = MAX(heurdata->nwaitingsols, heurdata->nfailures);

      if( nsols <= heurdata->nlastsols + i )
         return SCIP_OKAY;
   }
   else
   {
      i = MIN(heurdata->nwaitingsols, heurdata->nfailures);
      if( nsols <= heurdata->nlastsols + i
         && heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         return SCIP_OKAY;
   }

   /* get edge variables */
   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);
   assert(vars[0] != NULL);

   nedges = graph->edges;
   nnodes = graph->knots;
   results = NULL;
   randomize = TRUE;


   heurdata->ncalls++;
   *result = SCIP_DIDNOTRUN;

   if( probtype == STP_PRIZE_COLLECTING || probtype == STP_MAX_NODE_WEIGHT || probtype == STP_ROOTED_PRIZE_COLLECTING
      || probtype == STP_HOP_CONS )
      runs = 8;
   else
      runs = 8;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, runs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgresults, nedges) );

   for( v = 0; v < runs; v++ )
      perm[v] = v;

   if( probtype == STP_MAX_NODE_WEIGHT || probtype == STP_HOP_CONS || probtype == STP_DEG_CONS )
   {
      newsol = (SCIPgetSols(scip))[0];
   }
   /* first run? */
   else if( heurdata->lastsolindex == -1 )
   {
      newsol = (SCIPgetSols(scip))[SCIPgetRandomInt(0, heurdata->nusedsols - 1, &(heurdata->randseed))];
   }
   else
   {
      newsol = NULL;
   }

   /* save last solution index */
   solindex = 0;
   nsols = SCIPgetNSols(scip);
   assert(nsols > 0);
   sols = SCIPgetSols(scip);
   for( i = 1; i < nsols; i++ )
      if( SCIPisGT(scip, SCIPgetSolTime(scip, sols[i]), SCIPgetSolTime(scip, sols[solindex])) )
         solindex = i;

   lastsolindex = SCIPsolGetIndex(sols[solindex]);

   count = 0;

   if( SCIPgetRandomInt(0, 10, &(heurdata->randseed)) == 1 )
      modcost = FALSE;
   else
      modcost = TRUE;

   solfound = FALSE;
   for( v = 0; v < 8 * runs && !SCIPisStopped(scip); v++ )
   {
      if( count >= runs )
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
      if( perm[count] <= 1 )
      {
         heurdata->nusedsols = 2;
	 if( SCIPgetRandomInt(0, 1, &(heurdata->randseed)) == 1 )
            randomize = TRUE;
         else
	    randomize = FALSE;
      }
      else if( perm[count] <= 4 )
      {
	 if( SCIPgetRandomInt(0, 2, &(heurdata->randseed)) == 1 )
            randomize = TRUE;
         else
	    randomize = FALSE;
         heurdata->nusedsols = DEFAULT_NUSEDSOLS - 1 ;
      }
      else if( perm[count] <= 6 )
      {
         heurdata->nusedsols = DEFAULT_NUSEDSOLS;
      }
      else if( perm[count] <= 7 )
      {
         heurdata->nusedsols = DEFAULT_NUSEDSOLS + 1;
      }

      /* build up a new graph, consisting of several solutions */
      SCIP_CALL( buildsolgraph(scip, heurdata, &newsol, graph, &solgraph, &edgeancestor, &edgeweight, &success, randomize) );

      if( success )
      {
         assert(newsol != NULL);

         /* get TM heuristic data */
         assert(SCIPfindHeur(scip, "TM") != NULL);
         tmheurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

         /* presolve new graph */

         /* reduce new graph */
         if( probtype == STP_PRIZE_COLLECTING || probtype == STP_MAX_NODE_WEIGHT || probtype == STP_ROOTED_PRIZE_COLLECTING
            || probtype == STP_HOP_CONS || probtype == STP_DEG_CONS )
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 0, 2) );
         else
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 1, 2) );

         SCIP_CALL( graph_pack(scip, solgraph, &psolgraph, FALSE) );
	 solgraph = psolgraph;
	 assert(graph_valid(solgraph));
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
                     if(  SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
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
                  nodepriority[solgraph->head[e]] += avg - 1;
                  nodepriority[solgraph->tail[e]] += avg - 1;
                  if( modcost )
                  {
                     mult = costMultiplier(scip, heurdata, avg);
                     cost[e] = cost[e] * mult;
                  }
               }

               if( probtype == STP_HOP_CONS && SCIPisLT(scip, cost[e], BLOCKED) && SCIPisGT(scip, cost[e], maxcost) )
                  maxcost = cost[e];
            }
            for( e = 0; e < nsoledges; e++)
               costrev[e] = cost[flipedge(e)];

            /* init shortest path algorithm */
            SCIP_CALL( graph_path_init(scip, solgraph) );

            /* set (edge) result array to default */
            for( e = 0; e < nsoledges; e++ )
               results[e] = UNKNOWN;
            /* run TM heuristic */
            SCIP_CALL( SCIPheurComputeSteinerTree(scip, tmheurdata, solgraph, NULL, &best_start, results, heurdata->ntmruns,
                  solgraph->source[0], cost, costrev, &hopfactor, nodepriority, maxcost, &success) );

            if( !success )
            {
               graph_path_exit(scip, solgraph);
               SCIPfreeBufferArrayNull(scip, &nodepriority);
               SCIPfreeBufferArrayNull(scip, &costrev);
               SCIPfreeBufferArrayNull(scip, &cost);
               SCIPfreeBufferArrayNull(scip, &results);
               SCIPfreeMemoryArray(scip, &edgeancestor);
               SCIPfreeMemoryArray(scip, &edgeweight);
               graph_free(scip, solgraph, TRUE);
               continue;
            }

            assert(graph_valid(solgraph));
            assert(graph_sol_valid(scip, solgraph, results));

            /* run local heuristic */
            if( solgraph->stp_type != STP_HOP_CONS && solgraph->stp_type != STP_MAX_NODE_WEIGHT && probtype != STP_DEG_CONS )
               SCIP_CALL( SCIPheurImproveSteinerTree(scip, solgraph, cost, costrev, results) );

            if( probtype == STP_UNDIRECTED )
               assert(graph_sol_valid(scip, solgraph, results));

            graph_path_exit(scip, solgraph);
         }

         for( i = 0; i < nedges; i++ )
            orgresults[i] = UNKNOWN;


         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, nedges) );

         for( i = 0; i < nnodes; i++ )
            stnodes[i] = FALSE;

         /* retransform solution found by TM heuristic */
         if( solgraph->terms > 1 )
         {
	    assert(results != NULL);
            for( e = 0; e < nsoledges; e++ )
            {
               if( results[e] == CONNECT )
               {
                  /* iterate through list of ancestors */
                  curr = ancestors[e];
                  while( curr != NULL )
                  {
                     i = edgeancestor[curr->index];
                     if( probtype == STP_DEG_CONS )
                        orgresults[i] = CONNECT;
                     head = graph->head[i];
                     tail = graph->tail[i];
                     stnodes[tail] = TRUE;
                     stnodes[head] = TRUE;
                     curr = curr->parent;
                  }
               }
            }
         }

         /* retransform edges fixed during graph reduction */
         curr = solgraph->fixedges;
         while( curr != NULL )
         {
            i = edgeancestor[curr->index];
            if( probtype == STP_DEG_CONS )
               orgresults[i] = CONNECT;
            head = graph->head[i];
            tail = graph->tail[i];
            stnodes[tail] = TRUE;
            stnodes[head] = TRUE;
            curr = curr->parent;
         }
         /* prune solution (in the original graph) */
         if( probtype == STP_PRIZE_COLLECTING || probtype == STP_MAX_NODE_WEIGHT || probtype == STP_ROOTED_PRIZE_COLLECTING )
            SCIP_CALL( SCIPheurPrunePCSteinerTree(scip, graph, graph->cost, orgresults, stnodes) );
         else if( probtype == STP_DEG_CONS )
            SCIP_CALL( SCIPheurPruneDegConsSteinerTree(scip, graph, orgresults, stnodes) );
         else
            SCIP_CALL( SCIPheurPruneSteinerTree(scip, graph, graph->cost, 0, orgresults, stnodes) );

         SCIPfreeBufferArray(scip, &stnodes);
         pobj = 0.0;

         for( e = 0; e < graph->edges; e++ )
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
            sol = NULL;
            SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

            if( success )
            {
               *result = SCIP_FOUNDSOL;
               solfound = TRUE;
               nsols = SCIPgetNSols(scip);
               assert(nsols > 0);
               sols = SCIPgetSols(scip);
               solindex = 0;
               for( i = 1; i < nsols; i++ )
                  if( SCIPisGT(scip, SCIPgetSolTime(scip, sols[i]), SCIPgetSolTime(scip, sols[solindex])) )
                     solindex = i;
               newsol = sols[solindex];
               assert(graph_sol_valid(scip, graph, orgresults));
            }
            else
            {
               count++;
            }
         }
         else
         {
            count++;
         }

         /* free memory */
         SCIPfreeBufferArrayNull(scip, &nodepriority);
         SCIPfreeBufferArrayNull(scip, &costrev);
         SCIPfreeBufferArrayNull(scip, &cost);
         SCIPfreeBufferArrayNull(scip, &results);

         SCIPfreeMemoryArray(scip, &edgeancestor);
         SCIPfreeMemoryArray(scip, &edgeweight);
         graph_free(scip, solgraph, TRUE);

      }
      else
      {
         count++;
      }
   }

   if( *result == SCIP_FOUNDSOL )
      heurdata->nfailures = 0;
   else
      heurdata->nfailures++;

   heurdata->lastsolindex = lastsolindex;
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));
   heurdata->nlastsols = SCIPgetNSolsFound(scip);
   SCIPfreeBufferArray(scip, &orgresults);
   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the rec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRec(
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

   /* add rec primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nwaitingsols",
         "number of solution findings to be in abeyance",
         &heurdata->nwaitingsols, FALSE, DEFAULT_NWAITINGSOLS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/randseed",
         "random seed for heuristic",
         NULL, FALSE, DEFAULT_RANDSEED, 0, INT_MAX, paramChgdRandomseed, (SCIP_PARAMDATA*)heurdata) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnsols",
         "max size of solution pool for heuristic",
         &heurdata->maxnsols, FALSE, DEFAULT_MAXNSOLS, 5, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/ntmruns",
         "number of runs in TM",
         &heurdata->ntmruns, FALSE, DEFAULT_NTMRUNS, 1, INT_MAX, NULL, NULL) );

   heurdata->nusedsols = DEFAULT_NUSEDSOLS;

   return SCIP_OKAY;
}
