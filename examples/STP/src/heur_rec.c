/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_rec.c
 * @brief  rec primal heuristic
 * @author Daniel Rehfeldt
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
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNSOLS       50             /* maximum number of (good) solutions be regarded in the subproblem                   */
#define DEFAULT_NUSEDSOLS     4              /* number of solutions that will be taken into account                   */
#define DEFAULT_RANDSEED       0              /* random seed                                                               */
#define DEFAULT_NTMRUNS       50             /**< number of runs in TM heuristic        */
#define DEFAULT_NWAITINGSOLS  2              /* max number of new solutions to be available before executing the heuristic again  */
#define HASHSIZE_SOLS         11113          /* size of hash table for solution tuples in rec heuristic         */

/* event handler properties */
#define EVENTHDLR_NAME         "Rec"
#define EVENTHDLR_DESC         "LP event handler for rec heuristic"

#ifdef WITH_UG
extern
int getUgRank();
#endif

/*
 * Data structures
 */

typedef struct SolTuple SOLTUPLE;

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   maxnsols;           /**< maximum number of (good) solutions be regarded in the subproblem                   */
   SCIP_Longint          ncalls;             /**< number of calls                              */
   SCIP_Longint          nlastsols;          /**< number of solutions during the last run                 */
   int                   ntmruns;            /**< number of runs in TM heuristic        */
   int                   nusedsols;          /**< number of solutions that will be taken into account               */
   int                   nselectedsols;      /**< number of solutions actually selected */
   int                   nwaitingsols;       /**< number of new solutions before executing the heuristic again      */
   int                   nfailures;          /**< number of failures since last successful call                     */
   unsigned int          randseed;           /**< seed value for random number generator                            */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable used to store the solution tuples already used          */
   SOLTUPLE*             lasttuple;          /**< last tuple of solutions created by rec                      */
};

/** n-tuple of solutions and their hashkey */
struct SolTuple
{
   int*                  indices;            /**< sorted array of solution indices                                 */
   int                   size;               /**< size of the array (should be heurdata->nusedsols)                */
   unsigned int          key;                /**< hashkey of the tuple                                             */
   SOLTUPLE*             prev;               /**< previous solution tuple created                                  */
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


#if 0
/** for debug purposes only */
static
SCIP_RETCODE printGraph(
   SCIP* scip,
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   int*                  result
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;
   char* stnodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, graph->knots ) );

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->knots; e++ )
   {
      stnodes[e] = FALSE;
   }
   for( e = 0; e < graph->edges; e++ )
   {
      if( result[e] == CONNECT )
      {
	 stnodes[graph->tail[e]] = TRUE;
	 stnodes[graph->head[e]] = TRUE;
      }
   }

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
      if( stnodes[n] )
      {
         if( n == graph->source[0] )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
            m = 1;
         }
         else if( graph->term[n] == 0 )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
            e += 1;
         }
         else
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
         }

      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e ++ )
   {
      if( result[e] == CONNECT )
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);
	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      }
   }
   SCIPfreeBufferArray(scip, &stnodes);
   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

#endif

static
SCIP_Real costMultiplier(
   SCIP*              scip,
   SCIP_HEURDATA*     heurdata,           /**< primal heuristic data */
   SCIP_Real                avg
   )
{
   int factor = 1;
   int nusedsols = heurdata->nusedsols;
   SCIP_Real mult = 1;
   assert(SCIPisGE(scip, avg, 1));
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


static
SCIP_RETCODE selectdiffsols(
   SCIP*                 scip,
   GRAPH*                graph,
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_VAR**            vars,
   SCIP_SOL**            newsol,
   int*                  selection,
   SCIP_Bool*            success,
   SCIP_Bool             random
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
      //printf("(maxnsols: %d) real soltime 0:  %d \n", maxnsols, SCIPsolGetIndex(sols[perm[i]]));
      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         *newsol = sols[perm[i]];
      }
      else
      {
	 //printf("last solution has already been used \n");
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
      //  printf("soltime 0:  %f \n", SCIPgetSolTime(scip, sols[perm[i]]));
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

         if( diffnedges > 3 && eqnedges > 0  ) //&& strcmp(SCIPheurGetName(SCIPsolGetHeur(sols[perm[i]])), "trivial") != 0
	 {
	    selection[nselectedsols++] = k;
            solselected[k] = TRUE;
	    *success = TRUE;
	    //printf("success diff: %d, eq: %d \n", diffnedges, eqnedges);
            if( nselectedsols >= nusedsols )
               break;
	 }
	 else
	 {
            //printf("no success diff: %d, eq: %d \n", diffnedges, eqnedges);
	 }
      }
   }
   /*
     printf("newsols \n");
     for( i = 0; i < nselectedsols; i++ )
     printf("newsol: found by: %s %d\n", SCIPheurGetName(SCIPsolGetHeur(sols[selection[i]])), SCIPsolGetIndex(sols[selection[i]]));
   */
   assert(nselectedsols <= nusedsols);
   heurdata->nselectedsols = nselectedsols;
   SCIPfreeBufferArray(scip, &soltimes);
   SCIPfreeBufferArray(scip, &solselected);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &soledges);

   return SCIP_OKAY;
}



static
SCIP_RETCODE selectsols(
   SCIP*                 scip,
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_SOL**            newsol,
   int*                  selection,
   SCIP_Bool             random
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
   int nusedsols;                            /* number of solutions to use in rec     */
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
      //printf("(maxnsols: %d) real soltime 0:  %d \n", maxnsols, SCIPsolGetIndex(sols[perm[i]]));
      //printf(" lastsolindex:%d new: %d\n", heurdata->lastsolindex, SCIPsolGetIndex(sols[perm[i]]));
      /* has the latest solution already been tried? */
      if( heurdata->lastsolindex != SCIPsolGetIndex(sols[perm[i]]) )
      {
         *newsol = sols[perm[i]];
      }
      else
      {
	 //printf("last solution has already been used \n");
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
      //  printf("soltime 0:  %f \n", SCIPgetSolTime(scip, sols[perm[i]]));
   }

   //   printf("(maxnsols: %d) newsol: found by: %s %d\n",  maxnsols, SCIPheurGetName(SCIPsolGetHeur(*newsol)), SCIPsolGetIndex(*newsol));
   if( !random )
   {
      end = (int) ((SCIPgetRandomReal(1, nusedsols - 1, &(heurdata->randseed))) );

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
            //         printf("found by: %s %d\n", SCIPheurGetName(SCIPsolGetHeur(sols[perm[i]])), SCIPsolGetIndex(sols[perm[i]]));
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
            /*     if( SCIPsolGetHeur(sols[perm[i]]) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(sols[perm[i]])), "rec") == 0 )
                   continue;*/
            //  printf("random found by: %s %d\n", SCIPheurGetName(SCIPsolGetHeur(sols[perm[i]])), SCIPsolGetIndex(sols[perm[i]]));
            selection[nselectedsols++] = perm[i];
            if( nselectedsols >= nusedsols )
               break;
         }
      }
   }
   assert(nselectedsols <= nusedsols);
   heurdata->nselectedsols = nselectedsols;
   SCIPfreeBufferArray(scip, &soltimes);
   SCIPfreeBufferArray(scip, &solselected);
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/* merge selected solutions to a new graph */
static
SCIP_RETCODE buildsolgraph(
   SCIP*              scip,
   SCIP_HEURDATA*     heurdata,           /**< primal heuristic data */
   SCIP_SOL**         newsol,
   GRAPH*             graph,
   GRAPH**            solgraph,
   int**              edgeancestor,
   int**              edgeweight,
   SCIP_Bool*         success,
   SCIP_Bool          random
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
   int*   solselection;          /**< pool of solutions rec will use */
   char*  solnode;               /**< marks nodes contained in at least one solution */
   char*  soledge;               /**< marks edges contained in at least one solution */

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
      SCIP_CALL( selectdiffsols(scip, graph, heurdata, vars, newsol, solselection, success, random) );
   else
      SCIP_CALL( selectsols(scip, heurdata, newsol, solselection, random) );

   if( *success )
   {
      assert(heurdata->nselectedsols > 0);

      /* count and mark selected nodes and edges */
      for( i = 0; i < nedges; i += 2 )
      {
         /* */
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
      newgraph = graph_init(nsolnodes, 2 * nsoledges, 1, 0);
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
	 newgraph->maxdeg = malloc((size_t)nsolnodes * sizeof(int));
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
            graph_edge_add(newgraph, dnodemap[graph->tail[i]], dnodemap[graph->head[i]], graph->cost[i], graph->cost[i + 1]);

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
   SCIPfreeBufferArray(scip, &solnode);
   SCIPfreeBufferArray(scip, &soledge);
   SCIPfreeBufferArray(scip, &dnodemap);
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
   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#else
   heurdata->randseed = 0;
#endif

#if 0
   /* initialize hash table */
   SCIP_CALL( SCIPhashtableCreate(&heurdata->hashtable, SCIPblkmem(scip), HASHSIZE_SOLS,
         hashGetKeySols, hashKeyEqSols, hashKeyValSols, NULL) );
   assert(heurdata->hashtable != NULL );
#endif
   return SCIP_OKAY;
}



/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRec)
{  /*lint --e{715}*/
#if 0
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   soltuple = heurdata->lasttuple;

   /* free all soltuples iteratively */
   while( soltuple != NULL )
   {
      SOLTUPLE* tmp;
      tmp = soltuple->prev;
      SCIPfreeBlockMemoryArray(scip,&soltuple->indices,soltuple->size);
      SCIPfreeBlockMemory(scip,&soltuple);
      soltuple = tmp;
   }

   /* free hash table */
   assert(heurdata->hashtable != NULL );
   SCIPhashtableFree(&heurdata->hashtable);
#endif
   return SCIP_OKAY;
}



/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolRec)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of rec primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolRec NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolRec)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of rec primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolRec NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRec)
{
   SCIP_HEUR** heurs;
   SCIP_HEURDATA* heurdata;
   SCIP_HEURDATA* tmheurdata;
   SCIP_PROBDATA* probdata;
   SCIP_VAR** vars;
   SCIP_SOL** sols;
   GRAPH* graph;                             /* graph structure */
   GRAPH* solgraph;                             /* graph structure */
   SCIP_SOL* sol;                            /* new solution */
   SCIP_SOL*  newsol;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real* nval;
   SCIP_Real pobj;
   SCIP_Real avg;
   SCIP_Real maxcost = 0.0;
   SCIP_Bool success;
   SCIP_Bool fixed;
   SCIP_Bool random = TRUE;
   SCIP_Bool modcost;
   SCIP_Bool solfound;
   IDX* curr;
   IDX** ancestors;
   int i;
   int e;
   int v;
   int artroot;
   SCIP_Longint nsols;                                /* number of all solutions found so far */
   int nedges;
   int index;
   int nnodes;
   int count;
   int head;
   int tail;
   int nsoledges;
   int nheurs;
   int runs;
   int best_start;
   int lastsolindex;
   int* perm;
   int* results = NULL;
   int* orgresults;
   int* edgeancestor;
   int* edgeweight;
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

   /* get edge variables */
   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);
   assert(vars[0] != NULL);

   nedges = graph->edges;
   nnodes = graph->knots;
   nsols = SCIPgetNSolsFound(scip);
   /*
     if( graph->stp_type == STP_DEG_CONS )
     return SCIP_OKAY;
   */
   /* only call heuristic, if sufficiently many solutions are available */
   if( SCIPgetNSols(scip) < heurdata->nusedsols + 1 )
      return SCIP_OKAY;

   /* suspend heuristic? */
   if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING
      || graph->stp_type == STP_HOP_CONS || graph->stp_type == STP_DEG_CONS )
   {
      if( heurdata->ncalls == 0 )
	 i = 0;
      else if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING || graph->stp_type == STP_HOP_CONS || graph->stp_type == STP_DEG_CONS )
         i = MAX(heurdata->nwaitingsols, 2 * heurdata->nfailures);
      else
	 i = MAX(heurdata->nwaitingsols, heurdata->nfailures);

      if( SCIPisLE(scip, nsols, heurdata->nlastsols + i) )
         return SCIP_OKAY;
   }
   else
   {
      i = MIN(heurdata->nwaitingsols, heurdata->nfailures);
      if( SCIPisLE(scip, nsols, heurdata->nlastsols + i)
         && heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         return SCIP_OKAY;
   }

   //printf("nsols: %d nfails: %d \n", (int) nsols, heurdata->nfailures);

   heurdata->ncalls++;
   *result = SCIP_DIDNOTRUN;

   if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING
      || graph->stp_type == STP_HOP_CONS )
      runs = 8;
   else
      runs = 8;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, runs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgresults, nedges) );

   for( v = 0; v < runs; v++ )
      perm[v] = v;
   //SCIPpermuteIntArray(perm, 0, 8, &(heurdata->randseed));

   if( graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_HOP_CONS || graph->stp_type == STP_DEG_CONS )
   {
      newsol = (SCIPgetSols(scip))[0];//(SCIPgetSols(scip))[SCIPgetRandomInt(0, 1, &(heurdata->randseed))];
   }
   /* first run? */
   else if( heurdata->lastsolindex == -1 )
   {
      newsol = (SCIPgetSols(scip))[SCIPgetRandomInt(0, heurdata->nusedsols - 1, &(heurdata->randseed))];
      //printf("selecting random best sol! \n");
   }
   else
   {
      newsol = NULL;
   }

   /* save last solution index */
   index = 0;
   nsols = SCIPgetNSols(scip);
   assert(nsols > 0);
   sols = SCIPgetSols(scip);
   for( i = 1; i < nsols; i++ )
      if( SCIPisGT(scip, SCIPgetSolTime(scip, sols[i]), SCIPgetSolTime(scip, sols[index])) )
         index = i;

   lastsolindex = SCIPsolGetIndex(sols[index]);

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
            random = TRUE;
         else
	    random = FALSE;
      }
      else if( perm[count] <= 4 )
      {
	 if( SCIPgetRandomInt(0, 2, &(heurdata->randseed)) == 1 )
            random = TRUE;
         else
	    random = FALSE;
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
      SCIP_CALL( buildsolgraph(scip, heurdata, &newsol, graph, &solgraph, &edgeancestor, &edgeweight, &success, random) );

      if( success )
      {
         assert(newsol != NULL);

         /* get TM heuristic data */
         heurs = SCIPgetHeurs(scip);
         nheurs = SCIPgetNHeurs(scip);
         for( i = 0; i < nheurs; i++ )
            if( strcmp(SCIPheurGetName(heurs[i]), "TM") == 0 )
               break;

         assert(i < nheurs);
         tmheurdata = SCIPheurGetData(heurs[i]);

         /* presolve new graph */

         /* reduce new graph */
         if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING
            || graph->stp_type == STP_HOP_CONS || graph->stp_type == STP_DEG_CONS )
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 0, 2) );
         else
            SCIP_CALL( reduce(scip, &solgraph, &pobj, 4, 2) );

         solgraph = graph_pack(scip, solgraph, FALSE);
	 assert(graph_valid(solgraph));
         ancestors = solgraph->ancestors;
         nsoledges = solgraph->edges;
         /* if graph reduction solved the whole problem, solgraph has only one node */
         if( solgraph->terms > 1 )
         {
            SCIP_Real mult;
            /* allocate memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &results, nsoledges) );
            SCIP_CALL( SCIPallocBufferArray(scip, &cost, nsoledges) );
            SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nsoledges) );

            BMScopyMemoryArray(cost, solgraph->cost, nsoledges);

            /* hop constraint problem? */
            if( graph->stp_type == STP_HOP_CONS )
            {
               maxcost = 0.0;
               for( e = 0; e < nsoledges; e++)
               {
                  curr = ancestors[e];
                  avg = 0.0;
                  i = 0;
                  fixed = FALSE;
                  while( curr != NULL )
                  {
                     i++;
                     avg += edgeweight[curr->index];
                     if(  SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
                     {
                        fixed = TRUE;
                     }
                     curr = curr->parent;
                  }
                  avg = (double) avg / (double) i;
                  assert(avg >= 1);
                  if( fixed )
                  {
                     cost[e] = 1e+10;
                  }
                  else if( modcost )
                  {
                     mult = costMultiplier(scip, heurdata, avg);
                     cost[e] = cost[e] * mult;
                  }

                  if( SCIPisLT(scip, cost[e], 1e+8 ) && SCIPisGT(scip, cost[e], maxcost) )
                     maxcost = cost[e];
               }
               for( e = 0; e < nsoledges; e++)
                  costrev[e] = cost[flipedge(e)];
            }
            else
            {

               for( e = 0; e < nsoledges; e += 2)
               {
                  fixed = FALSE;
                  curr = ancestors[e + 1];
                  assert(curr != NULL);
                  avg = 0.0;
                  i = 0;
                  while( curr != NULL )
                  {
                     i++;
                     avg += edgeweight[curr->index];
                     if(  SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
                     {
                        fixed = TRUE;
                     }
                     curr = curr->parent;
                  }
                  avg = (double) avg / (double) i;
                  assert(avg >= 1);
                  if( fixed )
                  {
                     costrev[e] = 1e+10;
                     cost[e + 1] = 1e+10;
                  }
                  else
                  {
                     if( modcost )
                     {
                        mult = costMultiplier(scip, heurdata, avg);
                        cost[e + 1] = cost[e + 1] * mult;
                     }

                     costrev[e] = cost[e + 1];
                     costrev[e + 1] = cost[e];
                  }

                  fixed = FALSE;
                  curr = ancestors[e];
                  assert(curr != NULL);
                  avg = 0.0;
                  i = 0;
                  while( curr != NULL )
                  {
                     i++;
                     avg += edgeweight[curr->index];
                     if( SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
                     {
                        fixed = TRUE;
                     }
                     curr = curr->parent;
                  }
                  avg = (double) avg / (double) i;
                  assert(avg >= 1);
                  if( fixed )
                  {
                     costrev[e + 1] = 1e+10;
                     cost[e] = 1e+10;
                  }
                  else
                  {
                     if( modcost )
                     {
                        mult = costMultiplier(scip, heurdata, avg);
                        cost[e] = cost[e] * mult;
                     }

                     costrev[e] = cost[e + 1];
                     costrev[e + 1] = cost[e];
                  }
               }
            }
            /* init shortest path algorithm */
            graph_path_init(solgraph);

            /* set (edge) result array to default */
            for( e = 0; e < nsoledges; e++ )
               results[e] = UNKNOWN;
            /* run TM heuristic */
            SCIP_CALL( do_layer(scip, tmheurdata, solgraph, NULL, &best_start, results, heurdata->ntmruns,
                  solgraph->source[0], cost, costrev, maxcost, &success) );

            if( !success )
            {
                  graph_path_exit(solgraph);
                  SCIPfreeBufferArrayNull(scip, &results);
                  SCIPfreeMemoryArray(scip, &edgeancestor);
                  SCIPfreeMemoryArray(scip, &edgeweight);
                  SCIPfreeBufferArrayNull(scip, &cost);
                  SCIPfreeBufferArrayNull(scip, &costrev);
                  graph_free(scip, solgraph, TRUE);
                  continue;
	    }
	    assert(graph_valid(solgraph));
	    assert(graph_sol_valid(solgraph, results));

            /* run local heuristic */
            if( solgraph->stp_type != STP_HOP_CONS && solgraph->stp_type != STP_MAX_NODE_WEIGHT && graph->stp_type != STP_DEG_CONS )
               SCIP_CALL( do_local(scip, solgraph, cost, costrev, results) );

	    if( graph->stp_type == STP_UNDIRECTED )
               assert(graph_sol_valid(solgraph, results));

            graph_path_exit(solgraph);
         }

         for( i = 0; i < nedges; i++ )
            orgresults[i] = UNKNOWN;
         artroot = graph->source[0];

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, nedges) );

         for( i = 0; i < nnodes; i++ )
            stnodes[i] = FALSE;

         /* retransform solution found by TM heuristic */
         if( solgraph->terms > 1 )
         {
            for( e = 0; e < nsoledges; e++ )
            {
               if( results[e] == CONNECT )
               {
                  /* iterate through list of ancestors */
                  curr = ancestors[e];
                  while( curr != NULL )
                  {
                     i = edgeancestor[curr->index];
                     if( graph->stp_type == STP_DEG_CONS )
                        orgresults[i] = CONNECT;
                     head = graph->head[i];
                     tail = graph->tail[i];
                     stnodes[tail] = TRUE;
                     stnodes[head] = TRUE;
#if 0
                     if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
                     {
                        if( head == graph->source[0] && !Is_term(graph->term[tail]) )
                        {
                           artroot = tail;
                        }
                        else if( tail == graph->source[0] && !Is_term(graph->term[head]) )
                        {
                           artroot = head;
                        }
                     }
#endif
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
	    if( graph->stp_type == STP_DEG_CONS )
	       orgresults[i] = CONNECT;
            head = graph->head[i];
            tail = graph->tail[i];
            stnodes[tail] = TRUE;
            stnodes[head] = TRUE;
#if 0
            if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
            {
               if( head == graph->source[0] && !Is_term(graph->term[tail]) )
               {
                  artroot = tail;
               }
               else if( tail == graph->source[0] && !Is_term(graph->term[head]) )
               {
                  artroot = head;
               }
            }
#endif
            curr = curr->parent;
         }
         /* prune solution (in the original graph) */
         if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
            SCIP_CALL( do_pcprune(scip, graph, graph->cost, orgresults, artroot, stnodes) );
         else if( graph->stp_type == STP_DEG_CONS )
	    SCIP_CALL( do_degprune(scip, graph, orgresults, stnodes) );
	 else
            SCIP_CALL( do_prune(scip, graph, graph->cost, 0, orgresults, stnodes) );

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
         /*
           printf("newcost: %f \n", pobj + SCIPprobdataGetOffset(scip));
           printf("oldcost: %f index: %d \n ", SCIPgetSolOrigObj(scip, newsol), SCIPsolGetIndex(newsol));
         */
         if( SCIPisGT(scip, SCIPgetSolOrigObj(scip, newsol) - SCIPprobdataGetOffset(scip), pobj) )
         {

            sol = NULL;
            SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

            if( success )
            {
               *result = SCIP_FOUNDSOL;
               solfound = TRUE;
               //printf("success in REC!!, cost: %f count: %d\n", pobj + SCIPprobdataGetOffset(scip), count);

               nsols = SCIPgetNSols(scip);
               assert(nsols > 0);
               sols = SCIPgetSols(scip);
               index = 0;
               for( i = 1; i < nsols; i++ )
                  if( SCIPisGT(scip, SCIPgetSolTime(scip, sols[i]), SCIPgetSolTime(scip, sols[index])) )
                     index = i;
               newsol = sols[index];
	       assert(graph_sol_valid(graph, orgresults));
            }
            else
            {
               //printf("NO success in REC \n");
               count++;
            }
         }
         else
         {
            count++;
	    //printf("NO success in REC \n");
         }

         //assert(graph_sol_valid(graph, orgresults));
         /* free memory */
         SCIPfreeBufferArrayNull(scip, &results);
         SCIPfreeMemoryArray(scip, &edgeancestor);
         SCIPfreeMemoryArray(scip, &edgeweight);

         SCIPfreeBufferArrayNull(scip, &cost);
         SCIPfreeBufferArrayNull(scip, &costrev);
         graph_free(scip, solgraph, TRUE);

      }
      else
      {
         count++;
         //printf("no success in merge \n");
      }
   }

   if( *result == SCIP_FOUNDSOL )
      heurdata->nfailures = 0;
   else
      heurdata->nfailures++;

   heurdata->lastsolindex = lastsolindex;
   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));
   heurdata->nlastsols = SCIPgetNSolsFound(scip);
   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &orgresults);
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
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRec) );

   /* add rec primal heuristic parameters */
   /*   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
        "number of nodes added to the contingent of the total nodes",
        &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

        SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
        "maximum number of nodes to regard in the subproblem",
        &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

        SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
        "minimum number of nodes required to start the subproblem",
        &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nusedsols",
         "number of solutions to be taken into account",
         &heurdata->nusedsols, FALSE, DEFAULT_NUSEDSOLS, 2, INT_MAX, NULL, NULL) );

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
   /*
     SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
     "number of nodes without incumbent change that heuristic should wait",
     &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

     SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
     "contingent of sub problem nodes in relation to the number of nodes of the original problem",
     &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

     SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
     "minimum percentage of integer variables that have to be fixed",
     &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

     SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
     "factor by which Rec should at least improve the incumbent",
     &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

     SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/lplimfac",
     "factor by which the limit on the number of LP depends on the node limit",
     &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

     SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/timelimit",
     "time limit for the sub problem to be solved, problem specifically chosen on default",
     &heurdata->timelimit, TRUE, DEFAULT_TIMELIMIT, -1.0, SCIP_REAL_MAX, NULL, NULL) );


     SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/dontwaitatroot",
     "should the nwaitingnodes parameter be ignored at the root node?",
     &heurdata->dontwaitatroot, TRUE, DEFAULT_DONTWAITATROOT, NULL, NULL) );

     SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
     "should subproblem be created out of the rows in the LP rows?",
     &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

     SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
     "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
     &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

     SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/permute",
     "should the subproblem be permuted to increase diversification?",
     &heurdata->permute, TRUE, DEFAULT_PERMUTE, NULL, NULL) ); */
   return SCIP_OKAY;
}
