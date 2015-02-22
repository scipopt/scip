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

#define DEFAULT_MAXNODES      5000LL         /* maximum number of nodes to be regarded in the subproblem                   */
#define DEFAULT_MINIMPROVE    0.001          /* factor by which Rec should at least improve the incumbent       */
#define DEFAULT_MINNODES      50LL           /* minimum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINFIXINGRATE 0.6            /* minimum percentage of variables to be fixed                           */
#define DEFAULT_NODESOFS      500LL          /* number of nodes added to the contingent of the total nodes            */
#define DEFAULT_NODESQUOT     0.1            /* subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_LPLIMFAC      2.0            /* factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_NUSEDSOLS     4              /* number of solutions that will be taken into account                   */
#define DEFAULT_RANDSEED       0              /* random seed                                                               */
#define DEFAULT_NWAITINGNODES 200LL          /* number of nodes without incumbent change heuristic should wait        */
#define DEFAULT_NTMRUNS       100            /**< number of runs in TM heuristic        */
#define DEFAULT_TIMELIMIT     -1             /* time limit for the sub problem to be solved, if -1, running time till first call of rec is chosen */
#define DEFAULT_NWAITINGSOLS  2              /* minimum number of new solutions to be available before executing the heuristic again  */
#define DEFAULT_DONTWAITATROOT FALSE         /* should the nwaitingnodes parameter be ignored at the root node?       */
#define DEFAULT_USELPROWS     FALSE          /* should subproblem be created out of the rows in the LP rows,
                                              * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      FALSE           /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                               * cutpool of the original scip be copied to constraints of the subscip
                                               */
#define DEFAULT_PERMUTE       FALSE          /* should the subproblem be permuted to increase diversification?        */
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
   SCIP_SOL*             prevlastsol;        /**< worst solution taken into account during the previous run         */
   SCIP_SOL*             prevbestsol;        /**< best solution during the previous run                             */
   int                   prevnsols;          /**< number of all solutions during the previous run                   */

   SCIP_Longint          ncalls;
   SCIP_Longint          nlastsols;          /**< number of solutions available during the last run                 */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by rec in earlier calls                  */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */
   int                   ntmruns;            /**< number of runs in TM heuristic        */
   int                   nusedsols;          /**< number of solutions that will be taken into account               */
   int                   nselectedsols;      /**< number of solutions actually selected */
   int                   nwaitingsols;       /**< number of new solutions before executing the heuristic again      */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             minimprove;         /**< factor by which Rec should at least improve the incumbent   */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Real             timelimit;          /**< time limit for the subproblem                                     */
   SCIP_Bool             dontwaitatroot;     /**< should the nwaitingnodes parameter be ignored at the root node?   */
   unsigned int          randseed;           /**< seed value for random number generator                            */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable used to store the solution tuples already used          */
   SOLTUPLE*             lasttuple;          /**< last tuple of solutions created by rec                      */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?                                     */
   SCIP_Bool             permute;            /**< should the subproblem be permuted to increase diversification?    */
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
/* select solutions */
static
SCIP_RETCODE selectsols(
   SCIP*       scip,
   SCIP_HEURDATA*     heurdata,           /**< primal heuristic data */
   int*        selection
   )
{
   SCIP_SOL** sols;                          /* array of all solutions found so far         */
   SCIP_SOL** permsols;
   // SCIP_Real* soltimes;
   int i;
   //int j;
   // int min;
   int end;
   //  int start;
   int nselectedsols = 0;
   int nsols;                                /* number of all solutions found so far        */
   int shift;
   int nusedsols;                            /* number of solutions to use in rec     */
   int* perm;
   int* solselected;

   /* get solution data */
   nsols = SCIPgetNSols(scip);
   nusedsols = heurdata->nusedsols;
   assert(nusedsols > 1);
   assert(nsols >= nusedsols);
   sols = SCIPgetSols(scip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &permsols, SCIPgetSols(scip), nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsols) );

   for( i = 0; i < nsols; i++ )
   {
      perm[i] = i;
      solselected[i] = FALSE;
   }

   if( SCIPgetRandomInt(1, 2, &(heurdata->randseed)) != 1 )
   {
      solselected[0] = TRUE;
      selection[nselectedsols++] = 0;
      //printf("0select sol found by: %s \n", SCIPheurGetName(SCIPsolGetHeur(SCIPgetSols(scip)[0])) );
      //printf("==2\n");
   }

   end = (nusedsols - nusedsols / 2) + (int) ((SCIPgetRandomReal(0.8, 1.5, &(heurdata->randseed))) );

   shift = SCIPgetRandomInt(1, (nusedsols - nusedsols / 2) + (double) nsols / 5.0 , &(heurdata->randseed));

   if( end + shift > nsols )
      shift = nsols - end;
   //printf("shift: %d \n", shift);
   SCIPpermuteIntArray(perm, nselectedsols, end + shift, &(heurdata->randseed));

   for( i = nselectedsols; i < end; i++ )
   {
      selection[nselectedsols++] = perm[i];
      solselected[perm[i]] = TRUE;
      //printf("select sol %d found by: %s \n", perm[i], SCIPheurGetName(SCIPsolGetHeur(SCIPgetSols(scip)[perm[i]])) );
   }

   SCIPpermuteIntArray(perm, 0, nsols, &(heurdata->randseed));
   for( i = 0; i < nsols && nselectedsols < heurdata->nusedsols; i++ )
   {
      if( solselected[perm[i]] == FALSE )
      {
	 if( strcmp(SCIPheurGetName(SCIPsolGetHeur(sols[perm[i]])), "rec") == 0 )
            continue;
	 //printf("p select sol %d found by: %s \n", perm[i], SCIPheurGetName(SCIPsolGetHeur(SCIPgetSols(scip)[perm[i]])) );
         selection[nselectedsols++] = perm[i];
      }
   }

   heurdata->nselectedsols = nselectedsols;

   SCIPfreeBufferArray(scip, &solselected);
   SCIPfreeBufferArray(scip, &permsols);
   SCIPfreeBufferArray(scip, &perm);
   return SCIP_OKAY;

   return SCIP_OKAY;
}

/* merge selected solutions to a new graph */
static
SCIP_RETCODE buildsolgraph(
   SCIP*       scip,
   SCIP_HEURDATA*     heurdata,           /**< primal heuristic data */
   GRAPH*      graph,
   GRAPH**     solgraph,
   int**       edgeancestor,
   int**       edgeweight
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
   char*  solnode;            /**< marks nodes contained in at least one solution */
   char*  soledge;            /**< marks edges contained in at least one solution */

   assert(scip != NULL);
   assert(graph != NULL);

   sols = SCIPgetSols(scip);
   nedges = graph->edges;
   nnodes = graph->knots;
   nsoledges = 0;
   nsolnodes = 0;

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
   SCIP_CALL( selectsols(scip, heurdata, solselection) );

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
   if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
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

   /* initialize new graph */
   newgraph = graph_init(nsolnodes, 2 * nsoledges, 1, 0);
   if( graph->stp_type == STP_GRID )
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
            graph_knot_add(newgraph, 0, -1, -1);
         else
            graph_knot_add(newgraph, -1, -1, -1);
      }
   }
   /* set root */
   newgraph->source[0] = dnodemap[graph->source[0]];
   assert(newgraph->source[0] >= 0);

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
   heurdata->usednodes = 0;
   heurdata->prevlastsol = NULL;
   heurdata->prevbestsol = NULL;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#else
   heurdata->randseed = 0;
#endif

   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;
   heurdata->prevnsols = 0;
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
   GRAPH* graph;                             /* graph structure */
   GRAPH* solgraph;                             /* graph structure */
   SCIP_SOL* sol;                            /* new solution */
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real* nval;
   SCIP_Real pobj;
   SCIP_Real avg;
   SCIP_Real maxcost = 0.0;
   SCIP_Bool success;
   SCIP_Bool fixed;
   IDX* curr;
   IDX** ancestors;
   int i;
   int e;
   int nsols;                                /* number of all solutions found so far */
   int nedges;
   int nnodes;
   int nsoledges;
   int nheurs;
   int nusedsols;
   int best_start;
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
   nusedsols = heurdata->nusedsols;

   /* only call heuristic, if sufficiently many solutions are available */
   if( SCIPgetNSols(scip) < nusedsols )
   {
      return SCIP_OKAY;
   }
   // if( graph->stp_type != STP_UNDIRECTED )
   //  return SCIP_OKAY;

   /* suspend heuristic? */
   if( SCIPisLT(scip, nsols, heurdata->nlastsols + heurdata->nwaitingsols + heurdata->nfailures) && heurdata->ncalls > 0 )
   {
      // printf("SUSPENDING REC \n");
      return SCIP_OKAY;
   }
   heurdata->ncalls++;
   heurdata->nlastsols = SCIPgetNSolsFound(scip);

   *result = SCIP_DIDNOTRUN;

   /* build up a new graph, consisting of several solutions */
   SCIP_CALL( buildsolgraph(scip, heurdata, graph, &solgraph, &edgeancestor, &edgeweight) );

   assert(graph_valid(solgraph));

   /* get TM heuristic data */
   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);
   for( i = 0; i < nheurs; i++ )
      if( strcmp(SCIPheurGetName(heurs[i]), "TM") == 0 )
         break;

   assert(i < nheurs);
   tmheurdata = SCIPheurGetData(heurs[i]);

   /* presolve new graph */

   /* init shortest path algorithm */
   graph_path_init(solgraph);

   /* reduce new graph */
   (void) reduce(solgraph, 4, scip);

   graph_path_exit(solgraph);

   solgraph = graph_pack(solgraph);

   /* if graph reduction solved the whole problem, solgraph has only one node */
   if( solgraph->knots > 1 )
   {
      int factor1 = 1;
      int factor2 = 1;
      nsoledges = solgraph->edges;
      ancestors = solgraph->ancestors;

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &results, nsoledges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cost, nsoledges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nsoledges) );
      /*  for( i = 0; i < nsoledges; i++ )
          {
	  assert(
	  cost[i] = solgraph->cost[i] / (1 + edgeweight[i] - 1
          }*/
      BMScopyMemoryArray(cost, solgraph->cost, nsoledges);
      //factor = SCIPgetRandomInt(2, 3, &(heurdata->randseed));
      /* hop constraint problem? */
      if( graph->stp_type == STP_HOP_CONS )
      {
         maxcost = 0.0;
	 factor1 = SCIPgetRandomInt(140, 200, &(heurdata->randseed));
	 factor2 = SCIPgetRandomInt(30, 60, &(heurdata->randseed));
         for( e = 0; e < nsoledges; e++)
         {
            /* TODO */

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
                  //break;
               }
               curr = curr->parent;
            }
            avg = (double) avg / (double) i;
	    assert(avg >= 1);
            if( fixed )
               cost[e] = 1e+10;
            else if( SCIPisLT(scip, avg, 3) )
               cost[e] = cost[e] * (double) factor1 * (1.0 / avg);
	    else if( SCIPisLT(scip, avg, 4) )
               cost[e] = cost[e] * (double) factor2 * (1.0 / avg);

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
	    factor1 = SCIPgetRandomInt(140, 200, &(heurdata->randseed));
	    factor2 = SCIPgetRandomInt(30, 60, &(heurdata->randseed));
            fixed = FALSE;
            curr = ancestors[e + 1];
	    avg = 0.0;
	    i = 0;
            while( curr != NULL )
            {
	       i++;
	       avg += edgeweight[curr->index];
               if(  SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
               {
                  fixed = TRUE;
                  //break;
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
	       if( SCIPisLT(scip, avg, 2) )
		  cost[e] = cost[e] * (double) factor1 * 4 * (1.0 / avg);
	       else if( SCIPisLT(scip, avg, 3) )
		  cost[e + 1] = cost[e + 1] * (double) factor1 * (1.0 / avg);
	       else if( SCIPisLT(scip, avg, 4) )
		  cost[e + 1] = cost[e + 1] * (double) factor2 * (1.0 / avg);

               costrev[e] = cost[e + 1];
               costrev[e + 1] = cost[e];
            }

            fixed = FALSE;
            curr = ancestors[e];
	    avg = 0.0;
	    i = 0;
            while( curr != NULL )
            {
	       i++;
	       avg += edgeweight[curr->index];
               if( SCIPvarGetUbGlobal(vars[edgeancestor[curr->index]] ) < 0.5 )
               {
                  fixed = TRUE;
                  //break;
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
	       if( SCIPisLT(scip, avg, 2) )
		  cost[e] = cost[e] * (double) factor1 * 4 * (1.0 / avg);
               else if( SCIPisLT(scip, avg, 3) )
		  cost[e] = cost[e] * (double) factor1 * (1.0 / avg);
	       else if( SCIPisLT(scip, avg, 4) )
		  cost[e] = cost[e] * (double) factor2 * (1.0 / avg);

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
      SCIP_CALL( do_layer(scip, tmheurdata, solgraph, &best_start, results, heurdata->ntmruns, solgraph->source[0], cost, costrev, maxcost) );

      assert(graph_sol_valid(solgraph, results));

      /* run local heuristic */
      if( solgraph->stp_type == STP_UNDIRECTED )
         SCIP_CALL( do_local(scip, solgraph, cost, costrev, results) );

      assert(graph_sol_valid(solgraph, results));

      graph_path_exit(solgraph);
   }
   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgresults, nnodes) );

   for( i = 0; i < nedges; i++ )
      orgresults[i] = UNKNOWN;

   if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
   {
      if( solgraph->knots > 1 )
      {
         /* retransform solution found by TM heuristic */
         for( e = 0; e < nsoledges; e++ )
         {
            if( results[e] == CONNECT )
            {
               /* iterate through list of ancestors */
               curr = ancestors[e];
               while( curr != NULL )
               {
                  orgresults[edgeancestor[curr->index]] = CONNECT;
                  curr = curr->parent;
               }
            }
         }
      }

      /* retransform edges fixed during graph reduction */
      while( curr != NULL )
      {
         orgresults[edgeancestor[curr->index]] = CONNECT;
         curr = curr->parent;
      }
   }
   else
   {
      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, nedges) );

      for( i = 0; i < nnodes; i++ )
         stnodes[i] = FALSE;

      /* retransform solution found by TM heuristic */
      if( solgraph->knots > 1 )
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
                  stnodes[graph->tail[i]] = TRUE;
                  stnodes[graph->head[i]] = TRUE;
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
         stnodes[graph->tail[i]] = TRUE;
         stnodes[graph->head[i]] = TRUE;
         curr = curr->parent;
      }

      /* prune solution (in the original graph) */
      SCIP_CALL( do_prune(scip, graph, graph->cost, 0, orgresults, stnodes) );

      SCIPfreeBufferArray(scip, &stnodes);
   }
   pobj = 0.0;
   for( e = 0; e < graph->edges; e++ )
   {
      if( orgresults[e] == CONNECT )
      {
         nval[e] = 1.0;
         pobj += graph->cost[i];
      }
      else
      {
         nval[e] = 0.0;
      }
   }
   assert(graph_sol_valid(graph, orgresults));

   SCIPfreeBufferArray(scip, &orgresults);

   /* printf("[%d] = %d \n", e, edgeancestor[e] );
    */
   //printf("COST: %f\n", pobj);
   //SCIP_CALL( printGraph(scip, graph, "OrgRec.gml", orgresults) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   SCIPfreeSol(scip, &sol);
   /* for( e = 0; e < nedges; e++ )
      {
      printf("val: %d : %e \n", e, nval[e]);
      }*/

   SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

   if( success )//&& (SCIPgetNSols(scip) < 3)? TRUE : SCIPisGE(scip, SCIPsolGetOrigObj(SCIPgetSols(scip)[2]),  SCIPsolGetOrigObj(sol))  )
   {
      *result = SCIP_FOUNDSOL;
      //printf("success in REC!! \n\n");
      heurdata->nfailures = 0;
   }
   else
   {
      // printf("NO success \n");
      heurdata->nfailures++;
   }

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &results);
   SCIPfreeMemoryArray(scip, &edgeancestor);
   SCIPfreeMemoryArray(scip, &edgeweight);
   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArrayNull(scip, &cost);
   SCIPfreeBufferArrayNull(scip, &costrev);
   graph_free(solgraph, TRUE);

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
