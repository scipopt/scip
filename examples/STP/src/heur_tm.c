/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_tm.c
 * @brief  TM primal heuristic
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_tm.h"
#include "probdata_stp.h"
#include "grph.h"
#include "portab.h"


#define HEUR_NAME             "TM"
#define HEUR_DESC             "takahashi matsuyama primal heuristic for steiner trees"
#define HEUR_DISPCHAR         '+'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_EVALRUNS 10
#define DEFAULT_INITRUNS 10
#define DEFAULT_LEAFRUNS 10
#define DEFAULT_ROOTRUNS 50
#define DEFAULT_DURINGLPFREQ 10


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint ncalls;
   int evalruns;
   int initruns;
   int leafruns;
   int rootruns;
   int duringlpfreq;
};


/*
 * Local methods
 */

/* Die Heuristic stoert sich nicht dran, wenn sie einzelne Wege nicht
 * routen kann, sondern erklaert das jeweilige Netz einfach fuer fertig.
 * Die Loesung muss also ueberprueft werden.
 */
static
SCIP_RETCODE do_heuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  g,
   int           layer,
   int*          result,
   int           start,
   char*         connected,
   const SCIP_Real* cost,
   PATH**        path
   )
{
   PATH*  mst;
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    count;
   SCIP_Real min;
   int    i;
   int    j;
   int    old;
   int    newval;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(path      != NULL);
   assert(layer >= 0 && layer < g->layers);

   SCIPinfoMessage(scip, NULL, "Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, g->knots) );

   for( i = 0; i < g->knots; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;
   cluster[csize++] = start;

   /* CONSTCOND */
   for(;;)
   {
      /* Suche das Terminal, das am dichtesten dran ist
       */
      min = FARAWAY;
      old = -1;
      newval = -1;

      for(i = 0; i < g->knots; i++)
      {
         if (g->grad[i] == 0)
            continue;

         if (g->term[i] != layer)
            continue;

         if (connected[i])
            continue;

         /* Jetzt brauchen wir die Entfernungen.
          */
         if (path[i] == NULL)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(path[i]), g->knots) );

            assert(path[i] != NULL);

            /* ! Ob das die guenstiges Richtung ist die Wege zu berechnen, wenn
             * ! die Kosten fuer Hin und Rueckweg unterschiedlich sind ist doch
             * ! sehr fraglich.
             * ! Koennte aber sein, weil wir die Kosten unten umgedreht haben.
             */
            graph_path_exec(g, FSP_MODE, i, cost, path[i]);
         }
         for(k = 0; k < csize; k++)
         {
            j = cluster[k];

            assert(i != j);
            assert(connected[j]);

            if (LT(path[i][j].dist, min))
            {
               min = path[i][j].dist;
               newval = i;
               old = j;
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if (newval == -1)
         break;

      /* Weg setzten
       */
      assert((old > -1) && (newval > -1));
      assert(path[newval] != NULL);
      assert(path[newval][old].dist < FARAWAY);
      assert(g->term[newval] == layer);
      assert(!connected[newval]);
      assert(connected[old]);

      fputc('R', stdout);
      fflush(stdout);

      /*    printf("Connecting Knot %d-%d dist=%d\n", newval, old, path[newval][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;

      while(k != newval)
      {
         e = path[newval][k].edge;
         k = g->tail[e];

         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }
   SCIPfreeBufferArray(scip, &cluster);

   fputc('M', stdout);
   fflush(stdout);

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, g->knots) );

   /* MST berechnen
    */
   for(i = 0; i < g->knots; i++)
      g->mark[i] = connected[i];

   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  g->knots);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, mst);

   for(i = 0; i < g->knots; i++)
   {
      if (connected[i] && (mst[i].edge != -1))
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = layer;
      }
   }

   /* Baum beschneiden
    */
   do
   {
      fputc('C', stdout);
      fflush(stdout);

      count = 0;

      for(i = 0; i < g->knots; i++)
      {
         if (!g->mark[i])
            continue;

         if (g->term[i] == layer)
            continue;

         for(j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j])
            if (result[j] == layer)
               break;

         if (j == EAT_LAST)
         {
            /* Es muss genau eine eingehende Kante geben
             */
            for(j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j])
            {
               if (result[j] == layer)
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while(count > 0);

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}

static
SCIP_RETCODE do_layer(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  graph,
   int           layer,
   int*          best_result,
   int           runs,
   const SCIP_Real* cost
   )
{
   PATH** path;
   char* connected;
   int* result;
   int* start;
   SCIP_Real obj;
   SCIP_Real min = FARAWAY;
   int best = -1;
   int k;
   int r;
   int e;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(best_result != NULL);
   assert(cost != NULL);
   assert(layer >= 0 && layer < graph->layers);

   SCIP_CALL( SCIPallocBufferArray(scip, &path, graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &connected, graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, graph->edges) );
#if 1
   BMSclearMemoryArray(path, graph->knots);
#else
   for( k = 0; k < graph->knots; k++ )
      path[k] = NULL;
#endif

   /* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    * Patch um die heuristic nach einem restruct starten zu koennen
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    */
   if( best >= graph->knots )
      best = -1;

   if( graph->layers > 1 )
   {
      SCIP_CALL( do_heuristic(scip, graph, layer, best_result, graph->source[layer], connected, cost, path) );
   }
   else
   {
      runs = (runs > graph->knots) ? graph->knots         : runs;
      best = (best < 0)            ? graph->source[layer] : best;

      for( k = 0; k < graph->knots; k++ )
      {
         assert(graph->grad[k] > 0);

         start[k] = k;
      }

      /* if we run over all nodes, we do not need to do the following */
      if( runs < graph->knots )
      {
         int random;
         int tmp;

         /* swap the starting values randomly */
         for( r = 0; r < runs; r++ )
         {
            random   = rand() % graph->knots;
            tmp      = start[r];
            start[r] = start[random];
            start[random] = tmp;
         }

         /* check if we have a best starting value */
         for( r = 0; r < runs; r++ )
            if( start[r] == best )
               break;

         /* do we need to set the start by hand */
         if( r == runs )
            start[0] = best;
      }

      for( r = 0; r < runs; r++ )
      {
         /* incorrekt if layers > 1 ! */
         assert(graph->layers == 1);

         for( e = 0; e < graph->edges; e++ )
            result[e] = -1;

         SCIP_CALL( do_heuristic(scip, graph, layer, result, start[r], connected, cost, path) );

         obj = 0.0;

         /* here we take another measure than in do_heuristic() */
         for( e = 0; e < graph->edges; e++)
            obj += (result[e] > -1) ? graph->cost[e] : 0.0;

         SCIPinfoMessage(scip, NULL, " Obj=%.12e\n", obj);

         if( LT(obj, min) )
         {
            min = obj;

            for( e = 0; e < graph->edges; e++ )
               best_result[e] = result[e];

            best = start[r];
         }
      }
   }

   SCIPinfoMessage(scip, NULL, "Freeing Memory\n");

   for( k = 0; k < graph->knots; k++ )
   {
      assert(path[k] == NULL || graph->term[k] == layer);
      SCIPfreeBufferArrayNull(scip, &(path[k]));
   }

   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &connected);
   SCIPfreeBufferArray(scip, &path);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTM(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->ncalls = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitTM NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolTM NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolTM NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTM)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   GRAPH* graph;
   SCIP_Real* cost;
   SCIP_Real* nval;
   SCIP_Real* xval;
   int* results;
   SCIP_Bool solcreated;
   SCIP_Real pobj;
   int nvars;
   int layer;
   int runs;
   int e;
   int v;

   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;
   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   runs = 0;

   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
      runs = heurdata->initruns;
   else if( ((heurtiming & SCIP_HEURTIMING_DURINGLPLOOP) && (heurdata->ncalls % heurdata->duringlpfreq == 0)) || (heurtiming & SCIP_HEURTIMING_AFTERLPLOOP) )
      runs = heurdata->evalruns;
   else if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      if( SCIPgetDepth(scip) == 0 )
         runs = heurdata->rootruns;
      else
         runs = heurdata->leafruns;
   }

   heurdata->ncalls++;

   if( runs == 0 )
      return SCIP_OKAY;

   SCIPinfoMessage(scip, NULL, "Heuristic Start\n");

   nvars = SCIPprobdataGetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &results, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

   *result = SCIP_DIDNOTFIND;
   solcreated = FALSE;

   /* */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      sol = NULL;
   else
   {
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      solcreated = TRUE;

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );
   }

   xval = SCIPprobdataGetXval(scip, sol);

   if( solcreated )
   {
      SCIPfreeSol(scip, &sol);
   }

   for( e = 0; e < graph->edges; e++ )
      results[e] = -1;

   for( layer = 0; layer < graph->layers; layer++ )
   {
      if( xval == NULL )
      {
         BMScopyMemoryArray(cost, graph->cost, graph->edges);
      }
      else
      {
         /* swap costs */
         for( e = 0; e < graph->edges; e += 2)
         {
            cost[e]     = ((1.0 - xval[layer * graph->edges + e + 1]) * graph->cost[e + 1]);
            cost[e + 1] = ((1.0 - xval[layer * graph->edges + e    ]) * graph->cost[e]);
         }
      }
      /* can we connect the network */
      SCIP_CALL( do_layer(scip, graph, layer, results, runs, cost) );

      /* take the path */
      if( graph->layers > 1 )
      {
         for( e = 0; e < graph->edges; e += 2)
         {
            if( (results[e] == layer) || (results[e + 1] == layer) )
               graph_edge_hide(graph, e);
         }
      }
   }

   if( graph->layers > 1 )
      graph_uncover(graph);

   for( v = 0; v < nvars; v++ )
      nval[v] = (results[v % graph->edges] == (v / graph->edges)) ? 1.0 : 0.0;

   SCIPinfoMessage(scip, NULL, "Validation\n");

   if( validate(graph, nval) )
   {
      SCIP_VAR** vars;
      SCIP_Bool success;

      SCIPinfoMessage(scip, NULL, "Is Valid\n");

      pobj = 0.0;

      for( v = 0; v < nvars; v++ )
         pobj += graph->cost[v % graph->edges] * nval[v];

      SCIPinfoMessage(scip, NULL, "pobj=%.12e\n", pobj);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      vars = SCIPprobdataGetVars(scip);
      assert(vars != NULL);

      /* store new solution value and decrease fractionality counter */
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, nval) );

      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, &success) );

      if( success )
         *result = SCIP_FOUNDSOL;
   }

   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &results);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the TM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create TM primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTM) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTM) );
#if 0
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitTM) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolTM) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolTM) );
#endif

   /* add TM primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   /* add TM primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/evalruns",
         "number of runs for eval",
         &heurdata->evalruns, FALSE, DEFAULT_EVALRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/initruns",
         "number of runs for init",
         &heurdata->initruns, FALSE, DEFAULT_INITRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/leafruns",
         "number of runs for leaf",
         &heurdata->leafruns, FALSE, DEFAULT_LEAFRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/rootruns",
         "number of runs for root",
         &heurdata->rootruns, FALSE, DEFAULT_ROOTRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/duringlpfreq",
         "frequency for calling heuristic during LP loop",
         &heurdata->duringlpfreq, FALSE, DEFAULT_DURINGLPFREQ, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
