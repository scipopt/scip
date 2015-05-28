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

/**@file   heur_rs.c
 * @brief  RS primal heuristic
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_rs.h"
#include "probdata_stp.h"
#include "grph.h"
#include "portab.h"


#define HEUR_NAME             "RS"
#define HEUR_DESC             "V. J. Rayward-Smith primal heuristic"
#define HEUR_DISPCHAR         '#'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

typedef struct component
{
   int    start;
   int    knot;
   SCIP_Real dist;
} COMP;

#define OUTPUT

#define COMP_MINI     0
#define COMP_LAST    -1
#define COMP_NONE    -2


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint nruns;
};


/*
 * Local methods
 */

static
int comp_cmp1(
   const void* a,
   const void* b
   )
{
   const SCIP_Real da = ((const COMP*)a)->dist;
   const SCIP_Real db = ((const COMP*)b)->dist;

   if (EQ(da, db))
      return(0);

   return ((da < db) ? -1 : 1);
}

static
int comp_cmp2(
   const void* a,
   const void* b
   )
{
   const SCIP_Real da = ((const COMP*)a)->dist;
   const SCIP_Real db = ((const COMP*)b)->dist;

   if (EQ(da, db))
      return(0);

   return ((db < da) ? -1 : 1);
}

inline static
SCIP_Real evaluate_f(
   int         r,
   const COMP* tree
   )
{
   SCIP_Real ret = 0;
   int i;

   assert(r    >= 2);
   assert(tree != NULL);

   for( i = 0; i < r; i++ )
      ret += tree[i].dist;

   return (ret / (r - 1));
}

static
SCIP_Real distance(
   const GRAPH* g,
   const PATH*  path,
   int          k,
   int          trees,
   COMP*        tree,
   const int*   next
   )
{
   SCIP_Real fval;
   int    i;
   int    j;

   assert(g       != NULL);
   assert(path    != NULL);
   assert(trees   >= 2);
   assert(trees   <= g->terms);
   assert(tree    != NULL);
   assert(next    != NULL);

   for(i = 0; i < trees; i++)
   {
      tree[i].dist = FARAWAY;

      for(j = tree[i].start; j >= COMP_MINI; j = next[j])
         if (GT(tree[i].dist, path[j].dist))
            tree[i].dist = path[j].dist;
   }
   qsort(tree, (size_t)trees, sizeof(COMP), comp_cmp1);

   i = 2;

   do
   {
      fval = evaluate_f(i, tree);
      i++;
   }
   while((i <= trees) && LT(tree[i].dist, fval));

   return fval;
}

static
int find_best(
   const GRAPH* g,
   PATH*        path,
   int          trees,
   COMP*        tree,
   const int*   next
   )
{
   SCIP_Real fval;
   SCIP_Real min_fval;
   int    min_knot;
   int    i;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(trees  <= g->terms);
   assert(tree   != NULL);
   assert(next   != NULL);

   min_fval = FARAWAY;
   min_knot = -1;

   for(i = 0; i < g->knots; i++)
   {
      graph_path_exec(g, FSP_MODE, i, g->cost, path);

      fval = distance(g, path, i, trees, tree, next);

      if (LT(fval, min_fval))
      {
         min_fval = fval;
         min_knot = i;
      }
   }
   assert(LT(min_fval, FARAWAY));
   assert(min_knot > -1);
   assert(min_knot < g->knots);

   return min_knot;
}

static
int connect2(
   const GRAPH* g,
   PATH*        path,
   int          knot,
   int          trees,
   COMP*        tree,
   int*         next
   )
{
   int i;
   int j;
   int e;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(knot   >= 0);
   assert(knot   <  g->knots);
   assert(trees  >= 2);
   assert(trees  <= g->terms);
   assert(tree   != NULL);
   assert(next   != NULL);

   graph_path_exec(g, FSP_MODE, knot, g->cost, path);

   for(i = 0; i < trees; i++)
   {
      tree[i].dist = FARAWAY;
      tree[i].knot = COMP_NONE;

      for(j = tree[i].start; j >= COMP_MINI; j = next[j])
      {
         if (GT(tree[i].dist, path[j].dist))
         {
            tree[i].dist = path[j].dist;
            tree[i].knot = j;
         }
      }
      assert(LT(tree[i].dist, FARAWAY));
      assert(tree[i].knot       >= COMP_MINI);
      assert(next[tree[i].knot] != COMP_NONE);
   }
   qsort(tree, (size_t)trees, sizeof(COMP), comp_cmp2);

   /* Beiden Baeume zusammen haengen
    */
   trees--;

   for(i = tree[trees - 1].start; next[i] >= COMP_MINI; i = next[i]);

   assert(next[i] == COMP_LAST);

   next[i] = tree[trees].start;

   /* Den Weg zwischen den letzten beiden Komponenten via 'knot' ermitteln
    */
   for(j = trees; j > trees - 2; j--)
   {
      i = tree[j].knot;

      assert(next[i] != COMP_NONE);

      while(i != knot)
      {
         e = path[i].edge;

         assert(e >= 0);
         assert(e < g->edges);

         i = g->tail[e];

         assert(i == g->tail[e]);
         assert(i >= 0);
         assert(i <  g->knots);

         if (next[i] == COMP_NONE)
         {
            assert(next[i] == COMP_NONE);

            /* Knoten zur Komponente hinzufuegen
             */
            next[i]               = tree[trees - 1].start;
            tree[trees - 1].start = i;
         }
      }
   }

   return trees;
}



/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRS)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRS(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRS)
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
SCIP_DECL_HEURINIT(heurInitRS)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->nruns = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitRS)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of RS primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitRS NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolRS)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of RS primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolRS NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolRS)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of RS primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolRS NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRS)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   GRAPH* graph;
   PATH* path;
   COMP* tree;
   SCIP_Real* nval;
   int* results;
   int* next;
   SCIP_Real pobj = 0.0;
   int nvars;
   int layer = 0;
   int knot;
   int count;
   int trees  = 0;
   int i;
   int j;

   /* xval wird momentan noch nicht gebraucht. */

   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;
   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   assert(graph->layers  == 1);

   SCIPinfoMessage(scip, NULL, "Heuristic (Rayward-Smith) Start\n");

   nvars = SCIPprobdataGetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &path, graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &next, graph->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree, graph->terms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &results, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

   *result = SCIP_DIDNOTFIND;

   for( i = 0; i < graph->knots; i++ )
   {
      if( !Is_term(graph->term[i]) )
         next[i] = COMP_NONE;
      else
      {
         tree[trees].start = i;
         trees++;
         next[i] = COMP_LAST;
      }
   }

   for( i = 0; i < graph->knots; i++ )
      graph->mark[i] = TRUE;

   /* Hauptschleife...
    */
   while( trees > 1 )
   {
      SCIPinfoMessage(scip, NULL, "X");

      knot = find_best(graph, path, trees, tree, next);

      trees = connect2(graph, path, knot, trees, tree, next);

#ifdef OUTPUT
      for( i = 0; i < trees; i++ )
      {
         SCIPinfoMessage(scip, NULL, "Tree: %d = { ", i);
         for( j = tree[i].start; j >= COMP_MINI; j = next[j] )
            SCIPinfoMessage(scip, NULL, "%d ", j);
         SCIPinfoMessage(scip, NULL, "}\n");
      }
      SCIPinfoMessage(scip, NULL, "\n");
#endif
   }
   for( i = 0; i < graph->knots; i++ )
      graph->mark[i] = (next[i] != COMP_NONE);

   graph_path_exec(graph, MST_MODE, graph->source[layer], graph->cost, path);

   for( i = 0; i < graph->edges; i++ )
      results[i] = -1;

   for( i = 0; i < graph->knots; i++ )
   {
      if( graph->mark[i] && (path[i].edge != -1) )
      {
         assert(graph->head[path[i].edge] == i);
         assert(results[path[i].edge] == -1);

         results[path[i].edge] = layer;
      }
   }

   /* Baum beschneiden */
   do
   {
      SCIPinfoMessage(scip, NULL, "C");

      count = 0;

      for( i = 0; i < graph->knots; i++ )
      {
         if( !graph->mark[i] )
            continue;

         if( graph->term[i] == layer )
            continue;

         for( j = graph->outbeg[i]; j != EAT_LAST; j = graph->oeat[j] )
            if( results[j] == layer )
               break;

         if( j == EAT_LAST )
         {
            /* Es muss genau eine eingehende Kante geben */
            for( j = graph->inpbeg[i]; j != EAT_LAST; j = graph->ieat[j] )
            {
               if (results[j] == layer)
               {
                  results[j]    = -1;
                  graph->mark[i]   = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while(count > 0);

   for( i = 0; i < nvars; i++ )
      nval[i] = (results[i % graph->edges] == (i / graph->edges)) ? 1.0 : 0.0;

   SCIPinfoMessage(scip, NULL, "Validation\n");

   if( validate(graph, nval) )
   {
      SCIP_VAR** vars;
      SCIP_Bool success;

      SCIPinfoMessage(scip, NULL, "Is Valid\n");

      pobj = 0.0;

      for( i = 0; i < nvars; i++ )
         pobj += graph->cost[i % graph->edges] * nval[i];

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
   SCIPfreeBufferArray(scip, &tree);
   SCIPfreeBufferArray(scip, &next);
   SCIPfreeBufferArray(scip, &path);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the RS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create RS primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRS, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRS) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRS) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRS) );
#if 0
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRS) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRS) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRS) );
#endif

   /* add RS primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
